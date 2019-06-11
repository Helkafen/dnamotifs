{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}

module Lib where

import qualified Data.Map as M
import qualified Data.List
import qualified Data.Set as Set
import qualified Data.Vector as V
import           Data.Vector.Algorithms.Merge (sort)
import qualified Data.Vector.Storable as STO
import qualified Data.ByteString as B
import qualified Data.Text as T
import qualified Data.Text.IO as TIO
--import           System.IO (appendFile)
import           TextShow (showt)
import           Data.Time.Clock.POSIX (getPOSIXTime, POSIXTime)


import Types
import Bed (readPeaks)
import Fasta (loadFasta)
import Vcf (readVcfWithGenotypes)
import PatternFind (findPatternsInBlock, mkPatterns, mkNucleotideAndPositionBlock, Patterns)


-- Inclusive [Start, End] interval
applyVariants :: (Position0 -> Position0 -> BaseSequencePosition) -> Position0 -> Position0 -> [Diff] -> BaseSequencePosition
applyVariants takeReferenceGenome start end allDiffs =
    let chunks = go start (filter (\(Diff p _ _) -> p >= start && p <= end) allDiffs)
    in BaseSequencePosition (B.concat (map seqOf chunks)) (STO.concat (map posOf chunks))
  where go :: Position0 -> [Diff] -> [BaseSequencePosition]
        go refPosition [] = [takeReferenceGenome refPosition end]
        go refPosition diffs@(Diff pos ref alt:vs)
            | pos > refPosition = takeReferenceGenome refPosition (pos-1): go pos diffs
            | pos == refPosition && B.length ref == 1 = -- SNV, insertion
                BaseSequencePosition alt (STO.replicate (B.length alt) (fromIntegral refPosition)): go (refPosition + 1) vs
            | pos == refPosition && B.length alt == 1 = -- deletion
                BaseSequencePosition alt (STO.singleton (fromIntegral refPosition)) : go (offSet refPosition (B.length ref)) vs
            | pos == refPosition = error "Missing case in applyVariants (we assume that ref or alt have length 1)"
            | refPosition >= end = [takeReferenceGenome refPosition end]
            | otherwise = [BaseSequencePosition B.empty STO.empty]
        seqOf (BaseSequencePosition nuc _) = nuc
        posOf (BaseSequencePosition _ pos) = pos
        offSet (Position x) o = Position (x + o)


loadPatterns :: FilePath -> IO Patterns
loadPatterns _ = return $ mkPatterns []

hasVariantLeft :: Genotype -> Bool
hasVariantLeft x  | x == geno00 = False
                  | x == geno10 = True
                  | x == geno11 = True
                  | x == geno01 = False
                  | otherwise = error "Bad geno"

hasVariantRight :: Genotype -> Bool
hasVariantRight x | x == geno00 = False
                  | x == geno01 = True
                  | x == geno11 = True
                  | x == geno10 = False
                  | otherwise = error "Bad geno"


combineBySecond :: Ord b => [(a,b)] -> M.Map b [a]
combineBySecond = M.fromListWith (<>) . map (\(x,y) -> (y,[x]))

variantsToDiffs :: [Variant] -> M.Map (V.Vector Diff) [HaplotypeId]
variantsToDiffs [] = M.empty
variantsToDiffs variants@(v:_) = combineBySecond (M.toList haplotypeToDiff)
    where allDiffs = V.concatMap variantToDiffs (V.fromList variants) :: V.Vector (Int, Haplotype, Diff)
          haplotypeToDiff = M.mapKeys toHaplotypeId $ V.modify sort <$> M.fromListWith (<>) (V.toList $ V.map (\(i, h, d) -> ((i,h), V.singleton d)) allDiffs) :: M.Map HaplotypeId (V.Vector Diff)
          toHaplotypeId (i,h) = HaplotypeId (sampleIds v V.! i) h


variantToDiffs :: Variant -> V.Vector (Int, Haplotype, Diff)
variantToDiffs v = V.map (\i -> (i, HaploLeft, d)) (STO.convert (genotypesL v)) <> V.map (\i -> (i, HaploRight, d)) (STO.convert (genotypesR v))
    where d = Diff (position v) (reference v) (alternative v)


findPatterns :: Chromosome -> Patterns -> FilePath -> FilePath -> FilePath -> FilePath -> IO Bool
findPatterns chr patterns peakFile referenceGenomeFile vcfFile resultFile = do
    (takeReferenceGenome, referenceGenomeSize) <- loadFasta chr referenceGenomeFile
    TIO.putStrLn $ "Chromosome " <> unChr chr <> " : " <> T.pack (show referenceGenomeSize) <> " bases"
    --patterns <- loadPatterns ""
    bed <- readPeaks chr peakFile
    case bed of
        Left err -> print err >> return False
        Right peaks -> do
            vcf <- readVcfWithGenotypes vcfFile peaks
            case vcf of
                Left err -> print err >> return False
                Right [] -> print ("No variant loaded" :: String) >> return False
                Right variants@(x:_) -> do
                    let sampleIdList = sampleIds x
                    putStrLn $ "Population: " <> show (V.length sampleIdList)
                    let sampleIndexes = V.iterateN (V.length sampleIdList) (+1) 0 :: V.Vector Int
                    let samples = M.fromList (zip (V.toList sampleIndexes) (V.toList sampleIdList))
                    t0 <- getPOSIXTime
                    processPeaks t0 chr patterns takeReferenceGenome samples (zip [1..] peaks) variants
                    return True

processPeaks :: POSIXTime
             -> Chromosome
             -> Patterns
             -> (Position0 -> Position0 -> BaseSequencePosition)
             -> M.Map Int SampleId
             -> [(Int, (Position0, Position0))]
             -> [Variant]
             -> IO ()
processPeaks _ _ _ _ _ [] _ = pure ()
processPeaks t0 chr patterns takeReferenceGenome samples ((peakId, peak):xs) variants = do
    t1 <- getPOSIXTime
    (nextVariants, matches, numberOfHaplotypes, numberOfVariants, numberOfMatches) <- processPeak chr patterns takeReferenceGenome samples peak variants
    t2 <- getPOSIXTime
    putStrLn $ formatStatus t0 t1 t2 peakId (length xs + peakId) numberOfHaplotypes numberOfVariants numberOfMatches
    processPeaks t0 chr patterns takeReferenceGenome samples xs nextVariants

formatStatus :: POSIXTime -> POSIXTime -> POSIXTime -> Int -> Int -> Int -> Int -> Int -> String
formatStatus t0 t1 t2 peakId totalPeakNumber numberOfHaplotypes numberOfVariants numberOfMatches = 
    Data.List.intercalate " \t" [peakNumberString, timeString, haploString, variantsString, matchesString]
    where peakTime = round $ (t2 - t1) * 1000 :: Integer
          totalTime = round $ (t2 - t0) * 1000 :: Integer
          peakNumberString = "Peak " <> show peakId <> "/" <> show totalPeakNumber
          timeString = show peakTime <> " ms (" <> show totalTime <> " total)"
          haploString = show numberOfHaplotypes <> " haplotypes"
          variantsString = show numberOfVariants <> " variants"
          matchesString = show numberOfMatches <> " matches"

processPeak :: Chromosome
            -> Patterns
            -> (Position0 -> Position0 -> BaseSequencePosition)
            -> M.Map Int SampleId
            -> (Position0, Position0)
            -> [Variant]
            -> IO ([Variant], V.Vector (Match [HaplotypeId]), Int, Int, Int)
processPeak chr patterns takeReferenceGenome samples (peakStart, peakEnd) variants = do
    let (variantsInPeak, nextVariants) = Data.List.span (\v -> position v <= peakEnd) $ dropWhile (\v -> position v < peakStart) variants

    let haplotypes = buildAllHaplotypes takeReferenceGenome (peakStart, peakEnd) (M.elems samples) variantsInPeak
    let matches = findPatternsInBlock (mkNucleotideAndPositionBlock (V.map fst haplotypes)) patterns
    let matchesWithSampleIds = V.map (\(Match patId score pos sampleId matched) -> Match patId score pos (snd $ haplotypes V.! sampleId) matched) matches :: V.Vector (Match [HaplotypeId])

    return (nextVariants, matchesWithSampleIds, length haplotypes, length variantsInPeak, V.sum (V.map (length . mSampleId) matchesWithSampleIds))


buildAllHaplotypes :: (Position0 -> Position0 -> BaseSequencePosition) -> (Position0, Position0) -> [SampleId] -> [Variant] -> V.Vector (BaseSequencePosition, [HaplotypeId])
buildAllHaplotypes takeReferenceGenome (peakStart, peakEnd) sampleIds variants = V.fromList $ (takeReferenceGenome peakStart peakEnd, Set.toList samplesWithNoVariant) : haplotypes
    where haplotypes = M.toList $ M.mapKeysWith (<>) (applyVariants takeReferenceGenome peakStart peakEnd . V.toList) (variantsToDiffs variants) :: [(BaseSequencePosition, [HaplotypeId])]
          samplesWithNoVariant = Set.difference allHaplotypeIds (Set.fromList $ concatMap snd haplotypes) :: Set.Set HaplotypeId
          allHaplotypeIds = Set.fromList [HaplotypeId x y | x <- sampleIds, y <- [HaploLeft, HaploRight]]


formatMatch :: Chromosome -> Position0 -> Position0 -> M.Map Int SampleId -> Match Int -> (Int, Haplotype) -> T.Text
formatMatch (Chromosome chr) (Position peakStart) (Position peakStop) samples (Match patId score pos _ matched) (i, haplo) =
    T.intercalate "\t" [chr
                       ,showt pos
                       ,showt peakStart <> "-" <> showt peakStop
                       ,showt patId
                       ,showt $ M.findWithDefault (error "Coding error: shoud have sampleId") i samples
                       ,if haplo == HaploLeft then "Left" else "Right"
                       ,showt score
                       ,showt matched
                       ,"\n"]
