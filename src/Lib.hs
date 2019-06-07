{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}

module Lib where

import qualified Data.Map
import qualified Data.List
import qualified Data.Set as Set
import qualified Data.Vector as V
import           Data.Vector.Algorithms.Merge (sort)
import qualified Data.Vector.Storable as STO
import qualified Data.ByteString as B
import qualified Data.Text as T
import qualified Data.Text.IO as TIO
import           Control.Monad (forM_)
--import           System.IO (appendFile)
import           TextShow (showt)
import           Data.Time.Clock.POSIX (getPOSIXTime, POSIXTime)


import Types
import Bed (readPeaks)
import Fasta (loadFasta)
import Vcf (readVcfWithGenotypes)
import PatternFind (findPatternsInBlock, mkPatterns, mkNucleotideAndPositionBlock, Patterns)


-- Inclusive [Start, End] interval
applyVariants :: (Int -> Int -> BaseSequencePosition) -> Position0 -> Position0 -> [Diff] -> BaseSequencePosition
applyVariants takeReferenceGenome (Position start) (Position end) allDiffs =
    let chunks = go start (filter (\(Diff (Position p) _ _) -> p >= start && p <= end) allDiffs)
    in BaseSequencePosition (B.concat (map seqOf chunks)) (STO.concat (map posOf chunks))
  where go :: Int -> [Diff] -> [BaseSequencePosition]
        go refPosition [] = [takeReferenceGenome refPosition end]
        go refPosition diffs@(Diff (Position pos) ref alt:vs)
            | pos > refPosition = takeReferenceGenome refPosition (pos-1): go pos diffs
            | pos == refPosition && B.length ref == 1 = -- SNV, insertion
                BaseSequencePosition alt (STO.replicate (B.length alt) (fromIntegral refPosition)): go (refPosition + 1) vs
            | pos == refPosition && B.length alt == 1 = -- deletion
                BaseSequencePosition alt (STO.singleton (fromIntegral refPosition)) : go (refPosition + B.length ref) vs
            | pos == refPosition = error "Missing case in applyVariants (we assume that ref or alt have length 1)"
            | refPosition >= end = [takeReferenceGenome refPosition end]
            | otherwise = [BaseSequencePosition B.empty STO.empty]
        seqOf (BaseSequencePosition nuc _) = nuc
        posOf (BaseSequencePosition _ pos) = pos


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


variantsToDiffs :: [Variant] -> Data.Map.Map (Int, Haplotype) (V.Vector Diff)
variantsToDiffs variants = V.modify sort <$> Data.Map.fromListWith (<>) (V.toList $ V.map (\(i, h, d) -> ((i,h), V.singleton d)) allDiffs)
    where allDiffs = V.concatMap variantToDiffs (V.fromList variants) :: V.Vector (Int, Haplotype, Diff)


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
                    let samples = Data.Map.fromList (zip (V.toList sampleIndexes) (V.toList sampleIdList))
                    --let sampleIndexesTwoHaplotypes = sampleIndexes <> sampleIndexes
                    t0 <- getPOSIXTime
                    processPeaks t0 chr patterns takeReferenceGenome samples (zip [1..] peaks) variants
                    return True

processPeaks :: POSIXTime
             -> Chromosome
             -> Patterns
             -> (Int -> Int -> BaseSequencePosition)
             -> Data.Map.Map Int SampleId
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
            -> (Int -> Int -> BaseSequencePosition)
            -> Data.Map.Map Int SampleId
            -> (Position0, Position0)
            -> [Variant]
            -> IO ([Variant], V.Vector (Match [(SampleId, Haplotype)]), Int, Int, Int)
processPeak chr patterns takeReferenceGenome samples (peakStart, peakEnd) variants = do
    let (variantsInPeak, nextVariants) = Data.List.span (\v -> position v <= peakEnd) $ dropWhile (\v -> position v < peakStart) variants
    
    let uniqueDiffs = Data.Map.fromListWith (<>) $ map (\(x,y) -> (y,[x])) $ Data.Map.toList $ variantsToDiffs variantsInPeak :: Data.Map.Map (V.Vector Diff) [(Int, Haplotype)]
    
    let uniqueSequences = Data.Map.toList $ fmap (map (\(i,h) -> (toSampleId i,h))) $ Data.Map.mapKeysWith (<>) (applyVariants takeReferenceGenome peakStart peakEnd . V.toList) uniqueDiffs :: [(BaseSequencePosition, [(SampleId, Haplotype)])]

    matches <- findPatternsInBlock (mkNucleotideAndPositionBlock (map fst uniqueSequences)) patterns
    let matchesWithSampleIds = V.map (\(Match patId score pos sampleId matched) -> Match patId score pos (map snd uniqueSequences !! sampleId) matched) matches :: V.Vector (Match [(Int, Haplotype)])

    
    --matches2 <- findPatternsInBlock (mkNucleotideAndPositionBlock [applyVariants takeReferenceGenome peakStart peakEnd []]) patterns

    
    --V.map (\(Match patId score pos sampleId matched) -> Match patId score pos (map snd uniqueSequences) matched) matches1

    let matches = V.empty

    return (nextVariants, matchesWithSampleIds, length uniqueSequences, length variantsInPeak, V.sum (V.map (length . mSampleId) matchesWithSampleIds))
    where toSampleId = flip Data.Map.findWithDefault (error "Coding error: shoud have sampleId") samples
    -- numberOfHaplotypes, (length variantsInPeak), V.length matches
    
--    -- We don't want to generate and scan the same haplotypes n times in a large population, so
--    -- we put the unique haplotypes in an array, scan that array, then use the indices (mSampleId) 
--    -- of each Match to recover the [(SampleId, Haplotype)] of origin
--
--
--
--    -- Just some book keeping
--    
--
--    -- The actual scan
--    let noDiff = []
--    let block = mkNucleotideAndPositionBlock (map (applyVariants takeReferenceGenome peakStart peakEnd) ((V.toList <$> uniqueDiffs) ++ noDiff))
--    let numberOfHaplotypes = length uniqueDiffs + 1 -- 1 for the reference haplotype
--    matches <- findPatternsInBlock block patterns
--
--    -- Recover the [(SampleId, Haplotype)] of each match and print
--    forM_ (V.toList matches) $ \(Match patId score pos sampleId matched) -> do
--        let diffsOfMatch = uniqueDiffs !! sampleId :: V.Vector Diff
--        let haploIdsOfMatch = Data.Map.findWithDefault haplotypesSampleIdWithNoVariant diffsOfMatch m :: [(SampleId, Haplotype)]
--        let matches1 = map (\(i,h) -> Match patId score pos (Data.Map.findWithDefault (error "Coding error: shoud have sampleId") i samples,h) matched) haploIdsOfMatch :: [Match (SampleId, Haplotype)]
--
--        --let strings = map (formatMatch chr peakStart peakEnd samples match) haploIdsOfMatch :: [T.Text]
--        --forM_ (take 2 strings) $ TIO.putStrLn
--        --if (length strings == -1)
--        --    then print strings
--        --    else pure()
--        --appendFile resultFile (Data.List.concat strings)
--        --putStrLn ("No append" :: String)
--        pure ()
--    
--    pure (nextVariants, V.empty, numberOfHaplotypes, (length variantsInPeak), V.length matches)
--  where 
--        diffs = variantsToDiffs variantsInPeak :: Data.Map.Map (Int, Haplotype) (V.Vector Diff)
--        uniqueDiffs = Set.toList $ Set.fromList (Data.Map.elems diffs) :: [(V.Vector Diff)]
--        sampleIndexes = V.iterateN (Data.Map.size samples) (+1) 0 :: V.Vector Int
--        
--        haplotypeIds = V.map (,HaploLeft) sampleIndexes <> V.map (,HaploRight) sampleIndexes :: V.Vector (Int, Haplotype)
--        haplotypeSampleIds = V.map (\(i,h) -> (Data.Map.findWithDefault (error "Coding error: shoud have sampleId") i samples,h)) haplotypeIds :: V.Vector (SampleId, Haplotype)
--        
--        haplotypesWithNoVariant = Set.toList $ Set.difference (Set.fromList $ V.toList haplotypeIds) (Data.Map.keysSet diffs) :: [(Int, Haplotype)]
--        haplotypesSampleIdWithNoVariant = map (\(i,h) -> Data.Map.findWithDefault (error "Coding error: shoud have sampleId") i samples,h) haplotypesWithNoVariant
--        
--        m = Data.Map.fromListWith (<>) (map (\(x,y) -> (y,[x])) $ Data.Map.toList diffs) :: Data.Map.Map (V.Vector Diff) [(Int, Haplotype)]
--
--        hydrate (Match patId score pos sampleId matched) = map (\(i,h) -> Match patId score pos (Data.Map.findWithDefault (error "Coding error: shoud have sampleId") i samples,h) matched) haploIdsOfMatch


formatMatch :: Chromosome -> Position0 -> Position0 -> Data.Map.Map Int SampleId -> Match Int -> (Int, Haplotype) -> T.Text
formatMatch (Chromosome chr) (Position peakStart) (Position peakStop) samples (Match patId score pos _ matched) (i, haplo) =
    T.intercalate "\t" [chr
                       ,showt pos
                       ,showt peakStart <> "-" <> showt peakStop
                       ,showt patId
                       ,showt $ Data.Map.findWithDefault (error "Coding error: shoud have sampleId") i samples
                       ,if haplo == HaploLeft then "Left" else "Right"
                       ,showt score
                       ,showt matched
                       ,"\n"]
