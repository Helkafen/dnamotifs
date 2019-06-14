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
import qualified Data.Text.Encoding as TE
--import           System.IO (appendFile)
import           TextShow (showt)
import           Data.Time.Clock.POSIX (getPOSIXTime, POSIXTime)

--import qualified Data.ByteString.Lazy as BL
--import qualified Codec.Compression.GZip as GZip

import           Pipes (Producer, liftIO, yield, (>->), runEffect)
import           Pipes.GZip (compress, defaultCompression)
import qualified Pipes.ByteString as PBS
import           System.IO (withFile, IOMode(..))




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
                    withFile resultFile WriteMode $ \fh ->
                        runEffect $ compress defaultCompression (processPeaks t0 chr patterns takeReferenceGenome samples (zip [1..] peaks) variants) >-> PBS.toHandle fh
                    return True

processPeaks :: POSIXTime
             -> Chromosome
             -> Patterns
             -> (Position0 -> Position0 -> BaseSequencePosition)
             -> M.Map Int SampleId
             -> [(Int, (Position0, Position0))]
             -> [Variant]
             -> Producer B.ByteString IO ()
processPeaks _ _ _ _ _ [] _ = yield B.empty
processPeaks t0 chr patterns takeReferenceGenome samples ((peakId, peak):xs) variants = do
    t1 <- liftIO getPOSIXTime
    let (nextVariants, matches, numberOfHaplotypes, numberOfVariants, numberOfMatches) = processPeak patterns takeReferenceGenome samples peak variants
    --case variants of
    --    v:_ -> print (sampleIds v)
    --    _   -> print "no var"
    t2 <- liftIO getPOSIXTime
    liftIO $ putStrLn $ formatStatus t0 t1 t2 peakId (length xs + peakId) numberOfHaplotypes numberOfVariants numberOfMatches
    yield (reportAsByteString chr (countsInPeak peakId matches))
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

processPeak :: Patterns
            -> (Position0 -> Position0 -> BaseSequencePosition)
            -> M.Map Int SampleId
            -> (Position0, Position0)
            -> [Variant]
            -> ([Variant], V.Vector (Match [HaplotypeId]), Int, Int, Int)
processPeak patterns takeReferenceGenome samples (peakStart, peakEnd) variants = do
    let (variantsInPeak, nextVariants) = Data.List.span (\v -> position v <= peakEnd) $ dropWhile (\v -> position v < peakStart) variants

    let haplotypes = buildAllHaplotypes takeReferenceGenome (peakStart, peakEnd) (M.elems samples) variantsInPeak
    let matches = findPatternsInBlock (mkNucleotideAndPositionBlock (V.map fst haplotypes)) patterns                                                      :: V.Vector (Match Int)
    let matchesWithSampleIds = V.map (\(Match patId score pos sampleId matched) -> Match patId score pos (snd $ haplotypes V.! sampleId) matched) matches :: V.Vector (Match [HaplotypeId])

    (nextVariants, matchesWithSampleIds, length haplotypes, length variantsInPeak, V.sum (V.map (length . mSampleId) matchesWithSampleIds))


buildAllHaplotypes :: (Position0 -> Position0 -> BaseSequencePosition) -> (Position0, Position0) -> [SampleId] -> [Variant] -> V.Vector (BaseSequencePosition, [HaplotypeId])
buildAllHaplotypes takeReferenceGenome (peakStart, peakEnd) sampleIdList variants = V.fromList $ (takeReferenceGenome peakStart peakEnd, Set.toList samplesWithNoVariant) : haplotypes
    where haplotypes = M.toList $ M.mapKeysWith (<>) (applyVariants takeReferenceGenome peakStart peakEnd . V.toList) (variantsToDiffs variants) :: [(BaseSequencePosition, [HaplotypeId])]
          samplesWithNoVariant = Set.difference allHaplotypeIds (Set.fromList $ concatMap snd haplotypes) :: Set.Set HaplotypeId
          allHaplotypeIds = Set.fromList [HaplotypeId x y | x <- sampleIdList, y <- [HaploLeft, HaploRight]]


countsInPeak :: Int -> V.Vector (Match [HaplotypeId]) -> M.Map (Int, Int, SampleId) Count2
countsInPeak peakId matches = M.unionsWith (<>) (map (countM peakId) (V.toList matches))

countM :: Int -> Match [HaplotypeId] -> M.Map (Int, Int, SampleId) Count2
countM peakId (Match patternId _ _ haplotypeIds _) = fmap (\s -> (Count2 (count HaploLeft s) (count HaploRight s))) x
    where x = M.fromListWith (<>) $ map (\(HaplotypeId sampleId haplo) -> ((peakId, patternId, sampleId), [haplo])) haplotypeIds :: M.Map (Int, Int, SampleId) [Haplotype]

count :: Eq a => a -> [a] -> Int
count x xs = (length . filter (== x)) xs

data Count2 = Count2 Int Int
    deriving (Show)

instance Semigroup Count2 where
    (Count2 x y) <> (Count2 z v) = Count2 (x+z) (y+v)

reportAsByteString :: Chromosome -> M.Map (Int, Int, SampleId) Count2 -> B.ByteString
reportAsByteString (Chromosome chr) d = TE.encodeUtf8 $ T.concat $ map formatLine (M.toList d)
    where formatLine ((peakId, patternId, sampleId), Count2 c1 c2) = T.intercalate "\t" [showt chr, showt peakId, showt patternId, showt sampleId, showt c1 <> "/" <> showt c2 <> "\n"]
    
    