{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE BangPatterns #-}

module Run where

import qualified Data.Map as M
import qualified Data.Vector as V
import qualified Data.Vector.Storable as STO
import qualified Data.Vector.Storable.Mutable as STOM
import qualified Data.ByteString as B
import qualified Data.Text as T
import qualified Data.Text.Encoding as TE
import           TextShow (showt)
import qualified RIO.Set as Set
import           RIO.List (iterate, intercalate)
import qualified RIO.List.Partial as LP
import           Pipes (Producer, yield, (>->), runEffect)
import           Pipes.GZip (compress, defaultCompression)
import qualified Pipes.ByteString as PBS
import           System.IO (IOMode(..))
import           Text.Printf (printf)
import           Data.STRef (newSTRef, readSTRef, writeSTRef)

import Types
import Range
import Bed (readAllPeaks)
import Fasta (loadFasta)
import Vcf (readVcfWithGenotypes)
import PatternFind (findPatternsInBlock, mkNucleotideAndPositionBlock, Patterns)
import Haplotype (buildAllHaplotypes)
import Import


findPatterns :: HasLogFunc env => Chromosome -> Patterns -> [FilePath] -> FilePath -> FilePath -> FilePath -> RIO env Bool
findPatterns chr patterns peakFiles referenceGenomeFile vcfFile resultFile = do
    (takeReferenceGenome, referenceGenomeSize) <- loadFasta chr referenceGenomeFile
    logInfo $ display $ "Chromosome " <> unChr chr <> " : " <> T.pack (show referenceGenomeSize) <> " bases"
    (allPeaks, peaksByFile) <- readAllPeaks chr peakFiles
    (sampleIdList, variants) <- readVcfWithGenotypes vcfFile allPeaks
    logInfo $ display $ T.pack $ "Population: " <> show (length sampleIdList)
    let samples = M.fromList (zip (iterate (+1) 0) sampleIdList)
    t0 <- getCurrentTime
    Import.withFile resultFile WriteMode $ \fh -> do
<<<<<<< HEAD
        let header = "chromosome\tsource\tpeakStart\tpeakStop\tpatternId\tsampleId\tmatchCount\t" <> TE.encodeUtf8 (T.intercalate "\t" (map (\(SampleId s _) -> s) sampleIdList))  <> "\n"
=======
        let header = "chromosome\tsource\tpeakStart\tpeakStop\tpatternId\tsampleId\tdistinctMatchCounts\t" <> TE.encodeUtf8 (T.intercalate "\t" (map (\(SampleId s) -> s) sampleIdList))  <> "\n"
>>>>>>> add column distinctMatchCounts
        runEffect $ compress defaultCompression (yield header >> processPeaks t0 chr patterns takeReferenceGenome peaksByFile samples (zip [1..] variants)) >-> PBS.toHandle fh
        logSticky ""
    return True

processPeaks :: HasLogFunc env => UTCTime
             -> Chromosome
             -> Patterns
             -> (Range Position0 -> BaseSequencePosition)
             -> M.Map T.Text (Ranges Position0)
             -> M.Map Int SampleId
             -> [(Int, (Range Position0, [Variant]))]
             -> Producer B.ByteString (RIO env) ()
processPeaks _ _ _ _ _ _ [] = yield B.empty
processPeaks t0 chr@(Chromosome chro) patterns takeReferenceGenome peaksByFile samples xs@((peakId,variants@(peak,_)):nextVariants) = do
    t1 <- getCurrentTime
    let (matches, numberOfHaplotypes, numberOfVariants, numberOfMatches) = processPeak patterns takeReferenceGenome samples variants
    let rangesInThisPeak = M.map (overlaps peak . getRanges) peaksByFile :: M.Map T.Text [Range Position0]
    let countsAsGenotypes = countMatchesAndExpressAsHaplotypes samples matches rangesInThisPeak
    forM_ (M.assocs countsAsGenotypes) $ \((source, Range (Position start) (Position end), patternId), genotypes) ->
        yield $ TE.encodeUtf8 $ T.intercalate "\t" [chro, source, showt start, showt end, showt patternId, T.intercalate "\t" (map showt (STO.toList genotypes))] <> "\n"
    t2 <- getCurrentTime
    when (M.size countsAsGenotypes == 0) (logWarn "No match for this peak")
    logInfo $ display $ formatStatus t0 t1 t2 peakId (length xs + peakId) numberOfHaplotypes numberOfVariants numberOfMatches
    logSticky $ display $ formatStatusBar t0 t2 peakId (length xs + peakId)
    processPeaks t0 chr patterns takeReferenceGenome peaksByFile samples nextVariants


formatStatus :: UTCTime -> UTCTime -> UTCTime -> Int -> Int -> Int -> Int -> Int -> T.Text
formatStatus t0 t1 t2 peakId totalPeakNumber numberOfHaplotypes numberOfVariants numberOfMatches =
    T.pack $ intercalate " \t" [peakNumberString, timeString, haploString, variantsString, matchesString]
    where peakTime  = round $ toRational (diffUTCTime t2 t1) * 1000 :: Integer
          totalTime = round $ toRational (diffUTCTime t2 t0) * 1000 :: Integer
          peakNumberString = "Peak " <> show peakId <> "/" <> show totalPeakNumber
          timeString = show peakTime <> " ms (" <> show totalTime <> " total)"
          haploString = show numberOfHaplotypes <> " haplotypes"
          variantsString = show numberOfVariants <> " variants"
          matchesString = show numberOfMatches <> " matches"

countMatchesAndExpressAsHaplotypes :: M.Map Int SampleId -> Vector (Match [HaplotypeId]) -> M.Map T.Text [Range Position0] -> M.Map (T.Text, Range Position0, Int) (STO.Vector Genotype)
countMatchesAndExpressAsHaplotypes samples matches rangesInThisPeak =
    M.map (STO.fromList . encodeNumberOfMatches . STO.toList) $ STO.createT $ do
        mc <- newSTRef M.empty
        forM_ matches $ \(Match patternId _ pos haplos _) -> do
            let files = M.assocs $ M.filter (\ranges -> inRanges ranges (Position pos)) rangesInThisPeak :: [(T.Text, [Range Position0])]
            forM_ files $ \(source, peaks) ->
                forM_ peaks $ \peak -> do
                    m <- readSTRef mc
                    case M.lookup (source, peak, patternId) m of
                        Nothing -> do
                            array <- STOM.replicate (M.size samples) 0
                            forM_ haplos (\(HaplotypeId (SampleId _ i) _) -> STOM.modify array (\x -> x+1 :: Int) i)
                            writeSTRef mc (M.insert (source, peak, patternId) array m)
                        Just array -> forM_ haplos (\(HaplotypeId (SampleId _ i) _) -> STOM.modify array (\x -> x+1 :: Int) i)
        readSTRef mc


formatStatusBar :: UTCTime -> UTCTime -> Int -> Int -> T.Text
formatStatusBar t0 t2 peakId totalPeakNumber = progressText <> remainingTimeText
    where
        totalTime = fromRational $ toRational (diffUTCTime t2 t0)
        progressPercentage = fromIntegral peakId / (fromIntegral totalPeakNumber :: Double)
        progressText = T.pack $ printf "Progress: %.3f%%. " (100 * progressPercentage)
        remainingTime = round $ (totalTime / progressPercentage) - totalTime
        remainingDays = remainingTime `div` (3600 * 24) :: Integer
        remainingHours = (remainingTime - (remainingDays * 3600 * 24)) `div` 3600
        remainingMinutes = (remainingTime - (remainingDays * 3600 * 24) - (remainingHours * 3600)) `div` 60
        remainingSeconds = remainingTime - (remainingDays * 3600 * 24) - (remainingHours * 3600) - (remainingMinutes * 60)
        remainingTimeText =
            T.pack $ if remainingDays > 0
                        then printf "Time to completion: %.1d-%.2d:%.2d:%.2d" remainingDays remainingHours remainingMinutes remainingSeconds
                        else printf "Time to completion: %.2d:%.2d:%.2d" remainingHours remainingMinutes remainingSeconds

processPeak :: Patterns
            -> (Range Position0 -> BaseSequencePosition)
            -> M.Map Int SampleId
            -> (Range Position0, [Variant])
            -> (V.Vector (Match [HaplotypeId]), Int, Int, Int)
processPeak patterns takeReferenceGenome samples (peak, variants) = do
    let haplotypes = buildAllHaplotypes takeReferenceGenome peak (M.elems samples) variants
    let matches = findPatternsInBlock (mkNucleotideAndPositionBlock (V.map fst haplotypes)) patterns                                                      :: V.Vector (Match Int)
    let matchesWithSampleIds = V.map (\(Match patId score pos sampleId matched) -> Match patId score pos (snd $ haplotypes V.! sampleId) matched) matches :: V.Vector (Match [HaplotypeId])
    (matchesWithSampleIds, length haplotypes, length variants, V.sum (V.map (length . mSampleId) matchesWithSampleIds))


<<<<<<< HEAD
encodeNumberOfMatches :: [Int] -> [Genotype]
encodeNumberOfMatches matchNumbers = case map fst (sortBy (comparing snd) (count matchNumbers)) of
    []  -> []
    [_] -> replicate (length matchNumbers) geno00
    [x, y] -> map (\w -> if w == min x y then geno00 else geno11) matchNumbers
    [x, y, z] -> let lowest = min (min x y) z
                     highest = max (max x y) z
                 in map (\w -> if w == lowest then geno00 else if w == highest then geno11 else geno01) matchNumbers
    xs -> let lowest = LP.minimum xs
              intermediate = (highest + lowest) `div` 2
              highest = LP.maximum xs
              closest x
                | x - lowest < intermediate - x = geno00
                | highest - x < x - intermediate = geno11
                | otherwise = geno01
          in map (\w -> if w == lowest then geno00 else if w == highest then geno11 else closest w) matchNumbers
  where
    count :: Ord a => [a] -> [(a, Int)]
    count [] = []
    count (x:xs) = let (same, rest) = span (==x) xs
                   in (x,1+length same):count rest
=======

encodeNumberOfMatches :: [Int] -> ([Genotype], Set.Set Int)
encodeNumberOfMatches matchNumbers = (replaceByGenotypes (Set.toList distinctCounts), distinctCounts)
  where
    distinctCounts = Set.fromList matchNumbers
    replaceByGenotypes [] = []
    replaceByGenotypes [_] = replicate (length matchNumbers) geno00
    replaceByGenotypes [x, y] = map (\w -> if w == min x y then geno00 else geno11) matchNumbers
    replaceByGenotypes [x, y, z] = let lowest = min (min x y) z
                                       highest = max (max x y) z
                                   in map (\w -> if w == lowest then geno00 else if w == highest then geno11 else geno01) matchNumbers
    replaceByGenotypes xs = let lowest = LP.minimum xs
                                intermediate = (highest + lowest) `div` 2
                                highest = LP.maximum xs
                                closest x
                                    | x - lowest < intermediate - x = geno00
                                    | highest - x < x - intermediate = geno11
                                    | otherwise = geno01
                            in map (\w -> if w == lowest then geno00 else if w == highest then geno11 else closest w) matchNumbers


countsInPeaks :: M.Map Int SampleId -> V.Vector (Match [HaplotypeId]) -> [Range Position0] -> [(Range Position0, Int, [Count2])]
countsInPeaks samples matches = concatMap (countsInPeak samples matches)

countsInPeak :: M.Map Int SampleId -> V.Vector (Match [HaplotypeId]) -> Range Position0 -> [(Range Position0, Int, [Count2])]
countsInPeak samples matches peak = map (\(patternId,counts) -> (peak, patternId, counts)) (M.assocs byPatternId)
    where byPatternId = M.fromListWith (<>) $ V.toList $ V.map (\(Match patternId _ _ haplotypeIds _) -> (patternId, toCounts haplotypeIds)) inPeak :: M.Map Int [Count2]
          inPeak = V.filter (inRange peak . Position . mPosition) matches
          toCounts :: [HaplotypeId] -> [Count2]
          toCounts haplotypeIds = M.elems abc
            where abc = M.fromListWith (<>) (map (\(HaplotypeId s h) -> (s, hToCount h)) haplotypeIds) `M.union` M.fromList (map (,Count2 0 0) (M.elems samples)) :: M.Map SampleId Count2

hToCount :: Haplotype -> Count2
hToCount HaploLeft = Count2 1 0
hToCount HaploRight = Count2 0 1

data Count2 = Count2 Int Int
    deriving (Eq, Show)

instance Semigroup Count2 where
    (Count2 x y) <> (Count2 z v) = Count2 (x+z) (y+v)

instance Monoid Count2 where
    mempty = Count2 0 0
    mappend = (<>)

reportAsByteString :: Chromosome -> T.Text -> [(Range Position0, Int, ([Genotype], Set.Set Int))] -> B.ByteString
reportAsByteString (Chromosome chr) filename xs = TE.encodeUtf8 $ T.concat $ map formatLine xs
    where formatLine (Range s e, patternId, (counts, countSet)) = T.intercalate "\t" [chr, filename, showt s, showt e, showt patternId, countSetStr countSet, countsStr counts] <> "\n"
          countsStr counts = T.intercalate "\t" (map showt counts) :: T.Text
          countSetStr countSet = T.intercalate "\t" (map showt (Set.toList countSet)) :: T.Text
>>>>>>> add column distinctMatchCounts
