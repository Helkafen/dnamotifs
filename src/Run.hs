{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}

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
import           Data.STRef (newSTRef, readSTRef, writeSTRef, modifySTRef)
import           Data.Vector.Storable.ByteString (vectorToByteString)


import Types
import Range
import Bed (readAllPeaks)
import Fasta (loadFasta)
import Vcf (readVcfWithGenotypes)
import PatternFind (findPatternsInBlock, mkNucleotideAndPositionBlock, Patterns)
import Haplotype (buildAllHaplotypes)
import MotifDefinition (loadHocomocoPatternsAndScoreThresholds)
import Import


findPatterns :: HasLogFunc env => Chromosome -> [T.Text] -> [FilePath] -> FilePath -> FilePath -> FilePath -> RIO env Bool
findPatterns chr motifNames peakFiles referenceGenomeFile vcfFile resultFile = do
    patterns <- loadHocomocoPatternsAndScoreThresholds motifNames
    (takeReferenceGenome, referenceGenomeSize) <- loadFasta chr referenceGenomeFile
    logInfo $ display $ "Chromosome " <> unChr chr <> " : " <> T.pack (show referenceGenomeSize) <> " bases"
    (allPeaks, peaksByFile) <- readAllPeaks chr peakFiles
    (sampleIdList, variants) <- readVcfWithGenotypes vcfFile allPeaks
    logInfo $ display $ T.pack $ "Population: " <> show (length sampleIdList)
    let samples = M.fromList (zip (iterate (+1) 0) sampleIdList)
    t0 <- getCurrentTime
    fakePos <- newIORef (1 :: Int)
    Import.withFile resultFile WriteMode $ \fh -> do
        let header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" <> TE.encodeUtf8 (T.intercalate "\t" (map (\(SampleId s _) -> s) sampleIdList))  <> "\n"
        runEffect $ compress defaultCompression (yield header >> processPeaks t0 fakePos chr patterns takeReferenceGenome peaksByFile samples (zip [1..] variants)) >-> PBS.toHandle fh
        logSticky ""
    return True

processPeaks :: HasLogFunc env => UTCTime
             -> IORef Int
             -> Chromosome
             -> ([T.Text], Patterns)
             -> (Range Position0 -> BaseSequencePosition)
             -> M.Map T.Text (Ranges Position0)
             -> M.Map Int SampleId
             -> [(Int, (Range Position0, [Variant]))]
             -> Producer B.ByteString (RIO env) ()
processPeaks _ _ _ _ _ _ _ [] = yield B.empty
processPeaks t0 fakePos chr@(Chromosome chro) (patternNames, patterns) takeReferenceGenome peaksByFile samples xs@((peakId,variants@(peak,_)):nextVariants) = do
    t1 <- getCurrentTime
    let (matches, numberOfHaplotypes, numberOfVariants, numberOfMatches) = processPeak patterns takeReferenceGenome samples variants
    let rangesInThisPeak = M.map (overlaps peak . getRanges) peaksByFile :: M.Map T.Text [Range Position0]
    let countsAsGenotypes = countMatchesAndExpressAsHaplotypes samples matches rangesInThisPeak
    forM_ (M.assocs countsAsGenotypes) $ \((source, Range (Position start) (Position end), patternId), (genotypes, s)) -> do
        posStr <- showt <$> readIORef fakePos
        let idStr = source <> "," <> patternNames LP.!! patternId <> "," <> showt start <> "-" <> showt end
        let refStr = "."
        let altStr = "."
        let qualStr = "."
        let filterStr = "."
        let infoStr = "COUNTS=" <> T.intercalate "," (map showt (Set.toList s))
        yield $ TE.encodeUtf8 $ T.intercalate "\t" [chro, posStr, idStr, refStr, altStr, qualStr, filterStr, infoStr]
        yield (genoText genotypes)
        yield "\n"
        modifyIORef fakePos (+1)
    t2 <- getCurrentTime
    when (M.size countsAsGenotypes == 0) (logWarn "No match for this peak")
    logInfo $ display $ formatStatus t0 t1 t2 peakId (length xs + peakId) numberOfHaplotypes numberOfVariants numberOfMatches
    logSticky $ display $ formatStatusBar t0 t2 peakId (length xs + peakId)
    processPeaks t0 fakePos chr (patternNames, patterns) takeReferenceGenome peaksByFile samples nextVariants


genoText :: STO.Vector Genotype -> B.ByteString
genoText ge = vectorToByteString $ STO.create $ do
        x <- STOM.replicate (len * 4) (ord '=')
        mi <- newSTRef (0 :: Int)
        STO.forM_ ge $ \v -> do
            i <- readSTRef mi
            if      v == geno00 then STOM.unsafeWrite x (i*4) 9 >> STOM.unsafeWrite x (i*4+1) 48 >> STOM.unsafeWrite x (i*4+2) 124 >> STOM.unsafeWrite x (i*4+3) 48
            else if v == geno01 then STOM.unsafeWrite x (i*4) 9 >> STOM.unsafeWrite x (i*4+1) 48 >> STOM.unsafeWrite x (i*4+2) 124 >> STOM.unsafeWrite x (i*4+3) 49
            else                     STOM.unsafeWrite x (i*4) 9 >> STOM.unsafeWrite x (i*4+1) 49 >> STOM.unsafeWrite x (i*4+2) 124 >> STOM.unsafeWrite x (i*4+3) 49
            modifySTRef mi (+1)
        return x
    where len = STO.length ge



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

countMatchesAndExpressAsHaplotypes :: M.Map Int SampleId -> Vector (Match [HaplotypeId]) -> M.Map T.Text [Range Position0] -> M.Map (T.Text, Range Position0, Int) (STO.Vector Genotype, Set.Set Int)
countMatchesAndExpressAsHaplotypes samples matches rangesInThisPeak =
    M.map (\xs ->  let (ge, s) = encodeNumberOfMatchesAsGenotypes (STO.toList xs) in (STO.fromList ge, s)) $ STO.createT $ do
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
    let matches = findPatternsInBlock (mkNucleotideAndPositionBlock (V.map fst haplotypes)) patterns      :: V.Vector (Match Int)
    let matchesWithSampleIds = V.map (\m -> m { mSampleId = snd (haplotypes V.! mSampleId m) }) matches :: V.Vector (Match [HaplotypeId])
    (matchesWithSampleIds, length haplotypes, length variants, V.sum (V.map (length . mSampleId) matchesWithSampleIds))


encodeNumberOfMatchesAsGenotypes :: [Int] -> ([Genotype], Set.Set Int)
encodeNumberOfMatchesAsGenotypes matchNumbers = (replaceByGenotypes (Set.toList distinctCounts), distinctCounts)
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

