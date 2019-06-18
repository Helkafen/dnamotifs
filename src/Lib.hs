{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE BangPatterns #-}

module Lib where

import qualified Data.Map as M
import qualified Data.List
import qualified Data.Vector as V
import qualified Data.ByteString as B
import qualified Data.Text as T
import qualified Data.Text.IO as TIO
import qualified Data.Text.Encoding as TE
--import           System.IO (appendFile)
import           TextShow (showt)
import           Data.Time.Clock.POSIX (getPOSIXTime, POSIXTime)
import           Data.List (iterate)

--import qualified Data.ByteString.Lazy as BL
--import qualified Codec.Compression.GZip as GZip

import           Pipes (Producer, liftIO, yield, (>->), runEffect)
import           Pipes.GZip (compress, defaultCompression)
import qualified Pipes.ByteString as PBS
import           System.IO (withFile, IOMode(..))

import           Control.Monad.Trans.Except
import           Control.Monad.Except (lift)
import           Control.Monad (when)


import Types
import Range
import Bed (readAllPeaks)
import Fasta (loadFasta)
import Vcf (readVcfWithGenotypes)
import PatternFind (findPatternsInBlock, mkNucleotideAndPositionBlock, Patterns)
import Haplotype (buildAllHaplotypes)


findPatterns :: Chromosome -> Patterns -> [FilePath] -> FilePath -> FilePath -> FilePath -> ExceptT Error IO Bool
findPatterns chr patterns peakFiles referenceGenomeFile vcfFile resultFile = do
    (takeReferenceGenome, referenceGenomeSize) <- loadFasta chr referenceGenomeFile
    lift $ TIO.putStrLn $ "Chromosome " <> unChr chr <> " : " <> T.pack (show referenceGenomeSize) <> " bases"
    (allPeaks, peaksByFile) <- readAllPeaks chr peakFiles
    (sampleIdList, variants) <- readVcfWithGenotypes vcfFile allPeaks
    lift $ putStrLn $ "Population: " <> show (length sampleIdList)
    let samples = M.fromList (zip (iterate (+1) 0) sampleIdList)
    t0 <- lift $ getPOSIXTime
    lift $ withFile resultFile WriteMode $ \fh -> do
        let header = "chromosome\tsource\tpeakStart\tpeakStop\tpatternId\tsampleId\tmatchCount\t" <> TE.encodeUtf8 (T.intercalate "\t" (map (\(SampleId s) -> s) sampleIdList))  <> "\n"
        runEffect $ compress defaultCompression (yield header >> processPeaks t0 chr patterns takeReferenceGenome peaksByFile samples (zip [1..] variants)) >-> PBS.toHandle fh
    return True

processPeaks :: POSIXTime
             -> Chromosome
             -> Patterns
             -> (Range Position0 -> BaseSequencePosition)
             -> M.Map FilePath (Ranges Position0)
             -> M.Map Int SampleId
             -> [(Int, (Range Position0, [Variant]))]
             -> Producer B.ByteString IO ()
processPeaks _ _ _ _ _ _ [] = yield B.empty
processPeaks t0 chr patterns takeReferenceGenome peaksByFile samples xs@((peakId,variants@(peak,_)):nextVariants) = do
    t1 <- liftIO getPOSIXTime
    let (matches, numberOfHaplotypes, numberOfVariants, numberOfMatches) = processPeak patterns takeReferenceGenome samples variants
    let rangesInThisPeak = M.map (overlaps peak . getRanges) peaksByFile :: M.Map FilePath [Range Position0]
    let report = M.elems $ M.mapWithKey (\file ranges -> reportAsByteString chr file (countsInPeaks samples matches ranges)) rangesInThisPeak :: [B.ByteString]
    let !bs = B.concat report
    yield bs
    t2 <- liftIO getPOSIXTime
    when (B.length bs == 0) (lift $ putStrLn "No match for this peak")
    liftIO $ putStrLn $ formatStatus t0 t1 t2 peakId (length xs + peakId) numberOfHaplotypes numberOfVariants numberOfMatches
    processPeaks t0 chr patterns takeReferenceGenome peaksByFile samples nextVariants

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
            -> (Range Position0 -> BaseSequencePosition)
            -> M.Map Int SampleId
            -> (Range Position0, [Variant])
            -> (V.Vector (Match [HaplotypeId]), Int, Int, Int)
processPeak patterns takeReferenceGenome samples (peak, variants) = do
    let haplotypes = buildAllHaplotypes takeReferenceGenome peak (M.elems samples) variants
    let matches = findPatternsInBlock (mkNucleotideAndPositionBlock (V.map fst haplotypes)) patterns                                                      :: V.Vector (Match Int)
    let matchesWithSampleIds = V.map (\(Match patId score pos sampleId matched) -> Match patId score pos (snd $ haplotypes V.! sampleId) matched) matches :: V.Vector (Match [HaplotypeId])
    (matchesWithSampleIds, length haplotypes, length variants, V.sum (V.map (length . mSampleId) matchesWithSampleIds))




countsInPeaks :: M.Map Int SampleId -> V.Vector (Match [HaplotypeId]) -> [Range Position0] -> [(Range Position0, Int, [Count2])]
countsInPeaks samples matches peaks = concatMap (countsInPeak samples matches) peaks

countsInPeak :: M.Map Int SampleId -> V.Vector (Match [HaplotypeId]) -> Range Position0 -> [(Range Position0, Int, [Count2])]
countsInPeak samples matches peak = map (\(Match patternId _ _ haplotypeIds _) -> (peak, patternId, toCounts haplotypeIds)) matchesInPeak
    where matchesInPeak = filter (inRange peak . Position . mPosition) (V.toList matches)
          toCounts :: [HaplotypeId] -> [Count2]
          toCounts haplotypeIds = M.elems abc
            where abc = M.fromListWith (<>) (map (\(HaplotypeId s h) -> (s, hToCount h)) haplotypeIds) `M.union` (M.fromList (map (,Count2 0 0) (M.elems samples))) :: M.Map SampleId Count2

hToCount :: Haplotype -> Count2
hToCount HaploLeft = Count2 1 0
hToCount HaploRight = Count2 0 1

data Count2 = Count2 Int Int
    deriving (Show)

instance Semigroup Count2 where
    (Count2 x y) <> (Count2 z v) = Count2 (x+z) (y+v)

reportAsByteString :: Chromosome -> FilePath -> [(Range Position0, Int, [Count2])] -> B.ByteString
reportAsByteString (Chromosome chr) filename xs = TE.encodeUtf8 $ T.concat $ map formatLine xs
    where formatLine ((Range s e, patternId, counts)) = (T.intercalate "\t" [chr, path, showt s, showt e, showt patternId, countsStr counts]) <> "\n"
          formatCount (Count2 c1 c2) = showt c1 <> "|" <> showt c2
          countsStr counts = T.intercalate "\t" (map formatCount counts) :: T.Text
          path = T.pack filename
