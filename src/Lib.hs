{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}

module Lib where

import qualified Data.Map
import qualified Data.List
import qualified Data.Set as Set
import qualified Data.Vector as V
import qualified Data.Vector.Storable as STO
import qualified Data.ByteString as B
import qualified Data.Text as T
import qualified Data.Text.IO as TIO
import           Control.Monad (forM_)
import           System.IO (appendFile)
import qualified Codec.Compression.GZip as GZip
import           Debug.Trace (trace)

import Types
import Bed (readPeaks)
import Fasta (loadFasta)
import Vcf (readVcfWithGenotypes)
import PatternFind (findPatternsInBlock, mkPatterns, mkNucleotideAndPositionBlock, Patterns, blockInfo)


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

variantsToDiffs :: Haplotype -> [Variant] -> Int -> [Diff]
variantsToDiffs _ [] _ = []
variantsToDiffs haplo variants i = 
    let ge = case haplo of
                  HaploLeft -> [geno10, geno11]
                  HaploRight -> [geno01, geno11]
    in [Diff (position v) (reference v) (alternative v) | v <- variants, (STO.!) (genotypes v) i `elem` ge ]

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
                    let sampleIndexes = V.iterateN (V.length (sampleIds x)) (+1) 0 :: V.Vector Int
                    let sampleIndexesTwoHaplotypes = sampleIndexes <> sampleIndexes
                    forM_ peaks $ \(peakStart, peakEnd) -> do
                        let variantsInPeak = takeWhile (\v -> position v <= peakEnd) $ dropWhile (\v -> position v < peakStart) variants
                        let diffs = V.map (variantsToDiffs HaploLeft variantsInPeak <> variantsToDiffs HaploRight variantsInPeak) sampleIndexesTwoHaplotypes

                        -- We don't want to generate and scan the same haplotypes n times in a large population, so
                        -- we put the unique haplotypes in an array, scan that array, then use the indices (mSampleId) 
                        -- of each Match to recover the [(SampleId, Haplotype)] of origin
                        let uniqueDiffs = Set.toList $ Set.fromList (V.toList diffs)
                        let numberOfUniqueSequences = length uniqueDiffs
                        let haplotypeIds = V.map (,HaploLeft) (sampleIds x) <> V.map (,HaploRight) (sampleIds x) :: V.Vector (SampleId, Haplotype)

                        -- Just some book keeping
                        let m = Data.Map.fromListWith (++) (zip (V.toList diffs) (map (:[]) (V.toList haplotypeIds))) :: Data.Map.Map [Diff] [(SampleId, Haplotype)]

                        -- The actual scan
                        let block = mkNucleotideAndPositionBlock (map (applyVariants takeReferenceGenome peakStart peakEnd) uniqueDiffs)
                        putStr $ (blockInfo block <> (", for " <> show (V.length (sampleIds x)) <> " people"))
                        matches <- findPatternsInBlock block patterns
                        --print matches

                        -- Recover the [(SampleId, Haplotype)] of each match and print
                        forM_ matches $ \match -> do
                            let diffsOfMatch = uniqueDiffs !! mSampleId match :: [Diff]
                            let haploIdsOfMatch = Data.Map.findWithDefault (error "Coding error: should find indices") diffsOfMatch m :: [(SampleId, Haplotype)]
                            let strings = map (formatMatch chr peakStart peakEnd match) haploIdsOfMatch :: [String]
                            if (length strings == -1)
                                then print strings
                                else pure()
                            --appendFile resultFile (Data.List.concat strings)
                            --putStrLn ("No append" :: String)
                            pure ()

                    return True


formatMatch :: Chromosome -> Position0 -> Position0 -> Match -> (SampleId, Haplotype) -> String
formatMatch (Chromosome chr) (Position peakStart) (Position peakStop) (Match patId score pos _ matched) (SampleId sample, haplo) =
    Data.List.intercalate "\t" [T.unpack chr
                               ,show pos
                               ,show peakStart <> "-" <> show peakStop
                               ,show patId
                               ,show sample
                               ,show haplo
                               ,show score
                               ,show matched
                               ,"\n"]
