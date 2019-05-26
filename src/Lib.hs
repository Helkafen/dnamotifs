{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Lib where

import qualified Data.Set as Set
import qualified Data.Vector as V
import qualified Data.Vector.Storable as STO
import qualified Data.ByteString as B
import qualified Data.Text as T
import qualified Data.Text.IO as TIO
import           Control.Monad (forM_)
import           System.IO (appendFile)
import           Debug.Trace (trace)

import Types
import Bed (readPeaks)
import Fasta (loadFasta)
import Vcf (readVcfWithGenotypes)
import PatternFind (findPatternsInBlock, mkPatterns, mkNucleotideAndPositionBlock, Patterns, blockInfo)


-- Inclusive [Start, End] interval
applyVariants :: (Int -> Int -> BaseSequencePosition) -> Position ZeroBased -> Position ZeroBased -> [Diff ZeroBased] -> BaseSequencePosition
applyVariants takeReferenceGenome (Position start) (Position end) allDiffs =
    let chunks = go start (filter (\(Diff (Position p) _ _) -> p >= start && p <= end) allDiffs)
    in BaseSequencePosition (B.concat (map seqOf chunks)) (STO.concat (map posOf chunks))
  where go :: Int -> [Diff ZeroBased] -> [BaseSequencePosition]
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

variantsToDiffs :: Haplotype -> [Variant a] -> Int -> [Diff a]
variantsToDiffs _ [] _ = []
variantsToDiffs haplo variants i = 
    let ge = case haplo of
                  HaploLeft -> [Geno10, Geno11]
                  HaploRight -> [Geno01, Geno11]
    in trace "onediff" [Diff (position v) (reference v) (alternative v) | v <- variants, (V.!) (genotypes v) i `elem` ge ]

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
                        let numberOfUniqueSequences = length (Set.fromList (V.toList diffs))
                        print $ "Number of unique sequences :" ++ show numberOfUniqueSequences ++ "/" ++ show (V.length (sampleIds x))
                        let block = mkNucleotideAndPositionBlock $ map (applyVariants takeReferenceGenome peakStart peakEnd) (V.toList diffs)
                        print $ blockInfo block
                        matches <- findPatternsInBlock block patterns
                        print matches
                        appendFile resultFile (show matches)
                    return True


headOr :: a -> [a] -> a 
headOr def [] = def
headOr _ (x:_) = x