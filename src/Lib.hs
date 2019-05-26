{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Lib where

import           Data.Vector ((!), toList)
import qualified Data.Vector.Storable as STO
import qualified Data.ByteString as B
import           Data.List (elemIndex)
import           Control.Monad (forM_)
import           System.IO (appendFile)
import           Debug.Trace (trace)
--
import Types
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



--[(Vector Nucleotide, Vector Position)]
--make :: 
readPeaks :: Chromosome -> FilePath -> IO [(Position ZeroBased, Position ZeroBased)]
readPeaks _ _ = return [(Position 1000, Position 1010000)]

loadPatterns :: FilePath -> IO Patterns
loadPatterns _ = return $ mkPatterns []

data Haplotype = HaploLeft | HaploRight
    deriving (Eq)

variantsToDiffs :: Haplotype -> [Variant a] -> SampleId -> [Diff a]
variantsToDiffs _ [] _ = []
variantsToDiffs haplo variants@(f:_) sample = 
    let ge = case haplo of
                  HaploLeft -> [Geno10, Geno11]
                  HaploRight -> [Geno01, Geno11]
    in case elemIndex sample (toList $ sampleIds f) of
            Nothing -> []
            Just i -> trace "onediff" [Diff (position v) (reference v) (alternative v) | v <- variants, (!) (genotypes v) i `elem` ge ]

findPatterns :: Chromosome -> FilePath -> FilePath -> FilePath -> IO Bool
findPatterns chr referenceGenomeFile vcfFile resultFile = do
    takeReferenceGenome <- loadFasta chr referenceGenomeFile
    peaks <- readPeaks chr ""
    patterns <- loadPatterns ""
    vcf <- readVcfWithGenotypes vcfFile peaks
    case vcf of
        Left err -> print err >> return False
        Right [] -> print ("No variant loaded" :: String) >> return False
        Right variants@(x:_) -> do
            let samples = toList (sampleIds x) :: [SampleId]
            forM_ peaks $ \(peakStart, peakEnd) -> do
                let variantsInPeak = takeWhile (\v -> position v <= peakEnd) $ dropWhile (\v -> position v < peakStart) variants
                let diffs = map (variantsToDiffs HaploLeft variantsInPeak <> variantsToDiffs HaploRight variantsInPeak) samples :: [[Diff ZeroBased]]
                let block = mkNucleotideAndPositionBlock $ map (applyVariants takeReferenceGenome peakStart peakEnd) diffs
                print $ blockInfo block
                matches <- findPatternsInBlock block patterns
                print matches
                appendFile resultFile (show matches)
            return True


headOr :: a -> [a] -> a 
headOr def [] = def
headOr _ (x:_) = x