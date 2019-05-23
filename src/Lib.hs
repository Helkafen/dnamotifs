{-# LANGUAGE OverloadedStrings #-}

module Lib where

import qualified Data.DList as DList
import           Data.Vector ((!), toList)
import qualified Data.Vector.Storable as V
import qualified Data.Vector.Unboxed as U
--import qualified Data.Text as T
import           Data.List (elemIndex)
import           Control.Monad (forM_)
import           Data.List.Split (divvy)
import           System.IO (appendFile)
--
import Types
import Fasta (loadFasta)
import Vcf (readVcfWithGenotypes)
import PatternFind (findPatternsInBlock, mkPatterns, mkNucleotideAndPositionBlock, NucleotideAndPositionBlock, Patterns)

-- Inclusive [Start, End] interval
applyVariants :: V.Vector Nucleotide -> Position -> Position -> [Diff] -> [(Nucleotide, Position)]
applyVariants referenceGenome (Position start) (Position end) allDiffs =
    DList.toList $ go start (filter (\(Diff (Position p) _ _) -> p >= start && p <= end) allDiffs)
  where go refPosition [] = takeRef refPosition end
        go refPosition diffs@(Diff (Position pos) ref alt:vs)
            | pos > refPosition = takeRef refPosition (pos-1) <> go pos diffs
            | pos == refPosition && U.length ref == 1 = DList.fromList (zip (U.toList alt) (repeat (Position refPosition))) <> go (refPosition + 1) vs -- SNV, insertion
            | pos == refPosition && U.length alt == 1 = DList.fromList (zip (U.toList alt) (repeat (Position refPosition))) <> go (refPosition + U.length ref) vs -- deletion
            | pos == refPosition = error "Missing case in applyVariants (we assume that ref or alt have length 1)"
            | refPosition >= end = takeRef refPosition end
            | otherwise = DList.empty
        takeRef :: Int -> Int -> DList.DList (Nucleotide, Position)
        takeRef s e = DList.fromList $ zip (V.toList $ V.take (e-s+1) (V.drop s referenceGenome)) (map Position [s..e+1])


--[(Vector Nucleotide, Vector Position)]
--make :: 
readPeaks :: Chromosome -> FilePath -> IO [(Position, Position)]
readPeaks _ _ = return [(Position 10, Position 20)]

loadPatterns :: FilePath -> IO Patterns
loadPatterns _ = return $ mkPatterns []

data Haplotype = HaploLeft | HaploRight
    deriving (Eq)

variantsToDiffs :: Haplotype -> [Variant] -> SampleId -> [Diff]
variantsToDiffs _ [] _ = []
variantsToDiffs haplo variants@(f:_) sample = 
    let ge = case haplo of
                  HaploLeft -> [Geno10, Geno11]
                  HaploRight -> [Geno01, Geno11]
    in case elemIndex sample (toList $ sampleIds f) of
            Nothing -> []
            Just i -> [Diff (position v) (reference v) (alternative v) | v <- variants, (!) (genotypes v) i `elem` ge ]

overlappingchunksOf :: Int -> Int -> [(Nucleotide,Position)] -> [(V.Vector Nucleotide, V.Vector Position)]
overlappingchunksOf window overlap xs = map toVec (divvy window (window - overlap) xs)
        where toVec :: [(Nucleotide, Position)] -> (V.Vector Nucleotide, V.Vector Position)
              toVec np = (V.fromList (map fst np), V.fromList (map snd np))

findPatterns :: Chromosome -> FilePath -> FilePath -> FilePath -> IO Bool
findPatterns chr referenceGenomeFile vcfFile resultFile = do
    referenceGenome <- loadFasta chr referenceGenomeFile
    peaks <- readPeaks chr ""
    patterns <- loadPatterns ""
    vcf <- readVcfWithGenotypes vcfFile
    case vcf of
        Left err -> print err >> return False
        Right [] -> print ("No variant loaded" :: String) >> return False
        Right variants@(x:_) -> do
            let samples = toList (sampleIds x) :: [SampleId]
            forM_ peaks $ \(peakStart, peakEnd) -> do
                let diffs = map (variantsToDiffs HaploLeft variants <> variantsToDiffs HaploRight variants) samples :: [[Diff]]
                let seqs = map (applyVariants referenceGenome peakStart peakEnd) diffs :: [[(Nucleotide, Position)]]
                let blocks = map (mkNucleotideAndPositionBlock . overlappingchunksOf 100 10) seqs :: [NucleotideAndPositionBlock]
                forM_ blocks $ \block -> do
                    matches <- findPatternsInBlock block patterns
                    print matches
                    appendFile resultFile (show matches)
                -- overlappingchunksOf 100 10 . 
                undefined --patched = map (applyVariants referenceGenome peakStart peakEnd variants) [0..length samples]
            return True
