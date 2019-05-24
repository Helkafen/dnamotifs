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
import           Debug.Trace (trace, traceShow)
--
import Types
import Fasta (loadFasta)
import Vcf (readVcfWithGenotypes)
import PatternFind (findPatternsInBlock, mkPatterns, mkNucleotideAndPositionBlock, NucleotideAndPositionBlock, Patterns, blockInfo)

toNuc :: Nucleotide -> Nucleotide -- TODO ugly
toNuc 65 = a
toNuc 67 = c
toNuc 71 = g
toNuc 84 = t
toNuc 78 = n
toNuc 97 = a
toNuc 99 = c
toNuc 103 = g
toNuc 116 = t
toNuc 110 = n
toNuc other = error $ "Bad nucleotide " <> show other

-- Inclusive [Start, End] interval
applyVariants :: V.Vector Nucleotide -> Position ZeroBased -> Position ZeroBased -> [Diff ZeroBased] -> [(Nucleotide, Position ZeroBased)]
applyVariants referenceGenome (Position start) (Position end) allDiffs =
    DList.toList $ go start (filter (\(Diff (Position p) _ _) -> p >= start && p <= end) allDiffs)
  where go refPosition [] = takeRef referenceGenome refPosition end
        go refPosition diffs@(Diff (Position pos) ref alt:vs)
            | pos > refPosition = takeRef referenceGenome refPosition (pos-1) <> go pos diffs
            | pos == refPosition && U.length ref == 1 = DList.fromList (zip (U.toList alt) (repeat (Position refPosition))) <> go (refPosition + 1) vs -- SNV, insertion
            | pos == refPosition && U.length alt == 1 = DList.fromList (zip (U.toList alt) (repeat (Position refPosition))) <> go (refPosition + U.length ref) vs -- deletion
            | pos == refPosition = error "Missing case in applyVariants (we assume that ref or alt have length 1)"
            | refPosition >= end = takeRef referenceGenome refPosition end
            | otherwise = DList.empty

takeRef :: V.Vector Nucleotide ->  Int -> Int -> DList.DList (Nucleotide, Position ZeroBased)
takeRef referenceGenome s e = DList.fromList $ zip (V.toList $ V.map toNuc $ V.take (e-s+1) (V.drop s referenceGenome)) (map Position [s..e+1])


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
            Just i -> traceShow "onediff" [Diff (position v) (reference v) (alternative v) | v <- variants, (!) (genotypes v) i `elem` ge ]

overlappingchunksOf :: Int -> Int -> [(Nucleotide,Position a)] -> [(V.Vector Nucleotide, V.Vector (Position a))]
overlappingchunksOf window overlap xs = map toVec (divvy window (window - overlap) xs)
        where toVec :: [(Nucleotide, Position a)] -> (V.Vector Nucleotide, V.Vector (Position a))
              toVec np = (V.fromList (map fst np), V.fromList (map snd np))


--vectorsOfSample referenceGenome peakStart peakEnd haplo variants = overlappingchunksOf 100 10 . applyVariants referenceGenome peakStart peakEnd . variantsToDiffs haplo variants

findPatterns :: Chromosome -> FilePath -> FilePath -> FilePath -> IO Bool
findPatterns chr referenceGenomeFile vcfFile resultFile = do
    referenceGenome <- loadFasta chr referenceGenomeFile
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
                --let d = variantsToDiffs HaploLeft variants (head samples)
                --let s = applyVariants referenceGenome peakStart peakEnd d :: [(Nucleotide, Position)]
                --let cs = overlappingchunksOf 100 10 s :: [(V.Vector Nucleotide, V.Vector Position)]
                --let bs = map mkNucleotideAndPositionBlock cs
                --print (map (show . blockInfo) bs)
                --vectorsOfSample

                let diffs = map (variantsToDiffs HaploLeft variantsInPeak <> variantsToDiffs HaploRight variantsInPeak) samples :: [[Diff ZeroBased]]
                let seqs = map (applyVariants referenceGenome peakStart peakEnd) diffs :: [[(Nucleotide, Position ZeroBased)]]
                let chunks = map (overlappingchunksOf 100 10) seqs :: [[(V.Vector Nucleotide, V.Vector (Position ZeroBased))]]
                --forM_ chunks $ \chunk -> do
                --    print (chunk)
                go chunks
                --forM_ blocks $ \block -> do
                --    print (blockInfo block)
                    --matches <- findPatternsInBlock block patterns
                    --print matches
                    --appendFile resultFile (show matches)
                -- overlappingchunksOf 100 10 . 
                --undefined --patched = map (applyVariants referenceGenome peakStart peakEnd variants) [0..length samples]
            return True
    where go :: [[(V.Vector Nucleotide, V.Vector (Position ZeroBased))]] -> IO ()
          go ch = do
            let heads = map (headOr (V.empty, V.empty)) ch
            case all ((==0) . V.length . fst) heads of
                False -> do
                    let block = mkNucleotideAndPositionBlock heads
                    print $ blockInfo block
                    go (map tail ch)
                True -> do
                    print "End of peak"


headOr :: a -> [a] -> a 
headOr def [] = def
headOr _ (x:_) = x