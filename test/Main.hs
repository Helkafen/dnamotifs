{-# OPTIONS_GHC -F -pgmF htfpp #-}
{-# LANGUAGE OverloadedStrings #-}

import Test.Framework

--import           Foreign.C.Types                 (CChar)
import           Data.Vector.Storable            (Vector)
import qualified Data.Vector.Storable            as V
import qualified Data.Vector.Unboxed             as U
import qualified Data.Vector.Generic             as G
import           Data.Monoid                     ((<>))
import qualified Data.Text                       as T

import Types
import PatternFind
import Lib
import Vcf (parseVariant, filterOrderedIntervals, parseVcfContent)

inputDataSample0 :: Vector Nucleotide
inputDataSample0 = V.fromList [a,a,a,a,c,g,a] <> V.fromList (replicate 93 a)

inputDataSample1 :: Vector Nucleotide
inputDataSample1 = V.fromList [c,g,a,a,a,a,a] <> V.fromList (replicate 93 a) 

inputDataSample2 :: Vector Nucleotide
inputDataSample2 = V.fromList [a,a,a,a,a,a,a] <> V.fromList (replicate 93 a)

pattern_CG, pattern_cCGA, pattern_CGA, pattern_cccccccccc :: [Pweight]
pattern_CG = [Pweight 0 1 0 0, Pweight 0 0 1 0]
pattern_cCGA = [Pweight 0 0.1 0 0, Pweight 0 1 0 0, Pweight 0 0 1 0, Pweight 0.5 0 0 0]
pattern_CGA = [Pweight 0 1 0 0, Pweight 0 0 1 0]
pattern_cccccccccc = [cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern]
  where cPattern = Pweight 0 0 0.001 0

test_patterns_basic :: IO ()
test_patterns_basic = do
    matches <- findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual expected matches
  where inputData :: [(Vector Nucleotide, Vector Position)]
        inputData =  [(inputDataSample0, inputDataPositions)]

        inputDataPositions :: Vector Position
        inputDataPositions = V.fromList $ take (V.length inputDataSample0) [0..]

        patterns = [pattern_CG]

        expected = [Match {mPatternId = 0,  mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c,g]}]

test_patterns_basic2 :: IO ()
test_patterns_basic2 = do
    matches <- findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual expected matches
  where inputData :: [(Vector Nucleotide, Vector Position)]
        inputData =  [(inputDataSample0, inputDataPositions)]

        inputDataPositions :: Vector Position
        inputDataPositions = V.fromList $ take (V.length inputDataSample0) (map Position [0..])

        patterns = [pattern_CG, pattern_cCGA]

        expected = [Match {mPatternId = 0, mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c,g]}
                   ,Match {mPatternId = 1, mScore = 961,  mPosition = 3, mSampleId = 0, mMatched = [a, c, g, a]}]

test_patterns_1 :: IO ()
test_patterns_1 = do
    matches <- findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual expected matches
  where numberOfPeople = 1000 :: Int

        inputData :: [(Vector Nucleotide, Vector Position)]
        inputData =  [(inputDataSample0, inputDataPositions), (inputDataSample1, inputDataPositions)] ++ replicate (numberOfPeople - 2) (inputDataSample2, inputDataPositions)

        inputDataPositions :: Vector Position
        inputDataPositions = V.fromList $ take (V.length inputDataSample0) (map Position [0..])

        patterns = [pattern_CG, pattern_cCGA, pattern_CGA] ++ replicate 50 pattern_cccccccccc ++ [pattern_CG]

        expected = [Match {mPatternId = 53, mScore = 1000, mPosition = 0, mSampleId = 1, mMatched = [c, g]}
                   ,Match {mPatternId = 2,  mScore = 1000, mPosition = 0, mSampleId = 1, mMatched = [c, g]}
                   ,Match {mPatternId = 0,  mScore = 1000, mPosition = 0, mSampleId = 1, mMatched = [c, g]}
                   ,Match {mPatternId = 53, mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c, g]}
                   ,Match {mPatternId = 2,  mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c, g]}
                   ,Match {mPatternId = 0,  mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c, g]}
                   ,Match {mPatternId = 1,  mScore = 961,  mPosition = 3, mSampleId = 0, mMatched = [a, c, g, a]}]

-- Check that the Matches contain the reference position, not the real position
test_patterns_2 :: IO ()
test_patterns_2 = do
    matches <- findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual expected matches
  where numberOfPeople = 1000 :: Int

        inputData :: [(Vector Nucleotide, Vector Position)]
        inputData =  [(inputDataSample0, inputDataPositions), (inputDataSample1, V.map (+7) inputDataPositions)] ++ replicate (numberOfPeople - 2) (inputDataSample2, inputDataPositions)

        inputDataPositions :: Vector Position
        inputDataPositions = V.fromList $ take (V.length inputDataSample0) (map Position [10..])

        patterns = [pattern_CG, pattern_cCGA, pattern_CGA] ++ replicate 50 pattern_cccccccccc ++ [pattern_CG]

        expected = [Match {mPatternId = 53, mScore = 1000, mPosition = 10 + 7, mSampleId = 1, mMatched = [c, g]}
                   ,Match {mPatternId = 2,  mScore = 1000, mPosition = 10 + 7, mSampleId = 1, mMatched = [c, g]}
                   ,Match {mPatternId = 0,  mScore = 1000, mPosition = 10 + 7, mSampleId = 1, mMatched = [c, g]}
                   ,Match {mPatternId = 53, mScore = 1000, mPosition = 14,     mSampleId = 0, mMatched = [c, g]}
                   ,Match {mPatternId = 2,  mScore = 1000, mPosition = 14,     mSampleId = 0, mMatched = [c, g]}
                   ,Match {mPatternId = 0,  mScore = 1000, mPosition = 14,     mSampleId = 0, mMatched = [c, g]}
                   ,Match {mPatternId = 1,  mScore = 961,  mPosition = 13,     mSampleId = 0, mMatched = [a, c, g, a]}]

-- Check that the end of the chromosome is padded with N
test_patterns_padding :: IO ()
test_patterns_padding = do
    matches <- findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual expected matches
  where numberOfPeople = 1000 :: Int

        inputData :: [(Vector Nucleotide, Vector Position)]
        inputData =  [(V.take 5 inputDataSample0, V.take 5 inputDataPositions), (inputDataSample1, inputDataPositions)] ++ (take (numberOfPeople - 2) (repeat (V.take 10 inputDataSample2, V.take 10 inputDataPositions)))

        inputDataPositions :: Vector Position
        inputDataPositions = V.fromList $ take (V.length inputDataSample0) (map Position [10..])

        patterns = [pattern_CG, pattern_cCGA, pattern_CGA] ++ replicate 50 pattern_cccccccccc ++ [pattern_CG]

        expected = [Match{mPatternId = 53, mScore = 1000, mPosition = 10, mSampleId = 1, mMatched = [c, g]},
                    Match{mPatternId = 2,  mScore = 1000, mPosition = 10, mSampleId = 1, mMatched = [c, g]},
                    Match{mPatternId = 0,  mScore = 1000, mPosition = 10, mSampleId = 1, mMatched = [c, g]},
                    Match{mPatternId = 53, mScore = 500,  mPosition = 14, mSampleId = 0, mMatched = [c]},
                    Match{mPatternId = 2,  mScore = 500,  mPosition = 14, mSampleId = 0, mMatched = [c]},
                    Match{mPatternId = 0,  mScore = 500,  mPosition = 14, mSampleId = 0, mMatched = [c]}]

refGenome :: Vector Nucleotide
refGenome = V.fromList [65,67,71,84] -- ACGT

-- applyVariants :: V.Vector Nucleotide -> Position -> Position -> [Variant] -> [(Nucleotide, Position)]
-- applyVariants referenceGenome (Position start) (Position end) variants = 
test_apply_variant_1 :: IO ()
test_apply_variant_1 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) []
  assertEqual [(c, Position 1), (g, Position 2)] patched

test_apply_variant_2 :: IO ()
test_apply_variant_2 = do
  let patched = applyVariants refGenome (Position 0) (Position 2) []
  assertEqual [(a, Position 0), (c, Position 1), (g, Position 2)] patched

test_apply_variant_3 :: IO ()
test_apply_variant_3 = do
  let patched = applyVariants refGenome (Position 0) (Position 5) []
  assertEqual [(a, Position 0), (c, Position 1), (g, Position 2), (t, Position 3)] patched

test_apply_variant_4 :: IO ()
test_apply_variant_4 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 100 (U.fromList [a]) (U.fromList [c])]
  assertEqual [(c, Position 1), (g, Position 2)] patched
----prop_reverse :: [Int] -> Bool
----prop_reverse xs = xs == (myReverse (myReverse xs))
test_apply_variant_5 :: IO ()
test_apply_variant_5 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 1 (U.fromList [c]) (U.fromList [n])]
  assertEqual [(n, Position 1), (g, Position 2)] patched

test_apply_variant_6 :: IO ()
test_apply_variant_6 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 2 (U.fromList [g]) (U.fromList [a])]
  assertEqual [(c, Position 1), (a, Position 2)] patched

test_apply_variant_7 :: IO ()
test_apply_variant_7 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 1 (U.fromList [c]) (U.fromList [n]), Diff 2 (U.fromList [g]) (U.fromList [a])]
  assertEqual [(n, Position 1), (a, Position 2)] patched

test_apply_variant_8 :: IO ()
test_apply_variant_8 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 0 (U.fromList [c]) (U.fromList [n]), Diff 4 (U.fromList [g]) (U.fromList [a])]
  assertEqual [(c, Position 1), (g, Position 2)] patched

test_apply_variant_9 :: IO ()
test_apply_variant_9 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 1 (U.fromList [c]) (U.fromList [n,n])]
  assertEqual [(n, Position 1), (n, Position 1), (g, Position 2)] patched

test_apply_variant_10 :: IO ()
test_apply_variant_10 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 2 (U.fromList [c]) (U.fromList [n,n])]
  assertEqual [(c, Position 1), (n, Position 2), (n, Position 2)] patched

test_apply_variant_11 :: IO ()
test_apply_variant_11 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 3 (U.fromList [c]) (U.fromList [n,n])]
  assertEqual [(c, Position 1), (g, Position 2)] patched

test_apply_variant_12 :: IO ()
test_apply_variant_12 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 1 (U.fromList [c, g]) (U.fromList [n])]
  assertEqual [(n, Position 1)] patched

test_apply_variant_13 :: IO ()
test_apply_variant_13 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 2 (U.fromList [g, t]) (U.fromList [n])]
  assertEqual [(c, Position 1), (n, Position 2)] patched

test_apply_variant_14 :: IO ()
test_apply_variant_14 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 0 (U.fromList [a, c, g]) (U.fromList [n, n])]
  assertEqual [(c, Position 1), (g, Position 2)] patched -- For simplicity, do not apply a Diff that starts before the window we're observing

test_parse_1 :: IO ()
test_parse_1 = do
  let line = "chr1\t69081\t1:69081:G:C\tC\tG\t.\t.\t.\tGT\t1/1\t1/1"
  let sampleIdentifiers = G.fromList [SampleId "sample1", SampleId "sample2"]
  assertEqual (Right $ Variant (Chromosome "1") (Position 69080) (Just "1:69081:G:C") (G.fromList [c]) (G.fromList [g]) (G.fromList [Geno11, Geno11]) sampleIdentifiers) (parseVariant sampleIdentifiers line)

test_filterOrderedIntervals_1 :: IO ()
test_filterOrderedIntervals_1 = do
  let inf = [10..]
  let filtered = filterOrderedIntervals Position [(Position 15, Position 18), (Position 20, Position 23)] inf
  assertEqual [15, 16, 17, 18, 20, 21, 22, 23] filtered

test_parseVcfContent_1 :: IO ()
test_parseVcfContent_1 = do
  let vcf = ["#\t\t\t\t\t\t\t\t\tsample1\tsample2",
             "chr1\t4\tname4\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t5\tname5\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t12\tname12\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t15\tname15\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t16\tname16\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t21\tname21\tC\tA\t\t\t\t\t0/0\t0/1"
            ] :: [T.Text]
  let parsed = parseVcfContent [(Position 5, Position 15), (Position 18, Position 25)] vcf
  let var p = Variant (Chromosome "1") (Position (p-1)) (Just $ "name" <> T.pack (show p)) (U.fromList [c]) (U.fromList [a]) (G.fromList [Geno00, Geno01]) (G.fromList [SampleId "sample1", SampleId "sample2"])
  let expected = Right [var 5, var 12, var 15, var 21]
  assertEqual expected parsed

main :: IO ()
main = 
  htfMainWithArgs ["--quiet"] htf_thisModulesTests