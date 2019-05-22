{-# OPTIONS_GHC -F -pgmF htfpp #-}

import Test.Framework

--import           Foreign.C.Types                 (CChar)
import           Data.Vector.Storable            (Vector)
import qualified Data.Vector.Storable            as V
import qualified Data.Vector.Unboxed             as U
import           Data.Monoid                     ((<>))

import Types
import PatternFind
import Lib

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


-- applyVariants :: V.Vector Nucleotide -> Position -> Position -> [Variant] -> [(Nucleotide, Position)]
-- applyVariants referenceGenome (Position start) (Position end) variants = 
test_apply_variant_1 :: IO ()
test_apply_variant_1 = do
  let patched = applyVariants (V.fromList [a,c,g,t]) (Position 1) (Position 2) []
  assertEqual [(c, Position 1), (g, Position 2)] patched

test_apply_variant_2 :: IO ()
test_apply_variant_2 = do
  let patched = applyVariants (V.fromList [a,c,g,t]) (Position 0) (Position 2) []
  assertEqual [(a, Position 0), (c, Position 1), (g, Position 2)] patched

test_apply_variant_3 :: IO ()
test_apply_variant_3 = do
  let patched = applyVariants (V.fromList [a,c,g,t]) (Position 0) (Position 5) []
  assertEqual [(a, Position 0), (c, Position 1), (g, Position 2), (t, Position 3)] patched

test_apply_variant_4 :: IO ()
test_apply_variant_4 = do
  let patched = applyVariants (V.fromList [a,c,g,t]) (Position 1) (Position 2) [Diff 100 (U.fromList [a]) (U.fromList [c])]
  assertEqual [(c, Position 1), (g, Position 2)] patched
----prop_reverse :: [Int] -> Bool
----prop_reverse xs = xs == (myReverse (myReverse xs))
test_apply_variant_5 :: IO ()
test_apply_variant_5 = do
  let patched = applyVariants (V.fromList [a,c,g,t]) (Position 1) (Position 2) [Diff 1 (U.fromList [c]) (U.fromList [n])]
  assertEqual [(n, Position 1), (g, Position 2)] patched

test_apply_variant_6 :: IO ()
test_apply_variant_6 = do
  let patched = applyVariants (V.fromList [a,c,g,t]) (Position 1) (Position 2) [Diff 2 (U.fromList [g]) (U.fromList [a])]
  assertEqual [(c, Position 1), (a, Position 2)] patched

test_apply_variant_7 :: IO ()
test_apply_variant_7 = do
  let patched = applyVariants (V.fromList [a,c,g,t]) (Position 1) (Position 2) [Diff 1 (U.fromList [c]) (U.fromList [n]), Diff 2 (U.fromList [g]) (U.fromList [a])]
  assertEqual [(n, Position 1), (a, Position 2)] patched

test_apply_variant_8 :: IO ()
test_apply_variant_8 = do
  let patched = applyVariants (V.fromList [a,c,g,t]) (Position 1) (Position 2) [Diff 0 (U.fromList [c]) (U.fromList [n]), Diff 4 (U.fromList [g]) (U.fromList [a])]
  assertEqual [(c, Position 1), (g, Position 2)] patched

test_apply_variant_9 :: IO ()
test_apply_variant_9 = do
  let patched = applyVariants (V.fromList [a,c,g,t]) (Position 1) (Position 2) [Diff 1 (U.fromList [c]) (U.fromList [n,n])]
  assertEqual [(n, Position 1), (n, Position 1), (g, Position 2)] patched

test_apply_variant_10 :: IO ()
test_apply_variant_10 = do
  let patched = applyVariants (V.fromList [a,c,g,t]) (Position 1) (Position 2) [Diff 2 (U.fromList [c]) (U.fromList [n,n])]
  assertEqual [(c, Position 1), (n, Position 2), (n, Position 2)] patched

test_apply_variant_11 :: IO ()
test_apply_variant_11 = do
  let patched = applyVariants (V.fromList [a,c,g,t]) (Position 1) (Position 2) [Diff 3 (U.fromList [c]) (U.fromList [n,n])]
  assertEqual [(c, Position 1), (g, Position 2)] patched

test_apply_variant_12 :: IO ()
test_apply_variant_12 = do
  let patched = applyVariants (V.fromList [a,c,g,t]) (Position 1) (Position 2) [Diff 1 (U.fromList [c, g]) (U.fromList [n])]
  assertEqual [(n, Position 1)] patched

test_apply_variant_13 :: IO ()
test_apply_variant_13 = do
  let patched = applyVariants (V.fromList [a,c,g,t]) (Position 1) (Position 2) [Diff 2 (U.fromList [g, t]) (U.fromList [n])]
  assertEqual [(c, Position 1), (n, Position 2)] patched

test_apply_variant_14 :: IO ()
test_apply_variant_14 = do
  let patched = applyVariants (V.fromList [a,c,g,t]) (Position 1) (Position 2) [Diff 0 (U.fromList [a, c, g]) (U.fromList [n, n])]
  assertEqual [(c, Position 1), (g, Position 2)] patched -- For simplicity, do not apply a Diff that starts before the window we're observing

main :: IO ()
main = 
  htfMain htf_thisModulesTests