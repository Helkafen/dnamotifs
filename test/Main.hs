{-# OPTIONS_GHC -F -pgmF htfpp #-}

import Test.Framework

import           Foreign.C.Types                 (CInt, CChar)
import           Data.Vector.Storable            (Vector)
import qualified Data.Vector.Storable            as V
import           Data.Monoid                     ((<>))

import Types
import Lib

inputDataSample0 :: Vector CChar
inputDataSample0 = V.fromList [a,a,a,a,c,g,a] <> V.fromList (take 93 (repeat a))

inputDataSample1 :: Vector CChar
inputDataSample1 = V.fromList [c,g,a,a,a,a,a] <> V.fromList (take 93 (repeat a)) 

inputDataSample2 :: Vector CChar
inputDataSample2 = V.fromList [a,a,a,a,a,a,a] <> V.fromList (take 93 (repeat a))

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
  where numberOfPeople = 1000 :: Int

        inputData :: [(Vector CChar, Vector CInt)]
        inputData =  [(inputDataSample0, inputDataPositions)]

        inputDataPositions :: Vector CInt
        inputDataPositions = V.fromList $ take (V.length inputDataSample0) [0..]

        patterns = [pattern_CG]

        expected = [Match {mPatternId = 0,  mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c,g]}]

test_patterns_basic2 :: IO ()
test_patterns_basic2 = do
    matches <- findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual expected matches
  where numberOfPeople = 1000 :: Int

        inputData :: [(Vector CChar, Vector CInt)]
        inputData =  [(inputDataSample0, inputDataPositions)]

        inputDataPositions :: Vector CInt
        inputDataPositions = V.fromList $ take (V.length inputDataSample0) [0..]

        patterns = [pattern_CG, pattern_cCGA]

        expected = [Match {mPatternId = 0, mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c,g]}
                   ,Match {mPatternId = 1, mScore = 961,  mPosition = 3, mSampleId = 0, mMatched = [a, c, g, a]}]

test_patterns_1 :: IO ()
test_patterns_1 = do
    matches <- findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual expected matches
  where numberOfPeople = 1000 :: Int

        inputData :: [(Vector CChar, Vector CInt)]
        inputData =  [(inputDataSample0, inputDataPositions), (inputDataSample1, inputDataPositions)] ++ (take (numberOfPeople - 2) (repeat (inputDataSample2, inputDataPositions)))

        inputDataPositions :: Vector CInt
        inputDataPositions = V.fromList $ take (V.length inputDataSample0) [0..]

        patterns = [pattern_CG, pattern_cCGA, pattern_CGA] ++ take 50 (repeat pattern_cccccccccc) ++ [pattern_CG]

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

        inputData :: [(Vector CChar, Vector CInt)]
        inputData =  [(inputDataSample0, inputDataPositions), (inputDataSample1, V.map (+7) inputDataPositions)] ++ (take (numberOfPeople - 2) (repeat (inputDataSample2, inputDataPositions)))

        inputDataPositions :: Vector CInt
        inputDataPositions = V.fromList $ take (V.length inputDataSample0) [10..]

        patterns = [pattern_CG, pattern_cCGA, pattern_CGA] ++ take 50 (repeat pattern_cccccccccc) ++ [pattern_CG]

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

        inputData :: [(Vector CChar, Vector CInt)]
        inputData =  [(V.take 5 inputDataSample0, V.take 5 inputDataPositions), (inputDataSample1, inputDataPositions)] ++ (take (numberOfPeople - 2) (repeat (V.take 10 inputDataSample2, V.take 10 inputDataPositions)))

        inputDataPositions :: Vector CInt
        inputDataPositions = V.fromList $ take (V.length inputDataSample0) [10..]

        patterns = [pattern_CG, pattern_cCGA, pattern_CGA] ++ take 50 (repeat pattern_cccccccccc) ++ [pattern_CG]

        expected = [Match{mPatternId = 53, mScore = 1000, mPosition = 10, mSampleId = 1, mMatched = [c, g]},
                    Match{mPatternId = 2,  mScore = 1000, mPosition = 10, mSampleId = 1, mMatched = [c, g]},
                    Match{mPatternId = 0,  mScore = 1000, mPosition = 10, mSampleId = 1, mMatched = [c, g]},
                    Match{mPatternId = 53, mScore = 500,  mPosition = 14, mSampleId = 0, mMatched = [c]},
                    Match{mPatternId = 2,  mScore = 500,  mPosition = 14, mSampleId = 0, mMatched = [c]},
                    Match{mPatternId = 0,  mScore = 500,  mPosition = 14, mSampleId = 0, mMatched = [c]}]


----prop_reverse :: [Int] -> Bool
----prop_reverse xs = xs == (myReverse (myReverse xs))

main :: IO ()
main = htfMain htf_thisModulesTests