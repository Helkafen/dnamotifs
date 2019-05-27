{-# OPTIONS_GHC -F -pgmF htfpp #-}
{-# LANGUAGE OverloadedStrings #-}

import Test.Framework

--import           Foreign.C.Types                 (CChar)
import           Data.Vector.Storable            (Vector)
import qualified Data.Vector.Storable            as V
import qualified Data.Vector.Generic             as G
import           Data.Monoid                     ((<>))
import qualified Data.Text                       as T
import qualified Data.ByteString                 as B
import           Foreign.C.Types                 (CInt)

import Types
import PatternFind
import Lib
import Fasta (takeRef)
import Vcf (parseVariant, filterOrderedIntervals, parseVcfContent)
import Bed (parseBedContent)

inputDataSample0 :: B.ByteString
inputDataSample0 = B.pack [a,a,a,a,c,g,a] <> B.replicate 93 a

inputDataSample1 :: B.ByteString
inputDataSample1 = B.pack [c,g,a,a,a,a,a] <> B.replicate 93 a

inputDataSample2 :: B.ByteString
inputDataSample2 = B.pack [a,a,a,a,a,a,a] <> B.replicate 93 a

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
  where inputData =  [BaseSequencePosition inputDataSample0 inputDataPositions]

        inputDataPositions :: Vector CInt
        inputDataPositions = V.fromList $ map fromIntegral $ take (B.length inputDataSample0) [0..]

        patterns = [pattern_CG]

        expected = [Match {mPatternId = 0,  mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c,g]}]

test_patterns_basic2 :: IO ()
test_patterns_basic2 = do
    matches <- findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual expected matches
  where inputData =  [BaseSequencePosition inputDataSample0 inputDataPositions]
        inputDataPositions = V.fromList $ take (B.length inputDataSample0) (map fromIntegral [0..])

        patterns = [pattern_CG, pattern_cCGA]

        expected = [Match {mPatternId = 0, mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c,g]}
                   ,Match {mPatternId = 1, mScore = 961,  mPosition = 3, mSampleId = 0, mMatched = [a, c, g, a]}]

test_patterns_1 :: IO ()
test_patterns_1 = do
    matches <- findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual expected matches
  where numberOfPeople = 1000 :: Int

        inputData =  [BaseSequencePosition inputDataSample0 inputDataPositions, BaseSequencePosition inputDataSample1 inputDataPositions] ++ replicate (numberOfPeople - 2) (BaseSequencePosition inputDataSample2 inputDataPositions)

        inputDataPositions = V.fromList $ take (B.length inputDataSample0) (map fromIntegral [0..])

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

        inputData =  [BaseSequencePosition inputDataSample0 inputDataPositions, BaseSequencePosition inputDataSample1 (V.map (+7) inputDataPositions)] ++ replicate (numberOfPeople - 2) (BaseSequencePosition inputDataSample2 inputDataPositions)

        inputDataPositions = V.fromList $ take (B.length inputDataSample0) (map fromIntegral [10..])

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

        inputData =  [BaseSequencePosition (B.take 5 inputDataSample0) (V.take 5 inputDataPositions), BaseSequencePosition inputDataSample1 inputDataPositions] ++ replicate (numberOfPeople - 2) (BaseSequencePosition (B.take 10 inputDataSample2) (V.take 10 inputDataPositions))
        inputDataPositions = V.fromList $ take (B.length inputDataSample0) (map fromIntegral [10..])

        patterns = [pattern_CG, pattern_cCGA, pattern_CGA] ++ replicate 50 pattern_cccccccccc ++ [pattern_CG]

        expected = [Match{mPatternId = 53, mScore = 1000, mPosition = 10, mSampleId = 1, mMatched = [c, g]},
                    Match{mPatternId = 2,  mScore = 1000, mPosition = 10, mSampleId = 1, mMatched = [c, g]},
                    Match{mPatternId = 0,  mScore = 1000, mPosition = 10, mSampleId = 1, mMatched = [c, g]},
                    Match{mPatternId = 53, mScore = 500,  mPosition = 14, mSampleId = 0, mMatched = [c]},
                    Match{mPatternId = 2,  mScore = 500,  mPosition = 14, mSampleId = 0, mMatched = [c]},
                    Match{mPatternId = 0,  mScore = 500,  mPosition = 14, mSampleId = 0, mMatched = [c]}]

refGenome :: Int -> Int -> BaseSequencePosition
refGenome = takeRef (B.pack [65,67,71,84]) -- ACGT

-- applyVariants :: V.Vector Nucleotide -> Position -> Position -> [Variant] -> [(Nucleotide, Position)]
-- applyVariants referenceGenome (Position start) (Position end) variants = 
test_apply_variant_1 :: IO ()
test_apply_variant_1 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) []
  assertEqual (BaseSequencePosition (B.pack [c,g]) (V.fromList [1,2])) patched

test_apply_variant_2 :: IO ()
test_apply_variant_2 = do
  let patched = applyVariants refGenome (Position 0) (Position 2) []
  assertEqual (BaseSequencePosition (B.pack [a, c, g]) (V.fromList [0,1,2])) patched

test_apply_variant_3 :: IO ()
test_apply_variant_3 = do
  let patched = applyVariants refGenome (Position 0) (Position 5) []
  assertEqual (BaseSequencePosition (B.pack [a,c,g,t]) (V.fromList [0,1,2,3])) patched

test_apply_variant_4 :: IO ()
test_apply_variant_4 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 100 (B.pack [a]) (B.pack [c])]
  assertEqual (BaseSequencePosition (B.pack [c,g]) (V.fromList [1,2])) patched
----prop_reverse :: [Int] -> Bool
----prop_reverse xs = xs == (myReverse (myReverse xs))
test_apply_variant_5 :: IO ()
test_apply_variant_5 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 1 (B.pack [c]) (B.pack [n])]
  assertEqual (BaseSequencePosition (B.pack [n,g]) (V.fromList [1,2])) patched

test_apply_variant_6 :: IO ()
test_apply_variant_6 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 2 (B.pack [g]) (B.pack [a])]
  assertEqual (BaseSequencePosition (B.pack [c,a]) (V.fromList [1,2])) patched

test_apply_variant_7 :: IO ()
test_apply_variant_7 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 1 (B.pack [c]) (B.pack [n]), Diff 2 (B.pack [g]) (B.pack [a])]
  assertEqual (BaseSequencePosition (B.pack [n,a]) (V.fromList [1,2])) patched

test_apply_variant_8 :: IO ()
test_apply_variant_8 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 0 (B.pack [c]) (B.pack [n]), Diff 4 (B.pack [g]) (B.pack [a])]
  assertEqual (BaseSequencePosition (B.pack [c,g]) (V.fromList [1,2])) patched

test_apply_variant_9 :: IO ()
test_apply_variant_9 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 1 (B.pack [c]) (B.pack [n,n])]
  assertEqual (BaseSequencePosition (B.pack[n,n,g]) (V.fromList [1,1,2])) patched

test_apply_variant_10 :: IO ()
test_apply_variant_10 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 2 (B.pack [c]) (B.pack [n,n])]
  assertEqual (BaseSequencePosition (B.pack[c,n,n]) (V.fromList[1,2,2])) patched

test_apply_variant_11 :: IO ()
test_apply_variant_11 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 3 (B.pack [c]) (B.pack [n,n])]
  assertEqual (BaseSequencePosition (B.pack[c,g]) (V.fromList [1,2])) patched

test_apply_variant_12 :: IO ()
test_apply_variant_12 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 1 (B.pack [c, g]) (B.pack [n])]
  assertEqual (BaseSequencePosition (B.pack[n]) (V.fromList [1])) patched

test_apply_variant_13 :: IO ()
test_apply_variant_13 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 2 (B.pack [g, t]) (B.pack [n])]
  assertEqual (BaseSequencePosition (B.pack[c,n]) (V.fromList [1,2])) patched

test_apply_variant_14 :: IO ()
test_apply_variant_14 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 0 (B.pack [a, c, g]) (B.pack [n, n])]
  assertEqual (BaseSequencePosition (B.pack[c,g]) (V.fromList [1,2])) patched -- For simplicity, do not apply a Diff that starts before the window we're observing

test_parse_1 :: IO ()
test_parse_1 = do
  let line = "chr1\t69081\t1:69081:G:C\tC\tG\t.\t.\t.\tGT\t1/1\t1/1"
  let sampleIdentifiers = G.fromList [SampleId "sample1", SampleId "sample2"]
  assertEqual (Right $ Variant (Chromosome "1") (Position 69080) (Just "1:69081:G:C") (B.pack [c]) (B.pack [g]) (G.fromList [geno11, geno11]) sampleIdentifiers) (parseVariant sampleIdentifiers line)

test_filterOrderedIntervals_1 :: IO ()
test_filterOrderedIntervals_1 = do
  let inf = [10..]
  let filtered = filterOrderedIntervals Position [(Position 15, Position 18), (Position 20, Position 23)] inf
  assertEqual [15, 16, 17, 18, 20, 21, 22, 23] filtered

test_parseVcfContent_1 :: IO ()
test_parseVcfContent_1 = do
  let vcf = ["#\t\t\t\t\t\t\t\t\tsample1\tsample2",
             "chr1\t5\tname4\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t6\tname5\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t13\tname12\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t16\tname15\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t18\tname17\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t22\tname21\tC\tA\t\t\t\t\t0/0\t0/1"
            ] :: [T.Text]
  let parsed = parseVcfContent [(Position 5, Position 15), (Position 18, Position 25)] vcf
  let var p = Variant (Chromosome "1") (Position p) (Just $ "name" <> T.pack (show p)) (B.pack [c]) (B.pack [a]) (G.fromList [geno00, geno01]) (G.fromList [SampleId "sample1", SampleId "sample2"])
  let expected = Right [var 5, var 12, var 15, var 21]
  assertEqual expected parsed

test_parseBed_1 :: IO ()
test_parseBed_1 = do
  let bedContent = "chr1\t5\t6\nchr1\t8\t10\nchr2\t100\t110" :: T.Text
  let parsed = parseBedContent (Chromosome "1") bedContent
  let expected = Right [(Position 5, Position 6), (Position 8, Position 10)]
  assertEqual expected parsed

main :: IO ()
main = 
  htfMainWithArgs ["--quiet"] htf_thisModulesTests