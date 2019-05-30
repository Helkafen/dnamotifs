{-# OPTIONS_GHC -F -pgmF htfpp #-}
{-# LANGUAGE OverloadedStrings #-}

import Test.Framework

import qualified Data.Vector.Storable            as STO
import qualified Data.Vector                     as V
import qualified Data.Vector.Generic             as G
import           Data.Monoid                     ((<>))
import qualified Data.Text                       as T
import qualified Data.ByteString                 as B
import           Foreign.C.Types                 (CInt)
import           Data.Text.Encoding              (encodeUtf8)
import qualified Data.Map

import Types
import PatternFind
import Lib
import Fasta (takeRef)
import Vcf (parseVariant, filterOrderedIntervals, parseVcfContent)
import Bed (parseBedContent)
import MotifDefinition (parseMotifsContent)

mkSeq :: [Nucleotide] -> B.ByteString
mkSeq = B.pack . map unNuc

inputDataSample0 :: B.ByteString
inputDataSample0 = mkSeq [a,a,a,a,c,g,a] <> B.replicate 93 (unNuc a)

inputDataSample1 :: B.ByteString
inputDataSample1 = mkSeq [c,g,a,a,a,a,a] <> B.replicate 93 (unNuc a)

inputDataSample2 :: B.ByteString
inputDataSample2 = mkSeq [a,a,a,a,a,a,a] <> B.replicate 93 (unNuc a)

pattern_CG, pattern_cCGA, pattern_CGA, pattern_cccccccccc :: [Pweight]
pattern_CG = [Pweight 0 1 0 0, Pweight 0 0 1 0]
pattern_cCGA = [Pweight 0 0.1 0 0, Pweight 0 1 0 0, Pweight 0 0 1 0, Pweight 0.5 0 0 0]
pattern_CGA = [Pweight 0 1 0 0, Pweight 0 0 1 0]
pattern_cccccccccc = [cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern]
  where cPattern = Pweight 0 0 0.001 0

test_patterns_basic :: IO ()
test_patterns_basic = do
    matches <- findPatternsInBlock 500 (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual expected matches
  where inputData =  [BaseSequencePosition inputDataSample0 inputDataPositions]

        inputDataPositions :: STO.Vector CInt
        inputDataPositions = STO.fromList $ map fromIntegral $ take (B.length inputDataSample0) [0..]

        patterns = [pattern_CG]

        expected = V.fromList [Match {mPatternId = 0,  mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c,g]}]

test_patterns_basic2 :: IO ()
test_patterns_basic2 = do
    matches <- findPatternsInBlock 500 (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual expected matches
  where inputData =  [BaseSequencePosition inputDataSample0 inputDataPositions]
        inputDataPositions = STO.fromList $ take (B.length inputDataSample0) (map fromIntegral [0..])

        patterns = [pattern_CG, pattern_cCGA]

        expected = V.fromList [Match {mPatternId = 1, mScore = 961,  mPosition = 3, mSampleId = 0, mMatched = [a, c, g, a]},
                                    Match {mPatternId = 0, mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c,g]}]


test_patterns_1 :: IO ()
test_patterns_1 = do
    matches <- findPatternsInBlock 500 (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual (V.fromList  expected) matches
  where numberOfPeople = 1000 :: Int

        inputData =  [BaseSequencePosition inputDataSample0 inputDataPositions, BaseSequencePosition inputDataSample1 inputDataPositions] ++ replicate (numberOfPeople - 2) (BaseSequencePosition inputDataSample2 inputDataPositions)

        inputDataPositions = STO.fromList $ take (B.length inputDataSample0) (map fromIntegral [0..])

        patterns = [pattern_CG, pattern_cCGA, pattern_CGA] ++ replicate 50 pattern_cccccccccc ++ [pattern_CG]

        expected = [Match {mPatternId = 1,  mScore = 961,  mPosition = 3, mSampleId = 0, mMatched = [a, c, g, a]},
                    Match {mPatternId = 0,  mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c, g]},
                    Match {mPatternId = 2,  mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c, g]},
                    Match {mPatternId = 53, mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c, g]},
                    Match {mPatternId = 0,  mScore = 1000, mPosition = 0, mSampleId = 1, mMatched = [c, g]},
                    Match {mPatternId = 2,  mScore = 1000, mPosition = 0, mSampleId = 1, mMatched = [c, g]},
                    Match {mPatternId = 53, mScore = 1000, mPosition = 0, mSampleId = 1, mMatched = [c, g]}]


-- Check that the Matches contain the reference position, not the real position
test_patterns_2 :: IO ()
test_patterns_2 = do
    matches <- findPatternsInBlock 500 (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual (V.fromList  expected) matches
  where numberOfPeople = 1000 :: Int

        inputData =  [BaseSequencePosition inputDataSample0 inputDataPositions, BaseSequencePosition inputDataSample1 (STO.map (+7) inputDataPositions)] ++ replicate (numberOfPeople - 2) (BaseSequencePosition inputDataSample2 inputDataPositions)

        inputDataPositions = STO.fromList $ take (B.length inputDataSample0) (map fromIntegral [10..])

        patterns = [pattern_CG, pattern_cCGA, pattern_CGA] ++ replicate 50 pattern_cccccccccc ++ [pattern_CG]

        expected = [Match {mPatternId = 1,  mScore = 961,  mPosition = 13,     mSampleId = 0, mMatched = [a, c, g, a]},
                    Match {mPatternId = 0,  mScore = 1000, mPosition = 14,     mSampleId = 0, mMatched = [c, g]},
                    Match {mPatternId = 2,  mScore = 1000, mPosition = 14,     mSampleId = 0, mMatched = [c, g]},
                    Match {mPatternId = 53, mScore = 1000, mPosition = 14,     mSampleId = 0, mMatched = [c, g]},
                    Match {mPatternId = 0,  mScore = 1000, mPosition = 10 + 7, mSampleId = 1, mMatched = [c, g]},
                    Match {mPatternId = 2,  mScore = 1000, mPosition = 10 + 7, mSampleId = 1, mMatched = [c, g]},
                    Match {mPatternId = 53, mScore = 1000, mPosition = 10 + 7, mSampleId = 1, mMatched = [c, g]}]


-- Check that the end of the chromosome is padded with N
test_patterns_padding :: IO ()
test_patterns_padding = do
    matches <- findPatternsInBlock 500 (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual (V.fromList  expected) matches
  where numberOfPeople = 1000 :: Int

        inputData =  [BaseSequencePosition (B.take 5 inputDataSample0) (STO.take 5 inputDataPositions), BaseSequencePosition inputDataSample1 inputDataPositions] ++ replicate (numberOfPeople - 2) (BaseSequencePosition (B.take 10 inputDataSample2) (STO.take 10 inputDataPositions))
        inputDataPositions = STO.fromList $ take (B.length inputDataSample0) (map fromIntegral [10..])

        patterns = [pattern_CG, pattern_cCGA, pattern_CGA] ++ replicate 50 pattern_cccccccccc ++ [pattern_CG]

        expected = [Match{mPatternId = 0,  mScore = 500,  mPosition = 14, mSampleId = 0, mMatched = [c]},
                    Match{mPatternId = 2,  mScore = 500,  mPosition = 14, mSampleId = 0, mMatched = [c]},
                    Match{mPatternId = 53, mScore = 500,  mPosition = 14, mSampleId = 0, mMatched = [c]},
                    Match{mPatternId = 0,  mScore = 1000, mPosition = 10, mSampleId = 1, mMatched = [c, g]},
                    Match{mPatternId = 2,  mScore = 1000, mPosition = 10, mSampleId = 1, mMatched = [c, g]},
                    Match{mPatternId = 53, mScore = 1000, mPosition = 10, mSampleId = 1, mMatched = [c, g]}]


refGenome :: Int -> Int -> BaseSequencePosition
refGenome = takeRef (B.pack [65,67,71,84]) -- ACGT


test_apply_variant_1 :: IO ()
test_apply_variant_1 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) []
  assertEqual (BaseSequencePosition (mkSeq [c,g]) (STO.fromList [1,2])) patched

test_apply_variant_2 :: IO ()
test_apply_variant_2 = do
  let patched = applyVariants refGenome (Position 0) (Position 2) []
  assertEqual (BaseSequencePosition (mkSeq [a, c, g]) (STO.fromList [0,1,2])) patched

test_apply_variant_3 :: IO ()
test_apply_variant_3 = do
  let patched = applyVariants refGenome (Position 0) (Position 5) []
  assertEqual (BaseSequencePosition (mkSeq [a,c,g,t]) (STO.fromList [0,1,2,3])) patched

test_apply_variant_4 :: IO ()
test_apply_variant_4 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 100 (mkSeq [a]) (mkSeq [c])]
  assertEqual (BaseSequencePosition (mkSeq [c,g]) (STO.fromList [1,2])) patched
----prop_reverse :: [Int] -> Bool
----prop_reverse xs = xs == (myReverse (myReverse xs))
test_apply_variant_5 :: IO ()
test_apply_variant_5 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 1 (mkSeq [c]) (mkSeq [n])]
  assertEqual (BaseSequencePosition (mkSeq [n,g]) (STO.fromList [1,2])) patched

test_apply_variant_6 :: IO ()
test_apply_variant_6 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 2 (mkSeq [g]) (mkSeq [a])]
  assertEqual (BaseSequencePosition (mkSeq [c,a]) (STO.fromList [1,2])) patched

test_apply_variant_7 :: IO ()
test_apply_variant_7 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 1 (mkSeq [c]) (mkSeq [n]), Diff 2 (mkSeq [g]) (mkSeq [a])]
  assertEqual (BaseSequencePosition (mkSeq [n,a]) (STO.fromList [1,2])) patched

test_apply_variant_8 :: IO ()
test_apply_variant_8 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 0 (mkSeq [c]) (mkSeq [n]), Diff 4 (mkSeq [g]) (mkSeq [a])]
  assertEqual (BaseSequencePosition (mkSeq [c,g]) (STO.fromList [1,2])) patched

test_apply_variant_9 :: IO ()
test_apply_variant_9 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 1 (mkSeq [c]) (mkSeq [n,n])]
  assertEqual (BaseSequencePosition (mkSeq[n,n,g]) (STO.fromList [1,1,2])) patched

test_apply_variant_10 :: IO ()
test_apply_variant_10 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 2 (mkSeq [c]) (mkSeq [n,n])]
  assertEqual (BaseSequencePosition (mkSeq[c,n,n]) (STO.fromList[1,2,2])) patched

test_apply_variant_11 :: IO ()
test_apply_variant_11 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 3 (mkSeq [c]) (mkSeq [n,n])]
  assertEqual (BaseSequencePosition (mkSeq[c,g]) (STO.fromList [1,2])) patched

test_apply_variant_12 :: IO ()
test_apply_variant_12 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 1 (mkSeq [c, g]) (mkSeq [n])]
  assertEqual (BaseSequencePosition (mkSeq[n]) (STO.fromList [1])) patched

test_apply_variant_13 :: IO ()
test_apply_variant_13 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 2 (mkSeq [g, t]) (mkSeq [n])]
  assertEqual (BaseSequencePosition (mkSeq[c,n]) (STO.fromList [1,2])) patched

test_apply_variant_14 :: IO ()
test_apply_variant_14 = do
  let patched = applyVariants refGenome (Position 1) (Position 2) [Diff 0 (mkSeq [a, c, g]) (mkSeq [n, n])]
  assertEqual (BaseSequencePosition (mkSeq[c,g]) (STO.fromList [1,2])) patched -- For simplicity, do not apply a Diff that starts before the window we're observing

test_parse_1 :: IO ()
test_parse_1 = do
  let line = "chr1\t69081\t1:69081:G:C\tC\tG\t.\t.\t.\tGT\t1/1\t1/1"
  let sampleIdentifiers = G.fromList [SampleId "sample1", SampleId "sample2"]
  assertEqual (Right $ Variant (Chromosome "1") (Position 69080) (Just "1:69081:G:C") (mkSeq [c]) (mkSeq [g]) (G.fromList [geno11, geno11]) sampleIdentifiers) (parseVariant sampleIdentifiers line)

test_filterOrderedIntervals_1 :: IO ()
test_filterOrderedIntervals_1 = do
  let inf = [10..]
  let filtered = filterOrderedIntervals Position [(Position 15, Position 18), (Position 20, Position 23)] inf
  assertEqual [15, 16, 17, 18, 20, 21, 22, 23] filtered

test_parseVcfContent_1 :: IO ()
test_parseVcfContent_1 = do
  let vcf = map encodeUtf8 ["#\t\t\t\t\t\t\t\t\tsample1\tsample2",
             "chr1\t5\tname4\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t6\tname5\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t13\tname12\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t16\tname15\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t18\tname17\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t22\tname21\tC\tA\t\t\t\t\t0/0\t0/1"
            ] :: [B.ByteString]
  let parsed = parseVcfContent [(Position 5, Position 15), (Position 18, Position 25)] vcf
  let var p = Variant (Chromosome "1") (Position p) (Just $ "name" <> T.pack (show p)) (mkSeq [c]) (mkSeq [a]) (G.fromList [geno00, geno01]) (G.fromList [SampleId "sample1", SampleId "sample2"])
  let expected = Right [var 5, var 12, var 15, var 21]
  assertEqual expected parsed

test_parseBed_1 :: IO ()
test_parseBed_1 = do
  let bedContent = "chr1\t5\t6\nchr1\t8\t10\nchr2\t100\t110" :: T.Text
  let parsed = parseBedContent (Chromosome "1") bedContent
  let expected = Right [(Position 5, Position 6), (Position 8, Position 10)]
  assertEqual expected parsed

motifString :: T.Text
motifString = T.replace "   " "\t" (motifString1 <> motifString2)
  where motifString1 = ">ATGACTCATC\tAP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer        6.049537        -1.782996e+03   0       9805.3,5781.0,3085.1,5.0,-- 0.00e+00\n0.419\t0.275\t0.277\t0.028\n0.001\t0.001\t0.001\t0.997\n"
        motifString2 = ">SCCTSAGGSCAW\tAP-2gamma(AP2)/MCF7-TFAP2C-ChIP-Seq(GSE21234)/Homer     6.349794        -24627.169865   0       T:26194.0(44.86%),413.7-- (9.54%),P:1e-10695\n0.005\t0.431\t0.547\t0.017\n0.001\t0.997\t0.001\t0.001\n0.001\t0.947\t0.001\t0.051"


test_motif_parser_1 :: IO ()
test_motif_parser_1 = do
  let parsed = parseMotifsContent motifString
  let expected1 = ("AP-1", [Pweight{wa = 0.419, wc = 0.275, wg = 0.277, wt = 2.8e-2}, Pweight{wa = 1.0e-3, wc = 1.0e-3, wg = 1.0e-3, wt = 0.997}])
  let expected2 = ("AP-2gamma", [Pweight{wa = 5.0e-3, wc = 0.431, wg = 0.547, wt = 1.7e-2}
                               , Pweight{wa = 1.0e-3, wc = 0.997, wg = 1.0e-3, wt = 1.0e-3}
                               , Pweight{wa = 1.0e-3, wc = 0.947, wg = 1.0e-3, wt = 5.1e-2}])

  assertEqual (Right [expected1, expected2]) parsed

test_variantsToDiffs_1 :: IO ()
test_variantsToDiffs_1 = do
  let sampleIdentifiers = G.fromList [SampleId "sample1", SampleId "sample2"]
  let variant = Variant (Chromosome "1") (Position 69080) (Just "1:69081:G:C") (mkSeq [c]) (mkSeq [g]) (G.fromList [geno10, geno01]) sampleIdentifiers
  let diffs = variantsToDiffs [variant]
  let expected = Data.Map.fromList [((0,HaploLeft), V.fromList [Diff (Position 69080) (mkSeq [c]) (mkSeq [g])])
                                   ,((1,HaploRight),V.fromList [Diff (Position 69080) (mkSeq [c]) (mkSeq [g])])] :: Data.Map.Map (Int, Haplotype) (V.Vector Diff)
  assertEqual expected diffs

test_variantsToDiffs_2 :: IO ()
test_variantsToDiffs_2 = do
  let sampleIdentifiers = G.fromList [SampleId "sample1", SampleId "sample2"]
  let variant1 = Variant (Chromosome "1") (Position 69080) (Just "1:69081:G:C") (mkSeq [c]) (mkSeq [g]) (G.fromList [geno10, geno01]) sampleIdentifiers
  let variant2 = Variant (Chromosome "1") (Position 69079) (Just "1:69079:A:T") (mkSeq [a]) (mkSeq [t]) (G.fromList [geno10, geno00]) sampleIdentifiers
  let variant3 = Variant (Chromosome "1") (Position 69078) (Just "1:69078:T:G") (mkSeq [t]) (mkSeq [g]) (G.fromList [geno01, geno00]) sampleIdentifiers
  let diffs = variantsToDiffs [variant1, variant2, variant3]
  let diff1 = Diff (Position 69080) (mkSeq [c]) (mkSeq [g])
  let diff2 = Diff (Position 69079) (mkSeq [a]) (mkSeq [t])
  let diff3 = Diff (Position 69080) (mkSeq [c]) (mkSeq [g])
  let diff4 = Diff (Position 69078) (mkSeq [t]) (mkSeq [g])
  let expected = Data.Map.fromList [((0,HaploLeft), V.fromList [diff2, diff1]) -- Order of positions must be increasing
                                   ,((0,HaploRight), V.fromList [diff4])
                                   ,((1,HaploRight),V.fromList [diff3])] :: Data.Map.Map (Int, Haplotype) (V.Vector Diff)
  assertEqual expected diffs

main :: IO ()
main = 
  htfMainWithArgs ["--quiet"] htf_thisModulesTests