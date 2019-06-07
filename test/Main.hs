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
import           Data.Char                       (ord)

import Types
import PatternFind
import Lib
import Fasta (takeRef)
import Vcf (parseVariant, filterOrderedIntervals, parseVcfContent, fillVector)
import Bed (parseBedContent)
import MotifDefinition (parseHocomocoMotifsContent)

mkSeq :: [Nucleotide] -> B.ByteString
mkSeq = B.pack . map unNuc

inputDataSample0 :: B.ByteString
inputDataSample0 = mkSeq [a,a,a,a,c,g,a] <> B.replicate 93 (unNuc a)

inputDataSample1 :: B.ByteString
inputDataSample1 = mkSeq [c,g,a,a,a,a,a] <> B.replicate 93 (unNuc a)

inputDataSample2 :: B.ByteString
inputDataSample2 = mkSeq [a,a,a,a,a,a,a] <> B.replicate 93 (unNuc a)

pattern_CG, pattern_cCGA, pattern_cccccccccc :: [Pweight]
pattern_CG = [Pweight 0 1000 0 0, Pweight 0 0 1000 0]
pattern_cCGA = [Pweight 0 100 0 0, Pweight 0 1000 0 0, Pweight 0 0 1000 0, Pweight 500 0 0 0]
pattern_cccccccccc = [cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern]
  where cPattern = Pweight 0 0 1 0

test_patterns_basic :: IO ()
test_patterns_basic = do
    matches <- findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual expected matches
  where inputData =  [BaseSequencePosition inputDataSample0 inputDataPositions]

        inputDataPositions :: STO.Vector CInt
        inputDataPositions = STO.fromList $ map fromIntegral $ take (B.length inputDataSample0) [(0::CInt)..]

        patterns = [Pattern 2000 pattern_CG]

        expected = V.fromList [Match {mPatternId = 0,  mScore = 2000, mPosition = 4, mSampleId = 0, mMatched = [c,g]}]

test_patterns_basic2 :: IO ()
test_patterns_basic2 = do
    matches <- findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual expected matches
  where inputData =  [BaseSequencePosition inputDataSample0 inputDataPositions]
        inputDataPositions = STO.fromList $ take (B.length inputDataSample0) [(0::CInt)..]

        patterns = [Pattern 2000 pattern_CG, Pattern 900 pattern_cCGA]

        expected = V.fromList [Match {mPatternId = 1, mScore = 2500,  mPosition = 3, mSampleId = 0, mMatched = [a, c, g, a]},
                               Match {mPatternId = 0, mScore = 2000, mPosition = 4, mSampleId = 0, mMatched = [c,g]}]


test_patterns_1 :: IO ()
test_patterns_1 = do
    matches <- findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual (V.fromList  expected) matches
  where numberOfPeople = 1000 :: Int

        inputData =  [BaseSequencePosition inputDataSample0 inputDataPositions, BaseSequencePosition inputDataSample1 inputDataPositions] ++ replicate (numberOfPeople - 2) (BaseSequencePosition inputDataSample2 inputDataPositions)

        inputDataPositions = STO.fromList $ take (B.length inputDataSample0) [(0::CInt)..]

        patterns = [Pattern 2000 pattern_CG, Pattern 900 pattern_cCGA] ++ replicate 50 (Pattern 10000 pattern_cccccccccc) ++ [Pattern 2000 pattern_CG]

        expected = [Match {mPatternId = 1,  mScore = 2500, mPosition = 3, mSampleId = 0, mMatched = [a, c, g, a]},
                    Match {mPatternId = 0,  mScore = 2000, mPosition = 4, mSampleId = 0, mMatched = [c, g]},
                    Match {mPatternId = 52, mScore = 2000, mPosition = 4, mSampleId = 0, mMatched = [c, g]},
                    Match {mPatternId = 0,  mScore = 2000, mPosition = 0, mSampleId = 1, mMatched = [c, g]},
                    Match {mPatternId = 52, mScore = 2000, mPosition = 0, mSampleId = 1, mMatched = [c, g]}]


-- Check that the Matches contain the reference position, not the real position
test_patterns_2 :: IO ()
test_patterns_2 = do
    matches <- findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual (V.fromList  expected) matches
  where numberOfPeople = 1000 :: Int

        inputData =  [BaseSequencePosition inputDataSample0 inputDataPositions, BaseSequencePosition inputDataSample1 (STO.map (+7) inputDataPositions)] ++ replicate (numberOfPeople - 2) (BaseSequencePosition inputDataSample2 inputDataPositions)

        inputDataPositions = STO.fromList $ take (B.length inputDataSample0) [(10::CInt)..]

        patterns = [Pattern 2000 pattern_CG, Pattern 900 pattern_cCGA] ++ replicate 50 (Pattern 10000 pattern_cccccccccc) ++ [Pattern 2000 pattern_CG]

        expected = [Match {mPatternId = 1,  mScore = 2500, mPosition = 13,     mSampleId = 0, mMatched = [a, c, g, a]},
                    Match {mPatternId = 0,  mScore = 2000, mPosition = 14,     mSampleId = 0, mMatched = [c, g]},
                    Match {mPatternId = 52, mScore = 2000, mPosition = 14,     mSampleId = 0, mMatched = [c, g]},
                    Match {mPatternId = 0,  mScore = 2000, mPosition = 10 + 7, mSampleId = 1, mMatched = [c, g]},
                    Match {mPatternId = 52, mScore = 2000, mPosition = 10 + 7, mSampleId = 1, mMatched = [c, g]}]


-- Check that the end of the chromosome is padded with N
test_patterns_padding :: IO ()
test_patterns_padding = do
    matches <- findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    assertEqual (V.fromList  expected) matches
  where numberOfPeople = 1000 :: Int

        inputData =  [BaseSequencePosition (B.take 5 inputDataSample0) (STO.take 5 inputDataPositions), BaseSequencePosition inputDataSample1 inputDataPositions] ++ replicate (numberOfPeople - 2) (BaseSequencePosition (B.take 10 inputDataSample2) (STO.take 10 inputDataPositions))
        inputDataPositions = STO.fromList $ take (B.length inputDataSample0) [(10::CInt)..]

        patterns = [Pattern 1000 pattern_CG, Pattern 1000 pattern_cCGA, Pattern 1000 pattern_CG]

        expected = [Match{mPatternId = 1, mScore = 1000, mPosition = 13, mSampleId = 0, mMatched = [a, c]},
                    Match{mPatternId = 0, mScore = 1000, mPosition = 14, mSampleId = 0, mMatched = [c]},
                    Match{mPatternId = 2, mScore = 1000, mPosition = 14, mSampleId = 0, mMatched = [c]},
                    Match{mPatternId = 0, mScore = 2000, mPosition = 10, mSampleId = 1, mMatched = [c, g]},
                    Match{mPatternId = 2, mScore = 2000, mPosition = 10, mSampleId = 1, mMatched = [c, g]}]


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
  assertEqual (Right $ Variant (Chromosome "1") (Position 69080) (Just "1:69081:G:C") (mkSeq [c]) (mkSeq [g]) (G.fromList [0,1]) (G.fromList [0,1]) sampleIdentifiers) (parseVariant sampleIdentifiers line)

test_filterOrderedIntervals_1 :: IO ()
test_filterOrderedIntervals_1 = do
  let inf = [10..]
  let filtered = filterOrderedIntervals Position [(Position 15, Position 18), (Position 20, Position 23)] inf
  assertEqual [15, 16, 17, 18, 20, 21, 22, 23] filtered

test_fillVector1 :: IO ()
test_fillVector1 = do
  let (vecL, vecR) = fillVector (encodeUtf8 "0/0\t0/1")
  let expected = (STO.fromList [], STO.fromList [1])
  assertEqual expected (vecL, vecR)

test_fillVector2 :: IO ()
test_fillVector2 = do
  let (vecL, vecR) = fillVector (encodeUtf8 "0/1\t0/0")
  let expected = (STO.fromList [], STO.fromList [0])
  assertEqual expected (vecL, vecR)

test_fillVector3 :: IO ()
test_fillVector3 = do
  let (vecL, vecR) = fillVector (encodeUtf8 "1/0\t0/0")
  let expected = (STO.fromList [0], STO.fromList [])
  assertEqual expected (vecL, vecR)

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
  let var p = Variant (Chromosome "1") (Position p) (Just $ "name" <> T.pack (show p)) (mkSeq [c]) (mkSeq [a]) (G.fromList []) (G.fromList [1])  (G.fromList [SampleId "sample1", SampleId "sample2"])
  let expected = Right [var 5, var 12, var 15, var 21]
  assertEqual expected parsed

test_parseBed_1 :: IO ()
test_parseBed_1 = do
  let bedContent = "chr1\t5\t6\nchr1\t8\t10\nchr2\t100\t110" :: T.Text
  let parsed = parseBedContent (Chromosome "1") bedContent
  let expected = Right [(Position 5, Position 6), (Position 8, Position 10)]
  assertEqual expected parsed


motifString :: T.Text
motifString = T.replace "   " "\t" (m1 <> m2 <> m3)
  where m1 = "ALX1_HUMAN.H11MO.0.B\n"
        m2 ="0.5145001398573071   -0.15943627836250168   0.08001130176246185   -0.9383465141290595\n"
        m3 = "1.9857651245551615   2.9786476868327387   0.0   26.03558718861208"

test_motif_parser_1 :: IO ()
test_motif_parser_1 = do
  let parsed = parseHocomocoMotifsContent motifString
  let expected = ("ALX1_HUMAN.H11MO.0.B", [Pweight{wa = 515, wc = -159, wg = 80, wt = -938}
                                          ,Pweight{wa = 1986, wc = 2979, wg = 0, wt = 26036}])
  assertEqual (Right [expected]) parsed


test_variantsToDiffs_1 :: IO ()
test_variantsToDiffs_1 = do
  let sampleIdentifiers = G.fromList [SampleId "sample1", SampleId "sample2"]
  let var = Variant (Chromosome "1") (Position 69080) (Just "1:69081:G:C") (mkSeq [c]) (mkSeq [g]) (G.fromList [0]) (G.fromList [1]) sampleIdentifiers
  let diffs = variantsToDiffs [var]
  let expected = Data.Map.fromList [((0,HaploLeft), V.fromList [Diff (Position 69080) (mkSeq [c]) (mkSeq [g])])
                                   ,((1,HaploRight),V.fromList [Diff (Position 69080) (mkSeq [c]) (mkSeq [g])])] :: Data.Map.Map (Int, Haplotype) (V.Vector Diff)
  assertEqual expected diffs

test_variantsToDiffs_2 :: IO ()
test_variantsToDiffs_2 = do
  let sampleIdentifiers = G.fromList [SampleId "sample1", SampleId "sample2"]
  let variant1 = Variant (Chromosome "1") (Position 69080) (Just "1:69081:G:C") (mkSeq [c]) (mkSeq [g]) (G.fromList [0]) (G.fromList [1]) sampleIdentifiers
  let variant2 = Variant (Chromosome "1") (Position 69079) (Just "1:69079:A:T") (mkSeq [a]) (mkSeq [t]) (G.fromList [0]) (G.fromList [])  sampleIdentifiers
  let variant3 = Variant (Chromosome "1") (Position 69078) (Just "1:69078:T:G") (mkSeq [t]) (mkSeq [g]) (G.fromList []) (G.fromList [0])  sampleIdentifiers
  let diffs = variantsToDiffs [variant1, variant2, variant3]
  let diff1 = Diff (Position 69080) (mkSeq [c]) (mkSeq [g])
  let diff2 = Diff (Position 69079) (mkSeq [a]) (mkSeq [t])
  let diff3 = Diff (Position 69080) (mkSeq [c]) (mkSeq [g])
  let diff4 = Diff (Position 69078) (mkSeq [t]) (mkSeq [g])
  let expected = Data.Map.fromList [((0,HaploLeft), V.fromList [diff2, diff1]) -- Order of positions must be increasing
                                   ,((0,HaploRight), V.fromList [diff4])
                                   ,((1,HaploRight),V.fromList [diff3])] :: Data.Map.Map (Int, Haplotype) (V.Vector Diff)
  assertEqual expected diffs

test_processPeak_1 :: IO()
test_processPeak_1 = do
  let chr = Chromosome "1"
  let patterns = mkPatterns [Pattern 2000 pattern_CG, Pattern 900 pattern_cCGA]
  let ref = takeRef (B.pack (map (fromIntegral . ord) "AAAACCCGGGTTT"))
  --                                                       T            -- Sample 0, haplotype Left
  --                                                              C     -- Sample 0, haplotype Right
  --                                                                    -- Sample 1, haplotype Left
  --                                                          A         -- Sample 1, haplotype Right



  let samples = Data.Map.fromList [(0, SampleId "sample0"), (1, SampleId "sample1")]
  let sampleIdentifiers = G.fromList [SampleId "sample0", SampleId "sample1"]
  let (peakStart, peakEnd) = (3, 7)
  let variant1 = Variant chr (Position 4)  (Just "1:4:C:T")  (mkSeq [c]) (mkSeq [t]) (G.fromList [0]) (G.fromList []) sampleIdentifiers
  let variant2 = Variant chr (Position 7)  (Just "1:7:G:A")  (mkSeq [g]) (mkSeq [a]) (G.fromList []) (G.fromList [1]) sampleIdentifiers
  let variant3 = Variant chr (Position 11) (Just "1:11:T:C") (mkSeq [t]) (mkSeq [c]) (G.fromList []) (G.fromList [0]) sampleIdentifiers
  let variants = [variant1, variant2, variant3]
  (nextVariants, numberOfHaplotypes, numberOfVariants, numberOfMatches) <- processPeak chr patterns ref samples (peakStart, peakEnd) variants
  assertEqual nextVariants [variant3]
  assertEqual numberOfHaplotypes 3
  assertEqual numberOfVariants 2
  assertEqual numberOfMatches 6


main :: IO ()
main = 
  htfMainWithArgs ["--quiet"] htf_thisModulesTests