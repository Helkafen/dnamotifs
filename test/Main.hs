{-# OPTIONS_GHC -F -pgmF htfpp #-}
{-# LANGUAGE OverloadedStrings #-}

import Test.Framework

import qualified Data.Vector.Storable            as STO
import qualified Data.Vector                     as V
import           Data.Monoid                     ((<>))
import           Foreign.C.Types                 (CInt)
import qualified Data.Vector.Generic             as G
import qualified Data.Text                       as T
import qualified Data.ByteString                 as B
import qualified Data.Set                        as Set
import           Data.Text.Encoding              (encodeUtf8)
import qualified Data.Map
import           Data.Char                       (ord)
import           Data.Maybe                      (fromJust)
import           Data.List                       (sort)

import Types
import Range
import PatternFind
import Run
import Fasta (takeRef)
import Vcf (parseVariant, filterOrderedIntervals, parseVcfContent, fillVector)
import Bed (parseBedContent)
import MotifDefinition (parseHocomocoMotifsContent)
import Haplotype (applyVariants, variantsToDiffs)
import Prelude

mkSeq :: [Nucleotide] -> B.ByteString
mkSeq = B.pack . map unNuc

sample0, sample1, sample2 :: SampleId
sample0 = SampleId "sample0" 0
sample1 = SampleId "sample1" 1
sample2 = SampleId "sample2" 2

inputDataSample0 :: B.ByteString
inputDataSample0 = mkSeq [a,a,a,a,c,g,a] <> B.replicate 93 (unNuc a)

inputDataSample1 :: B.ByteString
inputDataSample1 = mkSeq [c,g,a,a,a,a,a] <> B.replicate 93 (unNuc a)

inputDataSample2 :: B.ByteString
inputDataSample2 = mkSeq [a,a,a,a,a,a,a] <> B.replicate 93 (unNuc a)

pattern_CG, pattern_cCGA, pattern_cccccccccc, pattern_AT, pattern_TC :: [Pweight]
pattern_CG = [Pweight 0 1000 0 0, Pweight 0 0 1000 0]
pattern_AT = [Pweight 1000 0 0 0, Pweight 0 0 0 1000]
pattern_TC = [Pweight 0 0 0 1000, Pweight 0 1000 0 0]
pattern_cCGA = [Pweight 0 100 0 0, Pweight 0 1000 0 0, Pweight 0 0 1000 0, Pweight 500 0 0 0]
pattern_cccccccccc = [cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern]
  where cPattern = Pweight 0 0 1 0

test_patterns_basic :: IO ()
test_patterns_basic =
    assertEqual expected (findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns))
  where inputData = V.fromList [BaseSequencePosition inputDataSample0 inputDataPositions]
        inputDataPositions :: STO.Vector CInt
        inputDataPositions = STO.fromList $ take (B.length inputDataSample0) [(0::CInt)..]
        patterns = [Pattern 2000 pattern_CG]
        expected = V.fromList [Match {mPatternId = 0,  mScore = 2000, mPosition = 4, mSampleId = 0, mMatched = [c,g]}]

test_patterns_basic2 :: IO ()
test_patterns_basic2 =
    assertEqual expected (findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns))
  where inputData = V.fromList [BaseSequencePosition inputDataSample0 inputDataPositions]
        inputDataPositions = STO.fromList $ take (B.length inputDataSample0) [(0::CInt)..]

        patterns = [Pattern 2000 pattern_CG, Pattern 900 pattern_cCGA]

        expected = V.fromList [Match {mPatternId = 1, mScore = 2500,  mPosition = 3, mSampleId = 0, mMatched = [a, c, g, a]},
                               Match {mPatternId = 0, mScore = 2000, mPosition = 4, mSampleId = 0, mMatched = [c,g]}]


test_patterns_1 :: IO ()
test_patterns_1 =
    assertEqual (V.fromList  expected) (findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns))
  where numberOfPeople = 1000 :: Int

        inputData = V.fromList [BaseSequencePosition inputDataSample0 inputDataPositions, BaseSequencePosition inputDataSample1 inputDataPositions] <> V.replicate (numberOfPeople - 2) (BaseSequencePosition inputDataSample2 inputDataPositions)

        inputDataPositions = STO.fromList $ take (B.length inputDataSample0) [(0::CInt)..]

        patterns = [Pattern 2000 pattern_CG, Pattern 900 pattern_cCGA] ++ replicate 50 (Pattern 10000 pattern_cccccccccc) ++ [Pattern 2000 pattern_CG]

        expected = [Match {mPatternId = 1,  mScore = 2500, mPosition = 3, mSampleId = 0, mMatched = [a, c, g, a]},
                    Match {mPatternId = 0,  mScore = 2000, mPosition = 4, mSampleId = 0, mMatched = [c, g]},
                    Match {mPatternId = 52, mScore = 2000, mPosition = 4, mSampleId = 0, mMatched = [c, g]},
                    Match {mPatternId = 0,  mScore = 2000, mPosition = 0, mSampleId = 1, mMatched = [c, g]},
                    Match {mPatternId = 52, mScore = 2000, mPosition = 0, mSampleId = 1, mMatched = [c, g]}]


-- Check that the Matches contain the reference position, not the real position
test_patterns_2 :: IO ()
test_patterns_2 =
    assertEqual (V.fromList  expected) (findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns))
  where numberOfPeople = 1000 :: Int

        inputData = V.fromList [BaseSequencePosition inputDataSample0 inputDataPositions, BaseSequencePosition inputDataSample1 (STO.map (+7) inputDataPositions)] <> V.replicate (numberOfPeople - 2) (BaseSequencePosition inputDataSample2 inputDataPositions)

        inputDataPositions = STO.fromList $ take (B.length inputDataSample0) [(10::CInt)..]

        patterns = [Pattern 2000 pattern_CG, Pattern 900 pattern_cCGA] ++ replicate 50 (Pattern 10000 pattern_cccccccccc) ++ [Pattern 2000 pattern_CG]

        expected = [Match {mPatternId = 1,  mScore = 2500, mPosition = 13,     mSampleId = 0, mMatched = [a, c, g, a]},
                    Match {mPatternId = 0,  mScore = 2000, mPosition = 14,     mSampleId = 0, mMatched = [c, g]},
                    Match {mPatternId = 52, mScore = 2000, mPosition = 14,     mSampleId = 0, mMatched = [c, g]},
                    Match {mPatternId = 0,  mScore = 2000, mPosition = 10 + 7, mSampleId = 1, mMatched = [c, g]},
                    Match {mPatternId = 52, mScore = 2000, mPosition = 10 + 7, mSampleId = 1, mMatched = [c, g]}]


-- Check that the end of the chromosome is padded with N
test_patterns_padding :: IO ()
test_patterns_padding =
    assertEqual (V.fromList  expected) (findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns))
  where numberOfPeople = 1000 :: Int

        inputData = V.fromList [BaseSequencePosition (B.take 5 inputDataSample0) (STO.take 5 inputDataPositions), BaseSequencePosition inputDataSample1 inputDataPositions] <> V.replicate (numberOfPeople - 2) (BaseSequencePosition (B.take 10 inputDataSample2) (STO.take 10 inputDataPositions))
        inputDataPositions = STO.fromList $ take (B.length inputDataSample0) [(10::CInt)..]

        patterns = [Pattern 1000 pattern_CG, Pattern 1000 pattern_cCGA, Pattern 1000 pattern_CG]

        expected = [Match{mPatternId = 1, mScore = 1000, mPosition = 13, mSampleId = 0, mMatched = [a, c]},
                    Match{mPatternId = 0, mScore = 1000, mPosition = 14, mSampleId = 0, mMatched = [c]},
                    Match{mPatternId = 2, mScore = 1000, mPosition = 14, mSampleId = 0, mMatched = [c]},
                    Match{mPatternId = 0, mScore = 2000, mPosition = 10, mSampleId = 1, mMatched = [c, g]},
                    Match{mPatternId = 2, mScore = 2000, mPosition = 10, mSampleId = 1, mMatched = [c, g]}]


refGenome :: Range Position0 -> BaseSequencePosition
refGenome = takeRef (B.pack [65,67,71,84]) -- ACGT


test_apply_variant_1 :: IO ()
test_apply_variant_1 = do
  let patched = applyVariants refGenome (Range (Position 1) (Position 2)) []
  assertEqual (BaseSequencePosition (mkSeq [c,g]) (STO.fromList [1,2])) patched

test_apply_variant_2 :: IO ()
test_apply_variant_2 = do
  let patched = applyVariants refGenome (Range (Position 0) (Position 2)) []
  assertEqual (BaseSequencePosition (mkSeq [a, c, g]) (STO.fromList [0,1,2])) patched

test_apply_variant_3 :: IO ()
test_apply_variant_3 = do
  let patched = applyVariants refGenome (Range (Position 0) (Position 5)) []
  assertEqual (BaseSequencePosition (mkSeq [a,c,g,t]) (STO.fromList [0,1,2,3])) patched

test_apply_variant_4 :: IO ()
test_apply_variant_4 = do
  let patched = applyVariants refGenome (Range (Position 1) (Position 2)) [Diff 100 (mkSeq [a]) (mkSeq [c])]
  assertEqual (BaseSequencePosition (mkSeq [c,g]) (STO.fromList [1,2])) patched

test_apply_variant_5 :: IO ()
test_apply_variant_5 = do
  let patched = applyVariants refGenome (Range (Position 1) (Position 2)) [Diff 1 (mkSeq [c]) (mkSeq [n])]
  assertEqual (BaseSequencePosition (mkSeq [n,g]) (STO.fromList [1,2])) patched

test_apply_variant_6 :: IO ()
test_apply_variant_6 = do
  let patched = applyVariants refGenome (Range (Position 1) (Position 2)) [Diff 2 (mkSeq [g]) (mkSeq [a])]
  assertEqual (BaseSequencePosition (mkSeq [c,a]) (STO.fromList [1,2])) patched

test_apply_variant_7 :: IO ()
test_apply_variant_7 = do
  let patched = applyVariants refGenome (Range (Position 1) (Position 2)) [Diff 1 (mkSeq [c]) (mkSeq [n]), Diff 2 (mkSeq [g]) (mkSeq [a])]
  assertEqual (BaseSequencePosition (mkSeq [n,a]) (STO.fromList [1,2])) patched

test_apply_variant_8 :: IO ()
test_apply_variant_8 = do
  let patched = applyVariants refGenome (Range (Position 1) (Position 2)) [Diff 0 (mkSeq [c]) (mkSeq [n]), Diff 4 (mkSeq [g]) (mkSeq [a])]
  assertEqual (BaseSequencePosition (mkSeq [c,g]) (STO.fromList [1,2])) patched

test_apply_variant_9 :: IO ()
test_apply_variant_9 = do
  let patched = applyVariants refGenome (Range (Position 1) (Position 2)) [Diff 1 (mkSeq [c]) (mkSeq [n,n])]
  assertEqual (BaseSequencePosition (mkSeq[n,n,g]) (STO.fromList [1,1,2])) patched

test_apply_variant_10 :: IO ()
test_apply_variant_10 = do
  let patched = applyVariants refGenome (Range (Position 1) (Position 2)) [Diff 2 (mkSeq [c]) (mkSeq [n,n])]
  assertEqual (BaseSequencePosition (mkSeq[c,n,n]) (STO.fromList[1,2,2])) patched

test_apply_variant_11 :: IO ()
test_apply_variant_11 = do
  let patched = applyVariants refGenome (Range (Position 1) (Position 2)) [Diff 3 (mkSeq [c]) (mkSeq [n,n])]
  assertEqual (BaseSequencePosition (mkSeq[c,g]) (STO.fromList [1,2])) patched

test_apply_variant_12 :: IO ()
test_apply_variant_12 = do
  let patched = applyVariants refGenome (Range (Position 1) (Position 2)) [Diff 1 (mkSeq [c, g]) (mkSeq [n])]
  assertEqual (BaseSequencePosition (mkSeq[n]) (STO.fromList [1])) patched

test_apply_variant_13 :: IO ()
test_apply_variant_13 = do
  let patched = applyVariants refGenome (Range (Position 1) (Position 2)) [Diff 2 (mkSeq [g, t]) (mkSeq [n])]
  assertEqual (BaseSequencePosition (mkSeq[c,n]) (STO.fromList [1,2])) patched

test_apply_variant_14 :: IO ()
test_apply_variant_14 = do
  let patched = applyVariants refGenome (Range (Position 1) (Position 2)) [Diff 0 (mkSeq [a, c, g]) (mkSeq [n, n])]
  assertEqual (BaseSequencePosition (mkSeq[c,g]) (STO.fromList [1,2])) patched -- For simplicity, do not apply a Diff that starts before the window we're observing

test_parse_1 :: IO ()
test_parse_1 = do
  let line = "chr1\t69081\t1:69081:G:C\tC\tG\t.\t.\t.\tGT\t1/1\t1/1"
  let sampleIdentifiers = G.fromList [sample1, sample2]
  assertEqual (Right $ Variant (Chromosome "1") (Position 69080) (Just "1:69081:G:C") (mkSeq [c]) (mkSeq [g]) (G.fromList [0,1]) (G.fromList [0,1]) sampleIdentifiers) (parseVariant sampleIdentifiers line)

test_filterOrderedIntervals_1 :: IO ()
test_filterOrderedIntervals_1 = do
  let inf = [10..]
  let range1 = Range (Position 15) (Position 18)
  let range2 = Range (Position 20) (Position 23)
  let ranges = fromJust $ mkRanges [range1, range2]
  --assertEqual (isRight ranges) True
  let filtered = filterOrderedIntervals Position ranges inf
  assertEqual [(range1, [15, 16, 17, 18]), (range2, [20, 21, 22, 23])] filtered

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
  let vcf = map encodeUtf8 ["#\t\t\t\t\t\t\t\t\tsample0\tsample1",
             "chr1\t5\tname4\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t6\tname5\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t13\tname12\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t16\tname15\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t18\tname17\tC\tA\t\t\t\t\t0/0\t0/1",
             "chr1\t22\tname21\tC\tA\t\t\t\t\t0/0\t0/1"
            ] :: [B.ByteString]
  let range1 = Range (Position 5) (Position 15)
  let range2 = Range (Position 18) (Position 25)
  let ranges = fromJust $ mkRanges [range1, range2]
  let parsed = parseVcfContent ranges vcf
  let var p = Variant (Chromosome "1") (Position p) (Just $ "name" <> T.pack (show p)) (mkSeq [c]) (mkSeq [a]) (G.fromList []) (G.fromList [1])  (G.fromList [sample0, sample1])
  let expected = Right [(range1, [var 5, var 12, var 15]), (range2, [var 21])]
  assertEqual expected parsed

test_parseBed_1 :: IO ()
test_parseBed_1 = do
  let bedContent = "chr1\t5\t6\nchr1\t8\t10\nchr2\t100\t110" :: T.Text
  let parsed = parseBedContent (Chromosome "1") bedContent
  let expected = fromJust $ mkRanges [Range (Position 5) (Position 6), Range (Position 8) (Position 10)]
  assertEqual (Right expected) parsed


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
  let sampleIdentifiers = G.fromList [sample1, sample2]
  let var = Variant (Chromosome "1") (Position 69080) (Just "1:69081:G:C") (mkSeq [c]) (mkSeq [g]) (G.fromList [0]) (G.fromList [1]) sampleIdentifiers
  let diffs = variantsToDiffs [var]
  let expected = Data.Map.fromList [(V.singleton (Diff (Position 69080) (mkSeq [c]) (mkSeq [g])),  [HaplotypeId sample2 HaploRight, HaplotypeId sample1 HaploLeft])]
  assertEqual expected diffs

test_variantsToDiffs_2 :: IO ()
test_variantsToDiffs_2 = do
  let sampleIdentifiers = G.fromList [sample1, sample2]
  let variant1 = Variant (Chromosome "1") (Position 69080) (Just "1:69081:G:C") (mkSeq [c]) (mkSeq [g]) (G.fromList [0]) (G.fromList [1]) sampleIdentifiers
  let variant2 = Variant (Chromosome "1") (Position 69079) (Just "1:69079:A:T") (mkSeq [a]) (mkSeq [t]) (G.fromList [0]) (G.fromList [])  sampleIdentifiers
  let variant3 = Variant (Chromosome "1") (Position 69078) (Just "1:69078:T:G") (mkSeq [t]) (mkSeq [g]) (G.fromList []) (G.fromList [0])  sampleIdentifiers
  let diffs = variantsToDiffs [variant1, variant2, variant3]
  let diff1 = Diff (Position 69080) (mkSeq [c]) (mkSeq [g])
  let diff2 = Diff (Position 69079) (mkSeq [a]) (mkSeq [t])
  let diff4 = Diff (Position 69078) (mkSeq [t]) (mkSeq [g])
  let expected = Data.Map.fromList [(V.fromList [diff4],       [HaplotypeId sample1 HaploRight])
                                   ,(V.fromList [diff2,diff1], [HaplotypeId sample1 HaploLeft])
                                   ,(V.fromList [diff1],       [HaplotypeId sample2 HaploRight])]
  assertEqual expected diffs

test_processPeak_1 :: IO()
test_processPeak_1 = do
  let chr = Chromosome "1"
  let pattern_CCCG = Pattern 4000 (replicate 3 (Pweight 0 1000 0 0) ++ [Pweight 0 0 1000 0])
  let patterns = mkPatterns [Pattern 2000 pattern_CG, Pattern 2000 pattern_AT, Pattern 2000 pattern_TC, pattern_CCCG]
  let ref = takeRef (B.pack (map (fromIntegral . ord) "AAAACCCGGGTTT"))
  --                                                       T            -- Sample 0, haplotype Left
  --                                                              C     -- Sample 0, haplotype Right
  --                                                                    -- Sample 1, haplotype Left
  --                                                          A         -- Sample 1, haplotype Right



  let samples = Data.Map.fromList [(0, sample0), (1, sample1)]
  let sampleIdentifiers = G.fromList [sample0, sample1]
  let range1 = Range (Position 3) (Position 7)
  let variant1 = Variant chr (Position 4)  (Just "1:4:C:T")  (mkSeq [c]) (mkSeq [t]) (G.fromList [0]) (G.fromList []) sampleIdentifiers
  let variant2 = Variant chr (Position 7)  (Just "1:7:G:A")  (mkSeq [g]) (mkSeq [a]) (G.fromList []) (G.fromList [1]) sampleIdentifiers
  --let variant3 = Variant chr (Position 11) (Just "1:11:T:C") (mkSeq [t]) (mkSeq [c]) (G.fromList []) (G.fromList [0]) sampleIdentifiers
  let variants = [variant1, variant2] -- No variant3. It must be filtered out before processPeak
  let (matches, numberOfHaplotypes, numberOfVariants, numberOfMatches) = processPeak patterns ref samples (range1, variants)
  assertEqual 3 numberOfHaplotypes -- Including the reference genome. Only on interval [3->7]
  assertEqual 2 numberOfVariants
  assertEqual 7 numberOfMatches
  let expectedMatches = [Match 3 4000 4 [HaplotypeId sample0 HaploRight, HaplotypeId sample1 HaploLeft] [c,c,c,g], -- This line matches the reference genome
                         Match 0 2000 6 [HaplotypeId sample0 HaploRight, HaplotypeId sample1 HaploLeft] [c,g],
                         Match 1 2000 3 [HaplotypeId sample0 HaploLeft]                                              [a,t],
                         Match 2 2000 4 [HaplotypeId sample0 HaploLeft]                                              [t,c],
                         Match 0 2000 6 [HaplotypeId sample0 HaploLeft]                                              [c,g]]
  assertEqual (V.fromList expectedMatches) matches

--test_countsInPeak_1 :: IO ()
--test_countsInPeak_1 = do
--  let samples = M.fromList [(0, sample0), (1, sample1)]
--  let matches = V.fromList [Match 0 2000 6 [HaplotypeId (sample0) HaploLeft] [c,g]]
--  let peak = Range (Position 0) (Position 10)
--  let counted = countsInPeak samples matches peak
--  let expected = [(peak, 0, [Count2 1 0, Count2 0 0])]
--  assertEqual expected counted
--
--
--test_countsInPeak_2 :: IO ()
--test_countsInPeak_2 = do
--  let samples = M.fromList [(0, sample0), (1, sample1)]
--  let matches = V.fromList [Match 0 2000 6 [HaplotypeId (sample0) HaploLeft] [c,g], Match 1 2000 6 [HaplotypeId (sample1) HaploRight] [c,g]]
--  let peak = Range (Position 0) (Position 10)
--  let counted = countsInPeak samples matches peak
--  let expected = [(peak, 0, [Count2 1 0, Count2 0 0]), (peak, 1, [Count2 0 0, Count2 0 1])]
--  assertEqual expected counted
--
--
---- One line per pattern found in the peak
--prop_countsInPeak_3 :: [Match [HaplotypeId]] -> Bool
--prop_countsInPeak_3 matches =
--  let fixedMatches = map (\m -> m { mSampleId = sort (Set.toList (Set.fromList (mSampleId m))) }) matches
--      samples = M.fromList [(0, sample0), (1, sample1)]
--      peak = Range (Position 10) (Position 90)
--      matchesInPeak = filter (inRange peak . Position . mPosition) fixedMatches
--      patternsFoundInPeak = Set.fromList (map mPatternId matchesInPeak)
--  in Set.size patternsFoundInPeak == length (countsInPeak samples (V.fromList fixedMatches) peak)
--
---- No Match is lost, and the peak coordinates are provided
--prop_countsInPeak_4 :: [Match [HaplotypeId]] -> Bool
--prop_countsInPeak_4 matches =
--  let fixedMatches = map (\m -> m { mSampleId = sort (Set.toList (Set.fromList (mSampleId m))) }) matches
--      samples = M.fromList [(0, sample0), (1, sample1)]
--      peak = Range (Position 10) (Position 90)
--      matchesInPeak = filter (inRange peak . Position . mPosition) fixedMatches
--      counts = countsInPeak samples (V.fromList fixedMatches) peak
--      Count2 l r = mconcat $ map (\(_, _, xs) -> mconcat xs) counts
--      outputPeaks = Set.fromList $ map (\(p, _, _) -> p) counts
--  in r + l == sum (map (length . mSampleId) matchesInPeak) && if not (null counts) then outputPeaks == Set.singleton peak else outputPeaks == Set.empty


test_encodeNumberOfMatchesAsGenotypes_1 :: IO ()
test_encodeNumberOfMatchesAsGenotypes_1 = assertEqual [geno00] (fst (encodeNumberOfMatchesAsGenotypes [3]))

test_encodeNumberOfMatchesAsGenotypes_2 :: IO ()
test_encodeNumberOfMatchesAsGenotypes_2 = assertEqual [] (fst (encodeNumberOfMatchesAsGenotypes []))

test_encodeNumberOfMatchesAsGenotypes_3 :: IO ()
test_encodeNumberOfMatchesAsGenotypes_3 = assertEqual [geno00, geno00, geno00] (fst (encodeNumberOfMatchesAsGenotypes [3, 3, 3]))

test_encodeNumberOfMatchesAsGenotypes_4 :: IO ()
test_encodeNumberOfMatchesAsGenotypes_4 = assertEqual [geno00, geno11, geno00] (fst (encodeNumberOfMatchesAsGenotypes [3, 4, 3]))

test_encodeNumberOfMatchesAsGenotypes_5 :: IO ()
test_encodeNumberOfMatchesAsGenotypes_5 = assertEqual [geno00, geno01, geno00, geno11, geno00, geno01] (fst (encodeNumberOfMatchesAsGenotypes [3, 4, 3, 5, 3, 4]))

test_encodeNumberOfMatchesAsGenotypes_6 :: IO ()
test_encodeNumberOfMatchesAsGenotypes_6 = assertEqual [geno00, geno00, geno01, geno01, geno11, geno11] (fst (encodeNumberOfMatchesAsGenotypes [2, 3, 6, 7, 9, 10]))

test_encodeNumberOfMatchesAsGenotypes_7 :: IO ()
test_encodeNumberOfMatchesAsGenotypes_7 = assertEqual [geno00, geno11] (fst (encodeNumberOfMatchesAsGenotypes [3, 4]))

test_encodeNumberOfMatchesAsGenotypes_8 :: IO ()
test_encodeNumberOfMatchesAsGenotypes_8 = assertEqual [geno00, geno00, geno01, geno11, geno11] (fst (encodeNumberOfMatchesAsGenotypes [0,0,3,5,6]))

prop_encodeNumberOfMatchesAsGenotypes_1 :: [Int] -> Bool
prop_encodeNumberOfMatchesAsGenotypes_1 xs = fst (encodeNumberOfMatchesAsGenotypes (sort positiveXs)) == sort (fst $ encodeNumberOfMatchesAsGenotypes positiveXs)
    where positive [] = []
          positive vs = map (+ abs (minimum vs)) vs
          positiveXs = positive xs

prop_encodeNumberOfMatchesAsGenotypes_2 :: [Int] -> Bool
prop_encodeNumberOfMatchesAsGenotypes_2 xs = snd (encodeNumberOfMatchesAsGenotypes positiveXs) == Set.fromList positiveXs
    where positive [] = []
          positive vs = map (+ abs (minimum vs)) vs
          positiveXs = positive xs

main :: IO ()
main =
  htfMainWithArgs ["--quiet"] htf_thisModulesTests