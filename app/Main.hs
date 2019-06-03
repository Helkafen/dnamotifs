{-# LANGUAGE OverloadedStrings #-}
module Main where

import qualified Data.Vector.Storable            as STO
import qualified Data.Vector                     as Vector
import qualified Data.ByteString                 as B
import qualified Data.Text                       as T
import           Data.Monoid                     ((<>))
import           Foreign.C.Types                 (CInt)
import           System.Environment              (getArgs)
import Data.List (elem)
import Types
import PatternFind
import Lib (findPatterns)
import MotifDefinition (loadMotifs)

inputDataSample0 :: B.ByteString
inputDataSample0 = B.pack (map unNuc [a,a,a,a,c,g,a]) <> B.replicate 93 (unNuc a)

inputDataSample1 :: B.ByteString
inputDataSample1 = B.pack (map unNuc [c,g,a,a,a,a,a]) <> B.replicate 93 (unNuc a)

inputDataSample2 :: B.ByteString
inputDataSample2 = B.pack (map unNuc [a,a,a,a,a,a,a]) <> B.replicate 93 (unNuc a)

pattern_CG, pattern_cCGA, pattern_CGA, pattern_cccccccccc :: [Pweight]
pattern_CG = [Pweight 0 1 0 0, Pweight 0 0 1 0]
pattern_cCGA = [Pweight 0 0.1 0 0, Pweight 0 1 0 0, Pweight 0 0 1 0, Pweight 0.5 0 0 0]
pattern_CGA = [Pweight 0 1 0 0, Pweight 0 0 1 0]
pattern_cccccccccc = [cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern]
  where cPattern = Pweight 0 0 0.001 0

patterns :: Patterns
patterns = mkPatterns $ [pattern_CG, pattern_cCGA, pattern_CGA] ++ replicate 50 pattern_cccccccccc ++ [pattern_CG]

bench :: IO ()
bench = do
    matches <- findPatternsInBlock 500 (mkNucleotideAndPositionBlock inputData) patterns
    print ((Vector.fromList  expected) == matches)
  where numberOfPeople = 10000 :: Int

        inputData :: [BaseSequencePosition]
        inputData =  [BaseSequencePosition inputDataSample0 inputDataPositions, BaseSequencePosition inputDataSample1 inputDataPositions] ++ (replicate (numberOfPeople - 2) (BaseSequencePosition inputDataSample2 inputDataPositions))

        inputDataPositions :: STO.Vector CInt
        inputDataPositions = STO.fromList $ take (B.length inputDataSample0) [0..]

        expected = [Match {mPatternId = 53, mScore = 1000, mPosition = 0, mSampleId = 1, mMatched = [c, g]}
                    ,Match {mPatternId = 2,  mScore = 1000, mPosition = 0, mSampleId = 1, mMatched = [c, g]}
                    ,Match {mPatternId = 0,  mScore = 1000, mPosition = 0, mSampleId = 1, mMatched = [c, g]}
                    ,Match {mPatternId = 53, mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c, g]}
                    ,Match {mPatternId = 2,  mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c, g]}
                    ,Match {mPatternId = 0,  mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c, g]}
                    ,Match {mPatternId = 1,  mScore = 961,  mPosition = 3, mSampleId = 0, mMatched = [a, c, g, a]}]


-- http://jaspar.genereg.net/api/v1/matrix/XXXXX.meme  http://jaspar.genereg.net/matrix/XXXXX/ http://jaspar.genereg.net/download/bed_files/XXXXX.bed
-- http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/pwm/XXX.pwm

data TranscriptionFactor = TranscriptionFactor {
  tfName :: String,
  tfJasparId :: Maybe String,
  tfHomerName :: Maybe String,
  tfHocomocoPwd :: Maybe String
}

-- DNA binding sites: representation and discovery: https://www.ncbi.nlm.nih.gov/pubmed/10812473
-- HOCOMOCO: https://repository.kaust.edu.sa/bitstream/handle/10754/613302/Nucl. Acids Res.-2016-Kulakovskiy-D116-25.pdf
-- Optimally choosing PWM motif databases and sequence scanning approaches based on ChIP-seq data: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0573-5
-- Algorithm of Match: https://academic.oup.com/nar/article/31/13/3576/2904207
-- HOCOMOCO: Explanation of PCM vs PWM files: http://www.cbrc.kaust.edu.sa/hocomoco/Details.php#400, https://academic.oup.com/nar/article/31/20/6016/1039515. And PWM cutoff: https://genome.cshlp.org/content/12/3/470.full
-- Not really useful to use Information content to account for GC vs AT content, since humans have 46% AT.

-- PCM files can be downloaded at http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/pcm/ + PCM filename
knownPatterns :: [TranscriptionFactor]
knownPatterns = [
  TranscriptionFactor "JUNB"    (Just "MA0490.1") Nothing (Just "JUNB_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "SMARCC1" Nothing Nothing Nothing,  -- Look for SMRC1
  TranscriptionFactor "FOS"     (Just "MA0476.1") Nothing (Just "JFOS_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "FOSB"    Nothing Nothing (Just "JFOSB_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "FOSL1"   (Just "MA0477.1") Nothing (Just "FOSL1_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "FOSL2"   (Just "MA0478.1") (Just "Fosl2") (Just "FOSL2_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "BCL11A"  Nothing Nothing Nothing,
  TranscriptionFactor "BCL11B"  Nothing Nothing Nothing,
  TranscriptionFactor "JDP2"    (Just "MA0655.1") Nothing (Just "JDP2_HUMAN.H11MO.0.D.pcm"),
  TranscriptionFactor "GATA1"   (Just "MA0140.2") (Just "Gata1") (Just "GATA1_HUMAN.H11MO.1.A.pcm"),
  TranscriptionFactor "GATA2"   (Just "MA0036.3") (Just "Gata2") (Just "GATA2_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "GATA3"   (Just "MA0037.3") (Just "Gata3") (Just "GATA3_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "GATA4"   Nothing (Just "Gata4") (Just "GATA4_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "GATA5"   (Just "MA0766.1") Nothing (Just "GATA5_HUMAN.H11MO.0.D.pcm"),
  TranscriptionFactor "GATA6"   (Just "MA1104.1") (Just "Gata6") (Just "GATA6_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "JUN"     Nothing Nothing (Just "JUN_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "JUND"    (Just "MA0491.1") (Just "JunD") (Just "JUND_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "BATF"    Nothing (Just "BATF") (Just "BATF_HUMAN.H11MO.1.A.pcm"),
  TranscriptionFactor "ATF3"    Nothing (Just "Atf3") (Just "ATF3_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "BACH1"   Nothing (Just "Bach1") (Just "BACH1_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "BACH2"   (Just "MA1101.1") (Just "Bach2") (Just "BACH2_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "NF-E2"   (Just "MA0841.1") (Just "NF-E2") (Just "NFE2_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "CEBPA"   (Just "MA0102.3") Nothing (Just "CEBPA_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "CEBPB"   (Just "MA0466.2") Nothing (Just "CEBPB_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "CEBPD"   (Just "MA0836.1") Nothing (Just "CEBPD_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "CEBPE"   (Just "MA0837.1") Nothing (Just "CEBPE_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "CEBPG"   (Just "MA0838.1") Nothing (Just "CEBPG_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "SPIB"    (Just "MA0081.1") (Just "SpiB") (Just "SPIB_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "IRF8"    (Just "MA0652.1") Nothing (Just "IRF8_HUMAN.H11MO.0.B.pcm"),
  TranscriptionFactor "SPI1"    (Just "MA0080.4") Nothing (Just "SPI1_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "SPI1D"   Nothing Nothing Nothing,
  TranscriptionFactor "LMO2"    Nothing Nothing Nothing,
  TranscriptionFactor "MESP1"   Nothing Nothing (Just "MESP1_HUMAN.H11MO.0.D.pcm"),
  TranscriptionFactor "MESP2"   Nothing Nothing Nothing,
  TranscriptionFactor "ID3"     Nothing Nothing Nothing,
  TranscriptionFactor "ID4"     (Just "MA0824.1") Nothing (Just "ID4_HUMAN.H11MO.0.D.pcm"),
  TranscriptionFactor "TCF12"   Nothing  (Just "Tcf12") (Just "HTF4_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "TCF4"    (Just "MA0830.1") (Just "Tcf4") (Just "ITF2_HUMAN.H11MO.0.C.pcm"),
  TranscriptionFactor "STAT1"   (Just "MA0137.3") (Just "STAT1") (Just "STAT1_HUMAN.H11MO.1.A.pcm"),
  TranscriptionFactor "STAT2"   Nothing Nothing (Just "STAT2_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "SPIC"    (Just "MA0687.1") Nothing (Just "SPIC_HUMAN.H11MO.0.D.pcm"),
  TranscriptionFactor "CTCF"    (Just "MA0139.1") (Just "CTCF") (Just "CTCF_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "IRF1"    (Just "MA0050.2") (Just "IRF1") (Just "IRF1_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "DBP"     (Just "MA0639.1") Nothing (Just "DBP_HUMAN.H11MO.0.B.pcm"),
  TranscriptionFactor "MAFK"    (Just "MA0496.2") (Just "MafK") (Just "MAFK_HUMAN.H11MO.1.A.pcm"),
  TranscriptionFactor "ATF4"    (Just "MA0833.1") (Just "Atf4") (Just "ATF4_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "ASCL1"   (Just "MA1100.1") (Just "Ascl1") (Just "ASCL1_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "ASCL2"   Nothing Nothing (Just "ASCL2_HUMAN.H11MO.0.D.pcm"),
  TranscriptionFactor "TCF3"    (Just "MA0522.2") (Just "Tcf3") (Just "TFE2_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "MYOD1"   Nothing Nothing (Just "MYOD1_HUMAN.H11MO.1.A.pcm"),
  TranscriptionFactor "ATOH8"   Nothing Nothing Nothing,
  TranscriptionFactor "MECOM"   Nothing Nothing (Just "EVI1_HUMAN.H11MO.0.B.pcm"),
  TranscriptionFactor "IRF3"    (Just "MA1418.1") Nothing (Just "IRF3_HUMAN.H11MO.0.B.pcm"),
  TranscriptionFactor "ZEB1"    (Just "MA0103.3") Nothing (Just "ZEB1_HUMAN.H11MO.0.A.pcm"),
  TranscriptionFactor "IRF9"    (Just "MA0653.1") Nothing (Just "IRF9_HUMAN.H11MO.0.C.pcm"),
  TranscriptionFactor "NHLH1"   (Just "MA0048.2") Nothing (Just "HEN1_HUMAN.H11MO.0.C.pcm"),
  TranscriptionFactor "LYL1"    Nothing Nothing (Just "LYL1_HUMAN.H11MO.0.A.pcm")]

main :: IO()
main = do
    args <- getArgs
    let patternNames = ["Fosl2", "JunD", "Gata1", "Jun-AP1", "Gata2", "Fosl2", "BATF", "Gata3", "Bach1", "Atf3", "Bach2", "NF-E2", "Gata4", "CEBP", "SpiB", "PU.1:IRF8",
                        "Gata6", "Tcf12", "Tcf4", "STAT1", "CTCF", "IRF1", "MafK", "Atf4", "Ascl1", "Tcf3"] :: [T.Text]
    
    case args of
      [] -> do
        --bench
        -- From http://schemer.buenrostrolab.com :
        --JUNB SMARCC1 FOSL2 FOSL1 JUND GATA1 JUN                  GATA2 FOS    BATF GATA3 BACH1 ATF3 BACH2 FOSB BCL11A BCL11B JDP2 GATA5 NFE2  SPI1D GATA4 CEBPB CEBPA SPIB IRF8      SPI1 CEBPD
        --x            Fosl2       JunD Gata1 Jun-AP1 or c-Jun-CRE Gata2 Fosl2? BATF Gata3 Bach1 Atf3 Bach2                               NF-E2       Gata4 CEBP        SpiB PU.1:IRF8           
        --
        --LMO2 GATA6 CEBPG MESP1 MESP2 ID3 ID4 TCF12 TCF4 STAT1 CEBPE SPIC CTCF IRF1 STAT2 DBP MAFK ATF4 ASCL1 TCF3 MYOD1 ATOH8 MECOM ASCL2 IRF3 ZEB1 IRF9 NHLH1 LYL1
        --x    Gata6                           Tcf12 Tcf4 STAT1            CTCF IRF1           MafK Atf4 Ascl1 Tcf3 
        _ <- findPatterns (Chromosome "1") patterns 900 "chr1.bed" "hg38.fa" "chr1.vcf.gz" "resultFile.tab"
        pure ()
      [chrom, referenceGenomeFastaFile, motifFile, peakBedFile, vcfFile, outputFile] -> do
        wantedPatterns <- (mkPatterns . map snd . filter ((`elem` patternNames) . fst)) <$> loadMotifs motifFile
        _ <- findPatterns (Chromosome $ T.pack chrom) wantedPatterns 900 peakBedFile referenceGenomeFastaFile vcfFile outputFile
        pure ()
      _ -> print ("Usage: xxx" :: String)
      