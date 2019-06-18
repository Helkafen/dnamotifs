{-# LANGUAGE OverloadedStrings #-}
module Main where

import qualified Data.Text                       as T
import qualified Data.Map                        as M
import           System.Environment              (getArgs)
import           Data.List                       (elem)
import           Data.Maybe                      (mapMaybe)
import           Control.Monad                   (mapM_)
import           Types
import           PatternFind
import           Lib                             (findPatterns)
import           MotifDefinition                 (loadHocomocoMotifs)
import           Data.List.Split                 (splitOn)
import           Control.Monad.Except            (runExceptT)


-- http://jaspar.genereg.net/api/v1/matrix/XXXXX.meme  http://jaspar.genereg.net/matrix/XXXXX/ http://jaspar.genereg.net/download/bed_files/XXXXX.bed
-- http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/pwm/XXX.pwm

data TranscriptionFactor = TranscriptionFactor {
  tfName :: T.Text,
  tfJasparId :: Maybe T.Text,
  tfHomerId :: Maybe T.Text,
  tfHocomocoId :: Maybe T.Text
}

-- Jaccard distance: https://arxiv.org/pdf/1811.00416.pdf
-- Explanation of PWMs: https://davetang.org/muse/2013/10/01/position-weight-matrix/. Paper: https://sci-hub.tw/10.1038/nrg1315 ("Applied bioinformatics for the identification of regulatory elements")
-- DNA binding sites: representation and discovery: https://www.ncbi.nlm.nih.gov/pubmed/10812473
-- HOCOMOCO: https://repository.kaust.edu.sa/bitstream/handle/10754/613302/Nucl. Acids Res.-2016-Kulakovskiy-D116-25.pdf
-- Optimally choosing PWM motif databases and sequence scanning approaches based on ChIP-seq data: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0573-5
-- Algorithm of Match: https://academic.oup.com/nar/article/31/13/3576/2904207
-- HOCOMOCO: Explanation of PCM vs PWM files: http://www.cbrc.kaust.edu.sa/hocomoco/Details.php#400, https://academic.oup.com/nar/article/31/20/6016/1039515. And PWM cutoff: https://genome.cshlp.org/content/12/3/470.full
-- Not really useful to use Information content to account for GC vs AT content, since humans have 46% AT.
-- Coverage of motifs: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0573-5

-- PCM files can be downloaded at http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/pcm/ + PCM filename
knownPatterns :: [TranscriptionFactor]
knownPatterns = [
  TranscriptionFactor "JUNB"    (Just "MA0490.1") Nothing (Just "JUNB_HUMAN.H11MO.0.A"),
  TranscriptionFactor "SMARCC1" Nothing Nothing Nothing,  -- Look for SMRC1
  TranscriptionFactor "FOS"     (Just "MA0476.1") Nothing (Just "JFOS_HUMAN.H11MO.0.A"),
  TranscriptionFactor "FOSB"    Nothing Nothing (Just "JFOSB_HUMAN.H11MO.0.A"),
  TranscriptionFactor "FOSL1"   (Just "MA0477.1") Nothing (Just "FOSL1_HUMAN.H11MO.0.A"),
  TranscriptionFactor "FOSL2"   (Just "MA0478.1") (Just "Fosl2") (Just "FOSL2_HUMAN.H11MO.0.A"),
  TranscriptionFactor "BCL11A"  Nothing Nothing Nothing,
  TranscriptionFactor "BCL11B"  Nothing Nothing Nothing,
  TranscriptionFactor "JDP2"    (Just "MA0655.1") Nothing (Just "JDP2_HUMAN.H11MO.0.D"),
  TranscriptionFactor "GATA1"   (Just "MA0140.2") (Just "Gata1") (Just "GATA1_HUMAN.H11MO.0.A"), -- Version 1 exists online, but not in the file that contains all PWMs
  TranscriptionFactor "GATA2"   (Just "MA0036.3") (Just "Gata2") (Just "GATA2_HUMAN.H11MO.0.A"),
  TranscriptionFactor "GATA3"   (Just "MA0037.3") (Just "Gata3") (Just "GATA3_HUMAN.H11MO.0.A"),
  TranscriptionFactor "GATA4"   Nothing (Just "Gata4") (Just "GATA4_HUMAN.H11MO.0.A"),
  TranscriptionFactor "GATA5"   (Just "MA0766.1") Nothing (Just "GATA5_HUMAN.H11MO.0.D"),
  TranscriptionFactor "GATA6"   (Just "MA1104.1") (Just "Gata6") (Just "GATA6_HUMAN.H11MO.0.A"),
  TranscriptionFactor "JUN"     Nothing Nothing (Just "JUN_HUMAN.H11MO.0.A"),
  TranscriptionFactor "JUND"    (Just "MA0491.1") (Just "JunD") (Just "JUND_HUMAN.H11MO.0.A"),
  TranscriptionFactor "BATF"    Nothing (Just "BATF") (Just "BATF_HUMAN.H11MO.1.A"),
  TranscriptionFactor "ATF3"    Nothing (Just "Atf3") (Just "ATF3_HUMAN.H11MO.0.A"),
  TranscriptionFactor "BACH1"   Nothing (Just "Bach1") (Just "BACH1_HUMAN.H11MO.0.A"),
  TranscriptionFactor "BACH2"   (Just "MA1101.1") (Just "Bach2") (Just "BACH2_HUMAN.H11MO.0.A"),
  TranscriptionFactor "NF-E2"   (Just "MA0841.1") (Just "NF-E2") (Just "NFE2_HUMAN.H11MO.0.A"),
  TranscriptionFactor "CEBPA"   (Just "MA0102.3") Nothing (Just "CEBPA_HUMAN.H11MO.0.A"),
  TranscriptionFactor "CEBPB"   (Just "MA0466.2") Nothing (Just "CEBPB_HUMAN.H11MO.0.A"),
  TranscriptionFactor "CEBPD"   (Just "MA0836.1") Nothing (Just "CEBPD_HUMAN.H11MO.0.A"),
  TranscriptionFactor "CEBPE"   (Just "MA0837.1") Nothing (Just "CEBPE_HUMAN.H11MO.0.A"),
  TranscriptionFactor "CEBPG"   (Just "MA0838.1") Nothing (Just "CEBPG_HUMAN.H11MO.0.A"),
  TranscriptionFactor "SPIB"    (Just "MA0081.1") (Just "SpiB") (Just "SPIB_HUMAN.H11MO.0.A"),
  TranscriptionFactor "IRF8"    (Just "MA0652.1") Nothing (Just "IRF8_HUMAN.H11MO.0.B"),
  TranscriptionFactor "SPI1"    (Just "MA0080.4") Nothing (Just "SPI1_HUMAN.H11MO.0.A"),
  TranscriptionFactor "SPI1D"   Nothing Nothing Nothing,
  TranscriptionFactor "LMO2"    Nothing Nothing Nothing,
  TranscriptionFactor "MESP1"   Nothing Nothing (Just "MESP1_HUMAN.H11MO.0.D"),
  TranscriptionFactor "MESP2"   Nothing Nothing Nothing,
  TranscriptionFactor "ID3"     Nothing Nothing Nothing,
  TranscriptionFactor "ID4"     (Just "MA0824.1") Nothing (Just "ID4_HUMAN.H11MO.0.D"),
  TranscriptionFactor "TCF12"   Nothing  (Just "Tcf12") (Just "HTF4_HUMAN.H11MO.0.A"),
  TranscriptionFactor "TCF4"    (Just "MA0830.1") (Just "Tcf4") (Just "ITF2_HUMAN.H11MO.0.C"),
  TranscriptionFactor "STAT1"   (Just "MA0137.3") (Just "STAT1") (Just "STAT1_HUMAN.H11MO.1.A"),
  TranscriptionFactor "STAT2"   Nothing Nothing (Just "STAT2_HUMAN.H11MO.0.A"),
  TranscriptionFactor "SPIC"    (Just "MA0687.1") Nothing (Just "SPIC_HUMAN.H11MO.0.D"),
  TranscriptionFactor "CTCF"    (Just "MA0139.1") (Just "CTCF") (Just "CTCF_HUMAN.H11MO.0.A"),
  TranscriptionFactor "IRF1"    (Just "MA0050.2") (Just "IRF1") (Just "IRF1_HUMAN.H11MO.0.A"),
  TranscriptionFactor "DBP"     (Just "MA0639.1") Nothing (Just "DBP_HUMAN.H11MO.0.B"),
  TranscriptionFactor "MAFK"    (Just "MA0496.2") (Just "MafK") (Just "MAFK_HUMAN.H11MO.1.A"),
  TranscriptionFactor "ATF4"    (Just "MA0833.1") (Just "Atf4") (Just "ATF4_HUMAN.H11MO.0.A"),
  TranscriptionFactor "ASCL1"   (Just "MA1100.1") (Just "Ascl1") (Just "ASCL1_HUMAN.H11MO.0.A"),
  TranscriptionFactor "ASCL2"   Nothing Nothing (Just "ASCL2_HUMAN.H11MO.0.D"),
  TranscriptionFactor "TCF3"    (Just "MA0522.2") (Just "Tcf3") (Just "TFE2_HUMAN.H11MO.0.A"),
  TranscriptionFactor "MYOD1"   Nothing Nothing (Just "MYOD1_HUMAN.H11MO.1.A"),
  TranscriptionFactor "ATOH8"   Nothing Nothing Nothing,
  TranscriptionFactor "MECOM"   Nothing Nothing (Just "EVI1_HUMAN.H11MO.0.B"),
  TranscriptionFactor "IRF3"    (Just "MA1418.1") Nothing (Just "IRF3_HUMAN.H11MO.0.B"),
  TranscriptionFactor "ZEB1"    (Just "MA0103.3") Nothing (Just "ZEB1_HUMAN.H11MO.0.A"),
  TranscriptionFactor "IRF9"    (Just "MA0653.1") Nothing (Just "IRF9_HUMAN.H11MO.0.C"),
  TranscriptionFactor "NHLH1"   (Just "MA0048.2") Nothing (Just "HEN1_HUMAN.H11MO.0.C"),
  TranscriptionFactor "LYL1"    Nothing Nothing (Just "LYL1_HUMAN.H11MO.0.A")]

loadHocomocoPatternsAndScoreThresholds :: FilePath -> FilePath -> [T.Text] -> IO Patterns
loadHocomocoPatternsAndScoreThresholds pwmFile thresholdsFile wantedHocomocoPatterns = do
   thresholds <- (M.fromList . map (\l -> (T.pack (takeWhile (/='\t') l), read (reverse (takeWhile (/='\t') (reverse l))))) . lines) <$> readFile thresholdsFile :: IO (M.Map T.Text Float)
   print (M.keys thresholds)
   patterns <- (M.fromList . filter ((`elem` wantedHocomocoPatterns) . fst)) <$> loadHocomocoMotifs pwmFile :: IO (M.Map T.Text [Pweight])
   print (M.keys patterns)
   let i = M.intersectionWith Pattern (fmap (floor . (*1000)) thresholds) patterns
   print i
   pure $ mkPatterns $ map snd $ M.toAscList i
   --pure $ mkPatterns $ map snd $ M.toAscList $ M.intersectionWith Pattern (fmap (floor . (*1000)) thresholds) patterns
   

main :: IO()
main = do
    args <- getArgs

    case args of
      [] -> do
        -- From http://schemer.buenrostrolab.com :
        --JUNB SMARCC1 FOSL2 FOSL1 JUND GATA1 JUN                  GATA2 FOS    BATF GATA3 BACH1 ATF3 BACH2 FOSB BCL11A BCL11B JDP2 GATA5 NFE2  SPI1D GATA4 CEBPB CEBPA SPIB IRF8      SPI1 CEBPD
        --x            Fosl2       JunD Gata1 Jun-AP1 or c-Jun-CRE Gata2 Fosl2? BATF Gata3 Bach1 Atf3 Bach2                               NF-E2       Gata4 CEBP        SpiB PU.1:IRF8           
        --
        --LMO2 GATA6 CEBPG MESP1 MESP2 ID3 ID4 TCF12 TCF4 STAT1 CEBPE SPIC CTCF IRF1 STAT2 DBP MAFK ATF4 ASCL1 TCF3 MYOD1 ATOH8 MECOM ASCL2 IRF3 ZEB1 IRF9 NHLH1 LYL1
        --x    Gata6                           Tcf12 Tcf4 STAT1            CTCF IRF1           MafK Atf4 Ascl1 Tcf3 
        --let wantedHocomocoPatterns = mapMaybe tfHocomocoId knownPatterns :: [T.Text]
        --patterns <- (mkPatterns . map snd . filter ((`elem` wantedHocomocoPatterns) . fst)) <$> loadHocomocoMotifs "HOCOMOCOv11_core_pwms_HUMAN_mono.txt"
    
        gata1 <- (map snd . filter ((== "GATA1_HUMAN.H11MO.0.A") . fst)) <$> loadHocomocoMotifs "HOCOMOCOv11_core_pwms_HUMAN_mono.txt"
        putStrLn "GATA1 motif:"
        mapM_ (mapM_ print) gata1

        patterns <- loadHocomocoPatternsAndScoreThresholds "HOCOMOCOv11_core_pwms_HUMAN_mono.txt" "hocomoco_thresholds.tab" ["GATA1_HUMAN.H11MO.0.A"]

        result <- runExceptT $ findPatterns (Chromosome "1") patterns ["chr1.bed"] "hg38.fa" "chr1.vcf.gz" "resultFile.tab.gz"
        print result
      [chrom, peakBedFiles, referenceGenomeFastaFile, vcfFile, motifsFile, score_thresholdsFile, outputFile] -> do
        let wantedHocomocoPatterns = mapMaybe tfHocomocoId knownPatterns :: [T.Text]
        --patterns <- (mkPatterns . map snd . filter ((`elem` wantedHocomocoPatterns) . fst)) <$> loadHocomocoMotifs motifsFile score_thresholdsFile
        patterns <- loadHocomocoPatternsAndScoreThresholds motifsFile score_thresholdsFile wantedHocomocoPatterns
        
        result <- runExceptT $ findPatterns (Chromosome $ T.pack chrom) patterns (splitOn "," peakBedFiles) referenceGenomeFastaFile vcfFile outputFile
        print result
      _ -> print ("Usage: xxx" :: String)
      
