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
      