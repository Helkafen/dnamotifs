{-# LANGUAGE OverloadedStrings #-}
module Main where

import qualified Data.Text                       as T
import           Types
import           Run                             (findPatterns)
import           MotifDefinition                 (loadHocomocoPatternsAndScoreThresholds)
import           Data.List.Split                 (splitOn)
import           RIO.Prelude.Simple              (runSimpleApp)
import           Options.Applicative


import Import

data Config = Config
  { cfgChromosome        :: String
  , cfgReferenceGenome   :: String
  , cfgMotifNames        :: String
  , cfgBedFiles          :: String
  , cfgInputVCF          :: String
  , cfgOutput            :: String
  , cfgHomocomo          :: String
  , cfgHomocomoThreshold :: String
  }

config :: Parser Config
config = Config
       <$> strOption (long "chromosome"          <> short 'c'  <> metavar "CHR"                 <> help "Ex: 1")
       <*> strOption (long "ref-genome"          <> short 'r'  <> metavar "FASTA"               <> help "Ex: hg38.fa")
       <*> strOption (long "motifs"              <> short 'm'  <> metavar "MOTIFS"              <> help "Ex: CTCF_HUMAN.H11MO.0.A,IRF1_HUMAN.H11MO.0.A")
       <*> strOption (long "bed-files"           <> short 'b'  <> metavar "BED"                 <> help "Ex: GATA1.bed,GATA2.bed")
       <*> strOption (long "vcf"                 <> short 'i'  <> metavar "VCF"                 <> help "Ex: chr1.vcf.gz")
       <*> strOption (long "output-file"         <> short 'o'  <> metavar "OUTPUT"              <> help "Ex: output.tab.gz")
       <*> strOption (long "hocomoco-file"       <> short 'h'  <> metavar "HOCOMOCO"            <> help "Ex: HOCOMOCOv11_core_pcms_HUMAN_mono.txt")
       <*> strOption (long "hocomoco-thresholds" <> short 't'  <> metavar "HOCOMOCO_THRESHOLD"  <> help "Ex: hocomoco_thresholds.tab")

main :: IO()
main = runSimpleApp $ do
        let p = prefs (disambiguate <> showHelpOnEmpty <> columns 100)
        cfg <- liftIO $ customExecParser p (info (config <**> helper) ( fullDesc <> progDesc "DNAMotif finds PWM motifs in a VCF files" ))
        patterns <- loadHocomocoPatternsAndScoreThresholds (cfgHomocomo cfg) (cfgHomocomoThreshold cfg) (map T.pack $ splitOn "," (cfgMotifNames cfg))
        success <- findPatterns (Chromosome $ T.pack (cfgChromosome cfg)) patterns (splitOn "," (cfgBedFiles cfg)) (cfgReferenceGenome cfg) (cfgInputVCF cfg) (cfgOutput cfg)
        if success
          then logInfo "Success"
          else logError "Error"

