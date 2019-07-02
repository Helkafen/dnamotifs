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
  }

config :: Parser Config
config = Config
       <$> strOption (long "chromosome"          <> short 'c'  <> metavar "CHR"                 <> help "Ex: 1")
       <*> strOption (long "ref-genome"          <> short 'r'  <> metavar "FASTA"               <> help "Ex: hg38.fa")
       <*> strOption (long "motifs"              <> short 'm'  <> metavar "MOTIFS"              <> help "Ex: CTCF_HUMAN.H11MO.0.A,IRF1_HUMAN.H11MO.0.A")
       <*> strOption (long "bed-files"           <> short 'b'  <> metavar "BED"                 <> help "Ex: GATA1.bed,GATA2.bed")
       <*> strOption (long "vcf"                 <> short 'i'  <> metavar "VCF"                 <> help "Ex: chr1.vcf.gz")
       <*> strOption (long "output-file"         <> short 'o'  <> metavar "OUTPUT"              <> help "Ex: output.tab.gz")

main :: IO()
main = runSimpleApp $ do
        let p = prefs (disambiguate <> showHelpOnEmpty <> columns 100)
        cfg <- liftIO $ customExecParser p (info (config <**> helper) ( fullDesc <> progDesc "DNAMotif finds PWM motifs in a VCF files" ))
        patterns <- loadHocomocoPatternsAndScoreThresholds (map T.pack $ splitOn "," (cfgMotifNames cfg))
        success <- findPatterns (Chromosome $ T.pack (cfgChromosome cfg)) patterns (splitOn "," (cfgBedFiles cfg)) (cfgReferenceGenome cfg) (cfgInputVCF cfg) (cfgOutput cfg)
        if success
          then logInfo "Success"
          else logError "Error"

