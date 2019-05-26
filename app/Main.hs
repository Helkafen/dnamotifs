{-# LANGUAGE OverloadedStrings #-}
module Main where

import           Data.Vector.Storable            (Vector)
import qualified Data.Vector.Storable            as V
import qualified Data.ByteString                 as B
import           Data.Monoid                     ((<>))
import           Data.Time.Clock.POSIX (getPOSIXTime)
import           Foreign.C.Types                 (CInt)
import Types
import PatternFind
import Lib (findPatterns)

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

bench :: IO ()
bench = do
    matches <- findPatternsInBlock (mkNucleotideAndPositionBlock inputData) (mkPatterns patterns)
    print (expected == matches)
  where numberOfPeople = 10000 :: Int

        inputData :: [BaseSequencePosition]
        inputData =  [BaseSequencePosition inputDataSample0 inputDataPositions, BaseSequencePosition inputDataSample1 inputDataPositions] ++ (replicate (numberOfPeople - 2) (BaseSequencePosition inputDataSample2 inputDataPositions))

        inputDataPositions :: Vector CInt
        inputDataPositions = V.fromList $ take (B.length inputDataSample0) [0..]

        patterns = [pattern_CG, pattern_cCGA, pattern_CGA] ++ replicate 50 pattern_cccccccccc ++ [pattern_CG]

        expected = [Match {mPatternId = 53, mScore = 1000, mPosition = 0, mSampleId = 1, mMatched = [c, g]}
                    ,Match {mPatternId = 2,  mScore = 1000, mPosition = 0, mSampleId = 1, mMatched = [c, g]}
                    ,Match {mPatternId = 0,  mScore = 1000, mPosition = 0, mSampleId = 1, mMatched = [c, g]}
                    ,Match {mPatternId = 53, mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c, g]}
                    ,Match {mPatternId = 2,  mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c, g]}
                    ,Match {mPatternId = 0,  mScore = 1000, mPosition = 4, mSampleId = 0, mMatched = [c, g]}
                    ,Match {mPatternId = 1,  mScore = 961,  mPosition = 3, mSampleId = 0, mMatched = [a, c, g, a]}]

main :: IO()
main = do
    print ("Hello world" :: String)
    t1 <- getPOSIXTime
    --bench
    _ <- findPatterns (Chromosome "1") "/media/seb/TERA/lab/hg38.fa" "/home/seb/masters/topmed/exome_modified_header_100000_chr_100s.vcf.gz" "resultFile.tab"
    t2 <- getPOSIXTime
    print (round $ (t2 - t1) * 1000 :: Integer) -- milliseconds
    pure ()
