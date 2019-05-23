{-# LANGUAGE OverloadedStrings #-}
module Main where

import           Data.Vector.Storable            (Vector)
import qualified Data.Vector.Storable            as V
import           Data.Monoid                     ((<>))
import           Data.Time.Clock.POSIX (getPOSIXTime)
import Types
import PatternFind
import Lib (findPatterns)

inputDataSample0 :: Vector Nucleotide
inputDataSample0 = V.fromList [a,a,a,a,c,g,a] <> V.fromList (replicate 93 a)

inputDataSample1 :: Vector Nucleotide
inputDataSample1 = V.fromList [c,g,a,a,a,a,a] <> V.fromList (replicate 93 a) 

inputDataSample2 :: Vector Nucleotide
inputDataSample2 = V.fromList [a,a,a,a,a,a,a] <> V.fromList (replicate 93 a)

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

        inputData :: [(Vector Nucleotide, Vector Position)]
        inputData =  [(inputDataSample0, inputDataPositions), (inputDataSample1, inputDataPositions)] ++ (take (numberOfPeople - 2) (repeat (inputDataSample2, inputDataPositions)))

        inputDataPositions :: Vector Position
        inputDataPositions = V.fromList $ take (V.length inputDataSample0) (map Position [0..])

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
    print "Hello world"
    t1 <- getPOSIXTime
    --bench
    _ <- findPatterns (Chromosome "1") "/media/seb/TERA/lab/hg38.fa" "/home/seb/masters/topmed/exome_modified_header_100000_chr_100s.vcf.gz" "resultFile.tab"
    t2 <- getPOSIXTime
    print (round $ (t2 - t1) * 1000 :: Integer) -- milliseconds
    pure ()
