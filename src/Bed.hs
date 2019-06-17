{-# LANGUAGE OverloadedStrings #-}

module Bed (parseBedContent, readPeaks) where

import           Data.Attoparsec.Text
import qualified Data.Text as T
import qualified Data.Text.IO as TIO
import           Control.Applicative ((<*))
import           Data.Functor (($>))
import           Control.Monad (guard)
import           Data.Range.Range (Range(..))

import Types

readPeaks :: Chromosome -> FilePath -> IO (Either String [Range Position0])
readPeaks chr path = do
    content <- TIO.readFile path
    return $ parseBedContent chr content --return [(Position 1000, Position 1010000)]

parseBedContent :: Chromosome -> T.Text -> Either String [Range Position0]
parseBedContent chr content = (map snd . filter ((== chr) . fst)) <$> parseOnly parser content

parser :: Parser [(Chromosome, Range Position0)]
parser = lineParser `sepBy` (char '\n')

lineParser :: Parser (Chromosome, Range Position0)
lineParser = do
     chr <- chromosomeParser <* char '\t'
     s <- decimal <* char '\t'
     e <- decimal
     guard (s < e)
     pure (chr, SpanRange (Position s) (Position e))

chromosomeParser :: Parser Chromosome
chromosomeParser = choice [autosomeParser, xParser, yParser]
    where autosomeParser = do
            _ <- option "" (string "chr")
            d <- decimal
            guard $ d > 0 && d < 23
            pure (Chromosome $ T.pack $ show (d :: Integer))
          xParser = string "X" $> Chromosome "X"
          yParser = string "Y" $> Chromosome "Y"