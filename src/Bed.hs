{-# LANGUAGE OverloadedStrings #-}

module Bed (parseBedContent, readPeaks, merge, readAllPeaks) where

import           Data.Attoparsec.Text
import qualified Data.Map as M
import qualified Data.Text as T
import qualified Data.Text.IO as TIO
import           Control.Applicative ((<*))
import           Data.Functor (($>))
import           Control.Monad (guard)
import           Data.Range.Range (Range(..), mergeRanges)
import           Control.Monad.Trans.Except
import           Control.Monad.Except (throwError, lift)
import           Data.Either (partitionEithers)

import Types

readAllPeaks :: Chromosome -> [FilePath] -> ExceptT Error IO ([Range Position0], M.Map FilePath [Range Position0])
readAllPeaks chr peakFiles = do
    eitherBed <- lift $ mapM (readPeaks chr) peakFiles
    case partitionEithers eitherBed of
        (errors@(_:_),_) -> throwError (BedLoadingError $ "Some bed files could not be read: " <> show errors)
        ([], xs) -> do
            _ <- lift $ mapM (\(peakFile, ranges) -> putStrLn ("Bed file: "<>peakFile <> ": " <> show (length ranges) <> " peaks")) (zip peakFiles xs)
            pure (merge (concat xs), M.fromList (zip peakFiles xs))

merge :: (Ord a, Enum a) => [Range a] -> [Range a]
merge = mergeRanges

readPeaks :: Chromosome -> FilePath -> IO (Either Error [Range Position0])
readPeaks chr path = do
    content <- TIO.readFile path
    case parseBedContent chr content of
        Left err -> return $ Left (BedLoadingError $ "File: " <> path <> ". Error: " <> err)
        Right parsed ->
            if (length parsed == 0)
                then return $ Left (BedLoadingError $ "File: " <> path <> ". No peak loaded")
                else return $ Right parsed


parseBedContent :: Chromosome -> T.Text -> Either String [Range Position0]
parseBedContent chr content = (merge . map snd . filter ((== chr) . fst)) <$> parseOnly parser content

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