{-# LANGUAGE OverloadedStrings #-}

module Bed (parseBedContent, readPeaks, readAllPeaks) where

import           Data.Attoparsec.Text
import qualified Data.Map as M
import qualified Data.Text as T
import           Data.List (sort)
import           System.FilePath.Posix (takeFileName)

import Import
import Range

readAllPeaks :: HasLogFunc env => Chromosome -> [FilePath] -> RIO env (Ranges Position0, M.Map T.Text (Ranges Position0))
readAllPeaks chr peakFiles = do
    xs <- mapM (readPeaks chr) peakFiles
    _ <- mapM (\(peakFile, ranges) -> logInfo (display $ "Bed file: " <> (T.pack peakFile) <> ": " <> (T.pack (show (rangesLength ranges))) <> " peaks")) (zip peakFiles xs)
    pure (mconcat xs, simplifyFilenames $ M.fromList (zip peakFiles xs))

simplifyFilenames :: M.Map FilePath a -> M.Map T.Text a
simplifyFilenames d = M.mapKeys T.pack (if M.size fileNamesOnly == M.size d then fileNamesOnly else d)
    where fileNamesOnly = M.mapKeys takeFileName d


readPeaks :: Chromosome -> FilePath -> RIO env (Ranges Position0)
readPeaks chr path = do
    content <- readFileUtf8 path
    case parseBedContent chr content of
        Left err -> throwM (BedLoadingError $ "File: " <> path <> ". Error: " <> err)
        Right parsed ->
            if (rangesLength parsed == 0)
                then throwM (BedLoadingError $ "File: " <> path <> ". No peak loaded")
                else pure parsed


parseBedContent :: Chromosome -> T.Text -> Either String (Ranges Position0)
parseBedContent chr content =
    case (sort . map snd . filter ((== chr) . fst)) <$> parseOnly parser content of
        Left err -> Left err
        Right intervals -> case mkRanges (map (\(s,e) -> Range s e) intervals) of
                             Just ranges -> Right ranges
                             Nothing -> Left "Peaks are not disjoint"



parser :: Parser [(Chromosome, (Position0, Position0))]
parser = lineParser `sepBy` (char '\n')

lineParser :: Parser (Chromosome, (Position0, Position0))
lineParser = do
     chr <- chromosomeParser <* char '\t'
     s <- decimal <* char '\t'
     e <- decimal
     guard (s < e)
     pure (chr, ((Position s), (Position e)))

chromosomeParser :: Parser Chromosome
chromosomeParser = choice [autosomeParser, xParser, yParser]
    where autosomeParser = do
            _ <- option "" (string "chr")
            d <- decimal
            guard $ d > 0 && d < 23
            pure (Chromosome $ T.pack $ show (d :: Integer))
          xParser = string "X" $> Chromosome "X"
          yParser = string "Y" $> Chromosome "Y"