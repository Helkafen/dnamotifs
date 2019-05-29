{-# LANGUAGE OverloadedStrings #-}
module MotifDefinition where

import           Data.Attoparsec.Text
import qualified Data.Text as T
import qualified Data.Text.IO as TIO
import           Control.Applicative ((<*))    
import           Data.Either (partitionEithers)

import Types

loadMotifs :: FilePath -> IO [(T.Text, Pattern)]
loadMotifs path = do putStr ("Load motif file " <> path <> "...")
                     content <- TIO.readFile path
                     case parseMotifsContent content of
                        Left errors -> print (show errors) >> error ("Motif loading error")
                        Right patterns -> putStrLn ("Done") >> return patterns


parseMotifsContent :: T.Text -> Either [String] [(T.Text, Pattern)]
parseMotifsContent content = 
    let paragraphs =  filter ((>0) . T.length) $ T.splitOn ">" content
        parsed = map (parseOnly patternParser) paragraphs
    in case partitionEithers parsed of
            ([], patterns) -> Right patterns
            (errors, _) -> Left errors

patternParser :: Parser (T.Text, Pattern)
patternParser = do
    _ <- takeWhile1 (/= '\t')
    name <- char '\t' *> takeWhile1 (/= '(')
    _ <- takeWhile1 (/= '\n') <* char '\n'
    weights <- weightsParser `sepBy` char '\n' 
    pure (name, weights)

weightsParser :: Parser Pweight
weightsParser = do
    aWeight <- realToFrac <$> (double <* char '\t')
    cWeight <- realToFrac <$> (double <* char '\t')
    gWeight <- realToFrac <$> (double <* char '\t')
    tWeight <- realToFrac <$> double
    pure (Pweight aWeight cWeight gWeight tWeight)

-- >ATGACTCATC     AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer        6.049537        -1.782996e+03   0       9805.3,5781.0,3085.1,2715.0,-- 0.00e+00
-- 0.419   0.275   0.277   0.028
-- 0.001   0.001   0.001   0.997
-- >SCCTSAGGSCAW   AP-2gamma(AP2)/MCF7-TFAP2C-ChIP-Seq(GSE21234)/Homer     6.349794        -24627.169865   0       T:26194.0(44.86%),B:5413.7-- (9.54%),P:1e-10695
-- 0.005   0.431   0.547   0.017
-- 0.001   0.997   0.001   0.001
-- 0.001   0.947   0.001   0.051

