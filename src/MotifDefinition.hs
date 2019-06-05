{-# LANGUAGE OverloadedStrings #-}
module MotifDefinition (loadHocomocoMotifs, loadHomerMotifs, parseHomerMotifsContent, parseHocomocoMotifsContent) where

import           Data.Attoparsec.Text
import qualified Data.Text as T
import qualified Data.Text.IO as TIO
import           Control.Applicative ((<*))    
import           Data.Either (partitionEithers)

import Types


loadHocomocoMotifs :: FilePath -> IO [(T.Text, Pattern)]
loadHocomocoMotifs path = do putStr ("Load motif file " <> path <> " ... ")
                             content <- TIO.readFile path
                             case parseHocomocoMotifsContent content of
                                Left errors -> print (show errors) >> error ("Motif loading error")
                                Right patterns -> putStrLn ("Done") >> return patterns

loadHomerMotifs :: FilePath -> IO [(T.Text, Pattern)]
loadHomerMotifs path = do putStr ("Load motif file " <> path <> " ... ")
                          content <- TIO.readFile path
                          case parseHomerMotifsContent content of
                            Left errors -> print (show errors) >> error ("Motif loading error")
                            Right patterns -> putStrLn ("Done") >> return patterns

parseHocomocoMotifsContent :: T.Text -> Either [String] [(T.Text, Pattern)]
parseHocomocoMotifsContent content =
    let paragraphs =  filter ((>0) . T.length) $ T.splitOn ">" content
        parsed = map (parseOnly hocomocoPatternParser) paragraphs
    in case partitionEithers parsed of
            ([], patterns) -> Right patterns
            (errors, _) -> Left errors

parseHomerMotifsContent :: T.Text -> Either [String] [(T.Text, Pattern)]
parseHomerMotifsContent content =
    let paragraphs =  filter ((>0) . T.length) $ T.splitOn ">" content
        parsed = map (parseOnly homerPatternParser) paragraphs
    in case partitionEithers parsed of
            ([], patterns) -> Right patterns
            (errors, _) -> Left errors

hocomocoPatternParser :: Parser (T.Text, Pattern)
hocomocoPatternParser = do
    name <- takeWhile1 (/= '\n')  <* char '\n'
    weights <- weightsParserHocomoco `sepBy` char '\n' 
    pure (name, weights)

homerPatternParser :: Parser (T.Text, Pattern)
homerPatternParser = do
    _ <- takeWhile1 (/= '\t')
    name <- char '\t' *> takeWhile1 (/= '(')
    _ <- takeWhile1 (/= '\n') <* char '\n'
    weights <- weightsParser `sepBy` char '\n'
    pure (name, weights)

weightsParserHocomoco :: Parser Pweight
weightsParserHocomoco = do
    aWeight <- realToFrac <$> (double <* char '\t')
    cWeight <- realToFrac <$> (double <* char '\t')
    gWeight <- realToFrac <$> (double <* char '\t')
    tWeight <- realToFrac <$> double
    let total = aWeight + cWeight + gWeight + tWeight
    pure (Pweight (aWeight / total) (cWeight / total) (gWeight / total) (tWeight / total))

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

