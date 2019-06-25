{-# LANGUAGE OverloadedStrings #-}
module MotifDefinition
  ( loadHocomocoMotifs
  , parseHocomocoMotifsContent
  )
where

import           Data.Attoparsec.Text
import qualified Data.Text                     as T
import           Import


loadHocomocoMotifs
  :: HasLogFunc env => FilePath -> RIO env [(T.Text, [Pweight])]
loadHocomocoMotifs path = do
  content <- readFileUtf8 path
  case parseHocomocoMotifsContent content of
    Left errors -> throwM
      (  MotifLoadingError
      $  "Load motif file "
      <> path
      <> ": Error. "
      <> show errors
      )
    Right patterns -> do logInfo (display $ "Loaded motif file: " <> (T.pack path))
                         return patterns

parseHocomocoMotifsContent :: T.Text -> Either [String] [(T.Text, [Pweight])]
parseHocomocoMotifsContent content =
  let paragraphs = filter ((> 0) . T.length) $ T.splitOn ">" content
      parsed     = map (parseOnly hocomocoPatternParser) paragraphs
  in  case partitionEithers parsed of
        ([]    , patterns) -> Right patterns
        (errors, _       ) -> Left errors

hocomocoPatternParser :: Parser (T.Text, [Pweight])
hocomocoPatternParser = do
  name    <- takeWhile1 (/= '\n') <* char '\n'
  weights <- weightsParser `sepBy` char '\n'
  pure (name, weights)

weightsParser :: Parser Pweight
weightsParser = do
  aWeight <- double <* char '\t'
  cWeight <- double <* char '\t'
  gWeight <- double <* char '\t'
  tWeight <- double
  pure
    (Pweight (round $ 1000 * (aWeight :: Double))
             (round $ 1000 * (cWeight :: Double))
             (round $ 1000 * (gWeight :: Double))
             (round $ 1000 * (tWeight :: Double))
    )


-- >ATGACTCATC     AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer        6.049537        -1.782996e+03   0       9805.3,5781.0,3085.1,2715.0,-- 0.00e+00
-- 0.419   0.275   0.277   0.028
-- 0.001   0.001   0.001   0.997
-- >SCCTSAGGSCAW   AP-2gamma(AP2)/MCF7-TFAP2C-ChIP-Seq(GSE21234)/Homer     6.349794        -24627.169865   0       T:26194.0(44.86%),B:5413.7-- (9.54%),P:1e-10695
-- 0.005   0.431   0.547   0.017
-- 0.001   0.997   0.001   0.001
-- 0.001   0.947   0.001   0.051

