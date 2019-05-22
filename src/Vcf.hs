{-# LANGUAGE OverloadedStrings #-}

module Vcf (readVcfWithGenotypes) where

--import Protolude hiding (decodeUtf8With)
import           Data.Text (Text)
import qualified Data.Text.Lazy as TL
import qualified Codec.Compression.GZip as GZip
import qualified Data.ByteString.Lazy as B
import           Data.Text.Lazy.Encoding (decodeUtf8With)
import           Data.Text.Encoding.Error (ignore)
import           Data.Either.Combinators (mapLeft)
import           Data.Char (ord)

import           Data.Attoparsec.Text
import qualified Data.Text as T
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import           Control.Monad (guard)
import           Data.Functor (($>))
import           Data.Either (rights)

import Types

readGzippedLines :: FilePath -> IO [Text]
readGzippedLines path = map TL.toStrict . TL.lines . decodeUtf8With ignore . GZip.decompress <$> B.readFile path

readVcfWithGenotypes :: FilePath -> IO (Either Error [Variant])
readVcfWithGenotypes path = (parseVcfContent . dropWhile ("##" `T.isPrefixOf`)) <$> readGzippedLines path

sampleIdsInHeader :: Text -> V.Vector SampleId
sampleIdsInHeader header = V.fromList $ map SampleId $ drop 9 (T.splitOn "\t" header)

-- Hand tested
parseVcfContent :: [Text] -> Either Error [Variant]
parseVcfContent vcfLines = case vcfLines of
    [] -> Left $ ParsingError "Empty vcf file"
    (header:rest) -> pure $ rights $ map (parseVariant sampleIdentifiers) rest -- TODO: exception for a left
        where sampleIdentifiers = sampleIdsInHeader header


chromosomeParser :: Parser Chromosome
chromosomeParser = choice [autosomeParser, xParser, yParser]
    where autosomeParser = do
            d <- decimal
            guard $ d > 0 && d < 23
            pure (Chromosome $ T.pack $ show (d :: Integer))
          xParser = string "X" $> Chromosome "X"
          yParser = string "Y" $> Chromosome "Y"

genoParser :: Parser Genotype
genoParser = choice [
    string "0/0" $> Geno00,
    string "0/1" $> Geno01,
    string "1/0" $> Geno10,
    string "1/1" $> Geno11,
    string "0|0" $> Geno00,
    string "0|1" $> Geno01,
    string "1|0" $> Geno10,
    string "1|1" $> Geno11]

variantIdParser :: Parser (Maybe Text)
variantIdParser = choice [none, rs] <?> "variant_name"
    where rs = (Just . T.pack) <$> many1 (notChar '\t')
          none = char '.' $> Nothing

variantParser :: V.Vector SampleId -> Parser Variant
variantParser sampleIdentifiers = do
    chr <- chromosomeParser <* tab
    pos <- (Position . (\x -> x - 1)) <$> decimal <* tab
    name <- variantIdParser <* tab
    ref <- (U.fromList . map (fromIntegral . ord)) <$> many1 letter <* tab
    alt <- (U.fromList . map (fromIntegral . ord))  <$> many1 letter <* tab
    skipField >> skipField >> skipField >> skipField
    geno <- V.fromList <$> genoParser `sepBy` char '\t'
    guard $ V.length geno == V.length sampleIdentifiers
    return $ Variant chr pos name ref alt geno sampleIdentifiers
  where
    skipField = skipWhile (/= '\t') >> skip (== '\t')
    tab = skip (== '\t')

parseVariant :: V.Vector SampleId -> Text -> Either Error Variant
parseVariant sampleIdentifiers s = mapLeft (ParsingError . (\e -> s <> " " <> T.pack e)) (parseOnly (variantParser sampleIdentifiers) s)