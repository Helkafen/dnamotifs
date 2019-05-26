{-# LANGUAGE OverloadedStrings #-}

module Vcf (readVcfWithGenotypes, parseVariant, filterOrderedIntervals, parseVcfContent) where

--import Protolude hiding (decodeUtf8With)
import           Data.Text (Text)
import qualified Data.Text.Lazy as TL
import qualified Codec.Compression.GZip as GZip
import qualified Data.ByteString as B
import qualified Data.ByteString.Lazy as BL
import           Data.Text.Lazy.Encoding (decodeUtf8With)
import           Data.Text.Encoding.Error (ignore)
import           Data.Either.Combinators (mapLeft)
import           Data.Char (ord)

import           Data.Attoparsec.Text
import qualified Data.Text as T
import qualified Data.Vector as V
import           Control.Monad (guard)
import           Data.Functor (($>))
import           Data.Either (rights, fromRight)
import qualified Data.DList as DList

import Types

readGzippedLines :: FilePath -> IO [Text]
readGzippedLines path = map TL.toStrict . TL.lines . decodeUtf8With ignore . GZip.decompress <$> BL.readFile path


filterOrderedIntervals :: Ord a => (b->a) -> [(a,a)] -> [b] -> [b]
filterOrderedIntervals pos r xs = DList.toList (go r xs)
    where go [] _ = DList.empty
          go ((start,end):rs) l = 
            let (vs, rest) =  span ((<=end) . pos) $ dropWhile ((<start) . pos) l
            in DList.fromList vs <> go rs rest


readVcfWithGenotypes :: FilePath -> [(Position ZeroBased, Position ZeroBased)] -> IO (Either Error [Variant ZeroBased])
readVcfWithGenotypes path regions = (parseVcfContent regions . dropWhile ("##" `T.isPrefixOf`)) <$> readGzippedLines path

sampleIdsInHeader :: Text -> V.Vector SampleId
sampleIdsInHeader header = V.fromList $ map SampleId $ drop 9 (T.splitOn "\t" header)

-- Hand tested
parseVcfContent :: [(Position ZeroBased, Position ZeroBased)] -> [Text] -> Either Error [Variant ZeroBased]
parseVcfContent regions vcfLines = case vcfLines of
    [] -> Left $ ParsingError "Empty vcf file"
    (header:rest) -> pure $ rights $ map (parseVariant sampleIdentifiers) (filterOrderedIntervals pos regions rest) -- TODO: exception for a left
        where sampleIdentifiers = sampleIdsInHeader header
              pos :: T.Text -> Position ZeroBased
              pos = Position . (\p -> p - 1) . fromRight (error "Bad position in vcf") . parseInt . T.takeWhile (/='\t') . T.drop 1  . T.dropWhile (/='\t')
              parseInt = parseOnly (signed decimal)


chromosomeParser :: Parser Chromosome
chromosomeParser = choice [autosomeParser, xParser, yParser]
    where autosomeParser = do
            _ <- option "" (string "chr")
            d <- decimal
            guard $ d > 0 && d < 23
            pure (Chromosome $ T.pack $ show (d :: Integer))
          xParser = string "X" $> Chromosome "X"
          yParser = string "Y" $> Chromosome "Y"

genoParser :: Parser Genotype
genoParser = choice [
    string "0/0" $> Geno00,
    string "0|0" $> Geno00,
    string "0/1" $> Geno01,
    string "1/0" $> Geno10,
    string "1/1" $> Geno11,
    string "0|1" $> Geno01,
    string "1|0" $> Geno10,
    string "1|1" $> Geno11]

variantIdParser :: Parser (Maybe Text)
variantIdParser = choice [none, rs] <?> "variant_name"
    where rs = (Just . T.pack) <$> many1 (notChar '\t')
          none = char '.' $> Nothing

variantParser :: V.Vector SampleId -> Parser (Variant ZeroBased)
variantParser sampleIdentifiers = do
    chr <- chromosomeParser <* tab
    pos <- (Position . (\x -> x - 1)) <$> decimal <* tab
    name <- variantIdParser <* tab
    ref <- (B.pack . map (toNuc . ord)) <$> many1 letter <* tab
    alt <- (B.pack . map (toNuc . ord))  <$> many1 letter <* tab
    skipField >> skipField >> skipField >> skipField
    geno <- V.fromList <$> genoParser `sepBy` char '\t'
    guard $ V.length geno == V.length sampleIdentifiers
    return $ Variant chr pos name ref alt geno sampleIdentifiers
  where
    skipField = skipWhile (/= '\t') >> skip (== '\t')
    tab = skip (== '\t')

parseVariant :: V.Vector SampleId -> Text -> Either Error (Variant ZeroBased)
parseVariant sampleIdentifiers s = mapLeft (ParsingError . (\e -> s <> " " <> T.pack e)) (parseOnly (variantParser sampleIdentifiers) s)

toNuc :: Int -> Nucleotide
toNuc 65 = a
toNuc 67 = c
toNuc 71 = g
toNuc 84 = t
toNuc 78 = n
toNuc 97 = a
toNuc 99 = c
toNuc 103 = g
toNuc 116 = t
toNuc 110 = n
toNuc other = error $ "Bad nucleotide " <> show other