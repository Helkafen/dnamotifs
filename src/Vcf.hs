{-# LANGUAGE OverloadedStrings #-}
{-# OPTIONS_GHC -fplugin Foreign.Storable.Generic.Plugin #-}

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
import qualified Data.Vector.Storable as STO

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


readVcfWithGenotypes :: FilePath -> [(Position0, Position0)] -> IO (Either Error [Variant])
readVcfWithGenotypes path regions = do
    putStrLn ("Loading VCF file " <> path)
    (parseVcfContent regions . dropWhile ("##" `T.isPrefixOf`)) <$> readGzippedLines path

sampleIdsInHeader :: Text -> V.Vector SampleId
sampleIdsInHeader header = V.fromList $ map SampleId $ drop 9 (T.splitOn "\t" header)

-- Hand tested
parseVcfContent :: [(Position0, Position0)] -> [Text] -> Either Error [Variant]
parseVcfContent regions vcfLines = case vcfLines of
    [] -> Left $ ParsingError "Empty vcf file"
    (header:rest) -> pure $ rights $ map (parseVariant sampleIdentifiers) (filterOrderedIntervals pos regions rest) -- TODO: exception for a left
        where sampleIdentifiers = sampleIdsInHeader header
              pos :: T.Text -> Position0
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

geno00_t1 = T.pack "0/0"
geno00_t2 = T.pack "0|0"
geno01_t1 = T.pack "0/1"
geno01_t2 = T.pack "0|1"
geno10_t1 = T.pack "1/0"
geno10_t2 = T.pack "1|0"
geno11_t1 = T.pack "1/1"
geno11_t2 = T.pack "1|1"

geno x | x == geno00_t1 = geno00
       | x == geno00_t2 = geno00
       | x == geno01_t1 = geno01
       | x == geno01_t2 = geno01
       | x == geno10_t1 = geno10
       | x == geno10_t2 = geno10
       | x == geno11_t1 = geno11
       | x == geno11_t2 = geno11
       | otherwise = error ("Unsupported genotype" <> show x)

variantIdParser :: Parser (Maybe Text)
variantIdParser = choice [none, rs] <?> "variant_name"
    where rs = (Just . T.pack) <$> many1 (notChar '\t')
          none = char '.' $> Nothing

variantParser :: V.Vector SampleId -> Parser Variant
variantParser sampleIdentifiers = do
    chr <- chromosomeParser <* tab
    pos <- (Position . (\x -> x - 1)) <$> decimal <* tab
    name <- variantIdParser <* tab
    ref <- (B.pack . map (unNuc . toNuc . fromIntegral . ord)) <$> many1 letter <* tab
    alt <- (B.pack . map (unNuc . toNuc . fromIntegral . ord))  <$> many1 letter <* tab
    skipField >> skipField >> skipField >> skipField
    geno <- (STO.fromListN (V.length sampleIdentifiers) . map geno . T.split (=='\t')) <$> takeWhile1 (/='\n')
    guard $ STO.length geno == V.length sampleIdentifiers
    _ <- option '0' (char '\n')
    return $ Variant chr pos name ref alt geno sampleIdentifiers
  where
    skipField = skipWhile (/= '\t') >> skip (== '\t')
    tab = skip (== '\t')

parseVariant :: V.Vector SampleId -> Text -> Either Error (Variant)
parseVariant sampleIdentifiers s = mapLeft (ParsingError . (\e -> s <> " " <> T.pack e)) (parseOnly (variantParser sampleIdentifiers) s)
