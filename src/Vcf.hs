{-# LANGUAGE OverloadedStrings #-}
{-# OPTIONS_GHC -fplugin Foreign.Storable.Generic.Plugin #-}

module Vcf (readVcfWithGenotypes, parseVariant, filterOrderedIntervals, parseVcfContent) where

import           Data.Text.Encoding (decodeUtf8With)
import           Data.Text.Encoding.Error (ignore)

import           Data.Text (Text)
import qualified Data.Text as T

import qualified Codec.Compression.GZip as GZip
import qualified Data.ByteString as B
import qualified Data.ByteString.Lazy as BL
import qualified Data.ByteString.Lazy.Char8 as BLC
import           Data.Either.Combinators (mapLeft)
import           Data.Char (ord)
import qualified Data.Vector.Storable as STO

import           Data.Attoparsec.ByteString
import           Data.Attoparsec.ByteString.Char8 (decimal, letter_ascii)

import qualified Data.Vector as V
import           Control.Monad (guard)
import           Data.Functor (($>))
import           Data.Either (rights, fromRight)
import qualified Data.DList as DList
import           Data.Word (Word8)

import Types

readGzippedLines :: FilePath -> IO [B.ByteString]
readGzippedLines path = map BL.toStrict . BLC.lines . GZip.decompress <$> BL.readFile path


filterOrderedIntervals :: Ord a => (b->a) -> [(a,a)] -> [b] -> [b]
filterOrderedIntervals pos r xs = DList.toList (go r xs)
    where go [] _ = DList.empty
          go ((start,end):rs) l = 
            let (vs, rest) =  span ((<=end) . pos) $ dropWhile ((<start) . pos) l
            in DList.fromList vs <> go rs rest


readVcfWithGenotypes :: FilePath -> [(Position0, Position0)] -> IO (Either Error [Variant])
readVcfWithGenotypes path regions = do
    putStrLn ("Loading VCF file " <> path)
    (parseVcfContent regions . dropWhile ("##" `B.isPrefixOf`)) <$> readGzippedLines path

sampleIdsInHeader :: B.ByteString -> V.Vector SampleId
sampleIdsInHeader header = V.fromList $ map (SampleId . decodeUtf8With ignore) $ drop 9 (B.split (fromIntegral $ ord '\t') header)

-- Hand tested
-- The first Variant is parsed by the slow parser in order to check more things. We assume that the other lines will have the same quirks
parseVcfContent :: [(Position0, Position0)] -> [B.ByteString] -> Either Error [Variant]
parseVcfContent regions vcfLines = case vcfLines of
    [] -> Left $ ParsingError "Empty vcf file"
    (header:rest) -> pure $ rights $ map (parseVariant sampleIdentifiers) (filterOrderedIntervals pos regions rest) -- TODO: exception for a left
        where sampleIdentifiers = sampleIdsInHeader header
              pos :: B.ByteString -> Position0
              pos = Position . (\p -> p - 1) . fromRight (error "Bad position in vcf") . parseInt . B.takeWhile (/= tab) . B.drop 1  . B.dropWhile (/= tab)
              parseInt = parseOnly decimal

--digit = satisfy isDigit
--    where isDigit w = w >= 48 && w <= 57

chromosomeParser :: Parser Chromosome
chromosomeParser = choice [autosomeParser, xParser, yParser]
    where autosomeParser = do
            _ <- option "" (string "chr")
            d <- decimal
            guard $ d > 0 && d < 23
            pure (Chromosome $ T.pack $ show (d :: Integer))
          xParser = string "X" $> Chromosome "X"
          yParser = string "Y" $> Chromosome "Y"

geno00_t1, geno00_t2 :: B.ByteString
geno00_t1 = B.pack (map (fromIntegral . ord) "0/0")
geno00_t2 = B.pack (map (fromIntegral . ord) "0|0")

geno01_t1, geno01_t2 :: B.ByteString
geno01_t1 = B.pack (map (fromIntegral . ord) "0/1")
geno01_t2 = B.pack (map (fromIntegral . ord) "0|1")

geno10_t1, geno10_t2 :: B.ByteString
geno10_t1 = B.pack (map (fromIntegral . ord) "1/0")
geno10_t2 = B.pack (map (fromIntegral . ord) "1|0")

geno11_t1, geno11_t2 :: B.ByteString
geno11_t1 = B.pack (map (fromIntegral . ord) "1/1")
geno11_t2 = B.pack (map (fromIntegral . ord) "1|1")

toGeno :: B.ByteString -> Genotype
toGeno x | x == geno00_t1 = geno00
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
    where rs = (Just . decodeUtf8With ignore . B.pack) <$> many1 (notWord8 tab)
          none = word8 dot $> Nothing

dot, tab, newline :: Word8
dot = (fromIntegral $ ord '.')
tab = (fromIntegral $ ord '\t')
newline = (fromIntegral $ ord '\n')

variantParser :: V.Vector SampleId -> Parser Variant
variantParser sampleIdentifiers = do
    chr <- chromosomeParser <* word8 tab
    pos <- (Position . (\x -> x - 1)) <$> decimal <* word8 tab
    name <- variantIdParser <* word8 tab
    ref <- (B.pack . map (unNuc . toNuc . fromIntegral . ord)) <$> many1 letter_ascii <* word8 tab
    alt <- (B.pack . map (unNuc . toNuc . fromIntegral . ord))  <$> many1 letter_ascii <* word8 tab
    skipField >> skipField >> skipField >> skipField
    geno <- STO.fromListN (V.length sampleIdentifiers) . map toGeno . B.split tab <$> takeWhile1 (/=newline)
    guard $ STO.length geno == V.length sampleIdentifiers
    _ <- option newline (word8 newline)
    return $ Variant chr pos name ref alt geno sampleIdentifiers
  where
    skipField = skipWhile (/= tab) >> skip (== tab)

parseVariant :: V.Vector SampleId -> B.ByteString -> Either Error (Variant)
parseVariant sampleIdentifiers s = 
    mapLeft (ParsingError . (\e -> (decodeUtf8With ignore s) <> " " <> T.pack e)) (parseOnly (variantParser sampleIdentifiers) s) 

