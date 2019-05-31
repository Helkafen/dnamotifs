{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}
-- {-# OPTIONS_GHC -fplugin Foreign.Storable.Generic.Plugin #-} -- Is supposed to make Storable instances faster. Doesn't seem to make a difference

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
import qualified Language.C.Inline               as C
import           System.IO.Unsafe (unsafePerformIO)
import           Foreign.C.Types                 (CInt)

import Types


C.context (C.baseCtx <> C.vecCtx <> C.fptrCtx <> C.bsCtx)

C.include "<stdio.h>"

parseV ::  B.ByteString -> CInt -> STO.Vector CInt ->  STO.Vector CInt -> IO CInt
parseV line line_length l_ptr r_ptr = [C.block| int {
          int count_l = 0;
          int count_r = 0;
          int* l = $vec-ptr:(int *l_ptr);
          int* r = $vec-ptr:(int *r_ptr);
          char* in = &($bs-ptr:line[0]);
          int length = $(int line_length);
          int res = 0;
          for(int i,j = 0; i<length - 2; i=i+4, j++) {
            if(in[i] == '1')
              l[count_l++] = j+1;
            if(in[i+2] == '1')
              r[count_r++] = j+1;
          }
          l[count_l] = 0;
          r[count_r] = 0;
          return 0;
          } |]

fillVector :: Int -> B.ByteString -> (STO.Vector CInt, STO.Vector CInt)
fillVector size s = let l = STO.replicate (size+1) 0 :: STO.Vector CInt
                        r = STO.replicate (size+1) 0 :: STO.Vector CInt
                        res = unsafePerformIO (parseV s (fromIntegral $ B.length s) l r)
                    in if res == res
                         then (STO.map (\x -> x - 1) $ STO.takeWhile (>0) l, STO.map (\x -> x - 1) $ STO.takeWhile (>0) r)
                         else error ("The impossible happened!") -- To be honest, this call to parseV is weird


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


parseVcfContent :: [(Position0, Position0)] -> [B.ByteString] -> Either Error [Variant]
parseVcfContent regions vcfLines = case vcfLines of
    [] -> Left $ ParsingError "Empty vcf file"
    (header:rest) -> pure $ rights $ map (parseVariant sampleIdentifiers) (filterOrderedIntervals pos regions rest) -- TODO: exception for a left
        where sampleIdentifiers = sampleIdsInHeader header
              pos :: B.ByteString -> Position0
              pos = Position . (\p -> p - 1) . fromRight (error "Bad position in vcf") . parseInt . B.takeWhile (/= tab) . B.drop 1  . B.dropWhile (/= tab)
              parseInt = parseOnly decimal


chromosomeParser :: Parser Chromosome
chromosomeParser = choice [autosomeParser, xParser, yParser]
    where autosomeParser = do
            _ <- option "" (string "chr")
            d <- decimal
            guard $ d > 0 && d < 23
            pure (Chromosome $ T.pack $ show (d :: Integer))
          xParser = string "X" $> Chromosome "X"
          yParser = string "Y" $> Chromosome "Y"


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
    (genoL, genoR) <- fillVector (V.length sampleIdentifiers) <$> takeWhile1 (/=newline)
    _ <- option newline (word8 newline)
    return $ Variant chr pos name ref alt genoL genoR sampleIdentifiers
  where
    skipField = skipWhile (/= tab) >> skip (== tab)


parseVariant :: V.Vector SampleId -> B.ByteString -> Either Error (Variant)
parseVariant sampleIdentifiers s = 
    mapLeft (ParsingError . (\e -> (decodeUtf8With ignore s) <> " " <> T.pack e)) (parseOnly (variantParser sampleIdentifiers) s) 
