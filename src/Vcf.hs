{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}
-- {-# OPTIONS_GHC -fplugin Foreign.Storable.Generic.Plugin #-} -- Is supposed to make Storable instances faster. Doesn't seem to make a difference

module Vcf (readVcfWithGenotypes, parseVariant, filterOrderedIntervals, parseVcfContent, fillVector) where

import qualified RIO.Text as T

import qualified Codec.Compression.GZip as GZip
import qualified RIO.ByteString as B
import qualified RIO.ByteString.Lazy as BL
import qualified Data.ByteString.Lazy.Char8 as BLC

import qualified RIO.Vector.Storable as STO
import qualified RIO.Vector.Storable.Partial as STO (head)
import qualified RIO.Vector.Storable.Unsafe as STO (unsafeFromForeignPtr0)

import           Data.Attoparsec.ByteString
import           Data.Attoparsec.ByteString.Char8 (decimal, letter_ascii)

import qualified Data.Vector as V

import qualified Language.C.Inline               as C
import           System.IO.Unsafe (unsafePerformIO)
import           Foreign.C.Types                 (CInt)
import           Foreign                         (Ptr, FunPtr)
import           Foreign.ForeignPtr              (newForeignPtr)


import Import
import Range

C.context (C.baseCtx <> C.vecCtx <> C.fptrCtx <> C.bsCtx)

C.include "<stdio.h>"
C.include "<stdlib.h>"

foreign import ccall "&free" freePtr :: FunPtr (Ptr CInt-> IO ())

parseV ::  B.ByteString -> CInt -> IO (Ptr CInt)
parseV line line_length = [C.block| int* {
          char* in = &($bs-ptr:line[0]);
          int length = $(int line_length);

          int* resultL = (int*)malloc(length * sizeof(int));
          int* resultR = (int*)malloc(length * sizeof(int));
          int count_l = 0;
          int count_r = 0;

          for(int i,j = 0; i<length - 2; i=i+4, j++) {
            if(in[i] == '1')
              resultL[count_l++] = j;
            if(in[i+2] == '1')
              resultR[count_r++] = j;
          }
          int* result = (int*)malloc((length+1) * sizeof(int) * 2);

          result[0] = count_l;
          result[1] = count_r;
          int i = 2;
          for(int j = 0;j<count_l;i++,j++) {
            result[i] = resultL[j];
          }

          for(int j = 0;j<count_r;i++,j++) {
            result[i] = resultR[j];
          }

          free(resultL);
          free(resultR);

          return result;
          } |]

fillVector :: B.ByteString -> (STO.Vector Int, STO.Vector Int)
fillVector s = unsafePerformIO go
  where go = do x <- parseV s len
                fptr <- newForeignPtr freePtr x
                let vec = STO.map fromIntegral $ STO.unsafeFromForeignPtr0 fptr (fromIntegral $ 2 + 2 * (len + 1))
                let lenL = STO.head vec
                let lenR = STO.head (STO.drop 1 vec)
                let left = STO.take lenL (STO.drop 2 vec)
                let right = STO.take lenR (STO.drop (2+lenL) vec)
                pure (left, right)
        len = fromIntegral (B.length s)


readGzippedLines :: FilePath -> RIO env [B.ByteString]
readGzippedLines path = map BL.toStrict . BLC.lines . GZip.decompress <$> BL.readFile path

filterOrderedIntervals :: Ord a => (b->a) -> Ranges a -> [b] -> [(Range a, [b])]
filterOrderedIntervals pos ranges xs = go (getRanges ranges) xs
    where go [] _ = []
          go (r@(Range start end):rs) l = 
            let (vs, rest) =  span ((<=end) . pos) $ dropWhile ((<start) . pos) l
            in (r, vs):go rs rest


readVcfWithGenotypes :: HasLogFunc env => FilePath -> Ranges Position0 -> RIO env ([SampleId], [(Range Position0, [Variant])])
readVcfWithGenotypes path regions = do
    logInfo (display $ T.pack ("Loading VCF file " <> path))
    vcf <- (parseVcfContent regions . dropWhile ("##" `B.isPrefixOf`)) <$> readGzippedLines path
    case vcf of
      Left err -> throwM err
      Right [] -> throwM (VcfLoadingError ("No variant loaded"))
      Right ((_, []):_) -> throwM (VcfLoadingError "No variant loaded in the first peak")
      Right variants@((_,x:_):_) -> pure (V.toList (sampleIds x), variants)


sampleIdsInHeader :: B.ByteString -> V.Vector SampleId
sampleIdsInHeader header = V.fromList $ map (\(i,te) -> SampleId (decodeUtf8Lenient te) i) $ zip [0..] $ Import.takeWhile ((>0) . B.length) $ drop 9 (B.split (fromIntegral $ ord '\t') header)


parseVcfContent :: Ranges Position0 -> [B.ByteString] -> Either Error [(Range Position0, [Variant])]
parseVcfContent regions vcfLines = case vcfLines of
    [] -> Left $ ParsingError "Empty vcf file"
    (header:rest) -> pure $ map (\(r,variants) -> (r, rights $ map (parseVariant sampleIdentifiers) variants)) (filterOrderedIntervals pos regions rest) -- TODO: exception for a left
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


variantIdParser :: Parser (Maybe T.Text)
variantIdParser = choice [none, rs] <?> "variant_name"
    where rs = (Just . decodeUtf8Lenient . B.pack) <$> many1 (notWord8 tab)
          none = word8 dot $> Nothing

dot, tab, newline :: Word8
dot = fromIntegral (ord '.')
tab = fromIntegral (ord '\t')
newline = fromIntegral (ord '\n')

variantParser :: V.Vector SampleId -> Parser Variant
variantParser sampleIdentifiers = do
    chr <- chromosomeParser <* word8 tab
    pos <- (Position . (\x -> x - 1)) <$> decimal <* word8 tab
    name <- variantIdParser <* word8 tab
    ref <- (B.pack . map (unNuc . toNuc . fromIntegral . ord)) <$> many1 letter_ascii <* word8 tab
    alt <- (B.pack . map (unNuc . toNuc . fromIntegral . ord))  <$> many1 letter_ascii <* word8 tab
    skipField >> skipField >> skipField >> skipField
    (genoL, genoR) <- fillVector <$> takeWhile1 (/=newline)
    _ <- option newline (word8 newline)
    return $ Variant chr pos name ref alt genoL genoR sampleIdentifiers
  where
    skipField = skipWhile (/= tab) >> skip (== tab)


parseVariant :: V.Vector SampleId -> B.ByteString -> Either Error Variant
parseVariant sampleIdentifiers s = 
    Import.mapLeft (ParsingError . (\e -> decodeUtf8Lenient s <> " " <> T.pack e)) (parseOnly (variantParser sampleIdentifiers) s) 
