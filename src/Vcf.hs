{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}
-- {-# OPTIONS_GHC -fplugin Foreign.Storable.Generic.Plugin #-} -- Is supposed to make Storable instances faster. Doesn't seem to make a difference

module Vcf (readVcfWithGenotypes, parseVariant, filterOrderedIntervals, parseVcfContent, fillVector) where

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
import           Data.Word (Word8)
import qualified Language.C.Inline               as C
import           System.IO.Unsafe (unsafePerformIO)
import           Foreign.C.Types                 (CInt)
import           Foreign                         (Ptr, FunPtr)
import           Foreign.ForeignPtr              (newForeignPtr)
import           Data.Range.Range                (Range(..))
import           Control.Monad.Trans.Except
import           Control.Monad.Except (throwError, lift)

import Types

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
                pure (STO.map fromIntegral left, STO.map fromIntegral right)
        len = fromIntegral (B.length s)


readGzippedLines :: FilePath -> IO [B.ByteString]
readGzippedLines path = map BL.toStrict . BLC.lines . GZip.decompress <$> BL.readFile path

-- xs needs to be ordered in ascending order
filterOrderedIntervals :: Ord a => (b->a) -> [Range a] -> [b] -> [(Range a, [b])]
filterOrderedIntervals pos ranges xs = go ranges xs
    where go [] _ = []
          go (r@(SpanRange start end):rs) l = 
            let (vs, rest) =  span ((<=end) . pos) $ dropWhile ((<start) . pos) l
            in (r, vs):go rs rest
          go ((SingletonRange x):rs) l = go ((SpanRange x x):rs) l
          go (r@(LowerBoundRange start):_) l = [(r, dropWhile ((<start) . pos) l)]
          go (r@(UpperBoundRange end):_) l = [(r, Prelude.takeWhile ((<=end) . pos) l)]
          go (r@InfiniteRange:_) l = [(r, l)]


readVcfWithGenotypes :: FilePath -> [Range Position0] -> ExceptT Error IO (V.Vector SampleId, [(Range Position0, [Variant])])
readVcfWithGenotypes path regions = do
    lift $ putStrLn ("Loading VCF file " <> path)
    vcf <- lift $ (parseVcfContent regions . dropWhile ("##" `B.isPrefixOf`)) <$> readGzippedLines path
    case vcf of
      Left err -> throwError err
      Right [] -> throwError (VcfLoadingError ("No variant loaded"))
      Right ((_, []):_) -> throwError (VcfLoadingError "No variant loaded in the first peak")
      Right variants@((_,x:_):_) -> pure (sampleIds x, variants)


sampleIdsInHeader :: B.ByteString -> V.Vector SampleId
sampleIdsInHeader header = V.fromList $ map (SampleId . decodeUtf8With ignore) $ Prelude.takeWhile ((>0) . B.length) $ drop 9 (B.split (fromIntegral $ ord '\t') header)


parseVcfContent :: [Range Position0] -> [B.ByteString] -> Either Error [(Range Position0, [Variant])]
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


variantIdParser :: Parser (Maybe Text)
variantIdParser = choice [none, rs] <?> "variant_name"
    where rs = (Just . decodeUtf8With ignore . B.pack) <$> many1 (notWord8 tab)
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
    mapLeft (ParsingError . (\e -> decodeUtf8With ignore s <> " " <> T.pack e)) (parseOnly (variantParser sampleIdentifiers) s) 
