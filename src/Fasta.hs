{-# LANGUAGE OverloadedStrings #-}

module Fasta (loadFasta, takeRef) where

import qualified Data.Vector.Storable as STO
import qualified Data.ByteString      as B
import qualified Data.ByteString.Lazy as BL
import           Data.Text.Encoding (decodeUtf8)
import           Data.Char (ord)
import           Data.Word (Word8)

import Types


-- Hand tested (not anymore)
loadFasta :: Chromosome -> FilePath -> IO (Int -> Int -> BaseSequencePosition, Int)
loadFasta (Chromosome chr) filename = do putStrLn ("Load reference genome " <> filename)
                                         alphaBaseSequence <- (BL.toStrict . transform . go) <$> BL.readFile filename
                                         pure (takeRef alphaBaseSequence, B.length alphaBaseSequence)
    where separator = fromIntegral (ord '>') :: Word8
          newLine = fromIntegral (ord '\n') :: Word8
          go content = let -- The '>chrXXX' line, and everything after
                          (nextStartName, restOfFile) = BL.break (== newLine) (BL.dropWhile (/= separator) content)
                       in if decodeUtf8 (BL.toStrict nextStartName) == (">chr" <> chr)
                            then BL.takeWhile (/= separator) (BL.tail restOfFile)
                            else go restOfFile
          transform = BL.filter (/= newLine)

takeRef :: B.ByteString -> Int -> Int -> BaseSequencePosition
takeRef referenceGenome s e = BaseSequencePosition bases positions
    where end = minimum [e, B.length referenceGenome - 1]
          bases = B.map (unNuc . toNuc) $ B.take (end-s+1) (B.drop s referenceGenome)
          positions = STO.fromListN (end-s+1) (map fromIntegral [s..end+1])

