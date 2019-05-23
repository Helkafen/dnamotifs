{-# LANGUAGE OverloadedStrings #-}

module Fasta (loadFasta) where

import qualified Data.Vector.Storable as V
import qualified Data.ByteString.Lazy as B
import           Data.Text.Encoding (decodeUtf8)
import           Data.Char (ord)
import           Data.Word (Word8)
import           Data.Vector.Storable.ByteString (byteStringToVector)

import Types


-- Hand tested
loadFasta :: Chromosome -> FilePath -> IO (V.Vector Nucleotide)
loadFasta (Chromosome chr) filename = (byteStringToVector . B.toStrict . transform . go) <$> B.readFile filename
 where separator = fromIntegral (ord '>') :: Word8
       newLine = fromIntegral (ord '\n') :: Word8
       go content = let -- The '>chrXXX' line, and everything after
                        (nextStartName, restOfFile) = B.break (== newLine) (B.dropWhile (/= separator) content) 
                    in if decodeUtf8 (B.toStrict nextStartName) == (">chr" <> chr)
                         then B.takeWhile (/= separator) (B.tail restOfFile)
                         else go restOfFile
       transform = B.filter (/= newLine)

