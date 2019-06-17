{-# LANGUAGE OverloadedStrings #-}

module Fasta (loadFasta, takeRef) where

import qualified Data.Vector.Storable as STO
import qualified Data.ByteString      as B
import qualified Data.ByteString.Lazy as BL
import qualified Data.ByteString.Lazy.Char8 as BLC
import qualified Data.Text as T
import           Data.Char (ord)
import           Data.Word (Word8)
import           Data.Range.Range (Range(..))

import Types


-- Hand tested (not anymore)
loadFasta :: Chromosome -> FilePath -> IO (Range Position0 -> BaseSequencePosition, Int)
loadFasta (Chromosome chr) filename = 
    do putStr ("Load reference genome " <> filename <> " ... ")
       -- BL.head is safe because we filtered out the empty lines earlier
       alphaBaseSequence <- (BL.toStrict . BL.concat . takeWhile (\l -> BL.head l /= separator) . tail . dropWhile (/=wantedHeader) . filter (not . BL.null) . BLC.lines) <$> BL.readFile filename
       if B.length alphaBaseSequence == 0
         then putStrLn "Error (sequence is empty)"
         else putStrLn "Done"
       pure (takeRef alphaBaseSequence, B.length alphaBaseSequence)
    where separator = fromIntegral (ord '>') :: Word8
          wantedHeader = BLC.pack (">chr"<>T.unpack chr) :: BL.ByteString

takeRef :: B.ByteString -> Range Position0 -> BaseSequencePosition
takeRef referenceGenome (SpanRange (Position s) (Position e)) = BaseSequencePosition bases positions
    where end = minimum [e, B.length referenceGenome - 1]
          bases = B.map (unNuc . toNuc) $ B.take (end-s+1) (B.drop s referenceGenome)
          positions = STO.fromListN (end-s+1) (map fromIntegral [s..end+1])
takeRef referenceGenome (SingletonRange x) = takeRef referenceGenome (SpanRange x x)
takeRef referenceGenome (LowerBoundRange s) = takeRef referenceGenome (SpanRange s (Position $ B.length referenceGenome - 1))
takeRef referenceGenome (UpperBoundRange e) = takeRef referenceGenome (SpanRange (Position 0) e)
takeRef referenceGenome InfiniteRange = takeRef referenceGenome (SpanRange (Position 0) (Position $ B.length referenceGenome - 1))

