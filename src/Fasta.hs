{-# LANGUAGE OverloadedStrings #-}

module Fasta (loadFasta, takeRef) where

import qualified RIO.Vector.Storable as STO
import qualified RIO.ByteString      as B
import qualified RIO.ByteString.Lazy as BL
import qualified RIO.ByteString.Lazy.Partial as BLP
import qualified Data.ByteString.Lazy.Char8 as BLC
import qualified RIO.Text as T

import Import
import Range
import RIO.List.Partial (tail)


-- Hand tested (not anymore)
loadFasta :: HasLogFunc env => Chromosome -> FilePath -> RIO env (Range Position0 -> BaseSequencePosition, Int)
loadFasta (Chromosome chr) filename = 
    do -- BL.head is safe because we filtered out the empty lines earlier
       alphaBaseSequence <- (BL.toStrict . BL.concat . takeWhile (\l -> BLP.head l /= separator) . tail . dropWhile (/=wantedHeader) . filter (not . BL.null) . BLC.lines) <$> BL.readFile filename
       if B.length alphaBaseSequence == 0
         then throwM (FastaLoadingError $ "Load reference genome " <> filename <> ": sequence is empty")
         else logInfo (display $ "Load reference genome " <> (T.pack filename))
       pure (takeRef alphaBaseSequence, B.length alphaBaseSequence)
    where separator = fromIntegral (ord '>') :: Word8
          wantedHeader = BLC.pack (">chr"<>T.unpack chr) :: BL.ByteString

takeRef :: B.ByteString -> Range Position0 -> BaseSequencePosition
takeRef referenceGenome (Range (Position s) (Position e)) = BaseSequencePosition bases positions
    where end = min e (B.length referenceGenome - 1)
          bases = B.map (unNuc . toNuc) $ B.take (end-s+1) (B.drop s referenceGenome)
          positions = STO.fromListN (end-s+1) (map fromIntegral [s..end+1])
