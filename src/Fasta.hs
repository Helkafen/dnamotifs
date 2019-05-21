{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE BangPatterns #-}

module Fasta (loadFasta) where

import qualified Data.Text.Lazy.IO as TIO (readFile)
import qualified Data.Text.Lazy as TL
import qualified Data.Text as T

import Types


-- Hand tested
loadFasta :: FilePath -> Chromosome -> IO (T.Text)
loadFasta fileName (Chromosome chr) = do
  chrLines <- (takeWhile (not . (">" `TL.isPrefixOf`)) . tail . dropWhile (/= (header)) . TL.lines) <$> TIO.readFile fileName :: IO ([TL.Text])
  let !content = T.toUpper (TL.toStrict $ TL.concat (chrLines :: [TL.Text]))
  pure content
  where header = (TL.pack ">chr") <> (TL.fromStrict chr)
