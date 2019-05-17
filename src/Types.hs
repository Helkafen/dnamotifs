module Types where

import           Foreign.C.Types                 (CChar, CInt)

data Pweight = Pweight {
    wa :: Float,
    wc :: Float,
    wg :: Float,
    wt :: Float
} deriving (Show)

type Pattern = [Pweight]

data Match = Match {
    mPatternId :: Int, -- 0 based
    mScore :: Int,     -- 0 - 1000
    mPosition :: Int,  -- 0 based
    mSampleId :: Int
} deriving (Eq, Show)

type Nucleotide = CChar
n, a, c, g, t :: Nucleotide
n = 0
a = 1
c = 2
g = 3
t = 4

type Position = CInt