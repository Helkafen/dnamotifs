{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Types where

import           Foreign.C.Types                 (CChar)
import           Data.Text (Text)
import qualified Data.Vector as V
import           Data.Vector.Storable


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
    mSampleId :: Int,
    mMatched :: [Nucleotide]
} deriving (Eq, Show)

type Nucleotide = CChar
n, a, c, g, t :: Nucleotide
n = 0
a = 1
c = 2
g = 3
t = 4

newtype Chromosome = Chromosome Text
  deriving (Eq, Show)

newtype Position = Position Int
  deriving (Eq, Ord, Show, Num, Enum, Real, Integral, Storable)

data Genotype = Geno00 | Geno01 | Geno10 | Geno11
  deriving (Eq, Show)

-- The ID of a person
newtype SampleId = SampleId Text
  deriving (Eq, Ord, Show)

newtype BaseSequence = BaseSequence Text
  deriving (Eq, Show)

data Variant = Variant {
    chromosome :: Chromosome,
    position :: Position,
    variantId :: Maybe Text,
    reference :: BaseSequence,
    alternative :: BaseSequence,
    genotypes :: V.Vector Genotype,
    sampleIds :: V.Vector SampleId
} deriving (Eq, Show)

data Error = ParsingError Text | MissingFile FilePath | NoVariantFound FilePath
  deriving (Show)