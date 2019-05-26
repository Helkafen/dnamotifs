{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Types where

import           Data.Text (Text)
import qualified Data.ByteString as B
import qualified Data.Vector as V
import qualified Data.Vector.Storable as STO
import           Foreign.C.Types                 (CInt)
import           Data.Word (Word8)

data Haplotype = HaploLeft | HaploRight
    deriving (Eq)

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

-- They need to keep theses values, because the C function uses them as memory offsets in a table
type Nucleotide = Word8
n, a, c, g, t :: Nucleotide
n = 0
a = 1
c = 2
g = 3
t = 4

type BaseSequence = B.ByteString

-- Base sequence + reference position
data BaseSequencePosition = BaseSequencePosition BaseSequence (STO.Vector CInt) -- TODO CInt16 for memory bandwidth
    deriving (Eq, Show)

newtype Chromosome = Chromosome { unChr :: Text }
  deriving (Eq, Show)


data ZeroBased
data OneBased
newtype Position a = Position Int
  deriving (Eq, Ord, Show, Num, Enum, Real, Integral, STO.Storable)

data Genotype = Geno00 | Geno01 | Geno10 | Geno11
  deriving (Eq, Show)
  
-- The ID of a person
newtype SampleId = SampleId Text
  deriving (Eq, Ord, Show)


-- 0-based (not like in a VCF, which is 1-based)
data Variant a = Variant {
    chromosome :: Chromosome,
    position :: Position a,
    variantId :: Maybe Text,
    reference :: BaseSequence,
    alternative :: BaseSequence,
    genotypes :: V.Vector Genotype,
    sampleIds :: V.Vector SampleId
} deriving (Eq, Show)

data Diff a = Diff (Position a) BaseSequence BaseSequence -- pos, ref, alt
  deriving (Eq, Ord, Show)

data Error = ParsingError Text | MissingFile FilePath | NoVariantFound FilePath
  deriving (Eq, Show)