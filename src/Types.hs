{-# LANGUAGE GeneralizedNewtypeDeriving, DeriveGeneric #-}

module Types (
  Haplotype(..),
  Pweight(..),
  Pattern,
  Match(..),
  Nucleotide(..), n, a, c, g, t, toNuc,
  Genotype, geno00, geno01, geno10, geno11, -- Do not export the data constructor to forbid the creation of illegal values
  BaseSequence,
  BaseSequencePosition(..),
  Chromosome(..),
  Position0(..),
  SampleId(..),
  Variant(..),
  Diff(..),
  Error(..),
) where

import           Data.Text (Text)
import qualified Data.ByteString as B
import qualified Data.Vector as V
import qualified Data.Vector.Storable as STO
import           Foreign.C.Types                 (CInt)
import           Data.Word (Word8)
import           Foreign.Storable.Generic
import           GHC.Generics

data Haplotype = HaploLeft | HaploRight
    deriving (Eq, Show)

data Pweight = Pweight {
    wa :: Float,
    wc :: Float,
    wg :: Float,
    wt :: Float
} deriving (Eq, Show)

type Pattern = [Pweight]

data Match = Match {
    mPatternId :: !Int, -- 0 based
    mScore :: !Int,     -- 0 - 1000
    mPosition :: !Int,  -- 0 based
    mSampleId :: !Int,
    mMatched :: ![Nucleotide]
} deriving (Eq, Show)


-- They need to keep theses values, because the C function uses them as memory offsets in a table
newtype Nucleotide = Nucleotide { unNuc :: Word8 }
  deriving (Eq, Show, Generic)

instance GStorable Nucleotide

n, a, c, g, t :: Nucleotide
n = Nucleotide 0
a = Nucleotide 1
c = Nucleotide 2
g = Nucleotide 3
t = Nucleotide 4

-- ACGTacgtNn -> 01234
{-# INLINE toNuc #-}
toNuc :: Word8 -> Nucleotide
toNuc 65 = a
toNuc 67 = c
toNuc 71 = g
toNuc 84 = t
toNuc 78 = n
toNuc 97 = a
toNuc 99 = c
toNuc 103 = g
toNuc 116 = t
toNuc 110 = n
toNuc other = error $ "Bad nucleotide " <> show other

newtype Genotype = Genotype Word8
    deriving (Eq, Show, Generic)

instance GStorable Genotype

geno00, geno01, geno10, geno11 :: Genotype
geno00 = Genotype 10
geno01 = Genotype 11
geno10 = Genotype 12
geno11 = Genotype 13

type BaseSequence = B.ByteString

-- Base sequence + reference position
data BaseSequencePosition = BaseSequencePosition {-# UNPACK #-} !BaseSequence {-# UNPACK #-} !(STO.Vector CInt) -- TODO CInt16 for memory bandwidth
    deriving (Eq, Show)

newtype Chromosome = Chromosome { unChr :: Text }
  deriving (Eq, Show)


-- Zero based position
newtype Position0 = Position Int
  deriving (Eq, Ord, Show, Num, Enum, Real, Integral, STO.Storable)


  
-- The ID of a person
newtype SampleId = SampleId Text
  deriving (Eq, Ord, Show)


-- 0-based (not like in a VCF, which is 1-based)
data Variant = Variant {
    chromosome :: !Chromosome,
    position :: !Position0,
    variantId :: !(Maybe Text),
    reference :: !BaseSequence,
    alternative :: !BaseSequence,
    genotypes :: !(STO.Vector Genotype),
    sampleIds :: !(V.Vector SampleId)
} deriving (Eq, Show)

data Diff = Diff !Position0 !BaseSequence !BaseSequence -- pos, ref, alt
  deriving (Eq, Ord, Show)

data Error = ParsingError Text | MissingFile FilePath | NoVariantFound FilePath
  deriving (Eq, Show)