{-# LANGUAGE GeneralizedNewtypeDeriving, DeriveGeneric, FlexibleInstances, OverloadedStrings #-}

module Types (
  Haplotype(..),
  Pweight(..),
  Pattern(..),
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
  HaplotypeId(..),
  Count2(..)
) where

import           RIO
import qualified RIO.Text as T
import qualified RIO.ByteString as B
import qualified RIO.Vector.Boxed as V
import qualified RIO.Vector.Storable as STO

import           Foreign.C.Types                 (CInt)
import           Foreign.Storable.Generic
import           GHC.Generics ()
import           TextShow
import           Control.DeepSeq ()
import           Test.QuickCheck
import           Control.Applicative (liftA2)

data Haplotype = HaploLeft | HaploRight
    deriving (Eq, Ord, Show, Generic)

instance NFData Haplotype


data Pweight = Pweight {
    wa :: CInt, -- 1000 times the score in the PWM file
    wc :: CInt,
    wg :: CInt,
    wt :: CInt
} deriving (Eq, Show)

data Pattern = Pattern {
  patternMinScore :: Int,     -- 1000 times the float score
  patternWeights :: [Pweight]
} deriving (Eq, Show)

data Match a = Match {
    mPatternId :: !Int, -- 0 based
    mScore :: !Int,   -- 1000 times the float score
    mPosition :: !Int,  -- 0 based
    mSampleId :: !a,
    mMatched :: ![Nucleotide]
} deriving (Eq, Show)

instance Arbitrary Nucleotide where
  arbitrary = elements [a, c, g, t, n]

instance Arbitrary SampleId where
  arbitrary = elements [SampleId (T.pack "sample0") 0, SampleId (T.pack "sample1") 1, SampleId (T.pack "sample2") 2]

instance Arbitrary Haplotype where
  arbitrary = elements [HaploLeft, HaploRight]

instance Arbitrary HaplotypeId where
  arbitrary = liftA2 HaplotypeId arbitrary arbitrary

instance Arbitrary (Match [HaplotypeId]) where
  arbitrary = do
    patId <- elements [0..10]
    score <- (`mod` 5000) <$> arbitrarySizedNatural
    pos <- (`mod` 100) <$> arbitrarySizedNatural
    sampleId <- listOf1 arbitrary
    matched <- listOf1 arbitrary
    return $ Match patId score pos sampleId matched


-- They need to keep theses values, because the C function uses them as memory offsets in a table
newtype Nucleotide = Nucleotide { unNuc :: Word8 }
  deriving (Eq, Generic)

instance GStorable Nucleotide

instance TextShow Nucleotide where
  showb (Nucleotide 0) = showb ("N" :: T.Text)
  showb (Nucleotide 1) = showb ("A" :: T.Text)
  showb (Nucleotide 2) = showb ("C" :: T.Text)
  showb (Nucleotide 3) = showb ("G" :: T.Text)
  showb (Nucleotide 4) = showb ("T" :: T.Text)
  showb (Nucleotide _) = showb ("_" :: T.Text)

n, a, c, g, t :: Nucleotide
n = Nucleotide 0
a = Nucleotide 1
c = Nucleotide 2
g = Nucleotide 3
t = Nucleotide 4

instance Show Nucleotide where
  show x | x == n = "N"
         | x == a = "A"
         | x == c = "C"
         | x == g = "G"
         | x == t = "T"
         | otherwise = error "Wrong nucleotide code in show instance"

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
    deriving (Eq, Ord, Generic)

instance GStorable Genotype

geno00, geno01, geno10, geno11 :: Genotype
geno00 = Genotype 10
geno01 = Genotype 11
geno10 = Genotype 12
geno11 = Genotype 13

instance Show Genotype where
  show x@(Genotype ge)
    | x == geno00 = "geno00"
    | x == geno01 = "geno01"
    | x == geno10 = "geno10"
    | x == geno11 = "geno11"
    | otherwise = error "bad genotype " <> show ge

instance TextShow Genotype where
  showb (Genotype 10) = fromText ("0|0" :: T.Text)
  showb (Genotype 11) = fromText ("0|1" :: T.Text)
  showb (Genotype 12) = fromText ("1|0" :: T.Text)
  showb (Genotype 13) = fromText ("1|1" :: T.Text)
  showb (Genotype _) = error "Bad genotype"

type BaseSequence = B.ByteString

-- Base sequence + reference position
data BaseSequencePosition = BaseSequencePosition {-# UNPACK #-} !BaseSequence {-# UNPACK #-} !(STO.Vector CInt) -- TODO CInt16 for memory bandwidth
    deriving (Eq, Ord)

instance Show BaseSequencePosition where
  show (BaseSequencePosition nuc pos) = show (map Nucleotide (B.unpack nuc), STO.toList pos)

newtype Chromosome = Chromosome { unChr :: T.Text }
  deriving (Eq, Show, Generic)

instance NFData Chromosome


-- Zero based position
newtype Position0 = Position Int
  deriving (Eq, Ord, Show, Num, Enum, Real, Integral, STO.Storable, Generic)

instance NFData Position0

instance TextShow Position0 where
  showb (Position p) = showb p

-- The ID of a person
-- Int is the position in the samples list
data SampleId = SampleId T.Text Int
  deriving (Eq, Ord, Show, Generic)

instance NFData SampleId


instance TextShow SampleId where
  showb (SampleId s _) = showb s


-- 0-based (not like in a VCF, which is 1-based)
data Variant = Variant {
    chromosome :: !Chromosome,
    position :: !Position0,
    variantId :: !(Maybe T.Text),
    reference :: !BaseSequence,
    alternative :: !BaseSequence,
    genotypesL :: !(STO.Vector Int),
    genotypesR :: !(STO.Vector Int),
    sampleIds :: !(V.Vector SampleId)
} deriving (Eq, Show, Generic)

instance NFData Variant

showBaseSequence :: BaseSequence -> String
showBaseSequence = map tr . B.unpack
  where tr 0 = 'N'
        tr 1 = 'A'
        tr 2 = 'C'
        tr 3 = 'G'
        tr 4 = 'T'
        tr _ = '_'

data Diff = Diff !Position0 !BaseSequence !BaseSequence -- pos, ref, alt
  deriving (Eq, Ord, Generic) -- Keep Position in the first place. We need it for the default ordering

instance NFData Diff


instance Show Diff where
  show (Diff (Position p) ref alt) = "Diff " <> show p <> ": " <> showBaseSequence ref <> "->" <> showBaseSequence alt

data Error = ParsingError T.Text | MissingFile FilePath | NoVariantFound FilePath | VcfLoadingError String | FastaLoadingError String | BedLoadingError String | MotifLoadingError String
  deriving (Eq, Show, Generic)

instance Exception Error

instance NFData Error

data HaplotypeId = HaplotypeId !SampleId !Haplotype
  deriving (Eq, Ord, Show)


data Count2 = Count2 Int Int
    deriving (Eq, Show)

instance Semigroup Count2 where
    (Count2 x y) <> (Count2 z v) = Count2 (x+z) (y+v)

instance Monoid Count2 where
    mempty = Count2 0 0
    mappend = (<>)