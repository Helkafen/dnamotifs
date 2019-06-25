module Haplotype (buildAllHaplotypes, applyVariants, variantsToDiffs) where

import qualified Data.Vector.Storable as STO
import qualified Data.ByteString as B
import qualified Data.Vector as V
import qualified Data.Map as M
import qualified Data.Set as Set
import           Data.Vector.Algorithms.Merge (sort)

import Import
import Range

-- Inclusive [Start, End] interval
applyVariants :: (Range Position0 -> BaseSequencePosition) -> Range Position0 -> [Diff] -> BaseSequencePosition
applyVariants takeReferenceGenome (Range start end) allDiffs =
    let chunks = go start (filter (\(Diff p _ _) -> p >= start && p <= end) allDiffs)
    in BaseSequencePosition (B.concat (map seqOf chunks)) (STO.concat (map posOf chunks))
  where go :: Position0 -> [Diff] -> [BaseSequencePosition]
        go refPosition [] = [takeReferenceGenome (Range refPosition end)]
        go refPosition diffs@(Diff pos ref alt:vs)
            | pos > refPosition = takeReferenceGenome (Range refPosition (pos-1)): go pos diffs
            | pos == refPosition && B.length ref == 1 = -- SNV, insertion
                BaseSequencePosition alt (STO.replicate (B.length alt) (fromIntegral refPosition)): go (refPosition + 1) vs
            | pos == refPosition && B.length alt == 1 = -- deletion
                BaseSequencePosition alt (STO.singleton (fromIntegral refPosition)) : go (offSet refPosition (B.length ref)) vs
            | pos == refPosition = error "Missing case in applyVariants (we assume that ref or alt have length 1)"
            | refPosition >= end = [takeReferenceGenome (Range refPosition end)]
            | otherwise = [BaseSequencePosition B.empty STO.empty]
        seqOf (BaseSequencePosition nuc _) = nuc
        posOf (BaseSequencePosition _ pos) = pos
        offSet (Position x) o = Position (x + o)


combineBySecond :: Ord b => [(a,b)] -> M.Map b [a]
combineBySecond = M.fromListWith (<>) . map (\(x,y) -> (y,[x]))

variantsToDiffs :: [Variant] -> M.Map (V.Vector Diff) [HaplotypeId]
variantsToDiffs [] = M.empty
variantsToDiffs variants@(v:_) = combineBySecond (M.toList haplotypeToDiff)
    where allDiffs = V.concatMap variantToDiffs (V.fromList variants) :: V.Vector (Int, Haplotype, Diff)
          haplotypeToDiff = M.mapKeys toHaplotypeId $ V.modify sort <$> M.fromListWith (<>) (V.toList $ V.map (\(i, h, d) -> ((i,h), V.singleton d)) allDiffs) :: M.Map HaplotypeId (V.Vector Diff)
          toHaplotypeId (i,h) = HaplotypeId (sampleIds v V.! i) h

variantToDiffs :: Variant -> V.Vector (Int, Haplotype, Diff)
variantToDiffs v = V.map (\i -> (i, HaploLeft, d)) (STO.convert (genotypesL v)) <> V.map (\i -> (i, HaploRight, d)) (STO.convert (genotypesR v))
    where d = Diff (position v) (reference v) (alternative v)


buildAllHaplotypes :: (Range Position0 -> BaseSequencePosition) -> Range Position0 -> [SampleId] -> [Variant] -> V.Vector (BaseSequencePosition, [HaplotypeId])
buildAllHaplotypes takeReferenceGenome peak sampleIdList variants = V.fromList $ (takeReferenceGenome peak, Set.toList samplesWithNoVariant) : haplotypes
    where haplotypes = M.toList $ M.mapKeysWith (<>) (applyVariants takeReferenceGenome peak . V.toList) (variantsToDiffs variants) :: [(BaseSequencePosition, [HaplotypeId])]
          samplesWithNoVariant = Set.difference allHaplotypeIds (Set.fromList $ concatMap snd haplotypes) :: Set.Set HaplotypeId
          allHaplotypeIds = Set.fromList [HaplotypeId x y | x <- sampleIdList, y <- [HaploLeft, HaploRight]]
