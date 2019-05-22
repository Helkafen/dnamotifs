module Lib where

import qualified Data.DList as DList
import qualified Data.Vector.Storable as V
import qualified Data.Vector.Unboxed as U
--import qualified Data.Text as T
--
import Types
--import Fasta (loadFasta)
--import Vcf (readVcfWithGenotypes)
--import PatternFind (findPatternsInBlock, mkPatterns, mkNucleotideAndPositionBlock)

-- Inclusive [Start, End] interval
applyVariants :: V.Vector Nucleotide -> Position -> Position -> [Diff] -> [(Nucleotide, Position)]
applyVariants referenceGenome (Position start) (Position end) variants = 
    DList.toList $ go start (filter (\(Diff (Position p) _ _) -> p >= start && p <= end) variants)
  where go refPosition [] = takeRef refPosition end
        go refPosition diffs@(Diff (Position pos) ref alt:vs)
            | pos > refPosition = takeRef refPosition (pos-1) <> go pos diffs
            | pos == refPosition && U.length ref == 1 = DList.fromList (zip (U.toList alt) (repeat (Position refPosition))) <> go (refPosition + 1) vs -- SNV, insertion
            | pos == refPosition && U.length alt == 1 = DList.fromList (zip (U.toList alt) (repeat (Position refPosition))) <> go (refPosition + U.length ref) vs -- deletion
            | pos == refPosition = error "Missing case in applyVariants (we assume that ref or alt have length 1)"
            | refPosition >= end = takeRef refPosition end
            | otherwise = DList.empty
        takeRef :: Int -> Int -> DList.DList (Nucleotide, Position)
        takeRef s e = DList.fromList $ zip (V.toList $ V.take (e-s+1) (V.drop s referenceGenome)) (map Position [s..e+1])

--findPatterns :: Chromosome -> FilePath -> FilePath -> FilePath -> IO Bool
--findPatterns chromosome referenceGenomeFile vcfFile resultFile = do
--    referenceGenome <- loadFasta chromosome referenceGenomeFile
--    vcf <- readVcfWithGenotypes vcfFile
--    case vcf of
--        Left error -> print error >> return False
--        Right [] -> print "No variant loaded" >> return False
--        Right variants@(x:xs) -> 