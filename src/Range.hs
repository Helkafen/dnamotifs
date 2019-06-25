module Range (
    Range(..),
    Ranges(..),
    mergeRanges,
    inRange,
    inRanges,
    overlaps,
    mkRanges,
    rangesLength
) where

import RIO
import RIO.List (sort)

data Range a = Range a a
    deriving (Eq, Ord, Show)

-- In ascending order
data Ranges a = Ranges {getRanges :: [Range a]}
    deriving (Eq, Show)

mkRanges :: Ord a => [Range a] -> Maybe (Ranges a)
mkRanges xs = if (isSortedAndDisjoint xs) then Just (Ranges sorted) else Nothing
    where sorted = sort xs

rangesLength :: Ranges a -> Int
rangesLength (Ranges rs) = length rs

isSortedAndDisjoint :: Ord a => [Range a] -> Bool
isSortedAndDisjoint [] = True
isSortedAndDisjoint [_] = True
isSortedAndDisjoint ((Range s e):(Range s2 e2):xs) = s <= e && s < s2 && e <= s2 && isSortedAndDisjoint ((Range s2 e2):xs)

mergeRanges :: Ord a => Ranges a -> Ranges a
mergeRanges (Ranges ranges) = Ranges (go sortedRanges)
    where sortedRanges = sort ranges
          go [] = []
          go [x] = [x]
          go (x:y:xs) = if (overlap x y) then go ((merge1 x y):xs) else x:go (y:xs)
          merge1 (Range x y) (Range z v) = Range (min x z) (max y v)


instance (Show a, Ord a) => Semigroup (Ranges a) where
    (<>) (Ranges r1) (Ranges r2) = mergeRanges (Ranges (r1 <> r2))

instance (Show a, Ord a) => Monoid (Ranges a) where
    mempty = Ranges []

--intersection1 :: Ord a => Range a -> Range a -> [Range a]
--intersection1 a@(Range x y) b@(Range z v) =
--    if (inRange a z || inRange a v)
--        then [Range (min x z) (max y v)]
--        else [a, b]

inRange :: Ord a => Range a -> a -> Bool
inRange (Range s e) x = s <= x && x <= e

inRanges :: Ord a => [Range a] -> a -> Bool
inRanges ranges x = any (\r -> inRange r x) ranges

overlap :: Ord a => Range a -> Range a -> Bool
overlap a (Range z v) = inRange a z || inRange a v

overlaps :: Ord a => Range a -> [Range a] -> [Range a]
overlaps r xs = takeWhile (overlap r) dropped
    where dropped = dropWhile (not . overlap r) xs