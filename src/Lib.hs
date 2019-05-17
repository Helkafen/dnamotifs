{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}

module Lib ( findPatternsInBlock, mkPatterns, mkNucleotideAndPositionBlock ) where

import           Foreign                         (Ptr, alloca, peek, FunPtr) 
import           Foreign.ForeignPtr              (newForeignPtr)
import           Foreign.C.Types                 (CInt, CChar)
import qualified Language.C.Inline               as C
import           Data.Vector.Storable            (Vector)
import qualified Data.Vector.Storable            as V
import           Data.Monoid                     ((<>))

import Types

newtype Patterns = Patterns (Vector CInt)

-- numberOfPeople, blockSize, nucleotides, positions
data NucleotideAndPositionBlock = NucleotideAndPositionBlock Int Int (Vector Nucleotide) (Vector Position)


vectorToMatches :: Vector CInt -> [Match]
vectorToMatches v = map toMatch (list4uple $ V.toList v)
    where toMatch (p,s,pos,sam) = Match (fromIntegral p) (fromIntegral s) (fromIntegral pos) (fromIntegral sam)

list4uple :: [a] -> [(a,a,a,a)]
list4uple (x1:x2:x3:x4:xs) = (x1,x2,x3,x4):(list4uple xs)
list4uple _ = []

patternToVector :: Pattern -> Vector CInt
patternToVector p = V.fromList [fromIntegral $ length p, sum (map (floor . (*1000) . m) p)] <> mconcat (map step p)
    where
    step x = V.map (floor . (*1000)) (V.fromList [0, wa x, wc x, wg x, wt x, m x]) -- 0 For "N"
    m x = maximum [wa x, wc x, wg x, wt x]

mkNucleotideAndPositionBlock :: [(Vector Nucleotide, Vector CInt)] -> NucleotideAndPositionBlock
mkNucleotideAndPositionBlock [] = NucleotideAndPositionBlock 0 0 V.empty V.empty
mkNucleotideAndPositionBlock xs = NucleotideAndPositionBlock numberOfPeople max_length (mconcat $ map (pad n . fst) xs) (mconcat $ map (pad 0 . snd) xs)  -- Pad with nucleotide N if sizes are different
    where numberOfPeople = length xs
          max_length = maximum (map (V.length . fst) xs)
          pad e v = v <> (V.fromList $ take (max_length - V.length v) (repeat e))


mkPatterns :: [Pattern] -> Patterns
mkPatterns patterns = Patterns $ mconcat (map patternToVector patterns) <> V.fromList [0]


C.context (C.baseCtx <> C.vecCtx <> C.fptrCtx)
C.include "<stdio.h>"
C.include "<math.h>"
C.include "<stdlib.h>"

foreign import ccall "&free" freePtr :: FunPtr (Ptr CInt-> IO ())

find_patterns :: CInt -> CInt -> Vector CChar -> Vector CInt -> Vector CInt -> Ptr CInt -> IO (Ptr CInt)
find_patterns sample_size block_size vec pos pat res_size = [C.block| int* {
    typedef struct Match {
        int score;
        int pattern;
        int sample;
        int position;
        struct Match* next;
    } Match;

    typedef struct MatchRecord {
        int score;
        int pattern;
        int sample;
        int position;
    } MatchRecord;

    int sample_size = $(int sample_size);
    int block_size = $(int block_size);
    int* patterns = $vec-ptr:(int *pat);
    char* in = $vec-ptr:(char *vec);
    int* positions = $vec-ptr:(int *pos);

    Match *match_head = (Match*)0;
    int match_number = 0;

    for(int sample = 0; sample < sample_size; sample++) {
        for(int position = 0; position < block_size; position++) {
            int i = sample*block_size + position;

            int k = 0;
            int j = 0;
            
            int pattern_id = 0;
            
            int pattern_length = -1;
            while (1) {
                pattern_length = patterns[k]; k = k + 1;
                int score_union = patterns[k]; k = k + 1;
                int score_inter = 0;
                
                // End of pattern list
                if (pattern_length <= 0) {
                    break;
                }
                                
                const int bound1 = ((sample + 1) * block_size) - i;
                const int max_index = bound1 < pattern_length ? bound1 : pattern_length;
            
                // Hotspot loop
                for (j = 0; j < max_index; k = k + 6, j++) {
                    score_inter += patterns[k + in[i + j]];
                }

                // Move to the next pattern anyway if j was restricted by the size of the block
                if(bound1 == max_index) {
                    k -= j*6;
                    k += pattern_length*6;
                }

                int score = (1000 * score_inter) / score_union;
                if(score >= 500) {
                    struct Match *match = malloc (sizeof (struct Match));
                    if (match == 0)
                        return ((int*)(-1));
                    match->score = score;
                    match->pattern = pattern_id;
                    match->sample = sample;
                    match->position = positions[i];
                    if(match_head == 0) {
                        match_head = match;
                        match->next = 0;
                    }
                    else {
                        match->next = match_head;
                        match_head = match;
                    }
                    match_number++;
                }
                pattern_id++;
            }
        }
    }

    int *matches = malloc (match_number * 4 * sizeof (int));
    Match* p = match_head;
    int n = 0;
    while(p) {
        printf("Match score=%d pattern_id=%d sample_id=%d position=%d \n", p->score, p->pattern, p->sample, p->position);
        matches[n*4+0] = p->pattern;
        matches[n*4+1] = p->score;
        matches[n*4+2] = p->position;
        matches[n*4+3] = p->sample;
        Match* to_be_freed = p;
        p = p->next;
        free(to_be_freed);
        n++;
    }
    printf("match_number in this block: %d\n", match_number);
    *($(int* res_size)) = match_number;

    return (int*)matches;
    } |]


findPatternsInBlock :: NucleotideAndPositionBlock -> Patterns -> IO [Match]
findPatternsInBlock (NucleotideAndPositionBlock numberOfPeople block_size inputData positionData) (Patterns patternData) = do
    result <- alloca $ \n_ptr -> do
        x <- find_patterns (fromIntegral numberOfPeople) (fromIntegral block_size) inputData positionData patternData (n_ptr :: Ptr CInt)
        result_size <- peek n_ptr
        fptr <- newForeignPtr freePtr x
        print ("result_size", result_size)
        return $ V.unsafeFromForeignPtr0 fptr (fromIntegral (4 * result_size))
    return $ vectorToMatches (result :: Vector CInt)