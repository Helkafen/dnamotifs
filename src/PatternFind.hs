{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}

module PatternFind ( findPatternsInBlock, mkPatterns, mkNucleotideAndPositionBlock, NucleotideAndPositionBlock, Patterns, blockInfo ) where

import           Foreign                         (Ptr, alloca, peek, FunPtr) 
import           Foreign.ForeignPtr              (newForeignPtr)
import           Foreign.C.Types                 (CInt, CChar)
import qualified Language.C.Inline               as C
import           Data.Vector.Storable            (Vector)
import qualified Data.Vector.Storable            as V
import           Data.Monoid                     ((<>))
import           Data.Maybe                      (mapMaybe)
--import           Data.Word                       (Word8)

import Types

newtype Patterns = Patterns (Vector CInt)

-- numberOfPeople, blockSize, nucleotides, positions
data NucleotideAndPositionBlock = NucleotideAndPositionBlock Int Int (Vector CChar) (Vector CInt)

blockInfo :: NucleotideAndPositionBlock -> String
blockInfo (NucleotideAndPositionBlock numberOfPeople blockSize nucleotides positions) = "block " <> (show numberOfPeople) <> " " <> show blockSize <> " " <> show (V.minimum positions) <> " " <> show (V.maximum positions)

vectorToMatches :: Vector CInt -> [Match]
vectorToMatches v = mapMaybe toMatch (list4uple $ V.toList v)
    where toMatch (p:s:pos:sam:matchedSequence) = Just $ Match (fromIntegral p) (fromIntegral s) (fromIntegral pos) (fromIntegral sam) (trimMatchedSequence $ map fromIntegral matchedSequence)
          toMatch [] = Nothing
          toMatch x = error ("Wrong size in Match vector: " ++ show (length x))
          trimMatchedSequence = reverse . dropWhile (==n) . reverse

list4uple :: [a] -> [[a]]
list4uple [] = []
list4uple xs = let (x,rest) = splitAt 34 xs in x:list4uple rest


patternToVector :: Pattern -> Vector CInt
patternToVector p = V.fromList [fromIntegral $ length p, sum (map (floor . (*1000) . m) p)] <> mconcat (map step p)
    where
    step x = V.map (floor . (*1000)) (V.fromList [0, wa x, wc x, wg x, wt x, m x]) -- 0 For "N"
    m x = maximum [wa x, wc x, wg x, wt x]

pad :: V.Storable a => a -> Int -> Vector a -> Vector a
pad e len v = v <> padding
    where padding = V.fromList $ replicate (len - V.length v) e

-- Pad with nucleotide N if sizes are different
mkNucleotideAndPositionBlock :: [(Vector Nucleotide, Vector (Position ZeroBased))] -> NucleotideAndPositionBlock
mkNucleotideAndPositionBlock [] = NucleotideAndPositionBlock 0 0 V.empty V.empty
mkNucleotideAndPositionBlock xs = NucleotideAndPositionBlock numberOfPeople max_length (V.map fromIntegral $ mconcat $ map (pad n max_length . fst) xs) (mconcat $ map (pad 0 max_length . V.map unPos . snd) xs)
    where numberOfPeople = length xs
          max_length = maximum (map (V.length . fst) xs)
          unPos (Position p) = fromIntegral p :: CInt


mkPatterns :: [Pattern] -> Patterns
mkPatterns patterns = Patterns $ mconcat (map patternToVector patterns) <> V.fromList [0]


C.context (C.baseCtx <> C.vecCtx <> C.fptrCtx)
C.include "<stdio.h>"
C.include "<math.h>"
C.include "<stdlib.h>"
C.include "<string.h>"

foreign import ccall "&free" freePtr :: FunPtr (Ptr CInt-> IO ())

findPatterns :: CInt -> CInt -> Vector CChar -> Vector CInt -> Vector CInt -> Ptr CInt -> IO (Ptr CInt)
findPatterns sample_size block_size vec pos pat res_size = [C.block| int* {
    typedef struct Match {
        int score;
        int pattern;
        int sample;
        int position;
        int matched[30];
        struct Match* next;
    } Match;

    typedef struct MatchRecord {
        int score;
        int pattern;
        int sample;
        int position;
        int matched[30];
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

                    for(int k = 0; k<30; k++) { match->matched[k] = 0; }
                    int nchar = 29;
                    if(pattern_length<nchar) nchar=pattern_length;
                    if((sample+1)*block_size -i < nchar) nchar=(sample+1)*block_size - i;

                    for(int k = 0; k<nchar; k++) {
                        match->matched[k] = (int)(in[i+k]);
                        if(!in[i+k]) break;
                    }

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

    int* matches = (int*)malloc (match_number * sizeof (MatchRecord));
    Match* p = match_head;
    int n = 0;
    while(p) {
        printf("Match score=%d pattern_id=%d sample_id=%d position=%d \n", p->score, p->pattern, p->sample, p->position);
        matches[n+0] = p->pattern;
        matches[n+1] = p->score;
        matches[n+2] = p->position;
        matches[n+3] = p->sample;
        
        for(int k = 0; k<30; k++) { matches[n*4+4+k] = 0; }
        for(int k = 0; k<29; k++) {
            matches[n+4+k] = (int)(p->matched[k]);
            //if(!p->matched[k]) break;
        }
        
        Match* to_be_freed = p;
        p = p->next;
        free(to_be_freed);
        n = n + (sizeof (MatchRecord) / 4);
    }
    printf("match_number in this block: %d\n", match_number);
    *($(int* res_size)) = match_number;

    return (int*)matches;
    } |]


findPatternsInBlock :: NucleotideAndPositionBlock -> Patterns -> IO [Match]
findPatternsInBlock (NucleotideAndPositionBlock numberOfPeople block_size inputData positionData) (Patterns patternData) = do
    result <- alloca $ \n_ptr -> do
        x <- findPatterns (fromIntegral numberOfPeople) (fromIntegral block_size) inputData positionData patternData (n_ptr :: Ptr CInt)
        result_size <- peek n_ptr
        fptr <- newForeignPtr freePtr x
        print ("result_size", result_size)
        return $ V.unsafeFromForeignPtr0 fptr (fromIntegral (34 * result_size))
    return $ vectorToMatches (result :: Vector CInt)