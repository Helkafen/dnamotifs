{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}

module PatternFind ( findPatternsInBlock, mkPatterns, mkNucleotideAndPositionBlock, NucleotideAndPositionBlock, Patterns, blockInfo ) where

import           Foreign                         (Ptr, alloca, peek, FunPtr) 
import           Foreign.ForeignPtr              (newForeignPtr)
import           Foreign.C.Types                 (CInt)
import qualified Language.C.Inline               as C
import qualified Data.Vector.Storable            as STO
import qualified Data.Vector                     as V
import qualified Data.ByteString                 as B
import           Data.Monoid                     ((<>))
import           Data.Maybe                      (mapMaybe)

import Types

newtype Patterns = Patterns (STO.Vector CInt)

-- numberOfPeople, blockSize, nucleotides, positions
data NucleotideAndPositionBlock = NucleotideAndPositionBlock {-# UNPACK #-}!Int {-# UNPACK #-}!Int {-# UNPACK #-}!B.ByteString {-# UNPACK #-}!(STO.Vector CInt)

blockInfo :: NucleotideAndPositionBlock -> String
blockInfo (NucleotideAndPositionBlock numberOfHaplotypes blockSize _ positions) = "block " <> show numberOfHaplotypes <> " " <> show blockSize <> " " <> show (STO.minimum positions) <> " " <> show (STO.maximum positions)

vectorToMatches :: STO.Vector CInt -> V.Vector Match
vectorToMatches v = V.fromList $ reverse $ mapMaybe toMatch (list4uple $ STO.toList v)
    where toMatch (p:s:pos:sam:matchedSequence) = Just $ Match (fromIntegral p) (fromIntegral s) (fromIntegral pos) (fromIntegral sam) (trimMatchedSequence $ map (Nucleotide . fromIntegral) matchedSequence)
          toMatch [] = Nothing
          toMatch x = error ("Wrong size in Match vector: " ++ show (length x))
          trimMatchedSequence = reverse . dropWhile (==n) . reverse

list4uple :: [a] -> [[a]]
list4uple [] = []
list4uple xs = let (x,rest) = splitAt 34 xs in x:list4uple rest


patternToVector :: Pattern -> STO.Vector CInt
patternToVector p = STO.fromList [fromIntegral $ length p, sum (map (floor . (*1000) . m) p)] <> mconcat (map step p)
    where
    step x = STO.map (floor . (*1000)) (STO.fromList [0, wa x, wc x, wg x, wt x, m x]) -- 0 For "N"
    m x = maximum [wa x, wc x, wg x, wt x]

pad :: STO.Storable a => a -> Int -> STO.Vector a -> STO.Vector a
pad e len v = v <> padding
    where padding = STO.fromList $ replicate (len - STO.length v) e

padBS :: Nucleotide -> Int -> B.ByteString -> B.ByteString
padBS (Nucleotide e) len v = B.append v (B.replicate (len - B.length v) e)

-- Pad with nucleotide N if sizes are different
mkNucleotideAndPositionBlock :: [BaseSequencePosition] -> NucleotideAndPositionBlock
mkNucleotideAndPositionBlock [] = NucleotideAndPositionBlock 0 0 B.empty STO.empty
mkNucleotideAndPositionBlock xs = NucleotideAndPositionBlock numberOfPeople max_length (B.concat $ map (padBS n max_length . seqOf) xs) (mconcat $ map (pad 0 max_length . posOf) xs)
    where numberOfPeople = length xs
          max_length = maximum (map (B.length . seqOf) xs)
          seqOf (BaseSequencePosition nuc _) = nuc
          posOf (BaseSequencePosition _ pos) = pos


mkPatterns :: [Pattern] -> Patterns
mkPatterns patterns = Patterns $ mconcat (map patternToVector patterns) <> STO.fromList [0]


C.context (C.baseCtx <> C.vecCtx <> C.fptrCtx <> C.bsCtx)
C.include "<stdio.h>"
C.include "<math.h>"
C.include "<stdlib.h>"
C.include "<string.h>"

foreign import ccall "&free" freePtr :: FunPtr (Ptr CInt-> IO ())

findPatterns :: CInt -> CInt -> CInt -> B.ByteString -> STO.Vector CInt -> STO.Vector CInt -> Ptr CInt -> IO (Ptr CInt)
findPatterns min_score sample_size block_size vec pos pat res_size = [C.block| int* {
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

    int min_score = $(int min_score);
    int sample_size = $(int sample_size);
    int block_size = $(int block_size);
    int* patterns = $vec-ptr:(int *pat);
    char* in = &($bs-ptr:vec[0]);
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
                if(score >= min_score) {
                    struct Match *match = malloc (sizeof (struct Match));
                    if (match == 0) {
                        printf("Failed malloc for a Match\n");
                        return ((int*)(-1));
                    }
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
    if(matches == 0) {
        printf("Failed malloc for a Match array of size %d\n", match_number);
        return ((int*)(-1));
    }
    Match* p = match_head;
    int n = 0;
    while(p) {
        //printf("Match score=%d pattern_id=%d sample_id=%d position=%d \n", p->score, p->pattern, p->sample, p->position);
        matches[n+0] = p->pattern;
        matches[n+1] = p->score;
        matches[n+2] = p->position;
        matches[n+3] = p->sample;
        
        for(int k = 0; k<30; k++) { matches[n+4+k] = 0; }
        for(int k = 0; k<29; k++) {
            matches[n+4+k] = (int)(p->matched[k]);
            if(!p->matched[k]) break;
        }
        
        Match* to_be_freed = p;
        p = p->next;
        free(to_be_freed);
        n = n + (sizeof (MatchRecord) / 4);
    }
    //printf("match_number in this block: %d\n", match_number);
    *($(int* res_size)) = match_number;

    return (int*)matches;
    } |]


findPatternsInBlock :: Int -> NucleotideAndPositionBlock -> Patterns -> IO (V.Vector Match)
findPatternsInBlock minScore (NucleotideAndPositionBlock numberOfPeople block_size inputData positionData) (Patterns patternData) = do
    result <- alloca $ \n_ptr -> do
        x <- findPatterns (fromIntegral minScore) (fromIntegral numberOfPeople) (fromIntegral block_size) inputData positionData patternData (n_ptr :: Ptr CInt)
        result_size <- peek n_ptr
        fptr <- newForeignPtr freePtr x
        --print ("result_size", result_size)
        return $ STO.unsafeFromForeignPtr0 fptr (fromIntegral (34 * result_size))
    return $ vectorToMatches (result :: STO.Vector CInt)