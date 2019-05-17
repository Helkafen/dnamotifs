{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}


{-|
Module      : Main
Description : OpenCL Hello World Example
Copyright   : (c) Jonathan Merritt, 2017
License     : BSD3
Maintainer  : j.s.merritt@gmail.com
Stability   : experimental
Portability : POSIX + OpenCL

This is a Hello World OpenCL example that stores buffer data in vectors.

More information in this blog post:
  https://lancelet.github.io/posts/2017-12-26-opencl-helloworld.html
-}

module Main where

import           CLUtil         (OpenCLState (OpenCLState),
                                                  clContext,
                                                  clDevice, clQueue)
import CLUtil.VectorBuffers (bufferToVector, writeVectorToBuffer)
import           Control.Parallel.OpenCL

import           Control.Monad                   (forM_)
import           Data.Vector.Storable            (Vector)
import qualified Data.Vector.Storable            as V
import qualified Data.Vector.Storable.Mutable    as VM
import           Foreign                         (nullPtr, sizeOf, Ptr, alloca, peek, ForeignPtr, FunPtr) 
import Foreign.ForeignPtr (newForeignPtr)
import           Foreign.C.Types                 (CFloat, CInt, CChar, CDouble)
import           Language.C.Quote.OpenCL         (cfun)
import           Text.PrettyPrint.Mainland       (prettyCompact)
import           Text.PrettyPrint.Mainland.Class (ppr)
import Data.Time.Clock.POSIX (getPOSIXTime)
import Data.Int
import Data.Monoid ((<>))
import Data.Maybe (mapMaybe)
import Text.RawString.QQ
import qualified Data.Bits as Bits
import qualified Language.C.Inline as C


kernelSource :: String
kernelSource = 
    [r|
    int add(int a) {
        return (a+1);
    }

    struct match {
        ushort score;
        uchar pattern_id;
        ushort position;
        int participant_id;
    };

    __kernel void jaccard_distance(
        int block_size,
        int max_matches,
        __global char *in,
        __global int *in_reference_position,
        __constant float *patterns,
        __global int *out_matches,
        __global int *out_matches2)
    {
    const int i = get_global_id(0);
    const int participant_id = i / block_size;
    int k = 0;
    
    ushort pattern_id = 0;
    int best_match_1 = 0;
    int best_match_2 = 0;
    
    while (1) {
        int pattern_length = (int) patterns[k];
        if (pattern_length <= 0) break;
        k = k + 1;
        
        float score_inter = 0;
        float score_union = 0;
        const int max_index = min(pattern_length, (participant_id + 1) * block_size - i);
        
        // 9 ops per base
        for (int j = 0; j < max_index; k = k + 5, j++) {
            score_inter += patterns[k + in[i + j]];
            score_union += patterns[k + 4]; // TODO that's actually a constant, except at the end of the chromosome. I should add the score_union right after the pattern length
        }
        
        unsigned short score = 0;

        // One best match per position.
        int match_record = (score << 16) + pattern_id;

        //if(match_record > best_match_1) {
        //    best_match_2 = best_match_1;
        //    best_match_1 = match_record;
        //}
        //else if(match_record > best_match_2) {vector4uples :: V.Storable a => Vector a -> [Vector a]
            vector4uples v = if (V.length a == 4) then a : vector4uples rest else []
                where (a, rest) = V.splitAt 4 v
        //    best_match_2 = match_record;
        //}
        //barrier(CLK_GLOBAL_MEM_FENCE); TODO?
        best_match_1 = max(best_match_1, match_record);

        pattern_id++;
    }

    out_matches[i] = best_match_1;
    out_matches2[i] = best_match_2;
    }|]


repeatIOAction :: Int -> (Int -> IO a) -> IO a     -- Inputs: integer and IO. Outputs: IO 
repeatIOAction 0 action = action 0             -- exit recursive loop here
repeatIOAction n action = do
    _ <- action n                                 -- action to perform
    repeatIOAction (n-1) action -- decrement n to make it recursive

inputData :: Vector CChar
inputData = mconcat $ [inputDataSample0, inputDataSample1] <> (take (numberOfPeople - 2) $ repeat inputDataSample2)

numberOfPeople = 100000


--inputDataNoMatch :: Vector CChar
--inputDataNoMatch = mconcat (take 100000 $ repeat inputDataSample2)

type Nucleotide = CChar
n = 0 :: Nucleotide
a = 1 :: Nucleotide
c = 2 :: Nucleotide
g = 3 :: Nucleotide
t = 4 :: Nucleotide



-- Must be 100 bases long
inputDataSample0 :: Vector CChar
inputDataSample0 = V.fromList [a,a,a,a,c,g,a] <> V.fromList (take 93 (repeat a))

inputDataSample1 :: Vector CChar
inputDataSample1 = V.fromList [c,g,a,a,a,a,a] <> V.fromList (take 93 (repeat a)) 

inputDataSample2 :: Vector CChar
inputDataSample2 = V.fromList [a,a,a,a,a,a,a] <> V.fromList (take 93 (repeat a))

inputDataPositions :: Vector CInt
inputDataPositions = V.fromList (mconcat $ take numberOfPeople $ repeat p)
    where p = take (V.length inputDataSample0) [0..]

data Pweight = Pweight {
    wa :: CFloat,
    wc :: CFloat,
    wg :: CFloat,
    wt :: CFloat
} deriving (Show)

type Pattern = [Pweight]

data Match = Match {
    mPatternId :: Int, -- 0 based
    mScore :: Int,     -- 0 - 1000
    mPosition :: Int,  -- 0 based
    mSampleId :: Int
} deriving (Eq, Show)

numberOfBases = 100

vectorToMatches :: Vector CInt -> [Match]
vectorToMatches v = map toMatch (list4uple $ V.toList v)
    where toMatch (p,s,pos,sam) = Match (fromIntegral p) (fromIntegral s) (fromIntegral pos) (fromIntegral sam)

list4uple :: [a] -> [(a,a,a,a)]
list4uple (a:b:c:d:xs) = (a,b,c,d):(list4uple xs)
list4uple _ = []

vector4uples :: V.Storable a => Vector a -> [Vector a]
vector4uples v = if (V.length a == 4) then a : vector4uples rest else []
    where (a, rest) = V.splitAt 4 v


patternData :: Vector CInt
patternData = patternsToVector $ [
    [Pweight 0 1 0 0, Pweight 0 0 1 0] -- CG
    ,[Pweight 0 0.1 0 0, Pweight 0 1 0 0, Pweight 0 0 1 0, Pweight 0.5 0 0 0] -- cCGA
    ,[Pweight 0 1 0 0, Pweight 0 0 1 0] -- CGA
  ] ++ take 50 (repeat [cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern, cPattern]) ++ [[Pweight 0 1 0 0, Pweight 0 0 1 0]]

cPattern = Pweight 0 0 0.001 0

nElem  = V.length inputData
nBytes = nElem * sizeOf (undefined :: CChar)
nBytesPositions = nElem * sizeOf (undefined :: CInt)

maxMatchesPerPosition = 3
--nBytesMatches = nElem * (maxMatchesPerPosition * 4 * sizeOf (undefined :: Int32))
nBytesMatches = nElem * sizeOf (undefined :: Int32)

patternToVector :: Pattern -> Vector CInt
patternToVector p = V.fromList [fromIntegral $ length p, sum (map (floor . (*1000) . max) p)] <> mconcat (map step p)
  where
    step x = V.map (floor . (*1000)) (V.fromList [0, wa x, wc x, wg x, wt x, max x]) -- 0 For "N"
    max x = maximum [wa x, wc x, wg x, wt x]

patternsToVector :: [Pattern] -> Vector CInt
patternsToVector patterns = mconcat (map patternToVector patterns) <> V.fromList [0]

nElemPattern  = V.length patternData
nBytesPattern = nElemPattern * sizeOf (undefined :: CFloat)

--loop :: OpenCLState -> CLKernel -> CLCommandQueue -> IO (Vector CFloat)
loop state kernel queue bufIn bufInPositions bufInPatterns bufOut_matches bufOut_matches2 bufOut_debug input = do

    -- Copy our input data Vector to the input buffer; blocks until complete
    --writeVectorToImage :: state bufIn (Int,Int,1) -> inputData -- dums: w h d
    writeVectorToBuffer state bufIn input
    writeVectorToBuffer state bufInPositions inputDataPositions
    writeVectorToBuffer state bufInPatterns patternData

    -- Run the kernel
    clSetKernelArgSto kernel 0 ((V.length inputDataSample0) :: Int)
    clSetKernelArgSto kernel 1 (maxMatchesPerPosition :: Int)
    clSetKernelArgSto kernel 2 bufIn
    clSetKernelArgSto kernel 3 bufInPositions
    clSetKernelArgSto kernel 4 bufInPatterns
    clSetKernelArgSto kernel 5 bufOut_matches
    clSetKernelArgSto kernel 6 bufOut_matches2
    --clSetKernelArgSto kernel 7 bufOut_debug

    execEvent <- clEnqueueNDRangeKernel queue kernel [nElem] [] []

    -- Get the result; blocks until complete  vectorToMatches <$> 
    outputMatches <- bufferToVector queue bufOut_matches (nElem) [execEvent] :: IO (Vector Int32)
    outputMatches2 <- bufferToVector queue bufOut_matches2 (nElem) [execEvent] :: IO (Vector Int32)
    --outputDebug <- bufferToVector queue bufOut_debug (nElem) [execEvent] :: IO (Vector Int32)

    clReleaseEvent execEvent

    --return (outputMatches, outputDebug)
    return [outputMatches, outputMatches2]



C.context (C.baseCtx <> C.vecCtx <> C.fptrCtx)
C.include "<stdio.h>"
C.include "<math.h>"
C.include "<stdlib.h>"

-- | @readAndSum n@ reads @n@ numbers from standard input and returns
-- their sum.
--readAndSum :: CInt -> IO CInt
--readAndSum n  = [C.block| int {
--    // Read and sum n integers
--    int i, sum = 0, tmp;
--    for (i = 0; i < $(int n); i++) {
--        scanf("%d", &tmp);
--        sum += tmp;
--    }
--    return sum;
--    } |]


--[2, 2000, 0,1000,0,0,1000,   0,0,1000,0,1000,
-- 12
-- 4, 2600, 0,100,0,0,100,     0,1000,0,0,1000,    0,0,1000,0,1000,    500,0,0,0,500,
-- 2, 2000, 0,1000,0,0,1000,   0,0,1000,0,1000,
-- 0]

foreign import ccall "&free" freePtr :: FunPtr (Ptr CInt-> IO ())

find_patterns :: CInt -> CInt -> Vector CChar -> Vector CInt -> Ptr CInt -> IO (Ptr CInt)
find_patterns sample_size block_size vec pat res_size = [C.block| int* {
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
                    match->position = position;
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



-- 100000 * (3000000000 * 0.003) * 50 * 10 * 10 / 1000000000 / 3600 / 24
-- 52.083333333333336
-- 100000. * 50 * 10 * 10 / 1000000000 = Gops for 1bp

-- To use the vector anti-quoters, we need the 'C.vecCtx' along with the 'C.baseCtx'.


main = do
    --x <- readAndSum 2
    --print x
    --pf <- VM.new 10000 :: VM.IOVector CInt
    let block_size = 100
    --print patternData
    t1 <- getPOSIXTime
    result <- alloca $ \n_ptr -> do
        x <- find_patterns (fromIntegral numberOfPeople) block_size (inputData :: Vector CChar) (patternData :: Vector CInt) (n_ptr :: Ptr CInt)  --(pf :: VM.IOVector CInt)
        result_size <- peek (n_ptr :: Ptr CInt)
        fptr <- newForeignPtr freePtr x
        print ("result_size", result_size)
        --return result_size
        return $ V.unsafeFromForeignPtr0 fptr (fromIntegral (4 * result_size))
    --print result
    t2 <- getPOSIXTime
    print "Time in ms"
    print $ round $ (t2 - t1) * 1000
    let a = vectorToMatches (result :: Vector CInt)
    print a
    print (a == [Match {mPatternId = 53, mScore = 1000, mPosition = 0, mSampleId = 1},Match {mPatternId = 2, mScore = 1000, mPosition = 0, mSampleId = 1},Match {mPatternId = 0, mScore = 1000, mPosition = 0, mSampleId = 1},Match {mPatternId = 53, mScore = 1000, mPosition = 4, mSampleId = 0},Match {mPatternId = 2, mScore = 1000, mPosition = 4, mSampleId = 0},Match {mPatternId = 0, mScore = 1000, mPosition = 4, mSampleId = 0},Match {mPatternId = 1, mScore = 961, mPosition = 3, mSampleId = 0}])
    --print $ patternToVector [Pweight 0 0.1 0 0, Pweight 0 1 0 0, Pweight 0 0 1 0, Pweight 0.5 0 0 0]
    


main2 :: IO ()
main2 = do
    putStrLn "* Hello World OpenCL Example *"

    -- Describe the OpenCL Environment
    putStrLn "\n* OpenCL Platform Environment *"
    describePlatforms

    -- Create a Context, Queue and a CLUtil OpenCLState
    context <- clCreateContextFromType [] [CL_DEVICE_TYPE_GPU] print
    device  <- head <$> clGetContextDevices context
    queue   <- clCreateCommandQueue context device []
    -- NB: OpenCLState is used by CLUtil when manipulating vector buffers
    let state = OpenCLState
                { clDevice  = device
                , clContext = context
                , clQueue   = queue
                }

    -- Create the Kernel
    program <- clCreateProgramWithSource context kernelSource
    clBuildProgram program [device] ""
    kernel <- clCreateKernel program "jaccard_distance"

    -- Buffers for input and output data.
    -- We request OpenCL to create a buffer on the host (CL_MEM_ALLOC_HOST_PTR)
    -- since we're using CPU. The performance here may not be ideal, because
    -- we're copying the buffer. However, it's safe, and not unduly nested.
    bufIn <- clCreateBuffer context [CL_MEM_READ_ONLY, CL_MEM_ALLOC_HOST_PTR] (nBytes, nullPtr)
    bufInPositions <- clCreateBuffer context [CL_MEM_READ_ONLY, CL_MEM_ALLOC_HOST_PTR] (nBytesPositions, nullPtr)
    bufIn2 <- clCreateBuffer context [CL_MEM_READ_ONLY, CL_MEM_ALLOC_HOST_PTR] (nBytesPattern, nullPtr)
    bufOut_matches <- clCreateBuffer context [CL_MEM_WRITE_ONLY, CL_MEM_ALLOC_HOST_PTR] (nBytesMatches, nullPtr)
    bufOut_matches2 <- clCreateBuffer context [CL_MEM_WRITE_ONLY, CL_MEM_ALLOC_HOST_PTR] (nBytesMatches, nullPtr)
    bufOut_debug <- clCreateBuffer context [CL_MEM_READ_ONLY, CL_MEM_ALLOC_HOST_PTR] (nBytes, nullPtr) -- One Int of debug per nucleotide

    --print inputData
    --print inputDataPositions
    --print patternData
    --t1 <- getPOSIXTime
    --(outputData, outputDebug) <- loop state kernel queue bufIn bufInPositions bufIn2 bufOut_matches bufOut_debug
    --outputData <- loop state kernel queue bufIn bufInPositions bufIn2 bufOut_matches bufOut_debug inputDataNoMatch
    --t2 <- getPOSIXTime
    --print "Time in ms (no match)"
    --print $ round $ (t2 - t1) * 1000

    t1 <- getPOSIXTime
    --(outputData, outputDebug) <- loop state kernel queue bufIn bufInPositions bufIn2 bufOut_matches bufOut_debug
    outputData <- loop state kernel queue bufIn bufInPositions bufIn2 bufOut_matches bufOut_matches2 bufOut_debug inputData
    t2 <- getPOSIXTime
    print "Time in ms"
    print $ round $ (t2 - t1) * 1000
    

    

    _ <- clReleaseMemObject bufIn
    _ <- clReleaseMemObject bufIn2
    _ <- clReleaseMemObject bufOut_matches
    _ <- clReleaseMemObject bufOut_matches2

    -- Clean up the Context
    _ <- clFlush queue
    _ <- clReleaseContext context

    -- Show our work
    putStrLn "\n* Results *"
    putStrLn $ "Input:  " ++ show (V.length inputData)
    --putStrLn $ "Output: " ++ show outputData
    --let matches = concatMap (vectorToMatches 500) outputData
    --putStrLn $ "Output: " ++ show (matches)
    --putStrLn $ "Non regression: " ++ show (matches == [Match {mPatternId = 1, mScore = 961, mPosition = 3, mSampleId = 0},Match {mPatternId = 0, mScore = 1000, mPosition = 4, mSampleId = 0},Match {mPatternId = 2, mScore = 1000, mPosition = 4, mSampleId = 0},Match {mPatternId = 53, mScore = 1000, mPosition = 4, mSampleId = 0},Match {mPatternId = 0, mScore = 1000, mPosition = 0, mSampleId = 1},Match {mPatternId = 2, mScore = 1000, mPosition = 0, mSampleId = 1},Match {mPatternId = 53, mScore = 1000, mPosition = 0, mSampleId = 1}])
    --putStrLn $ "Output: " ++ show (V.take 100 outputDebug)<



-- | Summarises the OpenCL Platforms and their Devices.
--
--   The names of Platforms and all Devices belonging to them are printed to
--   stdout.
describePlatforms :: IO ()
describePlatforms = do

    -- fetch the list of OpenCL Platforms
    platformList <- clGetPlatformIDs :: IO [CLPlatformID]

    -- for each platform,
    forM_ platformList $ \platform -> do

        -- fetch the list of OpenCL Devices
        devs <- clGetDeviceIDs platform CL_DEVICE_TYPE_ALL :: IO [CLDeviceID]

        -- print the Platform name and Device names
        pname platform
        forM_ devs dname

  where
      putPair name value = putStrLn (name ++ value)
      pname p = clGetPlatformInfo p CL_PLATFORM_NAME >>= putPair "Platform: "
      dname d = clGetDeviceName d                    >>= putPair "  Device: "

