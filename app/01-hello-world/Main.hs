{-# LANGUAGE QuasiQuotes #-}
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
import           Foreign                         (nullPtr, sizeOf)
import           Foreign.C.Types                 (CFloat, CInt, CChar)
import           Language.C.Quote.OpenCL         (cfun)
import           Text.PrettyPrint.Mainland       (prettyCompact)
import           Text.PrettyPrint.Mainland.Class (ppr)
import Data.Time.Clock.POSIX (getPOSIXTime)
import Data.Int
import Data.Monoid ((<>))
import Data.Maybe (mapMaybe)


kernelSource :: String
kernelSource = prettyCompact . ppr $ [cfun|
    kernel void jaccard_distance(
        int size_chromosome,
        int block_size,
        int max_matches,
        global char *in,
        global int *in_reference_position,
        constant float *patterns,
        global int *out_matches,
        global int *out_debug
    ) {
        // For a given position i in the chromosome
        const int i = get_global_id(0);
        const int participant_id = i / block_size;
        const int base_index = i * 3*sizeof(int);

        int k = 0; // Position in the patterns array
        int match_id = 0;
        for(int j = 0; j < max_matches + 1; j++) {
            out_matches[base_index + (j * 4) + 0] = 0;
            out_matches[base_index + (j * 4) + 1] = 0;
            out_matches[base_index + (j * 4) + 2] = 0;
            out_matches[base_index + (j * 4) + 3] = 0;
        }
        int pattern_id = 0;
        while(1) {
            int pattern_length = (int)(patterns[k]);
            if(pattern_length <= 0) break;
            k = k + 1;

            float score_inter = 0;
            float score_union = 0;
            // For each position in the pattern

            // Equivalent to j < pattern_length && i+j < (participant_id + 1) * block_size. Saves 5% of runtime
            const max_index = min(pattern_length, (participant_id + 1) * block_size - i);

            for(int j = 0; j < max_index; k = k + 5, j++) {
                score_inter += patterns[k+in[i+j]];
                // Equivalent to the following. Saves 33% of runtime

                //if(in[i+j] == 0) { score_inter += patterns[k+0]; }      // Nucleotide is A
                //else if(in[i+j] == 1) { score_inter += patterns[k+1]; } // Nucleotide is C
                //else if(in[i+j] == 2) { score_inter += patterns[k+2]; } // Nucleotide is G
                //else if(in[i+j] == 3) { score_inter += patterns[k+3]; } // Nucleotide is T

                score_union += patterns[k+4]; // max of potential scores for this position
            }
            float score = score_inter / score_union;

            if(score > 0.5 && match_id < max_matches) {
                out_matches[base_index + (match_id * 4) + 0] = pattern_id;               // Pattern ID (0 based)
                out_matches[base_index + (match_id * 4) + 1] = (int)(score * 1000);      // Matching score
                out_matches[base_index + (match_id * 4) + 2] = in_reference_position[i]; // Position in the genome
                out_matches[base_index + (match_id * 4) + 3] = participant_id;           // ID of the participant (0 based)
                match_id = match_id + 1;
            }
            pattern_id++;
        }
        out_debug[i] = 0;
    }
|]


repeatIOAction :: Int -> (Int -> IO a) -> IO a     -- Inputs: integer and IO. Outputs: IO 
repeatIOAction 0 action = action 0             -- exit recursive loop here
repeatIOAction n action = do
    _ <- action n                                 -- action to perform
    repeatIOAction (n-1) action -- decrement n to make it recursive

inputData :: Vector CChar
inputData = mconcat $ [inputDataSample0, inputDataSample1] <> (take (100000 - 2) $ repeat inputDataSample2)

inputDataSample0 :: Vector CChar
inputDataSample0 = V.fromList [0,0,0,0,1,2,0] <> V.fromList (take 93 (repeat 0)) -- AAAACGA

inputDataSample1 :: Vector CChar
inputDataSample1 = V.fromList [1,2,0,0,0,0,0] <> V.fromList (take 93 (repeat 0)) -- CGAAAAA

inputDataSample2 :: Vector CChar
inputDataSample2 = V.fromList [0,0,0,0,0,0,0] <> V.fromList (take 93 (repeat 0)) -- AAAAAAA

inputDataPositions :: Vector CInt
inputDataPositions = V.fromList (mconcat $ take 10000 $ repeat p)
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

vectorToMatches :: Vector Int32 -> [Match]
vectorToMatches v = mapMaybe toMatch (vector4uples v)
    where toMatch x = if x V.! 1 > 0 then Just (Match (fromIntegral $ x V.! 0) (fromIntegral $ x V.! 1) (fromIntegral $ x V.! 2) (fromIntegral $ x V.! 3)) else Nothing

vector4uples :: V.Storable a => Vector a -> [Vector a]
vector4uples v = if (V.length a == 4) then a : vector4uples rest else []
    where (a, rest) = V.splitAt 4 v


patternData :: Vector CFloat
patternData = patternsToVector $ [
    [Pweight 0 1 0 0, Pweight 0 0 1 0] -- CG
    ,[Pweight 0 0.1 0 0, Pweight 0 1 0 0, Pweight 0 0 1 0, Pweight 0.5 0 0 0] -- cCGA
    ,[Pweight 0 1 0 0, Pweight 0 0 1 0] -- CGA
  ] ++ take 50 (repeat [Pweight 0 0 0 0, Pweight 0 0 0 0, Pweight 0 0 0 0, Pweight 0 0 0 0, Pweight 0 0 0 0, Pweight 0 0 0 0, Pweight 0 0 0 0, Pweight 0 0 0 0, Pweight 0 0 0 0, Pweight 0 0 0 0]) ++ [[Pweight 0 1 0 0, Pweight 0 0 1 0]]

nElem  = V.length inputData
nBytes = nElem * sizeOf (undefined :: CChar)
nBytesPositions = nElem * sizeOf (undefined :: CInt)

maxMatchesPerPosition = 3
nBytesMatches = nElem * (maxMatchesPerPosition * 4 * sizeOf (undefined :: Int32))

patternToVector :: Pattern -> Vector CFloat
patternToVector p = V.fromList [fromIntegral $ length p] <> go p
  where
    go [] = V.empty
    go (x:xs) =  V.fromList [wa x, wc x, wg x, wt x, maximum [wa x, wc x, wg x, wt x]] <> go xs

patternsToVector :: [Pattern] -> Vector CFloat
patternsToVector patterns = mconcat (map patternToVector patterns) <> V.fromList [0]

nElemPattern  = V.length patternData
nBytesPattern = nElemPattern * sizeOf (undefined :: CFloat)

--loop :: OpenCLState -> CLKernel -> CLCommandQueue -> IO (Vector CFloat)
loop state kernel queue bufIn bufInPositions bufInPatterns bufOut_matches bufOut_debug = do

    -- Copy our input data Vector to the input buffer; blocks until complete
    --writeVectorToImage :: state bufIn (Int,Int,1) -> inputData -- dums: w h d
    writeVectorToBuffer state bufIn inputData
    writeVectorToBuffer state bufInPositions inputDataPositions
    writeVectorToBuffer state bufInPatterns patternData

    -- Run the kernel
    clSetKernelArgSto kernel 0 ((V.length patternData) :: Int)
    clSetKernelArgSto kernel 1 ((V.length inputDataSample0) :: Int)
    clSetKernelArgSto kernel 2 (maxMatchesPerPosition :: Int)
    clSetKernelArgSto kernel 3 bufIn
    clSetKernelArgSto kernel 4 bufInPositions
    clSetKernelArgSto kernel 5 bufInPatterns
    clSetKernelArgSto kernel 6 bufOut_matches
    clSetKernelArgSto kernel 7 bufOut_debug

    execEvent <- clEnqueueNDRangeKernel queue kernel [nElem] [] []

    -- Get the result; blocks until complete  vectorToMatches <$> 
    outputMatches <- bufferToVector queue bufOut_matches (nElem*maxMatchesPerPosition * 3) [execEvent] :: IO (Vector Int32)
    outputDebug <- bufferToVector queue bufOut_debug (nElem) [execEvent] :: IO (Vector Int32)

    clReleaseEvent execEvent

    return (outputMatches, outputDebug)

-- 100000 * (3000000000 * 0.003) * 50 * 10 * 10 / 1000000000 / 3600 / 24
-- 52.083333333333336
-- 100000. * 50 * 10 * 10 / 1000000000 = Gops for 1bp
main :: IO ()
main = do
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
    bufOut_debug <- clCreateBuffer context [CL_MEM_READ_ONLY, CL_MEM_ALLOC_HOST_PTR] (nBytes, nullPtr) -- One Int of debug per nucleotide

    --print inputData
    --print inputDataPositions
    --print patternData
    t1 <- getPOSIXTime
    (outputData, outputDebug) <- loop state kernel queue bufIn bufInPositions bufIn2 bufOut_matches bufOut_debug
    t2 <- getPOSIXTime
    print "Time in ms"
    print $ round $ (t2 - t1) * 1000

    _ <- clReleaseMemObject bufIn
    _ <- clReleaseMemObject bufIn2
    _ <- clReleaseMemObject bufOut_matches

    -- Clean up the Context
    _ <- clFlush queue
    _ <- clReleaseContext context

    -- Show our work
    putStrLn "\n* Results *"
    --putStrLn $ "Input:  " ++ show inputData
    --putStrLn $ "Output: " ++ show outputData
    putStrLn $ "Output: " ++ show (vectorToMatches outputData)
    putStrLn $ "Non regression: " ++ show (vectorToMatches outputData == [Match {mPatternId = 1, mScore = 961, mPosition = 3, mSampleId = 0},Match {mPatternId = 0, mScore = 1000, mPosition = 4, mSampleId = 0},Match {mPatternId = 2, mScore = 1000, mPosition = 4, mSampleId = 0},Match {mPatternId = 53, mScore = 1000, mPosition = 4, mSampleId = 0},Match {mPatternId = 0, mScore = 1000, mPosition = 0, mSampleId = 1},Match {mPatternId = 2, mScore = 1000, mPosition = 0, mSampleId = 1},Match {mPatternId = 53, mScore = 1000, mPosition = 0, mSampleId = 1}])
    putStrLn $ "Output: " ++ show (V.take 100 outputDebug)



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

