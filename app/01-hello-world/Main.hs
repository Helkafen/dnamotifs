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

import           Control.Monad                   (forM_, when)
import           Data.Vector.Storable            (Vector)
import qualified Data.Vector.Storable            as V
import           Foreign                         (nullPtr, sizeOf)
import           Foreign.C.Types                 (CFloat, CInt)
import           Language.C.Quote.OpenCL         (cfun)
import           Text.PrettyPrint.Mainland       (prettyCompact)
import           Text.PrettyPrint.Mainland.Class (ppr)
import Data.Time.Clock.POSIX (getPOSIXTime)
import Data.Int
import Data.Monoid ((<>))
import Data.Maybe (mapMaybe)



-- jaccardDistance (Motif elements length) nucleotides = Score (inter / union)
--    where inter = sum $ zipWith minNucleotideJaccardDistance elements nucleotides
--          union = sum $ zipWith maxNucleotideJaccardDistance elements nucleotides
--          minNucleotideJaccardDistance :: MotifElement -> Nucleotide -> Float
--          minNucleotideJaccardDistance (MotifElement (a, _, _, _)) A = a
--          minNucleotideJaccardDistance (MotifElement (_, c, _, _)) C = c
--          minNucleotideJaccardDistance (MotifElement (_, _, g, _)) G = g
--          minNucleotideJaccardDistance (MotifElement (_, _, _, t)) T = t
--          maxNucleotideJaccardDistance :: MotifElement -> Nucleotide -> Float
--          maxNucleotideJaccardDistance (MotifElement (a, c, g, t)) _ = max (max (max a c) g) t


--instance Storable


-- | The kernel to execute: the equivalient of 'map (*2)'.
kernelSource :: String
kernelSource = prettyCompact . ppr $ [cfun|
    /* This example kernel just does `map (*2)` */
    kernel void doubleArray(
        int size_chromosome,
        int max_matches,
        global float *in,
        global float *patterns,
        global int *out_matches
    ) {
        // For a given position i in the chromosome
        const int i = get_global_id(0);

        int k = 0; // Position in the patterns array
        int match_id = 0;
        for(int j = 0; j < max_matches + 1; j++) {
            out_matches[i * 3*sizeof(int) + (j * 3) + 0] = 0;
            out_matches[i * 3*sizeof(int) + (j * 3) + 1] = 0;
            out_matches[i * 3*sizeof(int) + (j * 3) + 2] = 0;
        }
        while(1) {
            int pattern_length = (int)(patterns[k]);
            if(pattern_length <= 0) break;
            k = k + 1;

            float score_inter = 0;
            float score_union = 0;
            // For each position in the pattern
            int j = 0;
            float m = -1;
            for(j = 0; j < pattern_length && i+j < size_chromosome; k = k + 5, j++) {
                float a = patterns[k+0];
                float c = patterns[k+1];
                float g = patterns[k+2];
                float t = patterns[k+3];
                m = patterns[k+4]; // max of potential scores for this position
                score_union += m;
                if(in[i+j] == 0) { score_inter += a; }
                if(in[i+j] == 1) { score_inter += c; }
                if(in[i+j] == 2) { score_inter += g; }
                if(in[i+j] == 3) { score_inter += t; }
            }
            float score = score_inter / score_union;
            if(score > 0.5 && match_id < max_matches) {
                //barrier(CLK_LOCAL_MEM_FENCE);
                //barrier(CLK_GLOBAL_MEM_FENCE);
                out_matches[i * 3*sizeof(int) + match_id * 3 + 0] = (int)((k-1) / 5) - 1; // Pattern ID (0 based)
                out_matches[i * 3*sizeof(int) + match_id * 3 + 1] = (int)(score * 1000);  // Matching score
                out_matches[i * 3*sizeof(int) + match_id * 3 + 2] = i;                    // Position in the genome (in the current batch?)
                match_id++;
            }
            //break;
        }

        
        barrier(CLK_LOCAL_MEM_FENCE);
    }
|]


repeatIOAction :: Int -> (Int -> IO a) -> IO a     -- Inputs: integer and IO. Outputs: IO 
repeatIOAction 0 action = action 0             -- exit recursive loop here
repeatIOAction n action = do
    _ <- action n                                 -- action to perform
    repeatIOAction (n-1) action -- decrement n to make it recursive

inputData :: Vector CFloat
inputData = V.fromList [0,0,0,0,1,2,0] -- AAAACGAAA

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
    mPosition :: Int   -- 0 based
} deriving (Show)

vectorToMatches :: Vector Int32 -> [Match]
vectorToMatches v = mapMaybe toMatch (vector3uples v)
    where toMatch x = if x V.! 1 > 0 then Just (Match (fromIntegral $ x V.! 0) (fromIntegral $ x V.! 1) (fromIntegral $ x V.! 2)) else Nothing

vector3uples :: V.Storable a => Vector a -> [Vector a]
vector3uples v = if (V.length a == 3) then a : vector3uples rest else []
    where (a, rest) = V.splitAt 3 v


patternData :: Vector CFloat
patternData = patternsToVector [[Pweight 0 1 0 0, Pweight 0 0 1 0]] -- CG --V.fromList [2, 100,101,102,103,103, 0]

nElem  = V.length inputData
nBytes = nElem * sizeOf (undefined :: CFloat)

maxMatchesPerPosition = 3
nBytesMatches = nElem * (maxMatchesPerPosition * 3 * sizeOf (undefined :: Int32))

patternToVector :: Pattern -> Vector CFloat
patternToVector [] = V.empty
patternToVector (x:xs) = V.fromList [wa x, wc x, wg x, wt x, maximum [wa x, wc x, wg x, wt x]] <> patternToVector xs

patternsToVector :: [Pattern] -> Vector CFloat
patternsToVector patterns = V.fromList [fromIntegral $ length patterns] <> mconcat (map patternToVector patterns) <> V.fromList [0]

nElemPattern  = V.length patternData
nBytesPattern = nElemPattern * sizeOf (undefined :: CFloat)

--loop :: CLContext -> OpenCLState -> CLKernel -> CLCommandQueue -> IO (Vector CFloat)
loop context state kernel queue bufIn bufIn2 bufOut_matches = do

    -- Copy our input data Vector to the input buffer; blocks until complete
    --writeVectorToBuffer state bufIn (10000 :: CInt)
    writeVectorToBuffer state bufIn inputData
    writeVectorToBuffer state bufIn2 patternData
    --print "loop"

    -- Run the kernel
    clSetKernelArgSto kernel 0 ((V.length patternData) :: Int)
    clSetKernelArgSto kernel 1 (maxMatchesPerPosition :: Int)
    clSetKernelArgSto kernel 2 bufIn
    clSetKernelArgSto kernel 3 bufIn2
    clSetKernelArgSto kernel 4 bufOut_matches

    execEvent <- clEnqueueNDRangeKernel queue kernel [nElem] [] []

    -- Get the result; blocks until complete  vectorToMatches <$> 
    outputMatches <- bufferToVector queue bufOut_matches (nElem*maxMatchesPerPosition * 3) [execEvent] :: IO (Vector Int32) -- :: IO [Match]

    clReleaseEvent execEvent

    return outputMatches

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
    kernel <- clCreateKernel program "doubleArray"

    -- Buffers for input and output data.
    -- We request OpenCL to create a buffer on the host (CL_MEM_ALLOC_HOST_PTR)
    -- since we're using CPU. The performance here may not be ideal, because
    -- we're copying the buffer. However, it's safe, and not unduly nested.
    bufIn <- clCreateBuffer context [CL_MEM_READ_ONLY, CL_MEM_ALLOC_HOST_PTR] (nBytes, nullPtr)
    bufIn2 <- clCreateBuffer context [CL_MEM_READ_ONLY, CL_MEM_ALLOC_HOST_PTR] (nBytesPattern, nullPtr)
    bufOut_matches <- clCreateBuffer context [CL_MEM_WRITE_ONLY, CL_MEM_ALLOC_HOST_PTR] (nBytesMatches, nullPtr)

    print ("maxMatchesPerPosition") >> print maxMatchesPerPosition
    print patternData
    print nBytesMatches
    print ((nElem + 1) * maxMatchesPerPosition)
    t1 <- (round . (* 1000000)) <$> getPOSIXTime
    outputData <- loop context state kernel queue bufIn bufIn2 bufOut_matches
    t2 <- (round . (* 1000000)) <$> getPOSIXTime
    print (t2-t1)

    _ <- clReleaseMemObject bufIn
    _ <- clReleaseMemObject bufIn2
    _ <- clReleaseMemObject bufOut_matches

    -- Clean up the Context
    clFlush queue
    _ <- clReleaseContext context

    -- Show our work
    putStrLn "\n* Results *"
    --putStrLn $ "Input:  " ++ show inputData
    putStrLn $ "Output: " ++ show outputData
    putStrLn $ "Output: " ++ show (vectorToMatches outputData)



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

