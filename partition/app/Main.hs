{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}

-- |
-- Module      : MCFISP2K
-- Description : Algorithms for Genome partition problems with flexible intergenic regions.
-- Copyright   : (c) Gabriel Siqueira, 2023
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
module Main (main) where

import Control.Concurrent.ParallelIO.Global (parallel, stopGlobalPool)
import Control.Exception (SomeException, displayException, evaluate, try)
import Control.Monad (unless)
import Control.DeepSeq (force)
import Balancer (balancedLowerBound, balancedPartitionWith)
import Data.ByteString.Builder (intDec, toLazyByteString)
import Data.ByteString.Char8 qualified as BS
import Data.ByteString.Lazy qualified as LBS
import Data.Char (isSpace)
import Data.Map qualified as Map
import Data.Time (diffUTCTime, getCurrentTime)
import Genomes (ChromList, GenesIRsF, GenesIRsR, Genome (allGIdx, getGene), Matcher, geneToInt, RigidFlexibleDirectMatcher (..), RigidFlexibleReverseMatcher (..), Sign (..), equalizeLinear, flexibilizeMC, readMCFGenome, readMCRGenome)
import LocalBase
import Options.Applicative
import PApprox (approxLowerBound, approxPartitionPhase2)
import PFpt (fptPartition)
import PGreedy (greedyPartition)
import PSOAR (soarPartition)
import Partition (CommonPartition, getBlocksCorrespondence, writePartition)
import System.Exit (exitFailure)
import System.IO (hPutStrLn, stderr)
import Text.Printf (printf)

data Args = Args
  { partAlg :: PartitionAlgorithm,
    input :: String,
    output :: String,
    noParallel :: Bool,
    signed :: Sign,
    direct :: Bool
  }

data PartitionAlgorithm = SOAR | SOARComb | Greedy | GreedySin | Approx | FPT | ApproxLB deriving (Read, Show)

argsParser :: Parser Args
argsParser =
  Args
    <$> option
      auto
      ( long "algorithm"
          <> short 'a'
          <> metavar "ALGO"
          <> help "Algorithm to use for the partition problem. The options are SOAR, SOARComb, Greedy, GreedySin, Approx, and FPT. There is also and ApproxLB option, which produces and integer corresponding to a lower bound of the number of breakpoints in the common partition (counting both partitions) instead of a partition on the output."
      )
    <*> strOption
      ( long "input"
          <> short 'i'
          <> metavar "IFILE"
          <> help "Input file. Each 4 lines of the input file correspond to a instance, each line has a list of comma or space separated values, and represent in order the origin string, the origin intergenic region list, the target string, and the target intergenic region list. The target intergenic regions may be flexible (intervals l:u) or rigid (integers). Chromosomes are separated by ';' and may start with L (linear, default) or C (circular)."
      )
    <*> strOption
      ( long "outfile"
          <> short 'o'
          <> metavar "OFILE"
          <> help "Output file. For each instance five lines are produces in the file, the first four lines correspond to two reduced genomes produced with the partition algorithm (each character of the string correspond to a block of the partition and each integer of the intergenic regions list correspond to a breakpoint). The last line shows the wall clock time required to produce the partition. With the Approx algorithm, the file OFILE_phase2 has one line per instance with the number of times the second stage of the algorithm broke a block (NA if the instance failed)."
      )
    <*> switch
      ( long "no-par"
          <> help "Do not process the genomes in parallel."
      )
    <*> flag
      Unsigned
      Signed
      ( long "signed"
          <> short 's'
          <> help "Whether the input Strings are signed."
      )
    <*> switch
      ( long "direct"
          <> short 'd'
          <> help "Find a direct common partition (PDCMG: blocks must be completely compatible, without reversals) instead of a reverse common partition (PRCMG: blocks may be compatible after a reversal). The direct problem is usually solved for unsigned genomes."
      )

opts :: ParserInfo Args
opts =
  info
    (argsParser <**> helper)
    ( fullDesc
        <> progDesc
          "Algorithm for genome partition problems with flexible intergenic regions."
    )

main :: IO ()
main = do
  args <- execParser opts
  contents <- BS.readFile (input args)
  let quadruples = zip [(1 :: Int) ..] . toQuadruples . filter isContent . BS.lines $ contents
  results <-
    if noParallel args
      then mapM (runOne args) quadruples
      else do
        results <- parallel $ map (runOne args) quadruples
        stopGlobalPool
        return results
  case partAlg args of
    Approx -> BS.writeFile (output args ++ "_phase2") . BS.unlines $ map (maybe "NA" (BS.pack . show) . snd) results
    _ -> return ()
  unless (all fst results) exitFailure
  where
    -- an error in one instance is written on its output file (and reported on
    -- stderr), the other instances are still processed
    runOne args (i, bstrs) = do
      start <- getCurrentTime
      let run m = produceBlockMatch m (direct args) (partAlg args) (signed args) bstrs
      result <- try (evaluate . force =<< if direct args then run RFDM else run RFRM)
      end <- getCurrentTime
      let time = BS.pack . (show :: Double -> String) . realToFrac $ diffUTCTime end start
          outFile = output args ++ "_" ++ printf "%04d" i
      case result of
        Right (bstrs', phase2) -> do
          BS.writeFile outFile . BS.unlines $ bstrs' ++ ["# Time: " <> time]
          return (True, phase2)
        Left e -> do
          let msg = takeWhile (/= '\n') $ displayException (e :: SomeException)
          hPutStrLn stderr $ "Instance " ++ show i ++ ": " ++ msg
          BS.writeFile outFile . BS.unlines $ ["# Error: " <> BS.pack msg, "# Time: " <> time]
          return (False, Nothing)

    -- skip blank lines and comments (lines starting with '#', possibly after spaces)
    isContent l = case BS.uncons (BS.dropWhile isSpace l) of
      Nothing -> False
      Just (c, _) -> c /= '#'

    toQuadruples (s1 : i1 : s2 : i2 : ss) = (s1, i1, s2, i2) : toQuadruples ss
    toQuadruples [] = []
    toQuadruples _ = error "Incorrect number of lines."

-- | Run the algorithm on one instance. The matcher gives the compatibility between
-- blocks: RFDM (direct common partition, PDCMG) or RFRM (reverse, PRCMG). Returns the
-- lines of the output and, for Approx, the number of times its second stage broke a block.
produceBlockMatch :: (Matcher m GenesIRsR GenesIRsF) => m GenesIRsR GenesIRsF -> Bool -> PartitionAlgorithm -> Sign -> (BS.ByteString, BS.ByteString, BS.ByteString, BS.ByteString) -> IO ([BS.ByteString], Maybe Int)
produceBlockMatch matcher isDirect alg sign (s1, i1, s2, i2)
  | isDirect && sign == Signed && not (directPossible g h) =
      error "ERROR: There is no direct common partition, a gene has more copies with one of the signs in one genome and more copies with the other sign in the other genome (use a reverse partition, without --direct)."
  | otherwise = case alg of
    ApproxLB -> return ([BS.pack . show $ balancedLowerBound approxLowerBound matcher g h], Nothing)
    _ -> do
      (maybe_part, fullComp, phase2) <- getPartition matcher alg g h
      case maybe_part of
        Nothing -> if fullComp then return (["# Error: Execution finish but no partition was found"], phase2) else return (["# No solution found in time"], phase2)
        Just part ->
          let (s1', i1', s2', i2') = writePartition part
              bc = writeBlocksCorrespondence $ getBlocksCorrespondence matcher part
           in return ([s1', i1', s2', i2', bc, if fullComp then "# Exact solution" else "# Partial solution"], phase2)
  where
    g0 = readMCRGenome sign s1 i1
    -- a target without flexible intervals (l:u) is read as rigid
    h0 = if BS.elem ':' i2 then readMCFGenome sign s2 i2 else flexibilizeMC 0 (readMCRGenome sign s2 i2)
    -- both genomes get the same number of linear chromosomes (with empty ones)
    (g, h) = equalizeLinear sign g0 h0

-- | A direct common partition only matches genes with the same sign, so for each
-- label the genome with more positive copies must also have at least as many
-- negative copies (the excess goes to exclusive blocks of only one genome).
directPossible :: ChromList GenesIRsR -> ChromList GenesIRsF -> Bool
directPossible g h = all ok (Map.keys (Map.union cg ch))
  where
    count q = Map.fromListWith (+) [(geneToInt (getGene i q), 1 :: Int) | i <- allGIdx q]
    cg = count g
    ch = count h
    diff x = Map.findWithDefault 0 x cg - Map.findWithDefault 0 x ch
    ok x = x <= 0 || diff x * diff (-x) >= 0

writeBlocksCorrespondence :: [[Int]] -> BS.ByteString
writeBlocksCorrespondence = BS.unwords . (\l -> interleavelists l (replicate (length l - 1) "|")) . map (BS.unwords . map (LBS.toStrict . toLazyByteString . intDec))

getPartition :: (Matcher m GenesIRsR GenesIRsF) => m GenesIRsR GenesIRsF -> PartitionAlgorithm -> ChromList GenesIRsR -> ChromList GenesIRsF -> IO (Maybe (CommonPartition GenesIRsR GenesIRsF), Bool, Maybe Int)
getPartition matcher FPT g h = (\(p, full) -> (p, full, Nothing)) <$> fptPartition 600000000 matcher g h
getPartition matcher SOAR g h = return (Just $ soarPartition False matcher g h, True, Nothing)
getPartition matcher SOARComb g h = return (Just $ soarPartition True matcher g h, True, Nothing)
getPartition matcher Greedy g h = return (Just $ greedyPartition False matcher g h, True, Nothing)
getPartition matcher GreedySin g h = return (Just $ greedyPartition True matcher g h, True, Nothing)
getPartition matcher Approx g h =
  let (part, phase2) = balancedPartitionWith 0 approxPartitionPhase2 matcher g h
   in return (Just part, True, Just phase2)
getPartition _ _ _ _ = error patternError
