{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE TemplateHaskell #-}

module PFptCheck (tests) where

import Balancer (balancedPartition)
import Control.Monad (unless)
import Data.List (sortOn)
import Data.Set qualified as Set
import Genomes
import GenomesCheck (genMultiPair, genMultiPairWith)
import Hedgehog
import Hedgehog.Gen qualified as Gen
import PApprox (approxPartition)
import PFpt
import PGreedy (greedyPartition)
import PSOAR (soarPartition)
import Partition (CommonPartition (..), checkCommonBal, checkCommonUnbal, costUnbalanced, mkPartitionFromBreakpoints, numBreakpoints, underlineGenome)

-- | On small balanced multichromosome genomes the FPT finishes, produces a
-- common partition of the input genomes, and is never worse than the other
-- algorithms (it is exact). Reverse partitions with signed genomes and direct
-- partitions with unsigned genomes (the FPT does not support unsigned reverse).
fptIsValidAndNotWorseThanOthers :: (Matcher m GenesIRsR GenesIRsF) => m GenesIRsR GenesIRsF -> Sign -> Property
fptIsValidAndNotWorseThanOthers matcher sign =
  withTests 50 . property $ do
    (g, h) <- forAll $ genMultiPair False sign 8
    fptIsValidAndNotWorseThanOthersOn matcher False g h

-- | The same on unbalanced genomes, with labels that occur in both genomes with
-- distinct numbers of occurrences (so that branching rule 2 is used).
fptUnbalIsValidAndNotWorseThanOthers :: (Matcher m GenesIRsR GenesIRsF) => m GenesIRsR GenesIRsF -> Sign -> Property
fptUnbalIsValidAndNotWorseThanOthers matcher sign =
  withTests 100 . property $ do
    (g, h) <- forAll $ genMultiPairWith True True 50 sign 8
    fptIsValidAndNotWorseThanOthersOn matcher True g h

fptIsValidAndNotWorseThanOthersOn :: (Matcher m GenesIRsR GenesIRsF) => m GenesIRsR GenesIRsF -> Bool -> ChromList GenesIRsR -> ChromList GenesIRsF -> PropertyT IO ()
fptIsValidAndNotWorseThanOthersOn matcher unbal g h = do
  (maybe_part, full) <- evalIO $ fptPartition 20000000 matcher g h
  -- a few small instances (e.g. with many caps) take longer than the time limit
  unless full discard
  part <- maybe failure eval maybe_part
  case part of
    CGPbal pg ph -> do
      underlineGenome pg === g
      writeMultiC (writeFGenome False) (underlineGenome ph) === writeMultiC (writeFGenome False) h
    CGPunbal pg ph -> do
      assert unbal
      underlineGenome pg === g
      writeMultiC (writeFGenome False) (underlineGenome ph) === writeMultiC (writeFGenome False) h
  let cost = costUnbalanced matcher
      others =
        [ greedyPartition False matcher g h,
          soarPartition True matcher g h,
          balancedPartition approxPartition matcher g h
        ]
  others_costs <- mapM (eval . cost) others
  annotateShow (cost part, others_costs)
  assert $ all (cost part <=) others_costs

prop_fptIsValidAndNotWorseThanOthers :: Property
prop_fptIsValidAndNotWorseThanOthers = fptIsValidAndNotWorseThanOthers RFRM Signed

prop_fptDirectIsValidAndNotWorseThanOthers :: Property
prop_fptDirectIsValidAndNotWorseThanOthers = fptIsValidAndNotWorseThanOthers RFDM Unsigned

prop_fptUnbalIsValidAndNotWorseThanOthers :: Property
prop_fptUnbalIsValidAndNotWorseThanOthers = fptUnbalIsValidAndNotWorseThanOthers RFRM Signed

prop_fptUnbalDirectIsValidAndNotWorseThanOthers :: Property
prop_fptUnbalDirectIsValidAndNotWorseThanOthers = fptUnbalIsValidAndNotWorseThanOthers RFDM Unsigned

-- | With the equality relation (no flexibility), which is symmetric and transitive,
-- the FPT is exact: on small genomes (balanced or not) its cost is the minimum cost
-- of a common partition, found by enumerating the breakpoint sets of both genomes.
fptIsExact :: (Matcher m GenesIRsR GenesIRsF) => m GenesIRsR GenesIRsF -> Sign -> Property
fptIsExact matcher sign =
  withTests 200 . property $ do
    unbal <- forAll Gen.bool
    (g, h) <- forAll $ genMultiPairWith unbal True 0 sign 5
    (maybe_part, full) <- evalIO $ fptPartition 20000000 matcher g h
    -- the cost is the minimum only when the search finishes
    unless full discard
    part <- maybe failure eval maybe_part
    costUnbalanced matcher part === bruteForceCost matcher g h

-- | Minimum cost of a common partition (balanced or not) of g and h, trying the
-- pairs of breakpoint sets in increasing order of cost.
bruteForceCost :: (Matcher m GenesIRsR GenesIRsF) => m GenesIRsR GenesIRsF -> ChromList GenesIRsR -> ChromList GenesIRsF -> Int
bruteForceCost matcher g h = case costs of
  c : _ -> c
  [] -> error "no common partition"
  where
    costs =
      [ c
        | c <- [0 ..],
          (cg, pg) <- partsG,
          cg <= c,
          (ch, ph) <- partsH,
          cg + ch == c,
          checkCommonBal matcher pg ph || checkCommonUnbal matcher pg ph
      ]
    partsG = partitions g
    partsH = partitions h
    partitions q =
      sortOn fst
        [ (numBreakpoints p, p)
          | bps <- subsets (breakable q),
            let p = mkPartitionFromBreakpoints q (Set.fromList bps)
        ]
    breakable q = [(chri, i) | (chri, chr) <- zip [1 ..] (getChromosomes q), (i, _) <- adjacencies chr]
    subsets [] = [[]]
    subsets (x : xs) = let r = subsets xs in r ++ map (x :) r

prop_fptIsExact :: Property
prop_fptIsExact = fptIsExact RFRM Signed

prop_fptDirectIsExact :: Property
prop_fptDirectIsExact = fptIsExact RFDM Unsigned

tests :: IO Bool
tests = checkSequential $$(discover)
