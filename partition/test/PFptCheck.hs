{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE TemplateHaskell #-}

module PFptCheck (tests) where

import Balancer (balancedPartition)
import Genomes
import GenomesCheck (genMultiPair)
import Hedgehog
import PApprox (approxPartition)
import PFpt
import PGreedy (greedyPartition)
import PSOAR (soarPartition)
import Partition (CommonPartition (..), costUnbalanced, underlineGenome)

-- | On small balanced multichromosome genomes the FPT finishes, produces a
-- common partition of the input genomes, and is never worse than the other
-- algorithms (it is exact). Reverse partitions with signed genomes and direct
-- partitions with unsigned genomes (the FPT does not support unsigned reverse).
fptIsValidAndNotWorseThanOthers :: (Matcher m GenesIRsR GenesIRsF) => m GenesIRsR GenesIRsF -> Sign -> Property
fptIsValidAndNotWorseThanOthers matcher sign =
  withTests 50 . property $ do
    (g, h) <- forAll $ genMultiPair False sign 8
    (maybe_part, full) <- evalIO $ fptPartition 20000000 matcher g h
    assert full
    part <- maybe failure eval maybe_part
    case part of
      CGPbal pg ph -> do
        underlineGenome pg === g
        writeMultiC (writeFGenome False) (underlineGenome ph) === writeMultiC (writeFGenome False) h
      CGPunbal _ _ -> failure
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

tests :: IO Bool
tests = checkSequential $$(discover)
