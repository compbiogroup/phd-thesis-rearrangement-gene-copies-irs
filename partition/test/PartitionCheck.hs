{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE TemplateHaskell #-}

module PartitionCheck (tests) where

import Data.IntMap qualified as IntMap
import Data.List qualified as List
import Data.Maybe (fromJust)
import Data.Set qualified as Set
import Genomes
import Balancer (balancedPartition)
import GenomesCheck (genMultiPair, genRGenome, rearrangeAndFlexibilizeGenome)
import Hedgehog
import Hedgehog.Gen qualified as Gen
import Hedgehog.Range qualified as Range
import PApprox
import PGreedy
import PSOAR
import Partition

-- | Genome with a single linear chromosome (its first and last genes are the caps)
toMCF :: GenesIRsF -> ChromList GenesIRsF
toMCF h = ChromList [getGene 1 h, getGene (mkIdx (size h)) h] [h]

-- | Pair of genomes: a random genome and k random reversals of it, made flexible
genPair :: Int -> PropertyT IO (ChromList GenesIRsR, ChromList GenesIRsF)
genPair size_lim = do
  g <- forAll (genRGenome size_lim)
  k <- forAll $ Gen.int (Range.linear 0 (size g))
  h <- forAll $ rearrangeAndFlexibilizeGenome k g
  return (toMC g, toMCF h)

-- | Block of a genome between two positions of the same chromosome
subChr :: (Chromosome c) => ChromList c -> GIdx (ChromList c) -> GIdx (ChromList c) -> c
subChr q (c, i) (_, j) = subGenome i j (getChromosome c q)

prop_partitionConstructionMethodsAreEquivalent :: Property
prop_partitionConstructionMethodsAreEquivalent =
  property $ do
    chr <- forAll (genRGenome 100)
    bps <- forAll (Gen.subsequence [1 .. mkIdx (size chr - 1)])
    let g = toMC chr
        bps' = 0 : bps ++ [mkIdx (size chr)]
        bs = zipWith (\i succi -> (1, subGenome (incIdx i) succi chr)) bps' $ tail bps'
    assert $ mkPartitionFromBreakpoints g (Set.fromList [(1, i) | i <- bps]) == mkPartitionFromBlocks g bs

prop_commonPrefixFromBothAreEqual :: Property
prop_commonPrefixFromBothAreEqual =
  property $ do
    (g, h) <- genPair 100
    idx_g <- forAll $ Gen.element (allGIdx g)
    idx_h <- forAll . Gen.element . fromJust $ geneMapLookup (getGene idx_g g) (positionMap h)
    let ((g_beg, g_end), (h_beg, h_end), _) = fromJust $ commonPrefix Nothing RFRM Set.empty Set.empty g h idx_g idx_h
    subG <- forAll . return $ subChr g g_beg g_end
    subH <- forAll . return $ subChr h h_beg h_end
    assert $ isMatch RFRM subG subH

prop_longestSubstringFromBothAreEqual :: Property
prop_longestSubstringFromBothAreEqual =
  property $ do
    (g, h) <- genPair 100
    let (((g_beg, g_end), (h_beg, h_end)), _, _) = fromJust $ longestSubstring Nothing RFRM Set.empty Set.empty g h
    subG <- forAll . return $ subChr g g_beg g_end
    subH <- forAll . return $ subChr h h_beg h_end
    assert $ isMatch RFRM subG subH

type Rule = RigidFlexibleReverseMatcher GenesIRsR GenesIRsF -> ChromList GenesIRsR -> ChromList GenesIRsF -> (ChromList GenesIRsR, ChromList GenesIRsF)

prop_suboptimalRuleIntervalTurnsReplicasIntoSingletons :: Property
prop_suboptimalRuleIntervalTurnsReplicasIntoSingletons = suboptimalRuleTurnsReplicasIntoSingletons suboptimalRuleInterval

prop_suboptimalRulePairsTurnsReplicasIntoSingletons :: Property
prop_suboptimalRulePairsTurnsReplicasIntoSingletons = suboptimalRuleTurnsReplicasIntoSingletons suboptimalRulePairs

suboptimalRuleTurnsReplicasIntoSingletons :: Rule -> Property
suboptimalRuleTurnsReplicasIntoSingletons rule =
  property $ do
    (g, h) <- genPair 100
    (g', h') <- forAll . return $ rule RFRM g h
    assert $ getSingletons g `List.isSubsequenceOf` getSingletons g'
    assert $ getSingletons h `List.isSubsequenceOf` getSingletons h'
  where
    -- positions of the singletons of a genome
    getSingletons :: (Genome q, Ord (GIdx q)) => q -> [GIdx q]
    getSingletons q = List.sort $ foldr (\pos acc -> case pos of [p] -> p : acc; _ -> acc) [] (positionMap q)

prop_suboptimalRuleIntervalKeepBalancedGenomes :: Property
prop_suboptimalRuleIntervalKeepBalancedGenomes = suboptimalRuleKeepBalancedGenomes suboptimalRuleInterval

prop_suboptimalRulePairsKeepBalancedGenomes :: Property
prop_suboptimalRulePairsKeepBalancedGenomes = suboptimalRuleKeepBalancedGenomes suboptimalRulePairs

suboptimalRuleKeepBalancedGenomes :: Rule -> Property
suboptimalRuleKeepBalancedGenomes rule =
  property $ do
    (g, h) <- genPair 100
    (g', h') <- forAll . return $ rule RFRM g h
    assert $ areBalanced RFRM (getChromosome 1 g') (getChromosome 1 h')

prop_equalGenomesHaveOnlyPerfectComponentsOnBMG :: Property
prop_equalGenomesHaveOnlyPerfectComponentsOnBMG =
  property $ do
    g <- toMC <$> forAll (genRGenome 100)
    let pg = trivialPartition g
    let (comps, _) = getConnectedComponents (mkBlockMatchGraph RRRM pg pg)
    assert $ countPerfect comps == IntMap.size comps

prop_greedyPartitionProduceValidCorrespondence :: Property
prop_greedyPartitionProduceValidCorrespondence =
  partitionProduceValidCorrespondence (greedyPartition False) 100

prop_greedyPartitionSingProduceValidCorrespondence :: Property
prop_greedyPartitionSingProduceValidCorrespondence =
  partitionProduceValidCorrespondence (greedyPartition True) 100

prop_soarPartitionProduceValidCorrespondence :: Property
prop_soarPartitionProduceValidCorrespondence =
  partitionProduceValidCorrespondence (soarPartition False) 100

prop_approxPartitionProduceValidCorrespondence :: Property
prop_approxPartitionProduceValidCorrespondence =
  partitionProduceValidCorrespondence approxPartition 50

partitionProduceValidCorrespondence :: (RigidFlexibleReverseMatcher GenesIRsR GenesIRsF -> ChromList GenesIRsR -> ChromList GenesIRsF -> CommonPartition GenesIRsR GenesIRsF) -> Int -> Property
partitionProduceValidCorrespondence partAlg size_lim =
  property $ do
    (g, h) <- genPair size_lim
    part <- forAll . return $ partAlg RFRM g h
    -- the genomes are balanced, so the partition (also after combine) must be balanced
    case Partition.combine RFRM part of
      CGPbal pg ph -> assert $ checkCommonBal RFRM pg ph
      CGPunbal _ _ -> failure

-- | combine with many repeated genes (alphabet of 4 genes) must keep the
-- partition common
prop_combineKeepsCommonPartitionWithReplicas :: Property
prop_combineKeepsCommonPartitionWithReplicas =
  property $ do
    n <- forAll $ Gen.int (Range.linear 3 30)
    genes <- forAll $ Gen.list (Range.singleton n) (Gen.element [-4, -3, -2, -1, 1, 2, 3, 4])
    irs <- forAll $ Gen.list (Range.singleton (n + 1)) (Gen.int (Range.linear 0 10))
    let g = mkRGenome Linear (Just (0, maxBound)) Signed genes irs
    k <- forAll $ Gen.int (Range.linear 0 (size g))
    h <- forAll $ rearrangeAndFlexibilizeGenome k g
    alg <- forAll $ Gen.element [False, True]
    let part = (if alg then soarPartition False else greedyPartition False) RFRM (toMC g) (toMCF h)
    case Partition.combine RFRM part of
      CGPbal pg ph -> assert $ checkCommonBal RFRM pg ph
      CGPunbal _ _ -> failure

-- | Every algorithm, as called by Main, with the given matcher
runAlg :: (Matcher m GenesIRsR GenesIRsF) => String -> m GenesIRsR GenesIRsF -> ChromList GenesIRsR -> ChromList GenesIRsF -> CommonPartition GenesIRsR GenesIRsF
runAlg name = case name of
  "Greedy" -> greedyPartition False
  "GreedySin" -> greedyPartition True
  "SOAR" -> soarPartition False
  "SOARComb" -> soarPartition True
  "Approx" -> balancedPartition approxPartition
  _ -> error ("unknown algorithm " ++ name)

-- | On multichromosome genomes (linear and circular chromosomes, many
-- replicas, balanced or not) every algorithm must produce a common partition
-- (mkCommonPartition fails otherwise) of the input genomes, which is balanced
-- when the genomes are balanced. Reverse partitions (PRCMG) are tested with
-- signed and unsigned genomes, direct partitions (PDCMG) with unsigned genomes.
multiPartitionIsValid :: Bool -> String -> Property
multiPartitionIsValid direct name =
  withTests 300 . property $ do
    unbal <- forAll Gen.bool
    sign <- if direct then return Unsigned else forAll $ Gen.element [Signed, Unsigned]
    (g, h) <- forAll $ genMultiPair unbal sign 12
    part <- eval $ if direct then runAlg name RFDM g h else runAlg name RFRM g h
    case part of
      CGPbal pg ph -> do
        underlineGenome pg === g
        writeMultiC (writeFGenome False) (underlineGenome ph) === writeMultiC (writeFGenome False) h
      CGPunbal pg ph -> do
        assert unbal
        underlineGenome pg === g
        writeMultiC (writeFGenome False) (underlineGenome ph) === writeMultiC (writeFGenome False) h

prop_greedyMultiPartitionIsValid :: Property
prop_greedyMultiPartitionIsValid = multiPartitionIsValid False "Greedy"

prop_greedySinMultiPartitionIsValid :: Property
prop_greedySinMultiPartitionIsValid = multiPartitionIsValid False "GreedySin"

prop_soarMultiPartitionIsValid :: Property
prop_soarMultiPartitionIsValid = multiPartitionIsValid False "SOAR"

prop_soarCombMultiPartitionIsValid :: Property
prop_soarCombMultiPartitionIsValid = multiPartitionIsValid False "SOARComb"

prop_approxMultiPartitionIsValid :: Property
prop_approxMultiPartitionIsValid = multiPartitionIsValid False "Approx"

prop_greedyDirectMultiPartitionIsValid :: Property
prop_greedyDirectMultiPartitionIsValid = multiPartitionIsValid True "Greedy"

prop_greedySinDirectMultiPartitionIsValid :: Property
prop_greedySinDirectMultiPartitionIsValid = multiPartitionIsValid True "GreedySin"

prop_soarDirectMultiPartitionIsValid :: Property
prop_soarDirectMultiPartitionIsValid = multiPartitionIsValid True "SOAR"

prop_soarCombDirectMultiPartitionIsValid :: Property
prop_soarCombDirectMultiPartitionIsValid = multiPartitionIsValid True "SOARComb"

prop_approxDirectMultiPartitionIsValid :: Property
prop_approxDirectMultiPartitionIsValid = multiPartitionIsValid True "Approx"

-- prop_fptPartitionProduceValidCorrespondence :: Property
-- prop_fptPartitionProduceValidCorrespondence =
--   property $ do
--     (g, h) <- genPair 20
--     (part, _) <- fmap (first fromJust) . evalIO $ fptPartition 100000000 RFRM g h
--     case Partition.combine RFRM part of
--       CGPbal pg ph -> assert $ checkCommonBal RFRM pg ph
--       CGPunbal _ _ -> failure

tests :: IO Bool
tests = checkSequential $$(discover)
