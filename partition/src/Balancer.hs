{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE TypeOperators #-}

-- |
-- Module      : Balancer
-- Description : Heuristic to turn a pair of unbalanced genomes into a pair of balanced
-- genomes, so that a partition algorithm for balanced genomes can be used. Each gene
-- that occurs more times in one genome than in the other is removed from it,
-- choosing the largest runs of consecutive removable genes first. The removed runs
-- are the exclusive blocks of the final partition.
-- Copyright   : (c) Gabriel Siqueira, 2023
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
module Balancer
  ( balancedPartition,
    balancedPartitionWith,
    balancedLowerBound,
    balanceGenomes,
    Balanced (..),
    markExclusive,
    occBalance,
  )
where

import Data.Function (on)
import Data.IntMap (IntMap)
import Data.IntMap qualified as IntMap
import Data.List qualified as List
import Data.Map qualified as Map
import Data.Maybe (listToMaybe, maybeToList)
import Data.Set (Set)
import Data.Set qualified as Set
import Genomes (CIdx, ChromList (..), Chromosome (isCircular, prevIdx), GenIdx, Gene, Genome (..), Matcher (allowsReverse), MultiChromosome (getChromosome, getChromosomes), geneToInt)
import LocalBase
import Partition (CommonPartition (..), breakpoints, mkCommonPartition2)

-- | Key used to count the copies of a gene: its label (ignoring orientation) or,
-- if signedKey is True, the gene with its sign. A direct common partition only
-- matches genes with the same sign, so its genes in excess are counted by signed
-- gene; in a reverse one a block may be reversed, so they are counted by label.
geneKey :: Bool -> Gene -> Int
geneKey signedKey = geneToInt . (if signedKey then id else abs)

-- | Occurrences of each gene in g minus its occurrences in h (ignoring orientation)
occBalance :: (Genome g1, Genome g2) => g1 -> g2 -> IntMap Int
occBalance = occBalanceBy False

-- | Occurrences of each key (see geneKey) in g minus its occurrences in h
occBalanceBy :: (Genome g1, Genome g2) => Bool -> g1 -> g2 -> IntMap Int
occBalanceBy signedKey g h = IntMap.unionWith (+) (count g) (negate <$> count h)
  where
    count q = IntMap.fromListWith (+) [(geneKey signedKey (getGene i q), 1) | i <- allGIdx q]

-- | Select the genes of q that must be removed to balance it. Genes in excess in q are
-- the ones with sign * balance > 0 (sign is 1 for the first genome and -1 for the
-- second). While there are genes in excess, the largest run of consecutive genes in
-- excess (that were not selected yet) is selected (the leftmost one on ties), each
-- selected gene decreases the excess of its value by one. Runs wrap around the origin
-- of circular chromosomes. Returns the selected positions and the updated balance.
markExclusive :: (Chromosome c, Ord (GIdx c)) => Bool -> Int -> ChromList c -> IntMap Int -> (Set (CIdx, GIdx c), IntMap Int)
markExclusive signedKey sign q = go Set.empty
  where
    key idx = geneKey signedKey (getGene idx q)
    inExcess marked occ idx = idx `Set.notMember` marked && sign * IntMap.findWithDefault 0 (key idx) occ > 0

    go marked occ =
      case listToMaybe . sortWith (negate . length) $ concatMap (chromRuns marked occ) (zip [1 ..] (getChromosomes q)) of
        Nothing -> (marked, occ)
        Just run -> go (foldr Set.insert marked run) (foldr (IntMap.adjust (subtract sign) . key) occ run)

    -- Runs of consecutive genes in excess of a chromosome. The excess is updated
    -- inside a run, so a run does not take more copies of a gene than its excess.
    chromRuns marked occ (c, chr) = finish $ foldl step ([], [], occ) order
      where
        idxs = (c,) <$> allGIdx chr
        -- on a circular chromosome, start after a gene that is not in excess (if there
        -- is one), so that a run may cross the origin
        order = case break (not . inExcess marked occ) idxs of
          (pre, post@(_ : _)) | isCircular chr -> post ++ pre
          _ -> idxs
        step (runs, cur, occ_cur) idx
          | inExcess marked occ_cur idx = (runs, idx : cur, IntMap.adjust (subtract sign) (key idx) occ_cur)
          | otherwise = (close runs cur, [], occ)
        close runs cur = if null cur then runs else reverse cur : runs
        finish (runs, cur, _) = reverse (close runs cur)

-- | Maximal segments of consecutive indices of chr with the same value of p (in circular
-- order on circular chromosomes), each with that value
segments :: (Chromosome c) => (GIdx c -> Bool) -> c -> [(Bool, [GIdx c])]
segments p chr = map (\s -> (p (localHead s), s)) groups
  where
    groups = case List.groupBy ((==) `on` p) (allGIdx chr) of
      (first : rest@(_ : _))
        | isCircular chr && p (localHead first) == p (localHead (last rest)) -> (last rest ++ first) : init rest
      gs -> gs

-- | Remove the marked positions of q. Each chromosome with marked positions is split
-- into its maximal unmarked segments (as linear chromosomes). Returns the reduced
-- genome, a map from positions of the reduced genome to positions of q, and the
-- breakpoints of q that separate the removed segments from the rest.
reduce :: (Chromosome c, Ord (GIdx c)) => Set (CIdx, GIdx c) -> ChromList c -> (ChromList c, Map.Map (CIdx, GIdx c) (CIdx, GIdx c), Set (CIdx, GIdx c))
reduce marks q = (ChromList (caps q) (map fst pieces), posMap, cuts)
  where
    perChrom = [(c, chr, segments (\i -> (c, i) `Set.member` marks) chr) | (c, chr) <- zip [1 ..] (getChromosomes q)]
    pieces =
      concat
        [ if not (any fst segs)
            then [(chr, (c,) <$> allGIdx chr)]
            else [(subGenome (localHead s) (last s) chr, (c,) <$> s) | (False, s) <- segs]
          | (c, chr, segs) <- perChrom
        ]
    posMap = Map.fromList [((r, i), orig) | (r, (piece, origs)) <- zip [1 ..] pieces, (i, orig) <- zip (allGIdx piece) origs]
    cuts =
      Set.fromList
        [ (c, i)
          | (c, chr, segs) <- perChrom,
            length segs > 1,
            (True, s) <- segs,
            i <- maybeToList (prevIdx chr (localHead s)) ++ [last s]
        ]

-- | Result of balancing two genomes: the reduced genomes, the maps from their positions
-- to positions of the original genomes, and the breakpoints of the original genomes
-- around the removed segments.
data Balanced c1 c2 = Balanced
  { bal_g :: ChromList c1,
    bal_h :: ChromList c2,
    bal_posMapG :: Map.Map (CIdx, GIdx c1) (CIdx, GIdx c1),
    bal_posMapH :: Map.Map (CIdx, GIdx c2) (CIdx, GIdx c2),
    bal_cutsG :: Set (CIdx, GIdx c1),
    bal_cutsH :: Set (CIdx, GIdx c2)
  }

-- | Remove the genes in excess of both genomes (selected with markExclusive, with
-- the genes counted as in geneKey). Returns Nothing if the genomes are already balanced.
balanceGenomes :: (Chromosome c1, Chromosome c2, Ord (GIdx c1), Ord (GIdx c2)) => Bool -> ChromList c1 -> ChromList c2 -> Maybe (Balanced c1 c2)
balanceGenomes signedKey g h
  | Set.null marksG && Set.null marksH = Nothing
  | otherwise = Just (Balanced g' h' posMapG posMapH cutsG cutsH)
  where
    (marksG, occG) = markExclusive signedKey 1 g (occBalanceBy signedKey g h)
    (marksH, _) = markExclusive signedKey (-1) h occG
    (g', posMapG, cutsG) = reduce marksG g
    (h', posMapH, cutsH) = reduce marksH h

-- | Apply a partition algorithm for balanced genomes to genomes that may be unbalanced.
-- The genes in excess are selected with markExclusive and removed, the algorithm is
-- applied to the remaining (balanced) genomes, and its partition is mapped back to
-- the original genomes, where the removed segments are the exclusive blocks.
balancedPartition ::
  (Matcher m c1 c2, Chromosome c1, Chromosome c2, GIdx c1 ~ GIdx c2, GenIdx (GIdx c1), Ord (GIdx c1)) =>
  (m c1 c2 -> ChromList c1 -> ChromList c2 -> CommonPartition c1 c2) ->
  m c1 c2 ->
  ChromList c1 ->
  ChromList c2 ->
  CommonPartition c1 c2
balancedPartition partAlg matcher g h = fst $ balancedPartitionWith () (\m g' h' -> (partAlg m g' h', ())) matcher g h

-- | As balancedPartition, for an algorithm that also returns a value (for example, a
-- count of some step). The value is the one of the algorithm on the balanced genomes, or
-- def if the algorithm is not run (the balanced genomes are empty).
balancedPartitionWith ::
  (Matcher m c1 c2, Chromosome c1, Chromosome c2, GIdx c1 ~ GIdx c2, GenIdx (GIdx c1), Ord (GIdx c1)) =>
  a ->
  (m c1 c2 -> ChromList c1 -> ChromList c2 -> (CommonPartition c1 c2, a)) ->
  m c1 c2 ->
  ChromList c1 ->
  ChromList c2 ->
  (CommonPartition c1 c2, a)
balancedPartitionWith def partAlg matcher g h =
  case balanceGenomes (not (allowsReverse matcher)) g h of
    Nothing -> partAlg matcher g h
    Just bal ->
      let -- after the removal both genomes have the same genes, so either both are empty or none is
          ((bpsG', bpsH'), x)
            | null (chromosomes (bal_g bal)) = ((Set.empty, Set.empty), def)
            | otherwise = case partAlg matcher (bal_g bal) (bal_h bal) of
                (CGPbal pg ph, x') -> ((breakpoints pg, breakpoints ph), x')
                (CGPunbal pg ph, x') -> ((breakpoints pg, breakpoints ph), x')
          bpsG = Set.union (bal_cutsG bal) (Set.map (bal_posMapG bal Map.!) bpsG')
          bpsH = Set.union (bal_cutsH bal) (Set.map (bal_posMapH bal Map.!) bpsH')
       in (mkCommonPartition2 matcher g bpsG h bpsH, x)

-- | Apply a lower bound for balanced genomes (on the number of breakpoints of both
-- partitions) to genomes that may be unbalanced. The bound is computed on the
-- balanced genomes and the breakpoints around the removed segments are added to it.
-- Note that this bounds the partitions whose exclusive blocks are the segments chosen
-- by markExclusive, not every partition of the unbalanced genomes.
balancedLowerBound ::
  (Matcher m c1 c2, Chromosome c1, Chromosome c2, Ord (GIdx c1), Ord (GIdx c2)) =>
  (m c1 c2 -> ChromList c1 -> ChromList c2 -> Int) ->
  m c1 c2 ->
  ChromList c1 ->
  ChromList c2 ->
  Int
balancedLowerBound lbAlg matcher g h =
  case balanceGenomes (not (allowsReverse matcher)) g h of
    Nothing -> lbAlg matcher g h
    Just bal ->
      let lb' = if null (chromosomes (bal_g bal)) then 0 else lbAlg matcher (bal_g bal) (bal_h bal)
       in lb' + realBreaks g (bal_cutsG bal) + realBreaks h (bal_cutsH bal)
  where
    -- the end of a linear chromosome is not a breakpoint
    realBreaks q = length . filter (\(c, i) -> let chr = getChromosome c q in isCircular chr || i /= lastGIdx chr) . Set.toList
