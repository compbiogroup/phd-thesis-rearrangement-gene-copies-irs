{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GeneralisedNewtypeDeriving #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE TypeOperators #-}

-- \|
-- Module      : PGreedy
-- Description : Implementation of the algorithm in the SOAR framework for string
-- partition adapted to include intergenic region information.
-- Copyright   : (c) Gabriel Siqueira, 2023
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
module PSOAR
  ( soarPartition,
    suboptimalRuleInterval,
    suboptimalRulePairs,
  )
where

import Control.Applicative ((<|>))
import Data.Bits (Bits (complement, setBit, testBit, zeroBits))
import Data.List qualified as List
import Data.Maybe (catMaybes, fromMaybe, listToMaybe, mapMaybe, maybeToList)
import Data.Set qualified as Set
import Genomes (ChromList, Chromosome (adjacencies, getDuo, nextIdx, prevIdx), GenIdx, Genome (..), Matcher (..), MultiChromosome (..), geneMapLookup, positionMap, singletonOnBoth)
import LocalBase
import Partition (Breakpoints, CommonPartition, combine, mkCommonPartition2)
import Text.Printf (PrintfArg, printf)

-- TODO: check all to consider multiple chromosomes and linear chromosomes
soarPartition :: (Chromosome c1, Chromosome c2, Matcher m c1 c2, GIdx c1 ~ GIdx c2, GenIdx (GIdx c1), Ord (GIdx c1)) => Bool -> m c1 c2 -> ChromList c1 -> ChromList c2 -> CommonPartition c1 c2
soarPartition use_combine matcher g_ h_ = (if use_combine then combine matcher else id) part
  where
    -- (g, h) = uncurry (suboptimalRulePairs matcher) (suboptimalRuleInterval matcher g_ h_)
    -- replace the previous line with this next one to disable the use of
    -- the suboptimal rules
    (g, h) = (g_, h_)
    part = mkCommonPartition2 matcher g_ breaksG h_ breaksH
    graph = makePMGraph4 matcher g h
    (breaksG, breaksH) = getBpsFromIS g h graph . independentSet $ graph

newtype BitMask = BM Integer deriving (Eq, Ord, Bits, PrintfArg)

instance Show BitMask where
  show = printf "%llb"

class PMGraph pmg where
  independentSet :: pmg -> BitMask
  vertexCover :: pmg -> BitMask

  independentSet = complement . vertexCover
  vertexCover = complement . independentSet

type Vertex = Int

type Edge = (Vertex, Vertex)

type PairMatch c1 c2 = (c1, c2, GIdx (ChromList c1), GIdx (ChromList c2), Vertex)

data PMGraph4 c1 c2 = PMGraph4 [Edge] [PairMatch c1 c2]

instance (Show c1, Show c2) => Show (PMGraph4 c1 c2) where
  show (PMGraph4 edges pairMatches) = show edges ++ "\n" ++ show (map (\(c1, c2, _, _, v) -> (c1, c2, v)) pairMatches)

-- Verify conflict between the match due1|due2 (in positions idx1 and idx2) and the match due1'|due2' (in positions idx1' and idx2').
isConflict :: (Matcher m c1 c2, Chromosome c1, Chromosome c2, Eq (GIdx c1), Eq (GIdx c2)) => m c1 c2 -> ChromList c1 -> ChromList c2 -> PairMatch c1 c2 -> PairMatch c1 c2 -> Bool
isConflict matcher g h (due1, due2, idx1, idx2, _) (due1', due2', idx1', idx2', _) =
  idx1 == idx1' && idx2 /= idx2'
    || not overlapG && overlapH
    || followsG idx1 idx1' && not (consistentH due1 due2 idx2 idx2')
    || followsG idx1' idx1 && not (consistentH due1' due2' idx2' idx2)
    -- the same in the other direction, otherwise two duos may form a chain
    -- in G and close a circular chromosome in H
    || followsH idx2 idx2' && not (consistentG due1 due2 idx1 idx1')
    || followsH idx2' idx2 && not (consistentG due1' due2' idx1' idx1)
    -- duos covering a whole circular chromosome of 2 genes in only one genome
    || (followsG idx1 idx1' && followsG idx1' idx1) /= (followsH idx2 idx2' && followsH idx2' idx2)
  where
    -- duo starting at b comes right after duo starting at a (wraps on circular chromosomes)
    followsG (chr_a, a) (chr_b, b) = chr_a == chr_b && nextIdx (getChromosome chr_a g) a == Just b
    followsH (chr_a, a) (chr_b, b) = chr_a == chr_b && nextIdx (getChromosome chr_a h) a == Just b
    overlapG = idx1 == idx1' || followsG idx1 idx1' || followsG idx1' idx1
    overlapH = idx2 == idx2' || followsH idx2 idx2' || followsH idx2' idx2
    -- the duo after a in G must go after (direct) or before (reverse) in H
    consistentH d1 d2 j j' =
      followsH j j' && isDirectMatch matcher d1 d2
        || followsH j' j && isReverseMatch matcher d1 d2
    -- the duo after b in H must go after (direct) or before (reverse) in G
    consistentG d1 d2 i i' =
      followsG i i' && isDirectMatch matcher d1 d2
        || followsG i' i && isReverseMatch matcher d1 d2

makePMGraph4 :: (Chromosome c1, Chromosome c2, Matcher m c1 c2, Ord (GIdx c1), Eq (GIdx c2)) => m c1 c2 -> ChromList c1 -> ChromList c2 -> PMGraph4 c1 c2
makePMGraph4 matcher g h = PMGraph4 edges pms
  where
    edges = catMaybes (testEdge <$> pms <*> pms)
    testEdge pm1@(_, _, _, _, v) pm2@(_, _, _, _, u) =
      if isConflict matcher g h pm1 pm2
        then Just (v, u)
        else Nothing
    pms = zipWith addVertex [1 ..] . catMaybes $ (check <$> duesG <*> duesH)
    addVertex v (due1, due2, idx1, idx2) = (due1, due2, idx1, idx2, v)
    check (idx1, due1) (idx2, due2) =
      if isMatch matcher due1 due2
        then Just (due1, due2, idx1, idx2)
        else Nothing
    duesG = getDuos g
    duesH = getDuos h
    getDuos g' = do
      idx <- allGIdx g'
      let c = getChromosome (getIdxOfChrom g' idx) g'
      maybeToList $ (idx,) <$> getDuo (getIdxInChrom g' idx) c

instance PMGraph (PMGraph4 c1 c2) where
  vertexCover (PMGraph4 edges _) = foldr coverEdges zeroBits edges
    where
      coverEdges :: Edge -> BitMask -> BitMask
      coverEdges (u, v) onCover
        | not (testBit onCover u) && not (testBit onCover v) = setBit (setBit onCover u) v
        | otherwise = onCover

getBpsFromIS :: (Chromosome c1, Chromosome c2, GIdx c1 ~ GIdx c2, Ord (GIdx c1)) => ChromList c1 -> ChromList c2 -> PMGraph4 c1 c2 -> BitMask -> (Breakpoints (ChromList c1), Breakpoints (ChromList c2))
getBpsFromIS g h (PMGraph4 _ pms) indSet = (bpsG, bpsH)
  where
    -- every adjacency is a breakpoint unless its duo was chosen
    -- (the end of a linear chromosome is not an adjacency)
    bpsG = Set.fromList $ adjIdxs g List.\\ getDuos False
    bpsH = Set.fromList $ adjIdxs h List.\\ getDuos True
    adjIdxs gx = [(chr_i, i) | (chr_i, chr) <- zip [1 ..] (getChromosomes gx), (i, _) <- adjacencies chr]
    getDuos isSecond =
      mapMaybe
        ( \(_, _, idx1, idx2, v) ->
            if testBit indSet v
              then Just (if isSecond then idx2 else idx1)
              else Nothing
        )
        pms

-- | Match two correspondent intervals a^l...a^r from g and b^l...b^r from h, where a^l, a^r, b^l and b^r are singletons.
suboptimalRuleInterval :: (Matcher m c1 c2, Chromosome c1, Chromosome c2, GIdx c1 ~ GIdx c2, Ord (GIdx c1)) => m c1 c2 -> ChromList c1 -> ChromList c2 -> (ChromList c1, ChromList c2)
suboptimalRuleInterval matcher g h = foldr checkAndReplace (g, h) is_js
  where
    posMapG = positionMap g
    posMapH = positionMap h
    originallySingleton gidx = singletonOnBoth g posMapG posMapH (getGene gidx g)

    -- Check possible interval matches and replace than with singletons
    checkAndReplace ((chr_i, i), (chr_j, j)) (g', h') =
      case checkDirMatch i j chrGi chrHj chr_i <|> checkRevMatch i j chrGi chrHj chr_i <|> checkMiddleMatch i j chrGi chrHj chr_i of
        Nothing -> (g', h')
        -- bl and br are the left and right ends of the interval in H
        Just (al, ar, bl, br, isRev) ->
          let (g'', singletons_g) = makeSingletons h' ((chr_i,) <$> inner chrGi al ar) g'
              singletons_h =
                if isRev
                  then reverse . map (invGene g') $ singletons_g
                  else singletons_g
              h'' = setGenes ((chr_j,) <$> inner chrHj bl br) singletons_h h'
           in (g'', h'')
      where
        chrGi = getChromosome chr_i g
        chrHj = getChromosome chr_j h

    -- indices strictly between l and r (wraps on circular chromosomes)
    inner chr l r = betweenGIdx chr (step nextIdx chr l) (step prevIdx chr r)
    step f chr = fromMaybe (error logicError) . f chr

    -- candidates for indices of a^l in a^l...a^r
    is = do
      gidx <- allGIdx g
      if originallySingleton gidx then return gidx else []
    -- combine indices in g with the correspondent indices in h
    is_js =
      map
        ( \i -> case geneMapLookup (getGene i g) posMapH of
            Just [j] -> (i, j)
            _ -> error patternError
        )
        is

    -- check if the interval is a direct match (a^l...a^r = b^l...b^r) and
    -- returns indices of (a^l,a^r,b^l,b^r) and the value False indicating that
    -- the match is not reverse.
    checkDirMatch al_i bl_i chr_g chr_h chr_g_i = do
      (ar_i, br_i) <- List.find (originallySingleton . (chr_g_i,) . fst) matches
      -- There is no reason to apply the rule if the singletons are adjacent
      if Just ar_i == nextIdx chr_g al_i
        then Nothing
        else Just (al_i, ar_i, bl_i, br_i, False)
      where
        ps = zip (stepsFrom nextIdx chr_g al_i) (stepsFrom nextIdx chr_h bl_i)
        matches = map snd . takeWhile testMatch $ zip ((al_i, bl_i) : ps) ps
        testMatch ((ap, bp), (ar, br)) = isDirectMatch matcher (subGenome ap ar chr_g) (subGenome bp br chr_h)

    -- check if the interval is a reverse match (a^l...a^r = -b^r...-b^l) and
    -- returns indices of (a^l,a^r,b^l,b^r) and the value True indicating that
    -- the match is reverse.
    checkRevMatch al_i bl_i chr_g chr_h chr_g_i = do
      (ar_i, br_i) <- List.find (originallySingleton . (chr_g_i,) . fst) matches
      -- There is no reason to apply the rule if the singletons are adjacent
      if Just ar_i == nextIdx chr_g al_i
        then Nothing
        else Just (al_i, ar_i, br_i, bl_i, True) -- H ends as (left, right)
      where
        ps = zip (stepsFrom nextIdx chr_g al_i) (stepsFrom prevIdx chr_h bl_i)
        matches = map snd . takeWhile testMatch $ zip ((al_i, bl_i) : ps) ps
        testMatch ((ap, bp), (ar, br)) = isReverseMatch matcher (subGenome ap ar chr_g) (subGenome br bp chr_h)

    -- check if the middle is a reverse match
    -- (a^la^{l+1}...a^{r-1}a^r = b^l-b^{r-1}...-b^{l+1}b^r) and
    -- returns indices of (a^l,a^r,b^l,b^r) and the value True indicating that
    -- the match is reverse.
    checkMiddleMatch al_i bl_i chr_g chr_h chr_g_i = listToMaybe $ do
      -- start at d = 2 because there is no reason to apply the rule if the singletons are adjacent
      (d, ar_i) <- filter (originallySingleton . (chr_g_i,) . snd) . drop 1 . zip [1 :: Int ..] $ ar_is
      br_i <- take 1 . drop (d - 1) $ br_is
      -- a^{l+1}...a^{r-1} and b^{r-1}...b^{l+1}
      let a_ints = take (d - 1) ar_is
          b_ints = reverse . take (d - 1) $ br_is
      [ (al_i, ar_i, bl_i, br_i, True)
        | isDirectMatch matcher (subGenome al_i al_i chr_g) (subGenome bl_i bl_i chr_h),
          isDirectMatch matcher (subGenome ar_i ar_i chr_g) (subGenome br_i br_i chr_h),
          and $ zipWith (\a b -> isReverseMatch matcher (subGenome a a chr_g) (subGenome b b chr_h)) (take 1 a_ints) (take 1 b_ints),
          and $ zipWith testMatch (lPairs a_ints) (lPairs b_ints)
        ]
      where
        ar_is = stepsFrom nextIdx chr_g al_i
        br_is = stepsFrom nextIdx chr_h bl_i
        testMatch (a1, a2) (b2, b1) = isReverseMatch matcher (subGenome a1 a2 chr_g) (subGenome b1 b2 chr_h)

-- | Match pair of a singleton and a replica (a^ha^t). The replica in mapped into a singleton.
suboptimalRulePairs :: (Matcher m c1 c2, Chromosome c1, Chromosome c2, GIdx c1 ~ GIdx c2) => m c1 c2 -> ChromList c1 -> ChromList c2 -> (ChromList c1, ChromList c2)
suboptimalRulePairs matcher g h = foldr checkAndReplace (g, h) (match_ah ++ match_at)
  where
    posMapG = positionMap g
    posMapH = positionMap h
    originallySingleton i = singletonOnBoth g posMapG posMapH (getGene i g)

    -- replace replica of each matched pair with a singleton
    checkAndReplace pairsToCheck (g', h') =
      case pairsToCheck g' h' of
        Nothing -> (g', h')
        Just (a, b, isRev) ->
          let (g'', singletons_g) = makeSingletons h' [a] g'
              singletons_h =
                if isRev
                  then map (invGene g') singletons_g
                  else singletons_g
              h'' = setGenes [b] singletons_h h'
           in (g'', h'')

    -- candidates for indices of singletons
    is = do
      gidx <- allGIdx g
      if originallySingleton gidx then return gidx else []

    -- a^h is a singleton in both genomes and a^t is a replica
    is_ah = filter (\(chr_i, i) -> maybe False (testReplicaG . (chr_i,)) (nextIdx (getChromosome chr_i g) i)) is

    -- find indices from G and H that must be turn into singletons in this case
    match_ah = map (\(gidxi, gidxj) g' h' -> checkPair_ah g' h' gidxi gidxj) $ findInH is_ah

    checkPair_ah g' h' (chr_i, i) (chr_j, j) =
      if
        | Just i' <- nextIdx chrGi i,
          Just j' <- nextIdx chrHj j,
          isDirectMatch matcher (subGenome i i' chrGi) (subGenome j j' chrHj) ->
            Just ((chr_i, i'), (chr_j, j'), False)
        | Just i' <- nextIdx chrGi i,
          Just j' <- prevIdx chrHj j,
          isReverseMatch matcher (subGenome i i' chrGi) (subGenome j' j chrHj) ->
            Just ((chr_i, i'), (chr_j, j'), True)
        | otherwise -> Nothing
      where
        chrGi = getChromosome chr_i g'
        chrHj = getChromosome chr_j h'

    -- a^t is a singleton in both genomes and a^h is a replica
    is_at = filter (\(chr_i, i) -> maybe False (testReplicaG . (chr_i,)) (prevIdx (getChromosome chr_i g) i)) is

    -- find indices from G and H that must be turn into singletons in this case
    match_at = map (\(gidxi, gidxj) g' h' -> checkPair_at g' h' gidxi gidxj) $ findInH is_at

    checkPair_at g' h' (chr_i, i) (chr_j, j) =
      if
        | Just i' <- prevIdx chrGi i,
          Just j' <- nextIdx chrHj j,
          isReverseMatch matcher (subGenome i' i chrGi) (subGenome j j' chrHj) ->
            Just ((chr_i, i'), (chr_j, j'), True)
        | Just i' <- prevIdx chrGi i,
          Just j' <- prevIdx chrHj j,
          isDirectMatch matcher (subGenome i' i chrGi) (subGenome j' j chrHj) ->
            Just ((chr_i, i'), (chr_j, j'), False)
        | otherwise -> Nothing
      where
        chrGi = getChromosome chr_i g'
        chrHj = getChromosome chr_j h'

    -- combine is with the correspondent indices in h
    findInH =
      map
        ( \i -> case geneMapLookup (getGene i g) posMapH of
            Just [j] -> (i, j)
            _ -> error patternError
        )

    -- test if gene at position i is a replica in the genome g
    testReplicaG i =
      case geneMapLookup (getGene i g) posMapG of
        Nothing -> False
        Just pos -> length pos > 1

-- | Indices visited walking from i with step (nextIdx or prevIdx), at most size - 1 steps
stepsFrom :: (Chromosome c) => (c -> GIdx c -> Maybe (GIdx c)) -> c -> GIdx c -> [GIdx c]
stepsFrom step chr i = List.unfoldr go (size chr - 1, i)
  where
    go (k, x) = if k <= 0 then Nothing else (\y -> (y, (k - 1, y))) <$> step chr x
