{-# LANGUAGE ExplicitForAll #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE TypeOperators #-}

-- \|
-- Module      : PGreedy
-- Description : Implementation of the greedy algorithm for string partition adapted
-- to include intergenic region information.
-- Copyright   : (c) Gabriel Siqueira, 2023
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
module PGreedy
  ( greedyPartition,
    longestSubstring,
    commonPrefix,
  )
where

import Data.Foldable (toList)
import Data.Maybe (fromMaybe, mapMaybe, maybeToList)
import Data.Set (Set)
import Data.Set qualified as Set
import Genomes (ChromList, Chromosome (isCircular, nextIdx, prevIdx), GenIdx, Gene, Genome (GIdx, allGIdx, alphabet, betweenGIdx, isGene, lastGIdx, size, subGenome), Matcher (isDirectMatch, isMatch, isReverseMatch), MultiChromosome (getChromosome, getIdxOfChrom, setIdx), getIdxInChrom, positionMap, singletonOnBoth)
import LocalBase
import Partition (CommonPartition, mkCommonPartition2)

-- | Greedy algorithm for string partition with intergenic regions or
-- flexible intergeic regions. If withSingleton is True it looks first
-- for matched blocks that contain a singleton.
greedyPartition :: (Chromosome c1, Chromosome c2, Matcher m c1 c2, GIdx c1 ~ GIdx c2, GenIdx (GIdx c1), Ord (GIdx c1)) => Bool -> m c1 c2 -> ChromList c1 -> ChromList c2 -> CommonPartition c1 c2
greedyPartition withSingleton matcher g h = mkCommonPartition2 matcher g bg h bh
  where
    bg = removeRedundantBreaks final_breaksG g
    bh = removeRedundantBreaks final_breaksH h

    removeRedundantBreaks breaks g' = Set.filter removeRedundantBreak breaks
      where
        removeRedundantBreak idx =
          let (chr_i, i) = (getIdxOfChrom g' idx, getIdxInChrom g' idx)
              chr = getChromosome chr_i g'
           in isCircular chr || i /= lastGIdx chr

    (final_breaksG, final_breaksH) = getAllLS all_singletons Set.empty Set.empty Set.empty Set.empty
    getAllLS singletons genesOnBlocksG genesOnBlocksH breaksG breaksH =
      case longestSubstring' singletons of
        (Nothing, Nothing) -> (breaksG, breaksH)
        (Just other_singletons, Nothing) -> getAllLS other_singletons genesOnBlocksG genesOnBlocksH breaksG breaksH
        (maybe_other_singletons, Just (((g_beg, g_end), (h_beg, h_end)), genesOnBlocksG', genesOnBlocksH', closed))
          -- a whole circular chromosome matched (closing adjacency included) has no breakpoint
          | closed -> getAllLS (fromMaybe [] maybe_other_singletons) genesOnBlocksG' genesOnBlocksH' breaksG breaksH
          | otherwise ->
          let breaksG' = maybe id Set.insert (prevChromIdx g g_beg) (Set.insert g_end breaksG)
              breaksH' = maybe id Set.insert (prevChromIdx h h_beg) (Set.insert h_end breaksH)
           in getAllLS (fromMaybe [] maybe_other_singletons) genesOnBlocksG' genesOnBlocksH' breaksG' breaksH'
      where
        longestSubstring' (singleton : other_singletons) = (Just other_singletons, longestSubstringC (Just singleton) matcher genesOnBlocksG genesOnBlocksH g h)
        longestSubstring' [] = (Nothing, longestSubstringC Nothing matcher genesOnBlocksG genesOnBlocksH g h)

    -- index before i in its chromosome (wraps on circular, Nothing at the start of a linear one)
    prevChromIdx gx i = setIdx gx chr_i <$> prevIdx (getChromosome chr_i gx) (getIdxInChrom gx i)
      where
        chr_i = getIdxOfChrom gx i
    all_singletons =
      if withSingleton
        then filter (singletonOnBoth g posMapG posMapH) . toList $ alphabet g
        else []
    posMapG = positionMap g
    posMapH = positionMap h

longestSubstring :: forall m c1 c2. (Chromosome c1, Chromosome c2, Matcher m c1 c2, Ord (GIdx c1), Ord (GIdx c2)) => Maybe Gene -> m c1 c2 -> Set (GIdx (ChromList c1)) -> Set (GIdx (ChromList c2)) -> ChromList c1 -> ChromList c2 -> Maybe (((GIdx (ChromList c1), GIdx (ChromList c1)), (GIdx (ChromList c2), GIdx (ChromList c2))), Set (GIdx (ChromList c1)), Set (GIdx (ChromList c2)))
-- Find longest substring of l1 and l2, ignore genes already in a block.
-- Returns:
-- - the indices of the correspondent subgenomes in g and h (chromosome, first and last gene)
-- -
-- If singleton has a value (Just x), only consider blocks containing that gene x.
longestSubstring maybe_singleton matcher genesOnBlocksG genesOnBlocksH g h =
  (\(inds, gsG, gsH, _) -> (inds, gsG, gsH)) <$> longestSubstringC maybe_singleton matcher genesOnBlocksG genesOnBlocksH g h

-- | Same as longestSubstring, but also returns whether the blocks are whole
-- circular chromosomes whose closing adjacencies also match.
longestSubstringC :: forall m c1 c2. (Chromosome c1, Chromosome c2, Matcher m c1 c2, Ord (GIdx c1), Ord (GIdx c2)) => Maybe Gene -> m c1 c2 -> Set (GIdx (ChromList c1)) -> Set (GIdx (ChromList c2)) -> ChromList c1 -> ChromList c2 -> Maybe (((GIdx (ChromList c1), GIdx (ChromList c1)), (GIdx (ChromList c2), GIdx (ChromList c2))), Set (GIdx (ChromList c1)), Set (GIdx (ChromList c2)), Bool)
longestSubstringC maybe_singleton matcher genesOnBlocksG genesOnBlocksH g h =
  case maybe_inds of
    Nothing -> Nothing
    Just ((g_beg, g_end), (h_beg, h_end), _, closed) ->
      Just
        ( ((g_beg, g_end), (h_beg, h_end)),
          foldr Set.insert genesOnBlocksG (betweenGIdx g g_beg g_end),
          foldr Set.insert genesOnBlocksH (betweenGIdx h h_beg h_end),
          closed
        )
  where
    -- on a tie in length, prefer a match that keeps a whole circular chromosome
    -- (closed), which saves the breakpoint of that chromosome
    maybe_inds =
      maxWith (\(_, _, d, c) -> (d, c))
        . mapMaybe longestMatch
        . filter (`Set.notMember` genesOnBlocksG)
        $ allGIdx g

    -- find longest match between a subgenome of h and a suffix of g,
    -- ignore genes already in a block,
    -- returns the indices of the correspondent subgenomes in g and h (first and last gene)
    longestMatch ig = maxWith (\(_, _, d, c) -> (d, c)) $ do
      ih <- filter (`Set.notMember` genesOnBlocksH) $ allGIdx h
      maybeToList $ commonPrefixC maybe_singleton matcher genesOnBlocksG genesOnBlocksH g h ig ih

commonPrefix :: (Matcher m c1 c2, Chromosome c1, Chromosome c2, Ord (GIdx c1), Ord (GIdx c2)) => Maybe Gene -> m c1 c2 -> Set (GIdx (ChromList c1)) -> Set (GIdx (ChromList c2)) -> ChromList c1 -> ChromList c2 -> GIdx (ChromList c1) -> GIdx (ChromList c2) -> Maybe ((GIdx (ChromList c1), GIdx (ChromList c1)), (GIdx (ChromList c2), GIdx (ChromList c2)), Int)
-- find common prefix of subgenome of g starting with gidx_g
-- and subgenome of h starting with gidx_h
-- or if the elements is gidx_g and gidx_h are a reversed match
-- find common prefix of subgenome of g starting with gidx_g
-- and subgenome of rev(h) starting with gidx_h
-- returns Nothing if the prefix is empty
-- returns indices for the ends of both subgenomes otherwise
commonPrefix maybe_singleton matcher genesOnBlocksG genesOnBlocksH g h gidx_g gidx_h =
  (\(inds_g, inds_h, d, _) -> (inds_g, inds_h, d)) <$> commonPrefixC maybe_singleton matcher genesOnBlocksG genesOnBlocksH g h gidx_g gidx_h

-- | Same as commonPrefix, but also returns whether the prefixes are whole
-- circular chromosomes whose closing adjacencies also match.
commonPrefixC :: (Matcher m c1 c2, Chromosome c1, Chromosome c2, Ord (GIdx c1), Ord (GIdx c2)) => Maybe Gene -> m c1 c2 -> Set (GIdx (ChromList c1)) -> Set (GIdx (ChromList c2)) -> ChromList c1 -> ChromList c2 -> GIdx (ChromList c1) -> GIdx (ChromList c2) -> Maybe ((GIdx (ChromList c1), GIdx (ChromList c1)), (GIdx (ChromList c2), GIdx (ChromList c2)), Int, Bool)
commonPrefixC maybe_singleton matcher genesOnBlocksG genesOnBlocksH g h gidx_g gidx_h = do
  -- when the first genes match in both orientations, both prefixes are computed
  -- and the longest is kept (on a tie, the direct one, as maximumBy keeps the last)
  (endCPG, endCPH, rev, d, closed) <-
    maxWith (\(_, _, _, d', c) -> (d', c)) $
      [commonPrefix' True 1 ig ih | isReverseMatch matcher gene_g gene_h]
        ++ [commonPrefix' False 1 ig ih | isDirectMatch matcher gene_g gene_h]
  case maybe_singleton of
    Just singleton ->
      if singleton `isGene` subGenome gidx_g endCPG g
        then return ((gidx_g, endCPG), if rev then (endCPH, gidx_h) else (gidx_h, endCPH), d, closed)
        else Nothing
    Nothing -> return ((gidx_g, endCPG), if rev then (endCPH, gidx_h) else (gidx_h, endCPH), d, closed)
  where
    gene_g = subGenome ig ig chr_g
    gene_h = subGenome ih ih chr_h
    (chr_ig, ig) = (getIdxOfChrom g gidx_g, getIdxInChrom g gidx_g)
    (chr_ih, ih) = (getIdxOfChrom h gidx_h, getIdxInChrom h gidx_h)
    chr_g = getChromosome chr_ig g
    chr_h = getChromosome chr_ih h
    -- k is the length of the prefix so far, it can not be longer than the chromosomes
    commonPrefix' rev k jg jh =
      case (nextIdx chr_g jg, if rev then prevIdx chr_h jh else nextIdx chr_h jh) of
        (Just jg', Just jh')
          | k < size chr_g
              && k < size chr_h
              && (chr_ig, jg') `Set.notMember` genesOnBlocksG
              && (chr_ih, jh') `Set.notMember` genesOnBlocksH
              && testMatch matcher (subGenome jg jg' chr_g) (if rev then subGenome jh' jh chr_h else subGenome jh jh' chr_h) ->
              commonPrefix' rev (k + 1) jg' jh'
        (Just jg', Just jh') ->
          (setIdx g chr_ig jg, setIdx h chr_ih jh, rev, k, k == size chr_g && k == size chr_h && closes jg' jh')
        _ -> (setIdx g chr_ig jg, setIdx h chr_ih jh, rev, k, False)
      where
        testMatch = if rev then isReverseMatch else isDirectMatch
        -- both circular chromosomes were traversed and they match as whole
        -- chromosomes (same test used for whole circular blocks in Partition)
        closes jg' jh' =
          isCircular chr_g
            && isCircular chr_h
            && jg' == ig
            && jh' == ih
            && isMatch matcher chr_g chr_h
