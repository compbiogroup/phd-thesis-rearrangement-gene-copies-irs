{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE UndecidableInstances #-}

module Partition
  ( writePartition,
    getBlocksCorrespondence,
    trivialPartition,
    Partition,
    Breakpoints,
    mkPartitionFromBreakpoints,
    mkPartitionFromBlocks,
    underlineGenome,
    breakpoints,
    blocks,
    cindexedBlocks,
    checkCommonBal,
    checkCommonUnbal,
    findUncommon,
    getBlocksFit,
    costBalanced,
    costUnbalanced,
    numBreakpoints,
    CommonPartition (CGPbal, CGPunbal),
    mkCommonPartition,
    mkCommonPartition2,
    combine,
    bpsToBlockDelims,
    BlockDelim (..),
    blockDelims,
    delimToBlock,
  )
where

-- \|
-- Module      : Genomes
-- Description : Genome Partition representation and functions to manipulate it.
-- Copyright   : (c) Gabriel Siqueira, 2023
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com

import Control.Arrow ((&&&))
import Data.ByteString qualified as BS
import Data.ByteString.Char8 qualified as BS
import Data.Foldable (foldrM, toList)
import Data.Function (on)
import Data.List qualified as List
import Data.Map qualified as Map
import Data.Maybe (fromMaybe, isNothing, listToMaybe)
import Data.Sequence (Seq (Empty, (:<|), (:|>)), (><))
import Data.Sequence qualified as Seq
import Data.Set (Set)
import Data.Set qualified as Set
import Genomes (CIdx, ChromList, Idx, Chromosome (isCircular, isLinear, nextIdx), GenIdx (idxToInt, incIdx, mkIdx), GenesIRsF, GenesIRsR, Genome (..), IntergenicChromosome (..), Matcher (..), MultiChromosome (..), emptyGeneMap, geneMapIntersection, geneMapLookup, geneMapUnion, positionMap, writeFGenome, writeIR, writeRGenome)
import LocalBase

type Breakpoints g = Set (GIdx g)

-- TODO: Verificar mudanças para genomas circulares

data Partition c where
  -- Underlining Genome, Blocks, Breakpoints, Number of Blocks.
  -- The separation between the blocks are the breakpoints and the natural separation of the chromosomes.
  -- The Int is the number of blocks in the partition.
  GP :: (Chromosome c) => ChromList c -> [(CIdx,c)] -> Breakpoints (ChromList c) -> Int -> Partition c

instance (Eq (ChromList c), Eq c, Eq (GIdx c)) => Eq (Partition c) where
  (GP g1 bs1 bps1 n1) == (GP g2 bs2 bps2 n2) = g1 == g2 && bs1 == bs2 && bps1 == bps2 && n1 == n2

data BlockDelim c = BlockDelim
  { b_chr :: CIdx,
    b_beg :: GIdx c,
    b_end :: GIdx c,
    -- True if the block is a whole circular chromosome without breakpoints
    b_whole :: Bool
  }

instance (Show (GIdx c)) => Show (BlockDelim c) where
  show (BlockDelim c b e _) = show c ++ " (" ++ show b ++ "," ++ show e ++ ")"

instance (Eq (GIdx c)) => Eq (BlockDelim c) where
  b1 == b2 = (b_beg b1 == b_beg b2) && (b_end b1 == b_end b2)

-- Get the blocks delimiters of each chromosome. On a circular chromosome
-- with breakpoints the first block may cross the origin, and a breakpoint
-- at its last gene is a real cut between genes n and 1.
bpsToBlockDelims :: (MultiChromosome g, Chromosome (Chrom g), GenIdx (GIdx (Chrom g)), Eq (GIdx (Chrom g))) => g -> Breakpoints g -> [BlockDelim (Chrom g)]
bpsToBlockDelims g bps = concatMap chromDelims [1 .. mkIdx (Genomes.numChromosomes g)]
  where
    bpsByChr =
      Map.fromListWith (flip (++)) $
        (\idx -> (getIdxOfChrom g idx, [getIdxInChrom g idx])) <$> Set.toAscList bps
    chromDelims c =
      let chr = getChromosome c g
          ks = Map.findWithDefault [] c bpsByChr
          mk b e = BlockDelim c b e False
       in if
            | isLinear chr ->
                let cuts = filter (/= lastGIdx chr) ks
                 in zipWith mk (firstGIdx chr : map incIdx cuts) (cuts ++ [lastGIdx chr])
            | null ks -> [BlockDelim c (firstGIdx chr) (lastGIdx chr) True]
            -- the first block contains gene 1 (it crosses the origin unless the last breakpoint is at n)
            | otherwise -> zipWith mk (map (wrapInc chr) (last ks : init ks)) ks
    wrapInc chr = fromMaybe (error logicError) . nextIdx chr

-- Get the block of the genome given by its delimiter
delimToBlock :: (MultiChromosome g, Chromosome (Chrom g)) => g -> BlockDelim (Chrom g) -> Chrom g
delimToBlock g b = if b_whole b then chr else subGenome (b_beg b) (b_end b) chr
  where
    chr = getChromosome (b_chr b) g

mkPartitionFromBreakpoints :: (Chromosome c, GenIdx (GIdx c), Eq (GIdx c)) => ChromList c -> Breakpoints (ChromList c) -> Partition c
mkPartitionFromBreakpoints g bps = GP g bs bps (length bs)
  where
    bs = map (\b -> (b_chr b, delimToBlock g b)) (bpsToBlockDelims g bps)

-- | Make a partition from a list of blocks.
-- The blocks of @bs@ must concatenate into the genome @g@ and @g@ cannot have size 0.
-- The blocks of each chromosome must start at its first gene (a block
-- cannot cross the origin of a circular chromosome).
mkPartitionFromBlocks :: forall c. (Chromosome c, GenIdx (GIdx c), Ord (GIdx c)) => ChromList c -> [(CIdx, c)] -> Partition c
mkPartitionFromBlocks g bs = GP g bs bps (length bs)
  where
    chrs_0 = map (Genomes.size &&& isCircular) . getChromosomes $ g
    bps = Set.fromList $ go chrs_0 (map snd bs) (mkIdx 1) (mkIdx 0) []

    go :: [(Int, Bool)] -> [c] -> CIdx -> GIdx c -> [(CIdx, GIdx c)] -> [(CIdx, GIdx c)]
    go [] [] _ _ acc = acc
    go ((chr_size, circ) : chrs) (b : bs') chri i acc =
      if
        | b_size < chr_size -> go ((chr_size - b_size, circ) : chrs) bs' chri i' ((chri, i') : acc)
        -- the end of a split circular chromosome is a cut between genes n and 1
        | b_size == chr_size -> go chrs bs' (incIdx chri) (mkIdx 0) (if circ && isLinear b then (chri, i') : acc else acc)
        | otherwise -> error logicError
      where
        b_size = Genomes.size b
        i' = mkIdx $ idxToInt i + b_size
    go _ _ _ _ _ = error logicError

-- | Make a partition from block delimiters (in the order of the genome).
mkPartitionFromDelims :: (Chromosome c, Ord (GIdx c)) => ChromList c -> [BlockDelim c] -> Partition c
mkPartitionFromDelims g delims = GP g bs bps (length bs)
  where
    bs = map (\b -> (b_chr b, delimToBlock g b)) delims
    -- every block ends at a breakpoint, except at the end of a linear chromosome
    bps =
      Set.fromList
        [ (b_chr b, b_end b)
          | b <- delims,
            not (b_whole b),
            let chr = getChromosome (b_chr b) g,
            isCircular chr || b_end b /= lastGIdx chr
        ]

underlineGenome :: Partition c -> ChromList c
underlineGenome (GP g _ _ _) = g

trivialPartition :: (Chromosome c) => ChromList c -> Partition c
trivialPartition g = GP g (zip [1..] (getChromosomes g)) Set.empty (numChromosomes g)

-- | Delimiters of the blocks of a partition (in the order of the genome)
blockDelims :: (GenIdx (GIdx c), Eq (GIdx c)) => Partition c -> [BlockDelim c]
blockDelims (GP g _ bps _) = bpsToBlockDelims g bps

breakpoints :: Partition c -> Breakpoints (ChromList c)
breakpoints (GP _ _ bps _) = bps

blocks :: Partition c -> [c]
blocks (GP _ bs _ _) = map snd bs

cindexedBlocks :: Partition c -> [(CIdx,c)]
cindexedBlocks (GP _ bs _ _) = bs

size :: Partition c -> Int
size (GP _ _ _ n) = n

data CommonPartition g1 g2 where
  CGPbal :: (Chromosome c1, Chromosome c2) => Partition c1 -> Partition c2 -> CommonPartition c1 c2
  CGPunbal :: (Chromosome c1, Chromosome c2) => Partition c1 -> Partition c2 -> CommonPartition c1 c2

mkCommonPartition :: (Matcher m c1 c2, Chromosome c1, Chromosome c2, GIdx c1 ~ GIdx c2, Eq (GIdx c1)) => m c1 c2 -> Partition c1 -> Partition c2 -> CommonPartition c1 c2
mkCommonPartition matcher pg ph =
  if
    | checkCommonBal matcher pg ph -> CGPbal pg ph
    | checkCommonUnbal matcher pg ph -> CGPunbal pg ph
    | otherwise -> error logicError

-- | Number of breakpoints of a partition: each block ends at a breakpoint,
-- except the last block of a linear chromosome and a whole circular
-- chromosome (a single circular block, which has no breakpoint).
numBreakpoints :: Partition c -> Int
numBreakpoints (GP g bs _ _) = length bs - length (filter isLinear (getChromosomes g)) - length (filter (isCircular . snd) bs)

-- | Cost of a balanced common partition: the number of breakpoints of both partitions
costBalanced :: CommonPartition c1 c2 -> Int
costBalanced (CGPbal pg ph) = numBreakpoints pg + numBreakpoints ph
costBalanced (CGPunbal _ _) = error patternError

-- | Cost of a common partition (balanced or not): the number of breakpoints of
-- both partitions (exclusive blocks included)
costUnbalanced :: m c1 c2 -> CommonPartition c1 c2 -> Int
costUnbalanced _ (CGPbal pg ph) = numBreakpoints pg + numBreakpoints ph
costUnbalanced _ (CGPunbal pg ph) = numBreakpoints pg + numBreakpoints ph

-- | Emparelhamento máximo: bloco de G (índice 0-based) -> bloco de H (índice 0-based)
getBlocksFit :: (Matcher m c1 c2) => m c1 c2 -> Partition c1 -> Partition c2 -> Map.Map Int Int
getBlocksFit matcher pg ph = flipMap $ foldl selectBlock Map.empty [0 .. length corresp - 1]
  where
    corresp = getBlocksCorrespondence_ matcher pg ph
    -- assigned: bloco de H -> bloco de G
    selectBlock assigned b_g = fromMaybe assigned . snd $ augment Set.empty b_g assigned

    -- Busca um caminho aumentante a partir de b_g. O conjunto seen (blocos de H
    -- já visitados) é devolvido também quando a busca falha, de modo que cada
    -- bloco de H é visitado no máximo uma vez por busca (algoritmo de Kuhn, polinomial).
    augment seen b_g assigned = go seen (corresp !! b_g)
      where
        go seen' [] = (seen', Nothing)
        go seen' (b_h : rest)
          | b_h `Set.member` seen' = go seen' rest
          | otherwise =
              let seen'' = Set.insert b_h seen'
               in case Map.lookup b_h assigned of
                    Nothing -> (seen'', Just (Map.insert b_h b_g assigned))
                    Just b_g' -> case augment seen'' b_g' assigned of
                      (s, Just assigned') -> (s, Just (Map.insert b_h b_g assigned'))
                      (s, Nothing) -> go s rest

-- | Verify if the elements of pg can be assign to elements of ph with a perfect
-- match.
checkCommonBal :: (Matcher m c1 c2, Chromosome c1) => m c1 c2 -> Partition c1 -> Partition c2 -> Bool
checkCommonBal matcher pg ph = Partition.size pg == Partition.size ph && isNothing (findUncommon matcher pg ph)

-- | Verify if the elements of pg (that are not in exclusive blocks) can be assign to elements of ph
-- (that are not in exclusive blocks) with a perfect match.
checkCommonUnbal :: (Matcher m c1 c2, Chromosome c1, Chromosome c2, GIdx c1 ~ GIdx c2, Eq (GIdx c1)) => m c1 c2 -> Partition c1 -> Partition c2 -> Bool
checkCommonUnbal matcher pg ph = geneMapIntersection gmap_g gmap_h == emptyGeneMap
  where
    -- the caps (telomeres) are labels like any other: they can only be in exclusive
    -- blocks of the genome where they are abundant (the one with more linear chromosomes)
    bs_g = blocks pg
    bs_h = blocks ph
    blocks_fit_g = getBlocksFit matcher pg ph
    blocks_fit_h = flipMap blocks_fit_g
    to_gmap blocks_fit bs = (foldr (geneMapUnion . positionMap . snd) emptyGeneMap . filter (\(i, _) -> i `Map.notMember` blocks_fit)) $ zip [0 ..] bs
    gmap_g = to_gmap blocks_fit_g bs_g
    gmap_h = to_gmap blocks_fit_h bs_h

-- | Find largest blocks that do not have sufficient correspondences between one partition and the other
-- match.
findUncommon :: (Matcher m c1 c2, Chromosome c1) => m c1 c2 -> Partition c1 -> Partition c2 -> Maybe [c1]
findUncommon matcher pg ph =
  case testPerfectMatch of
    Left (i, fixCorres) -> Just ((bsG !! i) : map ((bsG !!) . (fixCorres Map.!)) (cors !! i))
    Right _ -> Nothing
  where
    bsG = blocks pg
    cors = getBlocksCorrespondence_ matcher pg ph

    testPerfectMatch = foldrM (selectBlock Set.empty) Map.empty block_indices
      where
        -- Block indices sorted by block size in decreasing order
        block_indices = map snd . sortWith (negate . fst) $ zip (map Genomes.size bsG) [0 .. length cors - 1]

        -- Select block from ph that corresponds to block in position i of pg.
        -- fixCorres saves which blocks of pg are fixed in which blocks of ph (indices starting in 0).
        selectBlock seen i fixCorres =
          if done then Right fixCorres' else Left (i, fixCorres')
          where
            (_, fixCorres', done) = selectBlock' seen i fixCorres

        selectBlock' seen i = testCandidates i cands seen
          where
            cands = cors !! i

        testCandidates _ [] seen fixCorres = (seen, fixCorres, False)
        testCandidates i (u : us) seen fixCorres =
          if
            | u `Set.member` seen -> testCandidates i us seen fixCorres
            | u `Map.notMember` fixCorres -> (seen', fixCorres', True)
            | success -> (seen_rec, fixCorres_rec', True)
            | otherwise -> testCandidates i us seen_rec fixCorres_rec
          where
            seen' = Set.insert u seen
            fixCorres' = Map.insert u i fixCorres
            (seen_rec, fixCorres_rec, success) = selectBlock' seen' (fixCorres Map.! u) fixCorres
            fixCorres_rec' = Map.insert u i fixCorres_rec

mkCommonPartition2 :: (Matcher m c1 c2, Chromosome c1, Chromosome c2, GIdx c1 ~ GIdx c2, GenIdx (GIdx c1), Eq (GIdx c1)) => m c1 c2 -> ChromList c1 -> Breakpoints (ChromList c1) -> ChromList c2 -> Breakpoints (ChromList c2) -> CommonPartition c1 c2
mkCommonPartition2 matcher g bg h bh = mkCommonPartition matcher (mkPartitionFromBreakpoints g bg) (mkPartitionFromBreakpoints h bh)

mkCommonPartition3 :: (Matcher m c1 c2, Chromosome c1, Chromosome c2, GIdx c1 ~ GIdx c2, Ord (GIdx c1)) => m c1 c2 -> ChromList c1 -> [BlockDelim c1] -> ChromList c2 -> [BlockDelim c2] -> CommonPartition c1 c2
mkCommonPartition3 matcher g bsl_g h bsl_h = mkCommonPartition matcher (mkPartitionFromDelims g bsl_g) (mkPartitionFromDelims h bsl_h)

instance (Show c) => Show (Partition c) where
  show pg = combiStr subs_g
    where
      combiStr strs = unwords $ interleavelists strs (replicate (length strs - 1) "|")
      subs_g = map show $ blocks pg

instance (Show g1, Show g2) => Show (CommonPartition g1 g2) where
  show (CGPbal pg ph) = unlines [show pg, show ph]
  show (CGPunbal pg ph) = unlines [show pg, show ph]

-- | Write each side of the partition, one entry per chromosome (as in writeMultiC).
-- On circular chromosomes the last breakpoint IR is the one after the last block.
writePartition :: CommonPartition GenesIRsR GenesIRsF -> (BS.ByteString, BS.ByteString, BS.ByteString, BS.ByteString)
writePartition cp = (genes_bs_g, irs_bs_g, genes_bs_h, irs_bs_h)
  where
    (pg, ph) = case cp of
      CGPbal pg_ ph_ -> (pg_, ph_)
      CGPunbal pg_ ph_ -> (pg_, ph_)
    (genes_bs_g, irs_bs_g) = writeSide (writeRGenome False) pg
    (genes_bs_h, irs_bs_h) = writeSide (writeFGenome False) ph
    combiBS bss = BS.unwords $ interleavelists bss (replicate (length bss - 1) "|")

    writeSide :: (IntergenicChromosome c, GIdx c ~ Idx) => (c -> (BS.ByteString, BS.ByteString)) -> Partition c -> (BS.ByteString, BS.ByteString)
    writeSide writer (GP g bs bps _) = (BS.unwords (map fst chrs), BS.unwords (map snd chrs))
      where
        chrs = map writeChr [1 .. mkIdx (numChromosomes g)]
        writeChr c =
          let chr = getChromosome c g
              bss = map (writer . snd) . filter ((== c) . fst) $ bs
              bp_irs =
                [ writeIR (getIR i chr)
                  | (c', i) <- Set.toAscList bps,
                    c' == c,
                    isCircular chr || i /= lastGIdx chr
                ]
           in ( (if isLinear chr then "L " else "C ") <> combiBS (map fst bss) <> " ;",
                combiBS (interleavelists (map snd bss) bp_irs) <> " ;"
              )

-- | For each block in pg, the correspondent positions of blocks in pg, indices starting in 0
getBlocksCorrespondence :: (Matcher m c1 c2) => m c1 c2 -> CommonPartition c1 c2 -> [[Int]]
getBlocksCorrespondence matcher (CGPbal pg ph) = getBlocksCorrespondence_ matcher pg ph
getBlocksCorrespondence matcher (CGPunbal pg ph) = getBlocksCorrespondence_ matcher pg ph

getBlocksCorrespondence_ :: (Matcher m c1 c2) => m c1 c2 -> Partition c1 -> Partition c2 -> [[Int]]
getBlocksCorrespondence_ matcher pg ph =
  do
    sub_g <- blocks pg
    let sub_hs = blocks ph
    return . map fst . filter (isMatch matcher sub_g . snd) . zip [0 ..] $ sub_hs

data CombineWhat = CombineSinSin | CombineSinRep | CombineRepRep

combine :: forall m c1 c2. (Matcher m c1 c2, Chromosome c1, Chromosome c2, GIdx c1 ~ GIdx c2, Ord (GIdx c1), GenIdx (GIdx c1)) => m c1 c2 -> CommonPartition c1 c2 -> CommonPartition c1 c2
combine matcher cp = mkCommonPartition3 matcher g blo_g''' h blo_h'''
  where
    (pg, ph) = case cp of
      CGPbal pg_ ph_ -> (pg_, ph_)
      CGPunbal pg_ ph_ -> (pg_, ph_)
    g = underlineGenome pg
    h = underlineGenome ph
    bg = breakpoints pg
    bh = breakpoints ph
    (blo_g', blo_h') = combine_ CombineSinSin (Seq.fromList blo_g, Seq.fromList blo_h)
    (blo_g'', blo_h'') = mergeExclusive (toList blo_g') (toList blo_h')
    (blo_g''', blo_h''') = closeCirculars blo_g'' blo_h''
    blo_g = testSingletons g <$> bpsToBlockDelims g bg
    blo_h = testSingletons h <$> bpsToBlockDelims h bh

    testSingletons :: (Chromosome c) => ChromList c -> BlockDelim c -> (BlockDelim c, Bool)
    testSingletons gx bs = (bs, any testSingleton (betweenGIdx (getChromosome (b_chr bs) gx) (b_beg bs) (b_end bs)))
      where
        posMap = positionMap gx
        testSingleton i =
          case geneMapLookup (getGene (b_chr bs, i) gx) posMap of
            Nothing -> False
            Just pos -> length pos == 1

    combine_ :: CombineWhat -> (Seq (BlockDelim c1, Bool), Seq (BlockDelim c2, Bool)) -> (Seq (BlockDelim c1), Seq (BlockDelim c2))
    combine_ what (bs1, bs2) =
      let (bs1', bs2') = go1 (Empty, bs1, Empty, bs2)
       in if
            -- If there was a combination the length of bs1 changed
            | length bs1' /= length bs1 -> combine_ what (bs1', bs2')
            -- no consecutive blocks combined, try blocks across the origin
            | Just (bs1'', bs2'') <- wrapMerge -> combine_ what (bs1'', bs2'')
            | otherwise -> case what of
                CombineSinSin -> combine_ CombineSinRep (bs1, bs2)
                CombineSinRep -> combine_ CombineRepRep (bs1, bs2)
                CombineRepRep -> (fmap fst bs1, fmap fst bs2)
      where
        repsOk (ba, bb) = case what of
          CombineSinSin -> snd ba && snd bb
          CombineSinRep -> snd ba || snd bb
          CombineRepRep -> True

        wrapMerge =
          listToMaybe
            [ (mergeAt bs1 p1, mergeAt bs2 p2)
              | (p1, p2) <-
                  [(p1, p2) | p1 <- wraps1, p2 <- cons2 ++ wraps2]
                    ++ [(p1, p2) | p1 <- cons1, p2 <- wraps2],
                repsOk (Seq.index bs1 (fst p1), Seq.index bs1 (snd p1)),
                isMatch matcher (delimToBlock g . fst $ merged bs1 p1) (delimToBlock h . fst $ merged bs2 p2),
                isCommon (mergeAt bs1 p1) (mergeAt bs2 p2)
            ]
        (cons1, wraps1) = adjPairs g bs1
        (cons2, wraps2) = adjPairs h bs2

        -- We have to iterate through the elements of the two sequences,
        -- during the iteration they are moved from the sequence on the right
        -- to the one on the left.
        -- If the sequences are empty there is nothing to do
        go1 (Empty, Empty, Empty, Empty) = (Empty, Empty)
        -- We need at least one more element in the first and second sequence to combine
        go1 (pre, ba1 :<| Empty, Empty, blocks2) = (pre :|> ba1, blocks2)
        go1 (_, _ :<| Empty, _, _) = error patternError -- This should not happen, because the go1 is initially called with pre2 = Empty after each move in the first seq
        go1 (pre1, ba1 :<| bb1 :<| suf1, pre2, ba2 :<| Empty) =
          go1 (pre1 :|> ba1, bb1 :<| suf1, Empty, pre2 :|> ba2)
        -- Now we try to combine
        go1 (pre1, blocks1@(ba1 :<| bb1 :<| suf1), pre2, blocks2@(ba2 :<| bb2 :<| suf2))
          -- test if the two blocks fall into the correct case. If not go to next block in sequence 1.
          | not testReps = go1 (pre1 :|> ba1, bb1 :<| suf1, Empty, pre2 >< blocks2)
          -- test if the block can be combined
          | sameChr && isMatch matcher g_b1' h_b2' && isCommon (pre1 >< (b1' :<| suf1)) ((pre2 :|> b2') >< suf2) =
              go1 (pre1, b1' :<| suf1, Empty, (pre2 :|> b2') >< suf2)
          -- next block of sequence 2
          | otherwise = go1 (pre1, blocks1, pre2 :|> ba2, bb2 :<| suf2)
          where
            sameChr = b_chr (fst ba1) == b_chr (fst bb1) && b_chr (fst ba2) == b_chr (fst bb2) 
            testReps = repsOk (ba1, bb1)

            b1' = (BlockDelim (b_chr (fst ba1)) (b_beg (fst ba1)) (b_end (fst bb1)) False, snd ba1 || snd bb1)
            b2' = (BlockDelim (b_chr (fst ba2)) (b_beg (fst ba2)) (b_end (fst bb2)) False, snd ba2 || snd bb2)
            g_b1' = delimToBlock g (fst b1')
            h_b2' = delimToBlock h (fst b2')
        go1 (_, _, _, _) = error patternError

    -- The merged blocks may match each other while the blocks around them lose
    -- their only partner, so a merge is only accepted if the partition is still
    -- common (and of the same kind, balanced or not, as the original one)
    isCommon :: Seq (BlockDelim c1, Bool) -> Seq (BlockDelim c2, Bool) -> Bool
    isCommon bs1 bs2 = checkCommon matcher (mkPartitionFromDelims g (fst <$> toList bs1)) (mkPartitionFromDelims h (fst <$> toList bs2))
    checkCommon = case cp of
      CGPbal _ _ -> checkCommonBal
      CGPunbal _ _ -> checkCommonUnbal

    -- Exclusive blocks (blocks without a pair in a maximum matching) that are
    -- neighbours in the same chromosome are merged: the merged block is still
    -- exclusive and has the same labels, so the partition is still common.
    mergeExclusive :: [BlockDelim c1] -> [BlockDelim c2] -> ([BlockDelim c1], [BlockDelim c2])
    mergeExclusive ds1 ds2 =
      if isCommon (Seq.fromList (map (,False) ds1')) (Seq.fromList (map (,False) ds2'))
        then (ds1', ds2')
        else (ds1, ds2)
      where
        fit = getBlocksFit matcher (mkPartitionFromDelims g ds1) (mkPartitionFromDelims h ds2)
        exclusive1 = Set.fromList [k | k <- [0 .. length ds1 - 1], k `Map.notMember` fit]
        exclusive2 = Set.fromList [0 .. length ds2 - 1] `Set.difference` Set.fromList (Map.elems fit)
        ds1' = mergeRuns g exclusive1 ds1
        ds2' = mergeRuns h exclusive2 ds2

    -- Merge runs of consecutive blocks (in the order of the genome) whose
    -- indices are in @excl@, also across the origin of circular chromosomes.
    -- A circular chromosome left with a single block becomes a whole block.
    mergeRuns :: (Chromosome c) => ChromList c -> Set Int -> [BlockDelim c] -> [BlockDelim c]
    mergeRuns gx excl ds = concatMap mergeChr . List.groupBy ((==) `on` (b_chr . fst)) . linear $ zip ds [k `Set.member` excl | k <- [0 ..]]
      where
        joinD d1 d2 = BlockDelim (b_chr d1) (b_beg d1) (b_end d2) False
        linear ((d1, True) : (d2, True) : rest) | b_chr d1 == b_chr d2 = linear ((joinD d1 d2, True) : rest)
        linear (x : rest) = x : linear rest
        linear [] = []
        mergeChr bs@((d_first, True) : _)
          | (d_last, True) <- last bs,
            isCircular chr =
              case init bs of
                -- the whole circular chromosome is exclusive
                [] -> [whole]
                [_] -> [whole]
                _ : middle -> joinD d_last d_first : map fst middle
          where
            chr = getChromosome (b_chr d_first) gx
            whole = BlockDelim (b_chr d_first) (firstGIdx chr) (lastGIdx chr) True
        mergeChr bs = map fst bs

    -- A circular chromosome left with a single block that is not whole (cut only
    -- once, e.g. a circular chromosome of one gene, which SOAR always cuts) can not
    -- be fixed by wrapMerge, which needs two blocks. Such a chromosome of g is made
    -- whole together with one of h in the same situation, if the two chromosomes
    -- match as whole chromosomes and the partition is still common.
    closeCirculars :: [BlockDelim c1] -> [BlockDelim c2] -> ([BlockDelim c1], [BlockDelim c2])
    closeCirculars ds1 ds2 = foldl tryClose (ds1, ds2) (cutOnce g ds1)
      where
        tryClose (ds1', ds2') k1 =
          case [ (ds1'', ds2'')
                 | k2 <- cutOnce h ds2',
                   isMatch matcher (chrOf g ds1' k1) (chrOf h ds2' k2),
                   let ds1'' = closeAt g ds1' k1
                       ds2'' = closeAt h ds2' k2,
                   isCommon (Seq.fromList (map (,False) ds1'')) (Seq.fromList (map (,False) ds2''))
               ] of
            r : _ -> r
            [] -> (ds1', ds2')
        chrOf gx ds k = getChromosome (b_chr (ds !! k)) gx
        closeAt gx ds k =
          let d = ds !! k
              chr = getChromosome (b_chr d) gx
           in take k ds ++ [BlockDelim (b_chr d) (firstGIdx chr) (lastGIdx chr) True] ++ drop (k + 1) ds
    -- indices of the blocks that are the only block of a circular chromosome but not whole
    cutOnce :: (Chromosome c) => ChromList c -> [BlockDelim c] -> [Int]
    cutOnce gx ds =
      [ k
        | (k, d) <- zip [0 ..] ds,
          not (b_whole d),
          isCircular (getChromosome (b_chr d) gx),
          length (filter ((== b_chr d) . b_chr) ds) == 1
      ]

    -- Pairs of indices of neighbour blocks: consecutive in the same chromosome,
    -- and (last, first) of each circular chromosome with at least two blocks
    adjPairs :: (Chromosome c) => ChromList c -> Seq (BlockDelim c, Bool) -> ([(Int, Int)], [(Int, Int)])
    adjPairs gx bs = (consecutive, wraps)
      where
        chrOf = b_chr . fst . Seq.index bs
        consecutive = [(k, k + 1) | k <- [0 .. Seq.length bs - 2], chrOf k == chrOf (k + 1)]
        chrRanges = Map.fromListWith (\(f1, l1) (f2, l2) -> (min f1 f2, max l1 l2)) [(chrOf k, (k, k)) | k <- [0 .. Seq.length bs - 1]]
        wraps = [(l, f) | (chr, (f, l)) <- Map.toList chrRanges, l > f, isCircular (getChromosome chr gx)]

    -- Combined block of a pair, with the flag of either
    merged :: Seq (BlockDelim c, Bool) -> (Int, Int) -> (BlockDelim c, Bool)
    merged bs (a, b) = (BlockDelim (b_chr da) (b_beg da) (b_end db) False, fa || fb)
      where
        (da, fa) = Seq.index bs a
        (db, fb) = Seq.index bs b

    -- Put the combined block where the earlier block of the pair was and remove the other
    mergeAt :: Seq (BlockDelim c, Bool) -> (Int, Int) -> Seq (BlockDelim c, Bool)
    mergeAt bs (a, b) = Seq.deleteAt (max a b) (Seq.update (min a b) (merged bs (a, b)) bs)
