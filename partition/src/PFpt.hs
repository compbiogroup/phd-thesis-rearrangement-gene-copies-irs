{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE TypeOperators #-}

-- \|
-- Module      : Genomes
-- Description : Algorithm for genome partition with an polinomial complexity
-- parameterized by the partition size, maximum occurrence and/or number of distinct
-- elements between the alphabets. Genes that match both directly and in reverse
-- (every gene of an unsigned reverse instance, and the caps 0) get witnesses in both
-- directions: singletons of this kind are not fixed in the initial sample but added
-- by the branching rule, which tries both directions (findAllLimits). Unbalanced
-- genomes are handled by branching rule 2, which chooses the genes of the exclusive
-- blocks (the removed genes). The search starts with the best partition of the
-- heuristics Greedy and SOARComb as the best one found, so the pruning works from
-- the start and the result is never worse than theirs.
-- Copyright   : (c) Gabriel Siqueira, 2023
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
module PFpt
  ( fptPartition,
  )
where

import Control.Concurrent (newMVar, readMVar, swapMVar)
import Control.Concurrent.Async (race)
import Control.Concurrent.MVar (MVar)
import Control.Exception (assert, evaluate)
import Control.Monad (void, when)
import Data.Either (isRight)
import Data.Foldable (toList)
import Data.Map qualified as Map
import Data.List (group, minimumBy)
import Data.Maybe (catMaybes, fromMaybe, isNothing, listToMaybe)
import Data.MultiMap (MultiMap)
import Data.MultiMap qualified as MMap
import Data.Ord (comparing)
import Data.Set qualified as Set
import GHC.Conc (threadDelay)
import Genomes (CIdx, Gene, ChromList, Chromosome (isCircular, nextIdx, prevIdx), GIdx, GenIdx (idxToInt), GeneMap, Genome (alphabet, getGene, size, subGenome), Matcher (isDirectMatch, isMatch, isReverseMatch), MultiChromosome (getChromosome, getChromosomes, numChromosomes), adjacencies, allGIdx, geneMapLookup, occurrence, positionMap, singletonOnBoth)
import LocalBase
import PGreedy (greedyPartition)
import PSOAR (soarPartition)
import Partition (CommonPartition, costUnbalanced, mkCommonPartition2)

type FptCtx c1 c2 = (Chromosome c1, Chromosome c2, GIdx c1 ~ GIdx c2, Ord (GIdx c1), GenIdx (GIdx c1))

data Vertex i = Vtx
  { vInG :: Bool,
    vIdx :: i
  }
  deriving (Eq, Ord, Show)

-- In all Idx Pairs in Edge and Samples we assume that the first index is
-- from G and the second from H
data Edge i = Edge
  { eOther :: Vertex i,
    eVpair :: (i, i),
    eParent :: BlockWitness i,
    eBlack :: Bool
  }
  deriving (Eq, Ord, Show)

type Sample i = [BlockWitness i]

data BlockWitness i = BW
  { w_left_pair :: (i, i), -- leftmost pair (looking in G) parallel to black_pair
    w_black_pair :: (i, i),
    w_right_pair :: (i, i), -- rightmost pair (looking in G) parallel to black_pair
    w_inverted :: Bool
  }
  deriving (Eq, Ord, Show)

data SampleGraph m c1 c2 i = SG
  { sg_matcher :: m c1 c2,
    sg_g :: ChromList c1,
    sg_h :: ChromList c2,
    sg_posG :: GeneMap [i],
    sg_posH :: GeneMap [i],
    sample :: Sample i,
    edges :: MultiMap (Vertex i) (Edge i),
    -- genes chosen to be in exclusive blocks (the set D), they are not vertices of
    -- the graph (they have no edges)
    removed :: Set.Set (Vertex i)
  }

instance (Show g1, Show g2, Show i) => Show (SampleGraph m g1 g2 i) where
  show SG {..} = unlines [show sg_g, show sg_h, show sample, show (Set.toList removed), edges_str]
    where
      edges_str = unlines . map show . MMap.toList $ edges

fptPartition :: (FptCtx c1 c2, Matcher m c1 c2) => Int -> m c1 c2 -> ChromList c1 -> ChromList c2 -> IO (Maybe (CommonPartition c1 c2), Bool)
fptPartition time_limit matcher g h = do
  -- build the initial sample graph and the partitions of the heuristics before the
  -- time limit starts
  sg0' <- evaluate sg0
  cost0 <- evaluate (fst initial)
  best_sg <- newMVar (cost0, Nothing)
  full_comp <- isRight <$> race (threadDelay time_limit) (runFpt best_sg sg0')
  maybe_sg <- snd <$> readMVar best_sg
  return
    . (,full_comp) -- Add the bool indicating whether the computation finish
    . Just
    . maybe
      (snd initial) -- no sample graph is better than the heuristics
      ( sampleGraphToPartition
          . (\sg' -> assert (checkDegrees 1 1 sg') sg')
          . finalizeCorrespondence
          . (\sg' -> assert (checkDegrees 1 2 sg') sg')
      )
    $ maybe_sg
  where
    -- best partition of the heuristics (the first one on ties), with its cost
    initial =
      minimumBy
        (comparing fst)
        [(costUnbalanced matcher p, p) | p <- [greedyPartition False matcher g h, soarPartition True matcher g h]]
    sg0 = (\sg -> assert (checkDegrees 0 2 sg) sg) . removeAbsent $ createSampleGraph matcher g h sample0
    -- a gene whose label does not occur in the other genome is in an exclusive block
    -- in every partition (branching rule 2 has a single branch for it), so it is
    -- removed from the start, which also lets toBig count its block
    removeAbsent sg = foldr removeVertex sg (absent True g posMapH ++ absent False h posMapG)
    absent inG q posMapOther = [Vtx inG idx | idx <- allGIdx q, occurrence posMapOther (getGene idx q) == 0]
    -- singletons whose two occurrences match in exactly one direction. A singleton
    -- that matches both directly and in reverse (every gene of an unsigned reverse
    -- instance) is left out: its isolated vertex is a rare odd path, so branching
    -- rule 1 adds it with the witnesses of both directions (findAllLimits), instead
    -- of fixing the direct one here (findLimits)
    sample0 = filter uniqueWitness . filter (genesMatchOneWay matcher g h) . map getBlockWitnessIdxs . filter (singletonOnBoth g posMapG posMapH) . toList $ alphabet g

    -- a singleton whose witness is not unique (several windows in a circular
    -- chromosome) is also left to branching rule 1, which tries all of them
    uniqueWitness idxs = length (findAllLimits (createSampleGraph matcher g h []) idxs) == 1
    getBlockWitnessIdxs gene =
      case (geneMapLookup gene posMapG, geneMapLookup gene posMapH) of
        (Just [posG], Just [posH]) -> (posG, posH)
        _ -> error patternError
    -- not taken from sg0, which is built from sample0
    posMapG = positionMap g
    posMapH = positionMap h

-- The MVar keeps the cost of the best partition found so far and its sample graph
-- (Nothing while the best one is the initial partition of the heuristics).
runFpt :: (FptCtx c1 c2, Matcher m c1 c2) => MVar (Int, Maybe (SampleGraph m c1 c2 (GIdx (ChromList c1)))) -> SampleGraph m c1 c2 (GIdx (ChromList c1)) -> IO ()
runFpt best_sg sg0 = branch sg0
  where
    -- every label occurs as many times in G as in H (no gene is ever removed)
    balanced = fmap length (sg_posG sg0) == fmap length (sg_posH sg0)
    branch !sg = do
      test <- toBig sg
      if
        | test -> return ()
        | Just path <- getRareOddPath sg -> mapM_ branch' (rule1 sg path)
        | Just path <- getAbundantOddPath sg -> mapM_ branch' (rule2 sg path)
        | otherwise -> selectSG sg
    branch' sg = assert (checkDegrees 0 2 sg) $ branch sg

    -- Branching Rule 1: If the sample graph sg contains a connected component which is a rare odd path
    -- (u_1, v_1, ... ,v_{l−1}, u_l), then do the following for each vertex u_i, 1 <= i <= l: for each
    -- vertex v such that it has the same label as u_i, it is not in a match of the sample, and it is
    -- from a different string than u_i, branch into the case to add (u_i, v) to the sample, if (u_i, v)
    -- is not parallel to any match of the sample.
    rule1 sg path = [addBlackEdge bw sg | u <- largerSide path, bw <- candidateWitnesses sg u]

    -- Branching Rule 2: If the sample graph sg contains a connected component which is an abundant
    -- odd path (u_1, v_1, ... ,v_{l−1}, u_l), then do the following for each vertex u_i, 1 <= i <= l:
    -- branch into the case to remove u_i (put it in an exclusive block), if fewer than
    -- |occ(a, G) - occ(a, H)| genes with the label a of u_i were removed, and, as in rule 1, into
    -- the cases to add each candidate (u_i, v) to the sample.
    rule2 sg path =
      [ sg'
        | u <- largerSide path,
          sg' <- [removeVertex u sg | canRemove sg u] ++ [addBlackEdge bw sg | bw <- candidateWitnesses sg u]
      ]

    -- The sample and the removed genes only grow in a branch. Each witness is at least
    -- one common block on each side, and each run of removed genes that can not join
    -- another one is at least one exclusive block, so the cost (breakpoints of both
    -- partitions) of any partition from this branch is at least the number of these
    -- blocks minus one per chromosome (the last block of a linear chromosome or a
    -- whole circular chromosome has no breakpoint).
    toBig sg = do
      (best_cost, _) <- readMVar best_sg
      return $ 2 * length (sample sg) + (if balanced then 0 else exclusiveRuns sg) - numChrs sg >= best_cost
    numChrs sg = numChromosomes (sg_g sg) + numChromosomes (sg_h sg)
    selectSG sg = do
      (best_cost, _) <- readMVar best_sg
      let c = cost sg
      -- strict, so the previous best is not kept inside a thunk
      when (c < best_cost) $ void (swapMVar best_sg $! (c, Just sg))

-- | Vertices of an odd path in the genome with more vertices in it (the genome of
-- its ends).
largerSide :: [Vertex i] -> [Vertex i]
largerSide vertices = if length verticesG > length verticesH then verticesG else verticesH
  where
    verticesG = filter vInG vertices
    verticesH = filter (not . vInG) vertices

-- | Gene of a vertex.
geneOf :: (FptCtx c1 c2) => SampleGraph m c1 c2 (GIdx (ChromList c1)) -> Vertex (GIdx (ChromList c1)) -> Gene
geneOf SG {..} v = if vInG v then getGene (vIdx v) sg_g else getGene (vIdx v) sg_h

-- | Number of genes with the label of the gene that must be in exclusive blocks
-- (|occ(a, G) - occ(a, H)|).
excess :: SampleGraph m c1 c2 i -> Gene -> Int
excess SG {..} gene = abs (occurrence sg_posG gene - occurrence sg_posH gene)

-- | Number of removed genes with each label.
removedCount :: (FptCtx c1 c2) => SampleGraph m c1 c2 (GIdx (ChromList c1)) -> Map.Map Gene Int
removedCount sg = Map.fromListWith (+) [(canonicOri (geneOf sg v), 1) | v <- Set.toList (removed sg)]

-- | Whether u can be removed: fewer than |occ(a, G) - occ(a, H)| genes with its
-- label a were removed.
canRemove :: (FptCtx c1 c2) => SampleGraph m c1 c2 (GIdx (ChromList c1)) -> Vertex (GIdx (ChromList c1)) -> Bool
canRemove sg u = Map.findWithDefault 0 (canonicOri gene) (removedCount sg) < excess sg gene
  where
    gene = geneOf sg u

-- | Lower bound on the number of exclusive blocks of the partitions of a branch:
-- the maximal runs of removed genes of each chromosome, where two runs are counted
-- separately only if a gene between them can never be removed (a gene of a black
-- edge, or a gene whose label is not abundant in its genome or already has all
-- its genes in excess removed). Otherwise the genes between them may be removed
-- later and the runs become a single block. At the end of a branch each label a
-- has exactly |occ(a, G) - occ(a, H)| removed genes, so when the genes of a label
-- that can still be removed are exactly the ones it still needs, they are counted
-- as removed.
exclusiveRuns :: (FptCtx c1 c2) => SampleGraph m c1 c2 (GIdx (ChromList c1)) -> Int
exclusiveRuns sg@SG {..} = runsOf True sg_g + runsOf False sg_h
  where
    counts = removedCount sg
    runsOf inG q = sum [chrRuns chr (map (status . Vtx inG . (chri,)) (allGIdx chr)) | (chri, chr) <- zip [1 ..] (getChromosomes q)]
    allVertices = map (Vtx True) (allGIdx sg_g) ++ map (Vtx False) (allGIdx sg_h)
    label = canonicOri . geneOf sg
    missing v = excess sg (geneOf sg v) - Map.findWithDefault 0 (label v) counts
    -- number of genes of each label that can still be removed
    candidates = Map.fromListWith (+) [(label v, 1 :: Int) | v <- allVertices, isNothing (baseStatus v)]
    -- Just True for removed genes, Just False for genes that can never be removed
    status v = case baseStatus v of
      Nothing | Map.findWithDefault 0 (label v) candidates == missing v -> Just True
      st -> st
    baseStatus v
      | v `Set.member` removed = Just True
      | any eBlack (edges MMap.! v) = Just False
      | not (abundant v) || missing v <= 0 = Just False
      | otherwise = Nothing
    abundant v = (if vInG v then id else flip) (>) (occurrence sg_posG gene) (occurrence sg_posH gene)
      where
        gene = geneOf sg v
    chrRuns chr sts =
      let grps = map localHead . group . catMaybes $ sts
          rs = length (filter id grps)
          -- the first and the last runs of a circular chromosome are joined
          joined = isCircular chr && length grps > 1 && localHead grps && last grps
       in if joined then rs - 1 else rs

-- TODO: Change the implementation so that the partition do not need to be fully calculated
-- The cost of the partition correspondent to the sample graph
cost :: (FptCtx c1 c2, Matcher m c1 c2) => SampleGraph m c1 c2 (GIdx (ChromList c1)) -> Int
cost sg = costUnbalanced (sg_matcher sg) . sampleGraphToPartition . finalizeCorrespondence $ sg

sampleGraphToPartition :: (FptCtx c1 c2, Matcher m c1 c2) => SampleGraph m c1 c2 (GIdx (ChromList c1)) -> CommonPartition c1 c2
sampleGraphToPartition (SG {..}) = mkCommonPartition2 sg_matcher sg_g bg sg_h bh
  where
    bg = breaksOf True sg_g
    bh = breaksOf False sg_h
    -- an adjacency (i, i') is a breakpoint if its genes come from different witnesses,
    -- or if it closes a circular chromosome covered by a single witness (i' is the first
    -- gene of the witness), unless that witness also covers the closing adjacency. The
    -- removed genes form the exclusive blocks: an adjacency with one removed gene is a
    -- breakpoint, and one with two removed genes is not.
    breaksOf inG q =
      Set.fromList
        [ (chri, i)
          | (chri, chr) <- zip [1 ..] (getChromosomes q),
            (i, i') <- adjacencies chr,
            isBreak inG (chri, i) (chri, i')
        ]
    isBreak inG idx idx'
      | isRemoved idx || isRemoved idx' = isRemoved idx /= isRemoved idx'
      | otherwise = p /= parentOf inG idx' || idx' == startOf inG p && not (closesCircle p)
      where
        p = parentOf inG idx
        isRemoved = (`Set.member` removed) . Vtx inG
    parentOf inG idx = eParent . localHead . (\l -> if null l then error patternError else l) $ edges MMap.! Vtx inG idx

    -- first gene of the witness (in the order of the chromosome)
    startOf inG bw
      | inG = fst (w_left_pair bw)
      | otherwise = snd (if w_inverted bw then w_right_pair bw else w_left_pair bw)

    -- the witness covers two whole circular chromosomes that match (closing
    -- adjacency included); the whole chromosomes are compared, as in Partition,
    -- since for 1 gene subGenome gr gl does not include the closing IR
    closesCircle bw =
      isCircular chrG
        && isCircular chrH
        && size chrG == size chrH
        && nextIdx chrG gr == Just gl
        && (if w_inverted bw then nextIdx chrH hl == Just hr else nextIdx chrH hr == Just hl)
        && isMatch sg_matcher chrG chrH
      where
        -- hl and hr are the genes of H paired with the left (gl) and right (gr) ends in G
        ((cg, gl), (ch, hl)) = w_left_pair bw
        ((_, gr), (_, hr)) = w_right_pair bw
        chrG = getChromosome cg sg_g
        chrH = getChromosome ch sg_h

-- return vertices from a rare odd path if such path exists
getRareOddPath :: (FptCtx c1 c2) => SampleGraph m c1 c2 (GIdx (ChromList c1)) -> Maybe [Vertex (GIdx (ChromList c1))]
getRareOddPath = listToMaybe . getRareOddPaths

-- return vertices from a abundant odd path if such path exists
getAbundantOddPath :: (FptCtx c1 c2) => SampleGraph m c1 c2 (GIdx (ChromList c1)) -> Maybe [Vertex (GIdx (ChromList c1))]
getAbundantOddPath = listToMaybe . getAbundantOddPaths

-- return vertices from each rare odd path
getRareOddPaths :: (FptCtx c1 c2) => SampleGraph m c1 c2 (GIdx (ChromList c1)) -> [[Vertex (GIdx (ChromList c1))]]
getRareOddPaths sg = filter (isRare sg . localHead) (getOddPaths sg)

-- return vertices from each abundant odd path
getAbundantOddPaths :: (FptCtx c1 c2) => SampleGraph m c1 c2 (GIdx (ChromList c1)) -> [[Vertex (GIdx (ChromList c1))]]
getAbundantOddPaths sg = filter (not . isRare sg . localHead) (getOddPaths sg)

-- | The odd paths of the graph, each one starting at one of its ends (whose genome
-- has more vertices in the path, so the path is rare if the label of its first
-- vertex is rare). The paths are found from the vertices of G, plus the isolated
-- vertices of H (every other path has a vertex in G).
getOddPaths :: (FptCtx c1 c2) => SampleGraph m c1 c2 (GIdx (ChromList c1)) -> [[Vertex (GIdx (ChromList c1))]]
getOddPaths SG {..} = result_list ++ isolatedH
  where
    present = (`Set.notMember` removed)
    initialVertices = filter present . map (Vtx True) $ allGIdx sg_g
    isolatedH = [[v] | v <- filter present . map (Vtx False) $ allGIdx sg_h, null (edges MMap.! v)]
    result_list =
      snd $
        foldr getOddPathWithVertex (Set.empty, []) initialVertices
    getOddPathWithVertex v (seen, vs) =
      if
        | v `Set.member` seen -> (seen, vs)
        | isOdd -> (seen', vertices : vs)
        | otherwise -> (seen', vs)
      where
        isOdd = length vertices `mod` 2 == 1
        (vertices, seen') =
          case edges MMap.! v of
            [] -> ([v], aux_seen)
            [e] ->
              let (verticesRest, aux_seen') = getPath (eOther e) aux_seen
               in (v : verticesRest, aux_seen')
            [e1, e2] ->
              let (vertices1, aux_seen') = getPath (eOther e1) aux_seen
                  (vertices2, aux_seen'') = getPath (eOther e2) aux_seen'
               in (reverse vertices1 ++ (v : vertices2), aux_seen'')
            _ -> error patternError
          where
            aux_seen = Set.insert v seen

    getPath v aux_seen =
      if v `Set.member` aux_seen
        then ([], aux_seen)
        else case filter ((`Set.notMember` aux_seen) . eOther) (edges MMap.! v) of
          [] -> ([v], aux_seen')
          [e] ->
            let (verticesRest, aux_seen'') = getPath (eOther e) aux_seen'
             in (v : verticesRest, aux_seen'')
          _ -> error patternError
      where
        aux_seen' = Set.insert v aux_seen

-- check if element is rare
isRare :: (Chromosome c1, Chromosome c2, GIdx c1 ~ GIdx c2) => SampleGraph m c1 c2 (GIdx (ChromList c1)) -> Vertex (GIdx (ChromList c1)) -> Bool
isRare SG {..} v = if vInG v then occG <= occH else occH <= occG
  where
    occG = occurrence sg_posG gene
    occH = occurrence sg_posH gene
    gene =
      if vInG v
        then getGene (vIdx v) sg_g
        else getGene (vIdx v) sg_h

-- return true if the black edges from the witnesses share and edge
checkConflict :: (Eq i) => BlockWitness i -> BlockWitness i -> Bool
checkConflict bw bw' = idxG == idxG' || idxH == idxH'
  where
    (idxG, idxH) = w_black_pair bw
    (idxG', idxH') = w_black_pair bw'

-- | Steps from a forward to b in the chromosome (mod size on circular chromosomes)
fwdOffset :: (Chromosome c, GenIdx (GIdx c)) => c -> GIdx c -> GIdx c -> Int
fwdOffset chr a b = if isCircular chr then d `mod` size chr else d
  where
    d = idxToInt b - idxToInt a

-- return true if the black edges from the witnesses are parallel
checkParallel :: (FptCtx c1 c2) => ChromList c1 -> ChromList c2 -> BlockWitness (CIdx, GIdx c1) -> BlockWitness (CIdx, GIdx c1) -> Bool
checkParallel g h bw bw' = aligned && any matchH offsetsG
  where
    aligned = chrG == chrG' && chrH == chrH' && w_inverted bw == w_inverted bw'
    cG = getChromosome chrG g
    cH = getChromosome chrH h
    -- signed offsets of idxG' from idxG that are inside the witness bw
    ext_l = fwdOffset cG idxG_l idxG
    ext_r = fwdOffset cG idxG idxG_r
    dG = fwdOffset cG idxG idxG'
    offsetsG = filter (\t -> -ext_l <= t && t <= ext_r) (if isCircular cG then [dG, dG - size cG] else [dG])
    -- idxH' must be at the same offset from idxH (in the opposite direction if inverted)
    matchH t =
      let s = if w_inverted bw then -t else t
          dH = fwdOffset cH idxH idxH'
       in if isCircular cH then (dH - s) `mod` size cH == 0 else dH == s
    ((chrG, idxG), (chrH, idxH)) = w_black_pair bw
    ((_, idxG_l), _) = w_left_pair bw
    ((_, idxG_r), _) = w_right_pair bw
    ((chrG', idxG'), (chrH', idxH')) = w_black_pair bw'

-- check if every vertex has degree at least low and at most up, and it is connect only
-- with vertices of the other genome. In case of black edges, the degree must be
-- exactly 1. Removed genes have no edges.
checkDegrees :: (FptCtx c1 c2) => Int -> Int -> SampleGraph m c1 c2 (GIdx (ChromList c1)) -> Bool
checkDegrees low up SG {..} =
  all (checkDegree . Vtx True) (allGIdx sg_g)
    && all (checkDegree . Vtx False) (allGIdx sg_h)
  where
    checkDegree v =
      let l = edges MMap.! v
          (low', up')
            | v `Set.member` removed = (0, 0)
            | any eBlack l = (1, 1)
            | otherwise = (low, up)
       in low' <= length l
            && length l <= up'
            && not (any ((vInG v ==) . vInG . eOther) l)

-- | Candidate witnesses of a vertex u of an odd path: the witnesses (with each
-- possible window) of the pairs of u and a gene v of the other genome with the same
-- label, where v is not in a black edge (checkConflict) nor removed, that are not
-- parallel to a witness of the sample.
candidateWitnesses :: (FptCtx c1 c2, Matcher m c1 c2) => SampleGraph m c1 c2 (GIdx (ChromList c1)) -> Vertex (GIdx (ChromList c1)) -> [BlockWitness (GIdx (ChromList c1))]
candidateWitnesses sg@(SG {..}) v = filter (not . checkBWPairs) getNewEdges
  where
    checkBWPairs bw = any (\bw' -> checkParallel sg_g sg_h bw bw' || checkConflict bw bw') sample
    getNewEdges =
      concatMap (findAllLimits sg) $
        if vInG v
          then (vIdx v,) <$> getVertexInQ False sg_h (getGene (vIdx v) sg_g)
          else (,vIdx v) <$> getVertexInQ True sg_g (getGene (vIdx v) sg_h)
    -- return idx of Q including a vertex with label equal to the one of gene
    getVertexInQ inG q gene =
      filter
        ( \idx ->
            canonicOri (getGene idx q) == canonicOri gene
              && Vtx inG idx `Set.notMember` removed
        )
        (allGIdx q)

createSampleGraph :: (FptCtx c1 c2, Matcher m c1 c2) => m c1 c2 -> ChromList c1 -> ChromList c2 -> [(GIdx (ChromList c1), GIdx (ChromList c2))] -> SampleGraph m c1 c2 (GIdx (ChromList c1))
createSampleGraph matcher g h = foldr addPair sg0
  where
    -- a pair parallel to (or in conflict with) a witness of the sample is already
    -- covered by it, adding it would split that witness in two blocks
    addPair idxs sg =
      let bw = findLimits sg idxs
       in if any (\bw' -> checkParallel g h bw bw' || checkConflict bw bw') (sample sg)
            then sg
            else addBlackEdge bw sg
    sg0 = SG matcher g h posMapG posMapH [] MMap.empty Set.empty
    posMapG = positionMap g
    posMapH = positionMap h

-- | Whether the genes at the pair of positions match in exactly one direction
genesMatchOneWay :: (FptCtx c1 c2, Matcher m c1 c2) => m c1 c2 -> ChromList c1 -> ChromList c2 -> (GIdx (ChromList c1), GIdx (ChromList c2)) -> Bool
genesMatchOneWay matcher g h ((cg, i), (ch, j)) = isDirectMatch matcher a b /= isReverseMatch matcher a b
  where
    a = subGenome i i (getChromosome cg g)
    b = subGenome j j (getChromosome ch h)

-- Convert pair of indices into a block witness by finding the limits of
-- pairs parallel to it. If there are several witnesses (genes that match both directly
-- and in reverse, e.g. the caps 0, or several windows in a circular chromosome), the
-- first one is used.
findLimits :: (FptCtx c1 c2, Matcher m c1 c2) => SampleGraph m c1 c2 (GIdx (ChromList c1)) -> (GIdx (ChromList c1), GIdx (ChromList c2)) -> BlockWitness (GIdx (ChromList c1))
findLimits sg idxs = case findAllLimits sg idxs of
  bw : _ -> bw
  [] -> error patternError -- the pair must match, directly or in reverse

-- | Witnesses of a pair of indices, for each direction in which the genes match (both
-- directions for genes that match directly and in reverse, such as the caps 0), and
-- for each maximal window when a witness reaches the size of a circular chromosome
-- (see findLimitsDir), so that the branching considers all of them.
findAllLimits :: (FptCtx c1 c2, Matcher m c1 c2) => SampleGraph m c1 c2 (GIdx (ChromList c1)) -> (GIdx (ChromList c1), GIdx (ChromList c2)) -> [BlockWitness (GIdx (ChromList c1))]
findAllLimits sg ((chr_idx_g, i0), (chr_idx_h, j0)) =
  [bw | isDirectMatch (sg_matcher sg) a b, bw <- findLimitsDir sg True ((chr_idx_g, i0), (chr_idx_h, j0))]
    ++ [bw | isReverseMatch (sg_matcher sg) a b, bw <- findLimitsDir sg False ((chr_idx_g, i0), (chr_idx_h, j0))]
  where
    a = subGenome i0 i0 (getChromosome chr_idx_g (sg_g sg))
    b = subGenome j0 j0 (getChromosome chr_idx_h (sg_h sg))

-- | Witnesses of a pair of indices whose genes match directly (dirH0) or in reverse.
-- A witness extends the pair to the right and to the left while the duos match, but it
-- can not cover more than a whole circular chromosome (max_steps). If both extensions
-- together pass that limit, there are several maximal windows with max_steps steps
-- (rotations of the circular chromosome) and all of them are returned. Otherwise the
-- witness is unique.
findLimitsDir :: (FptCtx c1 c2, Matcher m c1 c2) => SampleGraph m c1 c2 (GIdx (ChromList c1)) -> Bool -> (GIdx (ChromList c1), GIdx (ChromList c2)) -> [BlockWitness (GIdx (ChromList c1))]
findLimitsDir SG {..} dirH0 ((chr_idx_g, i0), (chr_idx_h, j0)) =
  [ BW ((chr_idx_g, idxGL), (chr_idx_h, idxHL)) ((chr_idx_g, i0), (chr_idx_h, j0)) ((chr_idx_g, idxGR), (chr_idx_h, idxHR)) (not dirH0)
    | (r, l) <- windows,
      let (idxGR, idxHR, _) = findLimit True dirH0 r 0 (i0, j0)
          (idxGL, idxHL, _) = findLimit False (not dirH0) l 0 (i0, j0)
  ]
  where
    max_steps = min (size chrG) (size chrH) - 1
    (_, _, max_r) = findLimit True dirH0 max_steps 0 (i0, j0)
    (_, _, max_l) = findLimit False (not dirH0) max_steps 0 (i0, j0)
    -- steps to the right and to the left of each maximal window
    windows
      | max_r + max_l <= max_steps = [(max_r, max_l)]
      | otherwise = [(r, max_steps - r) | r <- [max_steps - max_l .. max_r]]
    chrG = getChromosome chr_idx_g sg_g
    chrH = getChromosome chr_idx_h sg_h

    -- The directions (dir) indicate what is the next gene to look at.
    -- If dir is False we are going from right to left if it is True we are going from
    -- left to right.
    -- k is the number of steps taken, at most budget
    findLimit dirG dirH budget k (i, j) =
      case (step dirG chrG i, step dirH chrH j) of
        (Just i', Just j')
          | k < budget && matchStep i' j' -> findLimit dirG dirH budget (k + 1) (i', j')
        _ -> (i, j, k)
      where
        step d chr = if d then nextIdx chr else prevIdx chr
        matchStep i' j' =
          let (a1, a2) = if dirG then (i, i') else (i', i)
              (b1, b2) = if dirH then (j, j') else (j', j)
              test = if dirG == dirH then isDirectMatch else isReverseMatch
           in test sg_matcher (subGenome a1 a2 chrG) (subGenome b1 b2 chrH)

-- Given a correspondence between a gene of G and a gene of H, include in the graph
-- the correspondent black edge and the edges parallel to that correspondence.
-- Also remove any edge that is no longer valid.
addBlackEdge :: (FptCtx c1 c2) => BlockWitness (GIdx (ChromList c1)) -> SampleGraph m c1 c2 (GIdx (ChromList c1)) -> SampleGraph m c1 c2 (GIdx (ChromList c1))
addBlackEdge bw sg = sg {sample = sample', edges = edges'''}
  where
    (idxG0, idxH0) = w_black_pair bw
    sample' = bw : sample sg
    vG0 = Vtx True idxG0
    vH0 = Vtx False idxH0
    edges_clean = edges (cleanEdges vH0 . cleanEdges vG0 $ sg)
    edges' =
      MMap.insert vH0 (Edge vG0 (idxG0, idxH0) bw True) $
        MMap.insert vG0 (Edge vH0 (idxG0, idxH0) bw True) edges_clean
    edges'' = addEdges True (not (w_inverted bw)) (idxG0, idxH0) edges'
    edges''' = addEdges False (w_inverted bw) (idxG0, idxH0) edges''

    -- The directions (dir) indicate what is the next gene to look at.
    -- If dir is False we are going from right to left if it is True we are going from
    -- left to right.
    addEdges dirG dirH (idxG, idxH) aux_edges =
      if not dirG && idxG == fst (w_left_pair bw) || dirG && idxG == fst (w_right_pair bw) || hasBlackEdge vG || hasBlackEdge vH || isRemoved vG || isRemoved vH
        then aux_edges
        else addEdges dirG dirH (idxG', idxH') aux_edges'
      where
        hasBlackEdge v = let l = aux_edges MMap.! v in not (null l) && eBlack (localHead l)
        isRemoved v = v `Set.member` removed sg
        -- the witness limits are reached before the end of a linear chromosome
        idxG' = fromMaybe (error logicError) $ stepIdx (sg_g sg) dirG idxG
        idxH' = fromMaybe (error logicError) $ stepIdx (sg_h sg) dirH idxH
        vG = Vtx True idxG'
        vH = Vtx False idxH'
        aux_edges' =
          MMap.insert vH (Edge vG (idxG', idxH') bw False) $
            MMap.insert vG (Edge vH (idxG', idxH') bw False) aux_edges

-- | Next (dir True) or previous index in the chromosome of idx
stepIdx :: (Chromosome c) => ChromList c -> Bool -> (CIdx, GIdx c) -> Maybe (CIdx, GIdx c)
stepIdx q dir (c, i) = (c,) <$> (if dir then nextIdx else prevIdx) (getChromosome c q) i

-- | Put a gene in an exclusive block: remove the vertex and the green edges that
-- are no longer valid (the ones with the gene between them and their parent).
removeVertex :: (FptCtx c1 c2) => Vertex (GIdx (ChromList c1)) -> SampleGraph m c1 c2 (GIdx (ChromList c1)) -> SampleGraph m c1 c2 (GIdx (ChromList c1))
removeVertex v sg = sg' {removed = Set.insert v (removed sg')}
  where
    sg' = cleanEdges v sg

-- Remove any green edge that is no longer valid after the inclusion of a new black
-- edge incident on v or after the deletion of v.
cleanEdges :: (FptCtx c1 c2) => Vertex (GIdx (ChromList c1)) -> SampleGraph m c1 c2 (GIdx (ChromList c1)) -> SampleGraph m c1 c2 (GIdx (ChromList c1))
cleanEdges v sg@(SG {..}) = sg {edges = edges'}
  where
    edges' = foldr delGreens edges (edges MMap.! v)
    delGreens edge = delFromVertex (eOther edge) . delFromVertex v
      where
        bw' = eParent edge
        (idxG, idxH) = w_black_pair bw'
        delFromVertex u_aux = go u_aux
          where
            -- move away from the black pair of bw'
            forward
              | vInG u_aux = isAhead sg_g idxG (fst (w_right_pair bw'))
              | otherwise = isAhead sg_h idxH (snd (if w_inverted bw' then w_left_pair bw' else w_right_pair bw'))
            -- u_aux is between the black pair (excluded) and the limit going forward
            isAhead q (c, i0) (_, lim) =
              let chr = getChromosome c q
                  d = fwdOffset chr i0 (snd (vIdx u_aux))
               in 1 <= d && d <= fwdOffset chr i0 lim
            -- limit of the witness in the direction of the walk
            limit
              | vInG u_aux = fst (if forward then w_right_pair bw' else w_left_pair bw')
              | otherwise = snd (if forward /= w_inverted bw' then w_right_pair bw' else w_left_pair bw')
            -- stop at the limit of the witness (when it covers a whole circular
            -- chromosome, the next gene is the other limit) and before the black pair
            -- (reachable only after going around a circle)
            move u
              | vIdx u == limit = Nothing
              | otherwise = do
                  u' <- Vtx (vInG u) <$> stepIdx' (vInG u) forward (vIdx u)
                  if vIdx u' == (if vInG u' then idxG else idxH) then Nothing else Just u'
            stepIdx' inG = if inG then stepIdx sg_g else stepIdx sg_h
            go u aux_edges =
              if length e_list == length e_list'
                then aux_edges
                else maybe aux_edges' (`go` aux_edges') (move u)
              where
                e_list = aux_edges MMap.! u
                e_list' = filter (\e -> eParent e /= bw') e_list
                aux_edges' = MMap.fromMap . Map.insert u e_list' . MMap.toMap $ aux_edges

-- Make sure that there is at most one edge incident in each pair of vertices
-- so that we have a correspondence between the rare genes of G and H.
-- At this point the graph should have only black edges and even green paths, so
-- we eliminate green edges to ensure the correspondence.
finalizeCorrespondence :: (FptCtx c1 c2) => SampleGraph m c1 c2 (GIdx (ChromList c1)) -> SampleGraph m c1 c2 (GIdx (ChromList c1))
finalizeCorrespondence sg = sg2
  where
    -- the removed genes have no edges
    verticesG = filter (`Set.notMember` removed sg) . map (Vtx True) . allGIdx . sg_g $ sg
    (seen1, sg1) = foldr changePaths (Set.empty, sg) verticesG
    sg2 = snd $ foldr changeCycle (seen1, sg1) verticesG

    changeCycle v (seen, sg_aux) =
      if v `Set.member` seen
        then (seen, sg_aux)
        else case edges sg_aux MMap.! v of
          [e1, e2] ->
            let e =
                  if (eOther e1 `Set.member` seen)
                    || (fst (w_black_pair . eParent $ e2) < fst (eVpair e2))
                      && (eOther e2 `Set.notMember` seen)
                    then e2
                    else e1
             in changePathsEven (eOther e) (seen, sg_aux)
          _ -> error patternError

    changePaths v (seen, sg_aux) =
      if v `Set.member` seen
        then (seen, sg_aux)
        else case edges sg_aux MMap.! v of
          [_, _] -> (seen, sg_aux) -- not an extreme or is part of a cycle, so we ignore for now
          [e] -> changePathsEven (eOther e) (aux_seen, sg_aux)
          _ -> error patternError
      where
        aux_seen = Set.insert v seen

    changePathsOdd v (seen, sg_aux) =
      if v `Set.member` seen
        then (seen, sg_aux)
        else case edges sg_aux MMap.! v of
          [e] -> changePathsEven (eOther e) (aux_seen, sg_aux)
          _ -> error patternError
      where
        aux_seen = Set.insert v seen

    changePathsEven v (seen, sg_aux) =
      case edges sg_aux MMap.! v of
        [_] -> (aux_seen, sg_aux) -- path has only one vertex remaining
        [e1, e2] ->
          let u =
                if eOther e1 `Set.member` seen
                  then eOther e2
                  else eOther e1
           in changePathsOdd u (aux_seen, delete_edge u v sg_aux)
        _ -> error patternError
      where
        aux_seen = Set.insert v seen
        delete_edge u1 u2 sg' = sg' {edges = updateEdgeInMap u1 u2 . updateEdgeInMap u2 u1 $ edges sg'}
        updateEdgeInMap u1 u2 m =
          MMap.fromMap . Map.insert u1 (delEdge [] (m MMap.! u1)) . MMap.toMap $ m
          where
            delEdge _ [] = error patternError
            delEdge acc (e : es) = if eOther e == u2 then reverse acc ++ es else delEdge (e : acc) es
