{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE TypeOperators #-}

-- \|
-- Module      : PApprox
-- Description : Approximation algorithm for partition problems. The approximation
-- factor depends on the maximum occurrence of a gene in the genomes.
-- Copyright   : (c) Gabriel Siqueira, 2023
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
module PApprox
  ( approxPartition,
    approxPartitionPhase2,
    approxLowerBound,
    countPerfect,
    getConnectedComponents,
    mkBlockMatchGraph,
  )
where

import Data.Foldable (foldrM)
import Data.IntMap (IntMap)
import Data.IntMap qualified as IntMap
import Data.IntSet (IntSet)
import Data.IntSet qualified as IntSet
import Data.Map (Map)
import Data.List qualified as List
import Data.Map qualified as Map
import Data.Maybe (fromMaybe, isJust, listToMaybe)
import Data.Set (Set)
import Data.Set qualified as Set
import Genomes (CIdx, ChromList, Chromosome (isCircular), FlipMatcher (FlipMatcher), GenIdx (idxToInt, mkIdx), Gene, Genome (GIdx, allGIdx, getGene, size, subGenome), Matcher (isDirectMatch, isMatch, isReverseMatch), MultiChromosome (getChromosome), isCompatibleWithSubGenome)
import LocalBase
import Partition (BlockDelim (..), CommonPartition, Partition, blockDelims, blocks, breakpoints, delimToBlock, findUncommon, getBlocksFit, mkCommonPartition, mkPartitionFromBreakpoints, trivialPartition, underlineGenome)

type ApproxCtx c1 c2 = (Chromosome c1, Chromosome c2, GIdx c1 ~ GIdx c2, Ord (GIdx c1), GenIdx (GIdx c1))

type Correspondents i = [Vertex i]

-- | Position of an adjacency: genome (True for G), chromosome and index of its first
-- gene (a breakpoint at this position separates that gene from the next one)
type AdjPos i = (Bool, CIdx, i)

-- | Classes of the compatibility graph G_0 of the initial partitions (each chromosome is
-- a block). The class of an adjacency is the component of G_0 that contains it (seen as
-- a substring with two genes).
data Classes i = Classes
  { -- | class of each adjacency
    clsOf :: Map (AdjPos i) Int,
    -- | component of G_0 of each vertex of a component of Tmin(G, H)
    clsTminComp :: Map (Vertex i) Int,
    -- | classes of the set B_C of each component C of Tmin(G, H)
    clsB :: IntMap [Int]
  }

-- | A substring of a block. Its genes are the one at vIdx and the following
-- vSize - 1 genes of the chromosome (wrapping around the origin of circular
-- chromosomes). If vWhole is True, the vertex is a whole circular chromosome
-- without breakpoints (so it also has the adjacency between its last and first genes).
data Vertex i = Vtx
  { vSize :: Int,
    vInG :: Bool,
    vChr :: CIdx,
    vIdx :: i,
    vWhole :: Bool
  }
  deriving (Eq, Ord, Show)

-- | Size and circularity of the chromosome of a vertex
geomIn :: (Chromosome c) => ChromList c -> Vertex i -> (Int, Bool)
geomIn q v = let chr = getChromosome (vChr v) q in (size chr, isCircular chr)

geomOf :: (ApproxCtx c1 c2) => ChromList c1 -> ChromList c2 -> Vertex (GIdx c1) -> (Int, Bool)
geomOf g h v = if vInG v then geomIn g v else geomIn h v

-- | Position p moved o genes forward (backward if o is negative) in a chromosome
shiftPos :: (Int, Bool) -> Int -> Int -> Int
shiftPos (n, circ) p o = if circ then (p - 1 + o) `mod` n + 1 else p + o

-- | Number of genes from a forward to b in a chromosome
offsetOf :: (Int, Bool) -> Int -> Int -> Int
offsetOf (n, circ) a b = if circ then (b - a) `mod` n else b - a

-- | Index of the gene in offset o of the vertex
posIn :: (Chromosome c, GenIdx (GIdx c)) => ChromList c -> Vertex (GIdx c) -> Int -> GIdx c
posIn q v o = mkIdx $ shiftPos (geomIn q v) (idxToInt (vIdx v)) o

-- | Subgenome correspondent to a vertex
pieceIn :: (Chromosome c, GenIdx (GIdx c)) => ChromList c -> Vertex (GIdx c) -> c
pieceIn q v = if vWhole v then chr else subGenome (vIdx v) (posIn q v (vSize v - 1)) chr
  where
    chr = getChromosome (vChr v) q

vSubGenome :: (ApproxCtx c1 c2) => ChromList c1 -> ChromList c2 -> Vertex (GIdx c1) -> Either c1 c2
vSubGenome g h v = if vInG v then Left (pieceIn g v) else Right (pieceIn h v)

-- | Position of the adjacency between the genes in offsets o and o + 1 of the vertex
adjPosAt :: (ApproxCtx c1 c2) => ChromList c1 -> ChromList c2 -> Vertex (GIdx c1) -> Int -> AdjPos (GIdx c1)
adjPosAt g h v o = (vInG v, vChr v, if vInG v then posIn g v o else posIn h v o)

-- | Vertex of the compatibility graph that is the adjacency in offset o of the vertex:
-- the substring with the two genes, or the chromosome itself for a whole circular
-- chromosome with one gene (whose adjacency joins the gene to itself)
adjVertexAt :: (ApproxCtx c1 c2) => ChromList c1 -> ChromList c2 -> Vertex (GIdx c1) -> Int -> Vertex (GIdx c1)
adjVertexAt g h v o
  | vWhole v && vSize v == 1 = v
  | otherwise = let (_, c, i) = adjPosAt g h v o in Vtx 2 (vInG v) c i False

-- | Whether a vertex of the compatibility graph is an adjacency (see adjVertexAt)
isAdjVertex :: Vertex i -> Bool
isAdjVertex v = if vWhole v then vSize v == 1 else vSize v == 2

-- | Offsets of the adjacencies of a vertex (a breakpoint in offset o separates o and o + 1)
adjOffsets :: Vertex i -> [Int]
adjOffsets v = [0 .. vSize v - (if vWhole v then 1 else 2)]

sameSeq :: Vertex i -> Vertex i -> Bool
sameSeq v u = vInG v == vInG u && vChr v == vChr u

-- | If u starts inside v (after the start of v) and continues after the end of v,
-- returns the size of their intersection (only if they share an adjacency)
rightOverlap :: (ApproxCtx c1 c2) => ChromList c1 -> ChromList c2 -> Vertex (GIdx c1) -> Vertex (GIdx c1) -> Maybe Int
rightOverlap g h v u
  | v == u || not (sameSeq v u) || vWhole v || vWhole u = Nothing
  | 0 < d && d + vSize u > vSize v && vSize v - d >= 2 = Just (vSize v - d)
  | otherwise = Nothing
  where
    d = offsetOf (geomOf g h v) (idxToInt (vIdx v)) (idxToInt (vIdx u))

-- | All vertices of the blocks of a partition (substrings with at least two genes)
partVertices :: (Chromosome c, GenIdx (GIdx c), Eq (GIdx c)) => Bool -> Partition c -> [Vertex (GIdx c)]
partVertices inG p = concatMap delimVertices (blockDelims p)
  where
    q = underlineGenome p
    delimVertices d =
      let chr = getChromosome (b_chr d) q
          geom = (size chr, isCircular chr)
          n = size chr
          beg = idxToInt (b_beg d)
          len = offsetOf geom beg (idxToInt (b_end d)) + 1
          mk s pos = Vtx s inG (b_chr d) (mkIdx pos) False
       in if b_whole d
            -- a whole circular chromosome is a vertex even with one gene (it may need a breakpoint)
            then Vtx n inG (b_chr d) (mkIdx 1) True : [mk s pos | pos <- [1 .. n], s <- [2 .. n]]
            else [mk s (shiftPos geom beg o) | o <- [0 .. len - 2], s <- [2 .. len - o]]

-- Graph where each vertex is a block and an edge connects two vertices if they are correspondents, also has a space to record the number of the connected components containing each graph
data BlockMatchGraph c1 c2 = BG (Partition c1) (Partition c2) (Map (Vertex (GIdx c1)) (Maybe Int, Correspondents (GIdx c1)))

data Component i = Component
  { cVtx :: Vertex i,
    cContG :: Int,
    cContH :: Int,
    cPerfect :: Bool
  }
  deriving (Show)

-- count the number of perfect components
countPerfect :: IntMap (Component i) -> Int
countPerfect = sum . map (fromEnum . cPerfect) . IntMap.elems

approxPartition :: (Orientable c1, Orientable c2, ApproxCtx c1 c2, Matcher m c1 c2) => m c1 c2 -> ChromList c1 -> ChromList c2 -> CommonPartition c1 c2
approxPartition matcher g h = fst (approxPartitionPhase2 matcher g h)

-- | The partition of the algorithm and the number of times its second stage broke a
-- block D (it is only used with relations that are not symmetric and transitive)
approxPartitionPhase2 :: (Orientable c1, Orientable c2, ApproxCtx c1 c2, Matcher m c1 c2) => m c1 c2 -> ChromList c1 -> ChromList c2 -> (CommonPartition c1 c2, Int)
approxPartitionPhase2 matcher g h = let (pg, ph, _, p2) = runApprox matcher g h in (mkCommonPartition matcher pg ph, p2)

-- | Lower bound on the cost of a minimum common partition: the number of components
-- whose sets B_C were activated by the algorithm (they do not share adjacencies, so a
-- common partition has a distinct breakpoint in each of them).
approxLowerBound :: (Orientable c1, Orientable c2, ApproxCtx c1 c2, Matcher m c1 c2) => m c1 c2 -> ChromList c1 -> ChromList c2 -> Int
approxLowerBound matcher g h = let (_, _, q, _) = runApprox matcher g h in q

-- | Runs the algorithm, returning the partitions, the number of components whose sets B_C
-- were activated, and the number of blocks D broken by the second stage
runApprox :: (Orientable c1, Orientable c2, ApproxCtx c1 c2, Matcher m c1 c2) => m c1 c2 -> ChromList c1 -> ChromList c2 -> (Partition c1, Partition c2, Int, Int)
runApprox matcher g h = (pg'', ph'', q'', p2'')
  where
    pg_0 = trivialPartition g
    ph_0 = trivialPartition h
    bmg_0@(BG _ _ corres_0) = mkBlockMatchGraph matcher pg_0 ph_0
    classes = uncurry mkClasses (getConnectedComponents bmg_0)
    -- The compatibility of two substrings does not depend on the partitions, so the graph of
    -- any partitions is the subgraph of the graph of the initial partitions (G_0) induced by
    -- the substrings without breakpoints
    mkGraph = subgraphOf corres_0

    (pg'', ph'', _, q'', p2'', _) = loop (pg_0, ph_0, IntSet.empty, 0, 0, False)
    -- repeat until the partition is common, stopping with an error if an iteration
    -- adds no breakpoint (otherwise the next iteration would repeat it forever)
    loop st
      | done2 st = st
      | noProgress st st' = error (logicError ++ " (approxPartition: no breakpoint can be added, but the partition is not common)")
      | otherwise = loop st'
      where
        st' = tryToBuildPartition st
    noProgress (pg, ph, _, _, _, _) st'@(pg', ph', _, _, _, _) =
      not (done2 st') && breakpoints pg == breakpoints pg' && breakpoints ph == breakpoints ph'
    -- first stage: add breakpoints until Tmin is empty; second stage (only with relations
    -- that are not symmetric and transitive): break one block without correspondent
    tryToBuildPartition (pg, ph, active, q, p2, _) = (pg2, ph2, active'', q', p2', isCommon)
      where
        (active', q', (_, BG pg' ph' _)) = until done1 (addBreaks classes) (active, q, getConnectedComponents (mkGraph pg ph))
        done1 (_, _, (comps', _)) = countPerfect comps' == IntMap.size comps'
        (pg2, ph2, active'', isCommon) = compatibilityCompletion matcher pg' ph' classes active'
        -- the second stage breaks at most one block D
        broke = Set.size (breakpoints pg2) + Set.size (breakpoints ph2) > Set.size (breakpoints pg') + Set.size (breakpoints ph')
        p2' = if broke then p2 + 1 else p2
    done2 (_, _, _, _, _, isCommon) = isCommon

-- | Add a breakpoint in a component C of Tmin. If a substring of C has an adjacency of an
-- active class, the breakpoint is added there. Otherwise, the classes of B_C are
-- activated and the breakpoint is added in an adjacency of B_C of a substring of C.
-- The state has the active classes and the number of activations of sets B_C.
addBreaks :: (ApproxCtx c1 c2) => Classes (GIdx c1) -> (IntSet, Int, (IntMap (Component (GIdx c1)), BlockMatchGraph c1 c2)) -> (IntSet, Int, (IntMap (Component (GIdx c1)), BlockMatchGraph c1 c2))
addBreaks classes (active, q, (comps, BG pg ph corres)) = (active', q', (comps', bmg'))
  where
    -- The substrings with fewest adjacencies in components that do not admit a perfect
    -- match are substrings of components of Tmin (a substring of a vertex has fewer
    -- adjacencies; the number of genes is not enough, as a whole circular chromosome with
    -- n genes contains substrings with n genes)
    numAdjs v_ = if vWhole v_ then vSize v_ else vSize v_ - 1
    minAdjs = minimum (map (numAdjs . fst) notPerfs)
    cand = filter ((== minAdjs) . numAdjs . fst) notPerfs
    notPerfs = filter isNotPerfect . Map.toAscList $ corres
    isNotPerfect (_, (Nothing, _)) = error patternError
    isNotPerfect (_, (Just c, _)) = not . cPerfect $ comps IntMap.! c
    sorted_cand = sortWith sortKey cand
    -- the component C is the one of the first candidate, and vsC are its substrings
    c_id = fst . snd $ localHead sorted_cand
    vsC = [v_ | (v_, (c, _)) <- sorted_cand, c == c_id]

    classAt v_ o = Map.lookup (adjPosAt g h v_ o) (clsOf classes)
    adjsIn ks = [(v_, o) | v_ <- vsC, o <- adjOffsets v_, Just k <- [classAt v_ o], k `elem` ks]
    -- classes of B_C, if C is a component of Tmin(G, H)
    bC = Map.lookup (localHead vsC) (clsTminComp classes) >>= (`IntMap.lookup` clsB classes)
    (v, break_off, active', q') =
      case listToMaybe [(v_, o) | v_ <- vsC, o <- adjOffsets v_, Just k <- [classAt v_ o], k `IntSet.member` active] of
        Just (v_, o) -> (v_, o, active, q)
        Nothing -> case bC of
          Just ks | (v_, o) : _ <- adjsIn ks -> (v_, o, foldr IntSet.insert active ks, q + 1)
          -- this does not happen while all the breakpoints are in active classes (C is
          -- then a component of Tmin(G, H)); the class of the breakpoint is activated so
          -- that this remains true
          _ -> let v_ = localHead vsC in (v_, 0, maybe active (`IntSet.insert` active) (classAt v_ 0), q)

    -- We give priority to cut a substring in the genome with more vertices of the
    -- component and break ties with the degree of the vertex.
    sortKey (v_, (Just c_id_, neigs)) = (posWeight, dg)
      where
        c = comps IntMap.! c_id_
        posWeight =
          if
            | cContG c == cContH c -> (1 :: Int)
            | (cContG c > cContH c) == vInG v_ -> 0
            | otherwise -> 2
        dg = length neigs
    sortKey (_, (_, _)) = error patternError

    g = underlineGenome pg
    h = underlineGenome ph
    (comps', bmg') = removeAdjacency pg' ph' (adjPosAt g h v break_off) (comps, BG pg ph corres)
    bg = breakpoints pg
    bh = breakpoints ph
    pg' = mkPartitionFromBreakpoints g bg'
    ph' = mkPartitionFromBreakpoints h bh'
    (bg', bh') =
      if vInG v
        then (Set.insert (vChr v, posIn g v break_off) bg, bh)
        else (bg, Set.insert (vChr v, posIn h v break_off) bh)

-- | Get a set with all connected componets of the Block Match Graph and
-- record the component indices in the Graph
getConnectedComponents :: (Ord (GIdx c1)) => BlockMatchGraph c1 c2 -> (IntMap (Component (GIdx c1)), BlockMatchGraph c1 c2)
getConnectedComponents (BG pg ph corres) = (comps, BG pg ph corres')
  where
    (comps, corres') = labelComponents 0 (Map.keys corres) IntMap.empty corres

-- | Remove from the graph (with its components) the vertices that contain the adjacency
-- in position p of the genome (True for G), which are the vertices with a breakpoint
-- after p is added to the partitions pg' and ph'. Only the components of these vertices
-- change, and only they are computed again.
removeAdjacency :: (ApproxCtx c1 c2) => Partition c1 -> Partition c2 -> AdjPos (GIdx c1) -> (IntMap (Component (GIdx c1)), BlockMatchGraph c1 c2) -> (IntMap (Component (GIdx c1)), BlockMatchGraph c1 c2)
removeAdjacency pg' ph' (inG, chr, i) (comps, BG pg ph corres) = (comps', BG pg' ph' corres')
  where
    g = underlineGenome pg
    h = underlineGenome ph
    hasAdj v =
      vInG v == inG
        && vChr v == chr
        && (vWhole v || let d = offsetOf (geomOf g h v) (idxToInt (vIdx v)) (idxToInt i) in 0 <= d && d <= vSize v - 2)
    removed = Map.filterWithKey (\v _ -> hasAdj v) corres
    affected = IntSet.fromList [c | (Just c, _) <- Map.elems removed]
    kept = corres `Map.difference` removed
    affectedVs = [v | (v, (Just c, _)) <- Map.toList kept, c `IntSet.member` affected]
    -- the vertices of the affected components lose their component and the removed neighbours
    corres_ =
      foldr (Map.adjust (\(_, ns) -> (Nothing, filter (`Map.notMember` removed) ns))) kept affectedVs
    next_id = maybe 0 ((+ 1) . fst) (IntMap.lookupMax comps)
    (comps', corres') = labelComponents next_id affectedVs (IntMap.withoutKeys comps affected) corres_

-- | Label the components (with identifiers from comp_id_0) of the given vertices, which
-- must not have a component yet, and record them
labelComponents :: (Ord i) => Int -> [Vertex i] -> IntMap (Component i) -> Map (Vertex i) (Maybe Int, Correspondents i) -> (IntMap (Component i), Map (Vertex i) (Maybe Int, Correspondents i))
labelComponents comp_id_0 vertices comps_0 corres = (comps, corres')
  where
    (comps, corres', _) = foldr checkComponent (comps_0, corres, comp_id_0) vertices

    checkComponent v (comps', corres'', comp_id) =
      case fst (corres'' Map.! v) of
        Just _ -> (comps', corres'', comp_id)
        Nothing ->
          let (corres''', vs, cont_G, cont_H) = visitComp comp_id v (corres'', [], 0, 0)
              perfect = testPerfectMatch corres''' vs cont_G cont_H
              comp = Component v cont_G cont_H perfect
              comps'' = IntMap.insert comp_id comp comps'
           in (comps'', corres''', comp_id + 1)

    visitComp comp_id v (aux_corres, vs, cont_G, cont_H) =
      let (maybe_comp, neigs) = aux_corres Map.! v
       in case maybe_comp of
            Just _ -> (aux_corres, vs, cont_G, cont_H)
            Nothing ->
              let aux_corres' = Map.insert v (Just comp_id, neigs) aux_corres
                  (cont_G', cont_H') =
                    if vInG v
                      then (cont_G + 1, cont_H)
                      else (cont_G, cont_H + 1)
               in foldr (visitComp comp_id) (aux_corres', v : vs, cont_G', cont_H') neigs

    -- Given a correspondence and a list of vertices from a component c,
    -- checks if c admits a perfect match.
    -- Note that we check both directions if the blocks in G can be assigned to blocks
    -- of H and if blocks of H can be assigned to blocks of G.
    testPerfectMatch aux_corres allVertices cont_G cont_H =
      cont_G == cont_H && not (null hVertices) && isJust (foldrM (selectBlock Set.empty) Map.empty hVertices)
      where
        -- select only vertices from H, to avoid double checking
        hVertices = filter (not . vInG) allVertices
        -- Select block from G that corresponds to block represented by v.
        -- fixCorres saves a map between blocks of G and their fix correspondent
        -- in H.
        selectBlock seen v fixCorres =
          if success then Just fixCorres' else Nothing
          where
            (_, fixCorres', success) = selectBlock' seen v fixCorres

        selectBlock' seen v = testCandidates v neigs seen
          where
            (_, neigs) = aux_corres Map.! v

        testCandidates _ [] seen fixCorres = (seen, fixCorres, False)
        testCandidates v (u : us) seen fixCorres =
          if
            | u `Set.member` seen -> testCandidates v us seen fixCorres
            | u `Map.notMember` fixCorres -> (seen', fixCorres', True)
            | success -> (seen_rec, fixCorres_rec', True)
            | otherwise -> testCandidates v us seen_rec fixCorres_rec
          where
            seen' = Set.insert u seen
            fixCorres' = Map.insert u v fixCorres
            (seen_rec, fixCorres_rec, success) = selectBlock' seen' (fixCorres Map.! u) fixCorres
            fixCorres_rec' = Map.insert u v fixCorres_rec

-- | Construct the block match graph, which is represented by an adjacency list decorated -- with information about the component (to be filled by getConnectedComponents later).
-- As we are considering balanced genomes, we ignore subgenomes of size 1 (except whole circular chromosomes).
mkBlockMatchGraph :: (ApproxCtx c1 c2, Matcher m c1 c2) => m c1 c2 -> Partition c1 -> Partition c2 -> BlockMatchGraph c1 c2
mkBlockMatchGraph matcher pg ph = BG pg ph corres
  where
    g = underlineGenome pg
    h = underlineGenome ph
    vsG = map (\v -> (v, pieceIn g v)) (partVertices True pg)
    vsH = map (\v -> (v, pieceIn h v)) (partVertices False ph)
    -- Compatible substrings have the same genes, up to orientation and rotation, so the
    -- match is only tested between substrings with the same key (same size, both whole
    -- circular chromosomes or not, and the same sorted list of genes without orientation)
    key v sub = (vSize v, vWhole v, List.sort [canonicOri (getGene i sub) :: Gene | i <- allGIdx sub])
    byKeyH = Map.fromListWith (flip (++)) [(key u sub, [(u, sub)]) | (u, sub) <- vsH]
    corresG = [(v, [u | (u, sub_h) <- Map.findWithDefault [] (key v sub_g) byKeyH, isMatch matcher sub_g sub_h]) | (v, sub_g) <- vsG]
    corresH = Map.fromListWith (flip (++)) [(u, [v]) | (v, us) <- corresG, u <- us]
    corres =
      Map.fromList $
        [(v, (Nothing, us)) | (v, us) <- corresG]
          ++ [(u, (Nothing, Map.findWithDefault [] u corresH)) | (u, _) <- vsH]

-- | The block match graph of the partitions pg and ph, given the neighbours of each vertex
-- in the graph of the initial partitions: it is the subgraph induced by the vertices of the
-- blocks of pg and ph, which are all vertices of the initial graph
subgraphOf :: (ApproxCtx c1 c2) => Map (Vertex (GIdx c1)) (Maybe Int, Correspondents (GIdx c1)) -> Partition c1 -> Partition c2 -> BlockMatchGraph c1 c2
subgraphOf corres_0 pg ph = BG pg ph corres
  where
    vs = partVertices True pg ++ partVertices False ph
    present = Set.fromList vs
    corres = Map.fromList [(v, (Nothing, filter (`Set.member` present) (snd (corres_0 Map.! v)))) | v <- vs]

-- | Classes of the compatibility graph of the initial partitions (with its components)
-- and the sets B_C of the components C of Tmin(G, H), defined by the cases of the lemma
-- tmin_comp of the text:
--   case i: C does not intersect other components of Tmin; B_C is the class of any
--     adjacency of C;
--   case ii: some adjacency of C is an intersection edge with respect to every component
--     that intersects C; B_C is the class of one of them;
--   case iii: B_C has the classes of the first and last adjacencies of C (every component
--     that intersects C contains one of them in the intersection).
mkClasses :: (Orientable c1, Orientable c2, ApproxCtx c1 c2) => IntMap (Component (GIdx c1)) -> BlockMatchGraph c1 c2 -> Classes (GIdx c1)
mkClasses comps (BG pg ph corres) = Classes classOf tminComp (IntMap.fromSet bOf tminSet)
  where
    g = underlineGenome pg
    h = underlineGenome ph

    -- class of the adjacency in offset o of v: the component of that adjacency
    bAt v o = case fst (corres Map.! adjVertexAt g h v o) of
      Just c -> c
      Nothing -> error patternError
    classOf =
      Map.fromList
        [ ((vInG u, vChr u, vIdx u), c)
          | (u, (Just c, _)) <- Map.toList corres,
            isAdjVertex u
        ]
    tminComp = Map.fromList tminVtxs

    -- Vertices from components that do not admit a perfect match
    notPerfs = filter (\(_, c) -> not . cPerfect $ comps IntMap.! c) . map getComp . Map.toAscList $ corres
    getComp (_, (Nothing, _)) = error patternError
    getComp (v, (Just c, _)) = (v, c)

    -- Vertices from components of Tmin (minimal components without perfect matching)
    tminVtxs = filter (\(_, c) -> c `IntSet.member` tminSet) notPerfs
    -- Pairs (v, u) where u starts inside v and continues after its end, with the size of the intersection
    -- (on a circular chromosome both (v, u) and (u, v) may appear)
    -- (only vertices of the same chromosome can overlap)
    vertexPairsTmin =
      [ ((v, cv), (u, cu), s)
        | vs <- Map.elems (Map.fromListWith (flip (++)) [((vInG v, vChr v), [(v, cv)]) | (v, cv) <- tminVtxs]),
          (v, cv) <- vs,
          (u, cu) <- vs,
          cv /= cu,
          Just s <- [rightOverlap g h v u]
      ]
    -- The components without perfect matching that do not contain (in a vertex) a
    -- smaller vertex of a component without perfect matching
    tminSet = IntSet.fromList (map snd notPerfs) IntSet.\\ IntSet.fromList [c | (v, c) <- notPerfs, containsNotPerf v]
    notPerfSet = Set.fromList (map fst notPerfs)
    -- a vertex contains another vertex (in the sense of contains) of a component without
    -- perfect matching if its prefix or suffix with one gene less is such a vertex or
    -- contains one; a whole circular chromosome contains its substrings with all genes
    containsNotPerf v
      | vWhole v = or [hasNotPerf (Vtx (vSize v) (vInG v) (vChr v) (mkIdx i) False) | vSize v >= 2, i <- [1 .. vSize v]]
      | vSize v <= 2 = False
      | otherwise =
          let (_, _, i1) = adjPosAt g h v 1
           in hasNotPerf (v {vSize = vSize v - 1}) || hasNotPerf (v {vSize = vSize v - 1, vIdx = i1})
    hasNotPerf v = hasNotPerfMap Map.! v
    hasNotPerfMap = Map.mapWithKey (\v _ -> v `Set.member` notPerfSet || containsNotPerf v) corres

    -- For each component, the intersection Map stores the smallest size of the intersection
    -- between components on the left (lv) and on the right (rv), and the set of components
    -- that intersect it.
    interMap = foldr updateInterMap interMap_0 vertexPairsTmin
    interMap_0 = IntMap.fromAscList . map (,((Nothing, Nothing), IntSet.empty)) . IntSet.toAscList $ tminSet

    -- v intersects u on the right of v (and on the left of u)
    updateInterMap ((v, cv), (u, cu), s) =
      updateEntry cu cv (vOri u) False s . updateEntry cv cu (vOri v) True s
      where
        vOri v' = either getOri getOri (vSubGenome g h v')

        -- Record in component c' an intersection of size s' with component c'', on the right
        -- (onRight) or on the left of the vertex with orientation ori. The side is flipped
        -- when ori differs from the orientation of the vertex representing c', since the
        -- blocks of a component may be inverted.
        updateEntry c' c'' ori onRight s' interMap_aux =
          IntMap.insert c' (sides', IntSet.insert c'' inter) interMap_aux
          where
            ((l, r), inter) = interMap_aux IntMap.! c'
            right = onRight == (ori == vOri (cVtx $ comps IntMap.! c'))
            select old = Just (maybe s' (min s') old)
            sides' = if right then (l, select r) else (select l, r)

    -- Intersection edges (as components of adjacencies) of c with respect to d: the
    -- components of the adjacencies shared by a vertex of c and a vertex of d (the
    -- text also includes the adjacencies of these components, which are the same sets)
    sharedMap :: Map (Int, Int) (Set Int)
    sharedMap =
      Map.fromListWith Set.union $
        concat
          [ [((cv, cu), ids), ((cu, cv), ids)]
            | ((v, cv), (_, cu), s) <- vertexPairsTmin,
              let ids = Set.fromList [bAt v o | o <- [vSize v - s .. vSize v - 2]]
          ]
    intersecting c = IntSet.toList (snd (interMap IntMap.! c))
    repOf c = cVtx (comps IntMap.! c)
    edgeIds c = Set.fromList [bAt (repOf c) o | o <- adjOffsets (repOf c)]
    -- Adjacencies of c (as components) that are intersection edges with respect to every
    -- component that intersects c (case ii of the lemma if not empty)
    case2Cands c = foldr (\d -> Set.intersection (Map.findWithDefault Set.empty (c, d) sharedMap)) (edgeIds c) (intersecting c)

    bOf c
      | null (intersecting c) = [bAt v 0]
      | not (Set.null cands) = [pick]
      | otherwise = IntSet.toList (IntSet.fromList [bAt v 0, bAt v (vSize v - 2)])
      where
        v = repOf c
        cands = case2Cands c
        -- the adjacency common (by position) to all the intersections is preferred
        positional = case fst (interMap IntMap.! c) of
          (Just l, Nothing) -> Just (bAt v (l - 2))
          (Nothing, Just r) -> Just (bAt v (vSize v - r))
          (Just l, Just r) | l - 1 > vSize v - r -> Just (bAt v (vSize v - r))
          _ -> Nothing
        pick = case positional of
          Just b | b `Set.member` cands -> b
          _ -> Set.findMin cands

-- | Given two partitions @pg@ and @ph@, with all components admitting a perfect matching,
-- if there are blocks without correspondents, take the largest of them (B, of either
-- partition) and break a block D of the other partition around a substring compatible
-- with B. Only one block is broken, so that the caller rebuilds the compatibility graph
-- and makes Tmin empty again before the next one (as in the algorithm of the text).
-- The classes of the new breakpoints are activated, so that every breakpoint is in an
-- active class (this stage does not add breakpoints with symmetric and transitive
-- relations, and the text does not guarantee an approximation factor otherwise).
-- The boolean in the result is True if the partitions are already common.
compatibilityCompletion :: (ApproxCtx c1 c2, Matcher m c1 c2) => m c1 c2 -> Partition c1 -> Partition c2 -> Classes (GIdx c1) -> IntSet -> (Partition c1, Partition c2, IntSet, Bool)
compatibilityCompletion matcher pg ph classes active
  | null unmatchedG && null unmatchedH = (pg, ph, active, True)
  | otherwise = case attempt of
      Just (pg', ph', cuts) -> (pg', ph', foldr IntSet.insert active [k | pos <- cuts, Just k <- [Map.lookup pos (clsOf classes)]], False)
      Nothing -> (pg, ph, active, False)
  where
    -- blocks without a pair in a maximum matching (as the genomes are balanced, the
    -- partitions are common if and only if there are none)
    fit = getBlocksFit matcher pg ph
    bsG = blocks pg
    bsH = blocks ph
    unmatchedG = [b | (k, b) <- zip [0 ..] bsG, k `Map.notMember` fit]
    unmatchedH = [b | (k, b) <- zip [0 :: Int ..] bsH, k `notElem` Map.elems fit]
    largest :: (Genome c) => [c] -> c
    largest = foldr1 (\x y -> if size x >= size y then x else y)
    -- B is the largest block without correspondent; the blocks of an alternating path
    -- that could not be matched (findUncommon) are tried next if B does not lead to a D
    attempt
      | not (null unmatchedG) && (null unmatchedH || size (largest unmatchedG) >= size (largest unmatchedH)) =
          findBreakpoint matcher False (largest unmatchedG : fromMaybe [] (findUncommon matcher pg ph)) pg ph
      | otherwise =
          (\(a, b, c) -> (b, a, c))
            <$> findBreakpoint (FlipMatcher matcher) True (largest unmatchedH : fromMaybe [] (findUncommon (FlipMatcher matcher) ph pg)) ph pg

-- | Find a block b2 of p2 that contains a substring compatible with a block of uncommon, break b2
-- around that substring, and break a block b1 of p1 that matches b2 at the correspondent positions.
-- Returns the new partitions and the positions of the new breakpoints (inG2 tells whether
-- p2 is a partition of G).
findBreakpoint :: (ApproxCtx c1 c2, Matcher m c1 c2) => m c1 c2 -> Bool -> [c1] -> Partition c1 -> Partition c2 -> Maybe (Partition c1, Partition c2, [AdjPos (GIdx c1)])
findBreakpoint matcher inG2 uncommon p1 p2 =
  -- B is the first block of uncommon; the others are only used if B gives no D
  breakAt <$> listToMaybe (candidatesFor (take 1 uncommon) ++ candidatesFor (drop 1 uncommon))
  where
    candidatesFor us =
      [ (d2, b2, idxToInt i - 1, size u)
        | (d2, b2) <- blocks2,
          u <- us,
          size u < size b2 || b_whole d2 && size u == size b2,
          Just i <- [isCompatibleWithSubGenome matcher u b2]
      ]
    -- b2 is broken into the offsets [o, o + k - 1] (of b2) and the rest
    breakAt (d2, b2, o, k) =
      let len = size b2
          (d1, t, dir) =
            fromMaybe (error logicError) . listToMaybe $
              [ (d1', t', dir')
                | (d1', b1) <- blocks1,
                  size b1 == len,
                  isMatch matcher b1 b2,
                  Just (t', dir') <- [alignment d1' b1 d2 b2]
              ]
          -- offset in b1 (from t) correspondent to offset o of b2
          a = if dir then o else len - o - k
          cuts1 = blockCuts q1 (b_chr d1) t (b_whole d1) len a k
          cuts2 = blockCuts q2 (b_chr d2) (b_beg d2) (b_whole d2) len o k
       in ( mkPartitionFromBreakpoints q1 (Set.union (breakpoints p1) (Set.fromList cuts1)),
            mkPartitionFromBreakpoints q2 (Set.union (breakpoints p2) (Set.fromList cuts2)),
            [(not inG2, c, i) | (c, i) <- cuts1] ++ [(inG2, c, i) | (c, i) <- cuts2]
          )
    q1 = underlineGenome p1
    q2 = underlineGenome p2
    blocks1 = [(d, delimToBlock q1 d) | d <- blockDelims p1]
    blocks2 = [(d, delimToBlock q2 d) | d <- blockDelims p2]

    -- Start of b1 (index in its chromosome) and direction such that b1 read from that start
    -- in that direction matches b2 read from its first gene
    alignment d1 b1 d2 b2
      | not (b_whole d1) && not (b_whole d2) = Just (b_beg d1, isDirectMatch matcher b1 b2)
      | b_whole d1 && b_whole d2 = listToMaybe [(t, dir) | t <- allGIdx b1, dir <- [True, False], testRotation t dir]
      | otherwise = Nothing
      where
        n = size b2
        testRotation t dir =
          let t_last = mkIdx $ shiftPos (n, True) (idxToInt t) (n - 1)
              test = if dir then isDirectMatch else isReverseMatch
           in test matcher (subGenome t t_last b1) (subGenome (mkIdx 1) (mkIdx n) b2)
                && test matcher (subGenome t_last t b1) (subGenome (mkIdx n) (mkIdx 1) b2)

-- | Breakpoints that separate the offsets [a, a + k - 1] of a block starting at index s
-- (of chromosome c) with the given length from the rest of the block
blockCuts :: (Chromosome c, GenIdx (GIdx c)) => ChromList c -> CIdx -> GIdx c -> Bool -> Int -> Int -> Int -> [(CIdx, GIdx c)]
blockCuts q c s whole len a k = before ++ after
  where
    chr = getChromosome c q
    pos o = mkIdx $ shiftPos (size chr, isCircular chr) (idxToInt s) o
    before = [(c, pos (a - 1)) | a > 0 || whole]
    after = [(c, pos (a + k - 1)) | if whole then k < len else a + k < len]
