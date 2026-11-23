{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE ExplicitForAll #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UndecidableInstances #-}

-- |
-- Module      : Genomes
-- Description : Representation of a Genome. The representation may include information
-- regarding intergenic regions and gene orientation.
-- Copyright   : (c) Gabriel Siqueira, 2023
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
module Genomes
  ( Genome (..),
    Chromosome (..),
    MultiChromosome (..),
    IntergenicChromosome (..),
    RigidIntergenicChromosome (..),
    RigidIntergenicMultiChromGenome (..),
    Gene,
    mkGene,
    geneToInt,
    IR,
    irToInt,
    Idx,
    CIdx,
    Sign (..),
    ChromType (..),
    GenIdx (..),
    occurrence,
    GeneMap,
    emptyGeneMap,
    positionMap,
    geneMapLookup,
    geneMapAdjust,
    geneMapUnion,
    geneMapDifference,
    geneMapDelete,
    geneMapIntersection,
    occurrenceMax,
    singletonOnBoth,
    ChromList (..),
    GenesIRsR,
    GenesIRsF,
    toMC,
    mkRGenome,
    mkFGenome,
    mkRGenome0,
    mkFGenome0,
    writeIR,
    readRGenome,
    readMCRGenome,
    writeRGenome,
    readFGenome,
    readMCFGenome,
    writeFGenome,
    writeMultiC,
    flexibilize,
    flexibilizeMC,
    equalizeLinear,
    Matcher (..),
    FlipMatcher (..),
    RigidRigidDirectMatcher (..),
    RigidFlexibleDirectMatcher (..),
    RigidRigidReverseMatcher (..),
    RigidFlexibleReverseMatcher (..),
    isCompatibleWithSubGenome,
    randomGenome,
    randomGenomeWithReplicas,
    shuffleGenome,
  )
where

import Control.Exception (assert)
import Control.Monad.Random (MonadRandom, getRandomRs, getRandoms)
import Control.Monad.ST (ST)
import Data.ByteString.Builder (intDec, toLazyByteString)
import Data.ByteString.Char8 qualified as BS
import Data.ByteString.Lazy qualified as LBS
import Data.Coerce (coerce)
import Data.Hashable (Hashable)
import Data.IntMap (IntMap)
import Data.IntMap qualified as IntMap
import Data.List (intersperse)
import Data.List qualified as List
import Data.Set (Set)
import Data.Set qualified as Set
import Data.Vector (Vector, (!))
import Data.Vector qualified as Vec
import Data.Vector.Mutable qualified as MVec
import LocalBase
import System.Random (Random)
import System.Random.Shuffle (shuffleM)

newtype Gene = Gene Int deriving newtype (Eq, Show, Read, Hashable, Ord, Num, Bounded, Enum, Random)

mkGene :: Int -> Gene
mkGene = Gene

geneToInt :: Gene -> Int
geneToInt (Gene i) = i

data IR = R Int | F Int Int deriving (Eq, Ord)

irToInt :: IR -> Int
irToInt ir = case ir of
  R i -> i
  F _ _ -> error patternError

sumIR :: IR -> IR -> IR
sumIR (R i) (R j) = R (i + j)
sumIR _ _ = error patternError

class GenIdx i where
  mkIdx :: Int -> i
  idxToInt :: i -> Int
  idxDist :: i -> i -> Int
  incIdx :: i -> i
  incKIdx :: Int -> i -> i
  decIdx :: i -> i
  decKIdx :: Int -> i -> i

  idxDist i1 i2 = idxToInt i2 - idxToInt i1 + 1
  incIdx = mkIdx . (+ 1) . idxToInt
  incKIdx k = fpow k incIdx
  decIdx = mkIdx . subtract 1 . idxToInt
  decKIdx k = fpow k decIdx

newtype Idx = Idx Int deriving newtype (Eq, Show, Read, Hashable, Ord, Num, Bounded, Enum, Integral, Real, Random)

newtype CIdx = CIdx Int deriving newtype (Eq, Show, Read, Hashable, Ord, Num, Bounded, Enum, Integral, Real, Random)

instance GenIdx Idx where
  mkIdx = Idx
  idxToInt (Idx i) = i
  incKIdx k = Idx . (+ k) . idxToInt
  decKIdx k = Idx . subtract k . idxToInt

instance GenIdx CIdx where
  mkIdx = CIdx
  idxToInt (CIdx i) = i
  incKIdx k = CIdx . (+ k) . idxToInt
  decKIdx k = CIdx . subtract k . idxToInt

instance Show IR where
  show (R x) = show x
  show (F x y) = show x ++ ":" ++ show y

instance Orientable Gene where
  getOri a = if a >= 0 then LR else RL
  invOri a = -a

data Sign = Signed | Unsigned deriving (Eq, Show, Enum)

data ChromType = Circular | Linear deriving (Eq, Show, Enum)

-- Representation with a list of genes and a list of intergenic regions
-- (only ints for rigid, or intervals for flexible)
-- We assume that linear chromosomes start and end with genes (which are
-- different from all other genes unless you use turnLinear) and the
-- intergenic region between the first and last gene in the list of a
-- circular chromosomes is the last one in the intergenic region list.
data GenesIRs = GenesIRs ChromType Sign (Vector Gene) (Vector IR) Gene

instance Show GenesIRs where
  show (GenesIRs Linear _ genes irs _) =
    unwords . (("(" ++ localHead str_s ++ ")") :) $
      zipWith (\ir a -> "- " ++ ir ++ " - (" ++ a ++ ")") str_i (skipFirst str_s)
    where
      str_s = Vec.toList $ (\i -> if i == maxBound then "inf" else show i) <$> genes
      str_i = Vec.toList $ show <$> irs
  show (GenesIRs Circular _ genes irs _) =
    unwords $ zipWith (\ir a -> "(" ++ a ++ ") - " ++ ir ++ " - ") str_i str_s
    where
      str_s = Vec.toList $ (\i -> if i == maxBound then "inf" else show i) <$> genes
      str_i = Vec.toList $ show <$> irs

-- Representation with a list of chromosomes.
-- A list of caps is required, so make sure the genes on the extremities of linear chromosomes are always in the caps list.
data ChromList c = ChromList {caps :: [Gene], chromosomes :: [c]}

instance (Show c) => Show (ChromList c) where
  show (ChromList _ chroms) = unwords . intersperse " | " . map show $ chroms

instance (Eq c) => Eq (ChromList c) where
  (ChromList _ g1s) == (ChromList _ g2s) = g1s == g2s

class Genome g where
  -- The indices must have an predefined order
  type GIdx g

  -- Get all indices of the genome in the order
  allGIdx :: g -> [GIdx g]

  -- Get all indices of the genome in the order between two indices
  betweenGIdx :: g -> GIdx g -> GIdx g -> [GIdx g]

  -- First index of the genome in the order. Does not work for empty genomes.
  firstGIdx :: g -> GIdx g
  firstGIdx = localHead . allGIdx

  -- Last index of the genome in the order. Does not work for empty genomes.
  lastGIdx :: g -> GIdx g
  lastGIdx = last . allGIdx

  isGene :: Gene -> g -> Bool

  -- Number of genes
  size :: g -> Int

  -- get gene at position
  getGene :: GIdx g -> g -> Gene

  -- invert gene according to the characteristics of the genome
  -- (takes into account if the genome has known orientation)
  invGene :: g -> Gene -> Gene

  -- set genes at indices.
  -- The signs must be specify on the new genes,
  -- because the original sings will be overwritten
  setGenes :: [GIdx g] -> [Gene] -> g -> g

  -- for indices i j, it gets the subgenome
  -- starting at gene position i
  -- and ending at gene position j
  -- requires that the first index must come before the second
  -- (on circular chromosomes the first index may come after the second,
  -- giving a block that crosses the origin)
  subGenome :: GIdx g -> GIdx g -> g -> g

  -- Return set with all the gene values
  alphabet :: g -> Set Gene

  -- Get a gene that does not appear on the genome (using Int as the
  -- underline representation, this value should be bigger than all
  -- other gene values)
  getNewGene :: g -> Gene

  -- Replace genes in indices with a new singleton gene
  -- (does not have any replica in the genome),
  -- Returns the genomes and the singleton that replace each index.
  -- Also ensures that they are singletons in g2.
  makeSingletons :: (Genome g2) => g2 -> [GIdx g] -> g -> (g, [Gene])

  -- Check if gene orientation is being included as signs
  hasSigns :: g -> Bool

class (Genome g) => Chromosome g where
  isLinear :: g -> Bool
  isCircular :: g -> Bool
  isCircular = not . isLinear
  turnCircular :: g -> g

  -- Turn into a linear chromossome by cutting after the position indicated by the index
  turnLinear :: GIdx g -> g -> g

  -- Combine two linear chromosomes into one putting a intergenic region between them.
  combine :: g -> IR -> g -> g

  -- Get a subgenome with the number of consecutive genes starting at the given index
  subGenomeBySize :: GIdx g -> Int -> g -> Maybe g
  getDuo :: GIdx g -> g -> Maybe g
  getDuo g = subGenomeBySize g 2

  -- Next/previous index in the chromosome, wrapping around
  -- the origin on circular chromosomes (Nothing past a linear end)
  nextIdx :: g -> GIdx g -> Maybe (GIdx g)
  prevIdx :: g -> GIdx g -> Maybe (GIdx g)

  -- All pairs of adjacent indices (includes (n, 1) on circular chromosomes)
  adjacencies :: g -> [(GIdx g, GIdx g)]
  adjacencies g = [(i, j) | i <- allGIdx g, Just j <- [nextIdx g i]]

class (Genome g) => MultiChromosome g where
  type Chrom g

  -- Get information from GIdx
  getIdxOfChrom :: g -> GIdx g -> CIdx
  getIdxInChrom :: g -> GIdx g -> GIdx (Chrom g)

  -- Make index
  setIdx :: g -> CIdx -> GIdx (Chrom g) -> GIdx g

  -- get chromosome at position i (index starts in 1), requires 1 <= i <= numChromosomes g
  getChromosome :: CIdx -> g -> Chrom g

  -- number of chromosomes in the genome
  numChromosomes :: g -> Int

  getChromosomes :: g -> [Chrom g]
  getChromosomes g = map (`getChromosome` g) [1 .. mkIdx (numChromosomes g)]

class (Chromosome g) => IntergenicChromosome g where
  intergenicFullReversal :: g -> g

  -- get the intergenic region after the gene at the given index
  getIR :: GIdx g -> g -> IR

class (IntergenicChromosome g) => RigidIntergenicChromosome g where
  intergenicReversal :: Idx -> Idx -> Int -> Int -> g -> g
  intergenicTransposition :: Idx -> Idx -> Idx -> Int -> Int -> Int -> g -> g
  intergenicInsertion :: Idx -> g -> g -> g
  intergenicDeletion :: Idx -> Idx -> Int -> g -> g

class (MultiChromosome g) => RigidIntergenicMultiChromGenome g where
  -- Cut between (idxi,idxi+1) = (a,b) and (idxj,idxj+1) = (c,d) followed by join of (a,c) and (b,d) or (a,d) and (b,c)
  intergenicDCJ :: GIdx g -> GIdx g -> Bool -> Int -> Int -> g -> g

  -- Cut the chromossome indicated by the first index after the position indicated
  -- by the second index, the cut is done after the given number of nucleotides.
  -- Add caps accordingly. Circular chromosomes are turned linear by this operation.
  cut :: GIdx g -> Int -> g -> g

  -- Join the chromosomes in the given indices (remove caps accordingly).
  -- The first bool indicates the direction (True to move the first and False to move the second).
  -- The second and third bools indicate whether the correspondent chromosome should be inverted.
  join :: CIdx -> CIdx -> Bool -> Bool -> Bool -> g -> g

  -- Turn chromosome at position into a circular chromosome
  turnCircularAtIndex :: CIdx -> g -> g

occurrence :: GeneMap [i] -> Gene -> Int
occurrence geneMap a = maybe 0 length (geneMapLookup a geneMap)

newtype GeneMap a = GM (IntMap a) deriving newtype (Functor, Foldable, Eq)

emptyGeneMap :: GeneMap a
emptyGeneMap = GM IntMap.empty

-- Get map from genes to the positions they appear in the genomes
positionMap :: (Genome g) => g -> GeneMap [GIdx g]
positionMap g = GM . IntMap.fromListWith (++) . map (\i -> (geneToInt . abs $ getGene i g, [i])) $ allGIdx g

geneMapLookup :: Gene -> GeneMap a -> Maybe a
geneMapLookup a (GM gm) = IntMap.lookup (geneToInt $ abs a) gm

geneMapAdjust :: (a -> a) -> Gene -> GeneMap a -> GeneMap a
geneMapAdjust f a (GM gm) = GM $ IntMap.adjust f (geneToInt $ abs a) gm

geneMapUnion :: GeneMap a -> GeneMap a -> GeneMap a
geneMapUnion (GM gm1) (GM gm2) = GM $ IntMap.union gm1 gm2

geneMapIntersection :: GeneMap a -> GeneMap a -> GeneMap a
geneMapIntersection (GM gm1) (GM gm2) = GM $ IntMap.intersection gm1 gm2

geneMapDifference :: GeneMap a -> GeneMap a -> GeneMap a
geneMapDifference (GM gm1) (GM gm2) = GM $ IntMap.difference gm1 gm2

geneMapDelete :: Gene -> GeneMap a -> GeneMap a
geneMapDelete a (GM gm) = GM $ IntMap.delete (geneToInt $ abs a) gm

occurrenceMax :: (Genome g) => g -> Int
occurrenceMax = maximum . fmap length . positionMap

singletonOnBoth :: g -> GeneMap [GIdx g] -> GeneMap [GIdx g] -> Gene -> Bool
-- ^ Check if the gene is a singleton in both genomes
singletonOnBoth _ posMapG posMapH gene =
  (case geneMapLookup gene posMapG of Nothing -> False; Just pos -> length pos == 1)
    && (case geneMapLookup gene posMapH of Nothing -> False; Just pos -> length pos == 1)

instance Genome GenesIRs where
  -- Indices must satisfy 1 <= i <= size g
  type GIdx GenesIRs = Idx

  firstGIdx _ = mkIdx 1
  lastGIdx = mkIdx . size
  allGIdx g = mkIdx <$> [1 .. size g]
  betweenGIdx g beg end
    | beg <= end = [beg .. end]
    | isCircular g = [beg .. lastGIdx g] ++ [firstGIdx g .. end]
    | otherwise = []

  isGene gene (GenesIRs _ sign genes _ _) =
    gene `elem` genes
      || case sign of
        Signed -> (-gene) `elem` genes
        Unsigned -> False
  size (GenesIRs _ _ genes _ _) = Vec.length genes
  getGene idx (GenesIRs _ _ genes _ _) = genes Vec.! (idxToInt idx - 1)

  invGene (GenesIRs _ sign _ _ _) a =
    case sign of Signed -> -a; Unsigned -> a

  setGenes idxs new_genes (GenesIRs chtype sign genes irs new) = GenesIRs chtype sign genes' irs new'
    where
      genes' = Vec.modify update genes
      update :: MVec.MVector s Gene -> ST s ()
      update v =
        mapM_ (\(idx, gene) -> MVec.write v (idxToInt idx - 1) gene) $
          zip idxs new_genes
      new' = foldr (\gene n -> if abs gene >= n then abs gene + 1 else n) new new_genes

  subGenome idxi idxj (GenesIRs chtype sign genes irs new)
    | idxi <= idxj = GenesIRs Linear sign (Vec.slice i n genes) (Vec.slice i (n - 1) irs) new
    -- block crossing the origin of a circular chromosome
    | chtype == Circular = GenesIRs Linear sign (Vec.drop i genes Vec.++ Vec.take j genes) (Vec.drop i irs Vec.++ Vec.take (j - 1) irs) new
    | otherwise = error indexError
    where
      i = idxToInt idxi - 1
      j = idxToInt idxj
      n = j - i

  alphabet (GenesIRs _ _ genes _ _) = Set.fromList . Vec.toList . fmap abs $ genes

  -- In this implementation the new gene is bigger than all other genes except for maxBound
  getNewGene (GenesIRs _ _ _ _ new) = new

  makeSingletons h idxs g@(GenesIRs _ _ _ _ new_g) = (g', new_genes)
    where
      g' = setGenes idxs new_genes g
      new_h = getNewGene h
      new_genes = take (length idxs) [max new_g new_h ..]

  hasSigns (GenesIRs _ sign _ _ _) = sign == Signed

instance Chromosome GenesIRs where
  -- turn into a circular chromosome. The first and last genes (caps) are deleted and the first
  -- and last IRs are combined.
  turnCircular (GenesIRs Linear sign genes irs new) = GenesIRs Circular sign (Vec.init . Vec.tail $ genes) (merge_ir irs) new
    where
      merge_ir v =
        let n = Vec.length v
         in Vec.tail $ v Vec.// [(n - 1, sumIR (v Vec.! 0) (v Vec.! (n - 1)))]
  turnCircular g = g
  turnLinear idx (GenesIRs Circular sign genes irs new) = GenesIRs Linear sign genes' irs' new
    where
      genes' = Vec.fromList . rotateL (idxToInt idx) . Vec.toList $ genes
      irs' = Vec.fromList . tail . rotateL (idxToInt idx - 1) . Vec.toList $ irs
  turnLinear _ (GenesIRs Linear _ _ _ _) = error patternError
  combine (GenesIRs Circular _ _ _ _) _ _ = error logicError
  combine _ _ (GenesIRs Circular _ _ _ _) = error logicError
  combine (GenesIRs Linear sign1 genes1 irs1 new1) ir (GenesIRs Linear sign2 genes2 irs2 new2) =
    if sign1 /= sign2 then error logicError else GenesIRs Linear sign1 genes' irs' new'
    where
      genes' = Vec.concat [genes1, genes2]
      irs' = Vec.concat [irs1, Vec.singleton ir, irs2]
      new' = max new1 new2
  isLinear (GenesIRs Linear _ _ _ _) = True
  isLinear (GenesIRs Circular _ _ _ _) = False

  nextIdx g i
    | i < lastGIdx g = Just (incIdx i)
    | isCircular g = Just (firstGIdx g)
    | otherwise = Nothing
  prevIdx g i
    | i > firstGIdx g = Just (decIdx i)
    | isCircular g = Just (lastGIdx g)
    | otherwise = Nothing

  subGenomeBySize idx_beg n g =
    let idx_end = fpow (n - 1) incIdx idx_beg
     in if
          | n > size g -> Nothing
          | idxToInt idx_end <= size g -> Just $ subGenome idx_beg idx_end g
          | isCircular g -> Just $ subGenome idx_beg (decKIdx (size g) idx_end) g
          | otherwise -> Nothing

instance (Genome c) => Genome (ChromList c) where
  type GIdx (ChromList c) = (CIdx, GIdx c)

  firstGIdx g = (mkIdx 1, firstGIdx . getChromosome (mkIdx 1) $ g)
  lastGIdx g = (mkIdx . numChromosomes $ g, lastGIdx . getChromosome (mkIdx . numChromosomes $ g) $ g)
  allGIdx g = [(i, j) | i <- mkIdx <$> [1 .. numChromosomes g], j <- allGIdx . getChromosome i $ g]
  betweenGIdx g (beg_c, beg_i) (end_c, end_i) =
    if beg_c == end_c
      then (beg_c,) <$> betweenGIdx (getChromosome beg_c g) beg_i end_i
      else
        concat
          [ (beg_c,) <$> betweenGIdx (getChromosome beg_c g) beg_i (lastGIdx . getChromosome beg_c $ g),
            do
              c <- [succ beg_c .. pred end_c]
              i <- allGIdx . getChromosome c $ g
              return (c, i),
            (end_c,) <$> betweenGIdx (getChromosome end_c g) (firstGIdx . getChromosome end_c $ g) end_i
          ]

  isGene gene = any (isGene gene) . chromosomes
  size = sum . map size . chromosomes
  getGene (chr_i, i) = getGene i . getChromosome chr_i
  invGene g = invGene (localHead . chromosomes $ g)

  -- TODO: Improve by combining the indices with the same CIdx
  setGenes idxs genes mc = foldr setGene mc (zip idxs genes)
    where
      setGene ((chr_i, i), gene) =
        let upChr = setGenes [i] [gene]
         in ChromList (caps mc) . replace (idxToInt chr_i - 1) upChr . chromosomes

  subGenome (chr_i, i) (chr_j, j) mc =
    ChromList
      (caps mc)
      ( if chr_i == chr_j
          then [subGenome i j . getChromosome chr_i $ mc]
          else chrI' : mid ++ [chrJ']
      )
    where
      mid = drop (idxToInt chr_i) . take (idxToInt chr_j - 1) . chromosomes $ mc
      chrI = getChromosome chr_i mc
      chrI' = subGenome (firstGIdx chrI) i chrI
      chrJ = getChromosome chr_j mc
      chrJ' = subGenome j (lastGIdx chrJ) chrJ

  alphabet = Set.unions . map alphabet . chromosomes
  getNewGene = maximum . map getNewGene . chromosomes
  makeSingletons h idxs g = (g', new_genes)
    where
      g' = setGenes idxs new_genes g
      new_g = getNewGene g
      new_h = getNewGene h
      new_genes = take (length idxs) [max new_g new_h ..]
  hasSigns = any hasSigns . chromosomes

instance IntergenicChromosome GenesIRs where
  getIR idx (GenesIRs _ _ _ irs _) = irs Vec.! (idxToInt idx - 1)
  intergenicFullReversal (GenesIRs chtype sign vs vi new) = GenesIRs chtype sign vs' vi' new
    where
      vs' = case sign of { Unsigned -> id; Signed -> fmap negate } $ Vec.reverse vs
      vi' = case chtype of
        Linear -> Vec.reverse vi
        -- the IR between the last and first genes stays at the end
        Circular -> if Vec.null vi then vi else Vec.reverse (Vec.init vi) `Vec.snoc` Vec.last vi

instance Orientable GenesIRs where
  getOri g_@(GenesIRs Circular _ _ _ _) = if (genes, irs) <= (genes', irs') then LR else RL
    where
      (GenesIRs _ _ genes irs _) = canonicPosition g_
      (GenesIRs _ _ genes' irs' _) = canonicPosition (invOri g_)
  getOri g@(GenesIRs Linear _ genes irs _) = if (genes, irs) <= (genes', irs') then LR else RL
    where
      (GenesIRs _ _ genes' irs' _) = invOri g
  invOri = intergenicFullReversal

-- Rotate GenesIRs such that after this rotation any equal pair of circular chromosomes will be on the same position
canonicPosition :: GenesIRs -> GenesIRs
canonicPosition g@(GenesIRs Linear _ _ _ _) = g
canonicPosition (GenesIRs Circular sign genes irs new) = GenesIRs Circular sign genes' irs' new
  where
    (genes', irs') = minimum . Vec.iterateN (Vec.length genes) rotate $ (genes, irs)
    rotate (gs, is) = (Vec.tail gs Vec.++ Vec.take 1 gs, Vec.tail is Vec.++ Vec.take 1 is)

newtype GenesIRsR = GLR GenesIRs deriving newtype (Show, Genome, Chromosome, IntergenicChromosome, Orientable)

newtype GenesIRsF = GLF GenesIRs deriving newtype (Show, Genome, Chromosome, IntergenicChromosome, Orientable)

instance Eq GenesIRsR where
  (GLR (GenesIRs Linear _ _ _ _)) == (GLR (GenesIRs Circular _ _ _ _)) = False
  (GLR (GenesIRs Circular _ _ _ _)) == (GLR (GenesIRs Linear _ _ _ _)) = False
  GLR g1 == GLR g2 =
    sign1 == sign2 && genes1 == genes2 && irs1 == irs2
    where
      (GenesIRs _ sign1 genes1 irs1 _) = canonicPosition g1
      (GenesIRs _ sign2 genes2 irs2 _) = canonicPosition g2

toMC :: GenesIRsR -> ChromList GenesIRsR
toMC g' = if isCircular g' then error undefinedError else ChromList [getGene 1 g', getGene (mkIdx (size g')) g'] [g']

-- Make an Rigid Genome with the given chromosome type, genes, and intergenic regions. The Maybe indicates if the genome have to be extended and which integers to use on the extension. The sign indicates if the genes are signed or unsigned.
mkRGenome :: ChromType -> Maybe (Int, Int) -> Sign -> [Int] -> [Int] -> GenesIRsR
mkRGenome chtype ext sign genes_ irs = GLR $ GenesIRs chtype sign (Vec.fromList . coerce $ genes) (Vec.fromList . map R $ irs) (Gene new)
  where
    genes = case ext of
      Just (b, e) | chtype == Linear -> b : (genes_ ++ [e])
      _ -> genes_
    -- bigger than every gene (caps included) except maxBound, which may be used as a cap
    new = maximum (0 : filter (/= maxBound) (map abs genes)) + 1

mkRGenome0 :: ChromType -> Sign -> [Int] -> GenesIRsR
mkRGenome0 chtype sign genes = mkRGenome chtype (Just (0, maxBound)) sign genes (replicate (numIRs chtype True (length genes)) 0)

-- Make an Flexible Genome (arguments are similar as mkRGenome)
mkFGenome :: ChromType -> Maybe (Int, Int) -> Sign -> [Int] -> [(Int, Int)] -> GenesIRsF
mkFGenome chtype ext sign genes_ irs = GLF $ GenesIRs chtype sign (Vec.fromList . coerce $ genes) (Vec.fromList . map (uncurry F) $ irs) (Gene new)
  where
    genes = case ext of
      Just (b, e) | chtype == Linear -> b : (genes_ ++ [e])
      _ -> genes_
    -- bigger than every gene (caps included) except maxBound, which may be used as a cap
    new = maximum (0 : filter (/= maxBound) (map abs genes)) + 1

mkFGenome0 :: ChromType -> Sign -> [Int] -> GenesIRsF
mkFGenome0 chtype sign genes = mkFGenome chtype (Just (0, maxBound)) sign genes (replicate (numIRs chtype True (length genes)) (0, 0))

-- | Number of intergenic regions of a chromosome with the given type and
-- number of genes, where the boolean indicates if a linear chromosome is
-- extended with caps (the caps are not included in the number of genes).
numIRs :: ChromType -> Bool -> Int -> Int
numIRs Circular _ n = n
numIRs Linear ext n = if ext then n + 1 else n - 1

writeIR :: IR -> BS.ByteString
writeIR ir = LBS.toStrict $ case ir of
  R i -> toLazyByteString . intDec $ i
  F l u -> (toLazyByteString . intDec $ l) <> ":" <> (toLazyByteString . intDec $ u)

splitElements :: BS.ByteString -> [BS.ByteString]
splitElements = filter (not . BS.null) . BS.splitWith (\x -> x == ',' || x == ' ')

readMCGenome :: (ChromType -> Maybe (Int, Int) -> Sign -> BS.ByteString -> BS.ByteString -> a) -> Sign -> BS.ByteString -> BS.ByteString -> ChromList a
readMCGenome reader sign bs_genes bs_irs =
  if length bs_genes_split /= length bs_irs_split
    then error inputError
    else ChromList ext chroms
  where
    ext = [0]
    chroms = zipWith3 (\chtype gs irs -> reader chtype (Just (0, 0)) sign gs irs) types bs_genes_split_clean bs_irs_split
    types = map (\bs -> if BS.head bs == 'C' then Circular else Linear) bs_genes_split
    -- writeMultiC ends each chromosome with ';', so skip blank segments
    nonBlank = filter (not . BS.all (== ' '))
    bs_genes_split = map (BS.dropWhile (== ' ')) . nonBlank $ BS.split ';' bs_genes
    bs_genes_split_clean = map (BS.dropWhile (\c -> c == 'C' || c == 'L')) bs_genes_split
    bs_irs_split = nonBlank $ BS.split ';' bs_irs

readMCRGenome :: Sign -> BS.ByteString -> BS.ByteString -> ChromList GenesIRsR
readMCRGenome = readMCGenome readRGenome

readMCFGenome :: Sign -> BS.ByteString -> BS.ByteString -> ChromList GenesIRsF
readMCFGenome = readMCGenome readFGenome

readRGenome :: ChromType -> Maybe (Int, Int) -> Sign -> BS.ByteString -> BS.ByteString -> GenesIRsR
readRGenome chtype extend sign bs_genes bs_irs =
  if invalidGenes sign genes
    || length irs /= numIRs chtype (extend /= Nothing) (length genes) || any (< 0) irs
    then error inputError
    else mkRGenome chtype extend sign genes irs
  where
    genes = map readInputInt . splitElements $ bs_genes
    irs = map readInputInt . splitElements $ bs_irs

writeRGenome :: Bool -> GenesIRsR -> (BS.ByteString, BS.ByteString)
writeRGenome rext g@(GLR (GenesIRs _ _ genes irs _)) =
  ( BS.unwords . fmap (LBS.toStrict . toLazyByteString . intDec . coerce) $ l_genes,
    BS.unwords . fmap (LBS.toStrict . toLazyByteString . intDec) $ l_irs
  )
  where
    l_genes = Vec.toList . (if rext && isLinear g then Vec.slice 1 (size g - 2) else id) $ genes
    l_irs = (\case R i -> i; F _ _ -> error patternError) <$> Vec.toList irs

readFGenome :: ChromType -> Maybe (Int, Int) -> Sign -> BS.ByteString -> BS.ByteString -> GenesIRsF
readFGenome chtype extend sign bs_genes bs_irs =
  if invalidGenes sign genes
    || length irs /= numIRs chtype (extend /= Nothing) (length genes) || any (\(l, u) -> l < 0 || l > u) irs
    then error inputError
    else mkFGenome chtype extend sign genes irs
  where
    genes = map readInputInt . splitElements $ bs_genes
    irs = map readF . splitElements $ bs_irs
    readF fir =
      case BS.splitWith (== ':') fir of
        [lir, uir] -> (readInputInt lir, readInputInt uir)
        _ -> error inputError

-- | Negative genes only make sense on signed genomes
invalidGenes :: Sign -> [Int] -> Bool
invalidGenes sign genes = sign == Unsigned && any (< 0) genes

-- | Read an integer that must take the whole string
readInputInt :: BS.ByteString -> Int
readInputInt bs = case BS.readInt bs of
  Just (x, rest) | BS.null rest -> x
  _ -> error inputError

writeFGenome :: Bool -> GenesIRsF -> (BS.ByteString, BS.ByteString)
writeFGenome rext g@(GLF (GenesIRs _ _ genes irs _)) =
  ( BS.unwords . fmap (LBS.toStrict . toLazyByteString . intDec . coerce) $ l_genes,
    BS.unwords . fmap LBS.toStrict $ l_irs
  )
  where
    l_genes = Vec.toList . (if rext && isLinear g then Vec.slice 1 (size g - 2) else id) $ genes
    l_irs = (\case R _ -> error patternError; F l u -> irToBS l u) <$> Vec.toList irs
    irToBS l u = (toLazyByteString . intDec $ l) <> ":" <> (toLazyByteString . intDec $ u)

writeMultiC :: (MultiChromosome g, Chromosome (Chrom g)) => (Chrom g -> (BS.ByteString, BS.ByteString)) -> g -> (BS.ByteString, BS.ByteString)
writeMultiC writer mc = foldr go ("", "") [1 .. mkIdx (numChromosomes mc)]
  where
    -- foldr: (bsg, bsi) already has chromosomes i + 1 .. n, so chromosome i goes first
    go i (bsg, bsi) =
      let chr = getChromosome i mc
          (chr_bsg, chr_bsi) = writer chr
       in ((if isLinear chr then "L " else "C ") <> chr_bsg <> " ; " <> bsg, chr_bsi <> " ; " <> bsi)

flexibilize :: Int -> GenesIRsR -> GenesIRsF
flexibilize l (GLR (GenesIRs chtype sign genes irs new)) = GLF $ GenesIRs chtype sign genes irs' new
  where
    irs' = (\i -> F (i - (l * i `div` 100)) (i + (l * i `div` 100))) . irToInt <$> irs

-- | Add empty linear chromosomes (only the two telomeres, with an intergenic region of
-- size zero between them) to the genome with fewer linear chromosomes, so that both
-- genomes have the same number of linear chromosomes, as in the definition of the text.
equalizeLinear :: Sign -> ChromList GenesIRsR -> ChromList GenesIRsF -> (ChromList GenesIRsR, ChromList GenesIRsF)
equalizeLinear sign g h = (addEmpty (nh - ng) emptyR g, addEmpty (ng - nh) emptyF h)
  where
    numLinear q = length . filter (not . isCircular) $ chromosomes q
    ng = numLinear g
    nh = numLinear h
    addEmpty k e (ChromList caps_ chrs) = ChromList caps_ (chrs ++ replicate k e)
    emptyR = mkRGenome Linear (Just (0, 0)) sign [] [0]
    emptyF = mkFGenome Linear (Just (0, 0)) sign [] [(0, 0)]

flexibilizeMC :: Int -> ChromList GenesIRsR -> ChromList GenesIRsF
flexibilizeMC l (ChromList caps_ chrs) = ChromList caps_ (map (flexibilize l) chrs)

-- TODO: Adapt operations to allow circular indices
instance RigidIntergenicChromosome GenesIRsR where
  intergenicReversal i j x y (GLR g@(GenesIRs chtype sign vs vi new)) =
    assert (case chtype of { Linear -> 2; Circular -> 1 } <= i)
      . assert (i < j)
      . assert (j <= coerce (case chtype of Linear -> size g - 1; Circular -> size g))
      -- on a circular chromosome, reversing genes 1..n would cut the same IR twice
      . assert (not (chtype == Circular && i == 1 && j == coerce (size g)))
      . assert (0 <= x && x <= ir_x)
      . assert (0 <= y && y <= ir_y)
      . GLR
      $ GenesIRs chtype sign vs' vi' new
    where
      vs' = Vec.modify updateG vs
      vi' = Vec.modify updateIR vi

      updateG :: MVec.MVector s Gene -> ST s ()
      updateG v = do
        mapM_ (\k -> MVec.swap v (coerce $ i + k - 1) (coerce $ j - k - 1)) [0 .. (j - i + 1) `div` 2 - 1]
        case sign of
          Unsigned -> pure ()
          Signed -> mapM_ (MVec.modify v invOri . coerce) [i - 1 .. j - 1]

      updateIR :: MVec.MVector s IR -> ST s ()
      updateIR v = do
        mapM_ (\k -> MVec.swap v (coerce i + coerce k - 1) (coerce j - coerce k - 2)) [0 .. (j - i + 1) `div` 2 - 1]
        MVec.write v ir_i_idx (R (x + y))
        MVec.write v (coerce j - 1) (R (x_rest + y_rest))

      -- IR before gene i (the last IR when i = 1 on a circular chromosome)
      ir_i_idx = (coerce i - 2) `mod` Vec.length vi
      ir_x = irToInt (vi ! ir_i_idx)
      x_rest = ir_x - x
      ir_y = irToInt (vi ! (coerce j - 1))
      y_rest = ir_y - y

  intergenicTransposition i j k x y z (GLR g@(GenesIRs chtype sign vs vi new)) =
    assert (case chtype of { Linear -> 2; Circular -> 1 } <= i)
      . assert (i < j)
      . assert (j < k)
      . assert (k <= coerce (case chtype of Linear -> size g; Circular -> size g + 1))
      -- on a circular chromosome, i = 1 and k = n + 1 would cut the same IR twice
      . assert (not (chtype == Circular && i == 1 && k == coerce (size g + 1)))
      . assert (0 <= x && x <= ir_x)
      . assert (0 <= y && y <= ir_y)
      . assert (0 <= z && z <= ir_z)
      . GLR
      $ GenesIRs chtype sign vs' vi' new
    where
      vs' = Vec.modify updateG vs
      vi' = Vec.modify updateIR vi

      updateG :: MVec.MVector s Gene -> ST s ()
      updateG v =
        do
          aux1 <- MVec.clone . MVec.slice (coerce i - 1) (coerce $ j - i) $ v
          aux2 <- MVec.clone . MVec.slice (coerce j - 1) (coerce $ k - j) $ v
          MVec.move (MVec.slice (coerce i - 1) (coerce $ k - j) v) aux2
          MVec.move (MVec.slice (coerce $ i + k - j - 1) (coerce $ j - i) v) aux1

      updateIR :: MVec.MVector s IR -> ST s ()
      updateIR v = do
        do
          aux1 <- MVec.clone . MVec.slice (coerce i - 1) (coerce $ j - i) $ v
          aux2 <- MVec.clone . MVec.slice (coerce j - 1) (coerce $ k - j) $ v
          MVec.move (MVec.slice (coerce i - 1) (coerce $ k - j) v) aux2
          MVec.move (MVec.slice (coerce $ i + k - j - 1) (coerce $ j - i) v) aux1
        MVec.write v ir_i_idx (R (x + y_rest))
        MVec.write v (coerce $ i + k - j - 2) (R (z + x_rest))
        MVec.write v (coerce k - 1 - 1) (R (y + z_rest))

      -- IR before gene i (the last IR when i = 1 on a circular chromosome)
      ir_i_idx = (coerce i - 2) `mod` Vec.length vi
      ir_x = irToInt (vi ! ir_i_idx)
      x_rest = ir_x - x
      ir_y = irToInt (vi ! (coerce j - 2))
      y_rest = ir_y - y
      ir_z = irToInt (vi ! (coerce k - 2))
      z_rest = ir_z - z

  intergenicDeletion i j x (GLR g@(GenesIRs chtype sign vs vi new)) =
    assert (case chtype of { Linear -> 2; Circular -> 1 } <= i)
      . assert (i < j)
      . assert (j <= coerce (case chtype of Linear -> size g; Circular -> size g + 1))
      -- a circular chromosome can not be deleted entirely
      . assert (not (chtype == Circular && i == 1 && j == coerce (size g + 1)))
      . assert (0 <= x && x <= ir_i + ir_j)
      . GLR
      $ GenesIRs chtype sign vs' vi' new
    where
      vs' =
        Vec.slice 0 (coerce i - 1) vs
          Vec.++ Vec.slice (coerce j - 1) (Vec.length vs - coerce j + 1) vs
      vi'
        -- circular deleting from gene 1: the merged IR is the last one
        | i == 1 = Vec.slice (coerce j - 1) (Vec.length vi - coerce j) vi `Vec.snoc` R x
        | otherwise =
            Vec.slice 0 (coerce i - 2) vi
              Vec.++ Vec.fromList [R x]
              -- IRs after gene j (on a circular chromosome this includes the closing IR)
              Vec.++ Vec.slice (coerce j - 1) (Vec.length vi - coerce j + 1) vi

      ir_i = irToInt (vi ! ((coerce i - 2) `mod` Vec.length vi))
      ir_j = irToInt (vi ! (coerce j - 2))

  -- Genome to be inserted must be open (n genes and n+1 itergenic regions)
  intergenicInsertion i (GLR (GenesIRs _ _ vsx vix _)) (GLR g@(GenesIRs chtype sign vs vi new)) =
    assert (1 <= i)
      . assert (i <= coerce (case chtype of Linear -> size g - 1; Circular -> size g))
      . GLR
      $ GenesIRs chtype sign vs' vi' new'
    where
      vs' =
        Vec.slice 0 (coerce i) vs
          Vec.++ vsx
          Vec.++ Vec.slice (coerce i) (Vec.length vs - coerce i) vs
      new' = Vec.foldr (\gene n -> if abs gene >= n then abs gene + 1 else n) new vsx
      vi' =
        ( if i == 1
            then Vec.empty
            else Vec.slice 0 (coerce i - 1) vi
        )
          Vec.++ vix
          Vec.++ Vec.slice (coerce i) (Vec.length vi - coerce i) vi

instance (Genome g) => MultiChromosome (ChromList g) where
  type Chrom (ChromList g) = g
  getChromosome idx (ChromList _ chroms) = chroms !! (idxToInt idx - 1)
  numChromosomes = length . chromosomes
  getIdxOfChrom _ (idx_c, _) = idx_c
  getIdxInChrom _ (_, idx_i) = idx_i
  setIdx _ idx_c idx_i = (idx_c, idx_i)

instance RigidIntergenicMultiChromGenome (ChromList GenesIRsR) where
  intergenicDCJ (chri, i) (chrj, j) connectAC x y mc =
    assert (1 <= chri)
      . assert (chri <= chrj)
      . assert (chrj <= mkIdx (numChromosomes mc))
      . assert (1 <= i)
      . assert ((if chri == chrj then i else 0) < j)
      . assert (i <= mkIdx (size chr_i))
      . assert (j <= mkIdx (size chr_j))
      . assert (0 <= x && x <= ir_x)
      . assert (0 <= y && y <= ir_y)
      $ if
        -- ...------i--------...------j--------...
        -- ...------i i+1----...------j j+1----...
        -- ...------i j------ i+1----...j+1----...
        -- ...--------------- i+1----...j+1----...
        -- ...---------------...----i+1 j+1----...
        -- ...---------------...---------------...
        | chri /= chrj && connectAC && isLinear chr_i && isLinear chr_j -> join (chri + 1) (chrj + 1) True True False . join chri chrj' False False True $ mc_cutted
        -- ...(-----i-------)...------j--------...
        -- ...(-----i i+1---)...------j j+1----...
        -- ...i+1-----------i...------j j+1----...
        -- ...i+1-----------i j------...j+1----...
        -- ...i+1--------------------...j+1----...
        -- ...----------------------i+1 j+1----...
        -- ...---------------------------------...
        | chri /= chrj && connectAC && isCircular chr_i && isLinear chr_j -> join chri chrj True True False . join chri chrj' False False True $ mc_cutted
        -- ------i--------...(-----j-------)...
        -- ------i i+1----...(-----j j+1---)...
        -- ------i--------...(-----j-------)...
        -- ------i i+1----...(-----j j+1---)...
        -- ------i i+1----...j+1-----------j...
        -- ------i j-----------j+1 i+1------...
        -- --------------------j+1 i+1------...
        -- ---------------------------------...
        | chri /= chrj && connectAC && isLinear chr_i && isCircular chr_j -> join chri (chri + 1) False False False . join chri chrj' False False True $ mc_cutted
        -- (-----i-------)...(-----j-------)...
        -- (-----i i+1---)...(-----j j+1---)...
        -- i+1-----------i...j+1-----------j...
        -- i+1-----------i j-------------j+1...
        -- i+1---------------------------j+1...
        -- (-------------------------------)...
        | chri /= chrj && connectAC && isCircular chr_i && isCircular chr_j -> turnCircularAtIndex chri . join chri chrj' False False True $ mc_cutted
        -- ------i--------...------j--------...
        -- ------i i+1----...------j j+1----...
        -- ------i j+1---- i+1----...------j...
        -- --------------- i+1----...------j...
        -- --------------- ----i+1 j--------...
        -- --------------- -----------------...
        | chri /= chrj && not connectAC && isLinear chr_i && isLinear chr_j -> join (chri + 1) chrj' False True True . join chri (chrj' + 1) False False False $ mc_cutted
        -- (-----i-------)...------j--------...
        -- (-----i i+1---)...------j j+1----...
        -- i+1-----------i...------j j+1----...
        -- i+1-----------i j+1----...------j...
        -- i+1--------------------...------j...
        -- --------------------i+1 j--------...
        -- ---------------------------------...
        | chri /= chrj && not connectAC && isCircular chr_i && isLinear chr_j -> join chri chrj' False True True . join chri (chrj' + 1) False False False $ mc_cutted
        -- ------i--------...(-----j-------)...
        -- ------i i+1----...(-----j j+1---)...
        -- ------i i+1----...j+1-----------j...
        -- ------i j+1-----------j i+1------...
        -- ----------------------j i+1------...
        -- ---------------------------------...
        | chri /= chrj && not connectAC && isLinear chr_i && isCircular chr_j -> join chri (chri + 1) False False False . join chri chrj' False False False $ mc_cutted
        -- (-----i-------)...(-----j-------)...
        -- (-----i i+1---)...(-----j j+1---)...
        -- i+1-----------i...j+1-----------j...
        -- i+1-----------i j+1-------------j...
        -- i+1-----------------------------j...
        -- (-------------------------------)...
        | chri /= chrj && not connectAC && isCircular chr_i && isCircular chr_j -> turnCircularAtIndex chri . join chri chrj' False False False $ mc_cutted
        -- ------i-----------------j--------...
        -- ------i i+1-------------j j+1----...
        -- ------i j-------------i+1 j+1----...
        -- ----------------------i+1 j+1----...
        -- ---------------------------------...
        | chri == chrj && connectAC && isLinear chr_i -> join chri (chri + 1) False False False . join chri (chri + 1) False False True $ mc_cutted
        -- (-----i-----------------j-------)...
        -- (-----i i+1-------------j-------)...
        -- i+1--------------j--------------i...
        -- i+1--------------j j+1----------i...
        -- j--------------i+1 j+1----------i...
        -- (-------------------------------)...
        | chri == chrj && connectAC && isCircular chr_i -> turnCircularAtIndex chri . join chri (chri + 1) True True False $ mc_cutted
        -- ------i-----------------j--------...
        -- ------i i+1-------------j j+1----...
        -- ------i j+1---- i+1-------------j...
        -- --------------- i+1-------------j...
        -- --------------- (i+1-----------j)...
        | chri == chrj && not connectAC && isLinear chr_i -> turnCircularAtIndex (chri + 1) . join chri (chri + 2) False False False $ mc_cutted
        -- (-----i-----------------j-------)...
        -- (-----i i+1-------------j-------)...
        -- i+1--------------j--------------i...
        -- i+1--------------j j+1----------i...
        -- i+1--------------j (j+1--------i)...
        -- (i+1------------j) (j+1--------i)...
        | chri == chrj && not connectAC && isCircular chr_i -> turnCircularAtIndex chri . turnCircularAtIndex (chri + 1) $ mc_cutted
        | otherwise -> error logicError
    where
      chr_i = getChromosome chri mc
      chr_j = getChromosome chrj mc
      mc_cutted = cut (chrj', j') y . cut (chri, i) x $ mc
      j' = if chri == chrj then j - i + 1 else j
      chrj' = if isLinear chr_i then chrj + 1 else chrj
      ir_x = irToInt (getIR i chr_i)
      ir_y = irToInt (getIR j chr_j)

  cut (chr, i) x mc =
    ChromList (caps mc) $
      case splitAt (idxToInt chr - 1) . chromosomes $ mc of
        (_, []) -> error indexError
        (hgs, g : tgs) -> if isLinear g then hgs ++ (g1 : g2 : tgs) else hgs ++ (g' : tgs)
          where
            sig = if hasSigns g then Signed else Unsigned
            cap = mkRGenome Linear Nothing sig [geneToInt . localHead . caps $ mc] []
            ir_idx = getIR i g
            g1 = combine (subGenome 1 i g) (R x) cap
            g2 = combine cap (sumIR ir_idx (R (-x))) (subGenome (i + 1) (mkIdx (size g)) g)
            g' = combine cap (sumIR ir_idx (R (-x))) (combine (turnLinear i g) (R x) cap)
  join chri chrj move_first invi invj mc =
    if chri > chrj
      then join chrj chri move_first invj invi mc
      else
        ChromList (caps mc) $
          let (gs1, gs23) = splitAt (idxToInt chri - 1) . chromosomes $ mc
           in case splitAt (idxToInt (chrj - chri)) gs23 of
                (_, []) -> error indexError
                ([], _) -> error indexError
                (gi : gs2, gj : gs3) -> gs1 ++ (if move_first then gs2 ++ g' : gs3 else g' : gs2 ++ gs3)
                  where
                    g' = combine (subGenome 1 (mkIdx $ size gi' - 1) gi') (R $ irToInt (getIR (mkIdx $ size gi' - 1) gi') + irToInt (getIR 1 gj')) (subGenome 2 (mkIdx $ size gj') gj')
                    gi' = if invi then intergenicFullReversal gi else gi
                    gj' = if invj then intergenicFullReversal gj else gj

  turnCircularAtIndex chri mc = ChromList (caps mc) . replace (idxToInt chri - 1) turnCircular . chromosomes $ mc

class Matcher m g1 g2 where
  isMatch :: m g1 g2 -> g1 -> g2 -> Bool
  isDirectMatch :: m g1 g2 -> g1 -> g2 -> Bool
  isReverseMatch :: m g1 g2 -> g1 -> g2 -> Bool
  areBalanced :: m g1 g2 -> g1 -> g2 -> Bool
  -- | Whether blocks may also be matched after a reversal (reverse common
  -- partition) or only directly (direct common partition)
  allowsReverse :: m g1 g2 -> Bool

newtype FlipMatcher m g1 g2 = FlipMatcher m

data RigidRigidDirectMatcher g1 g2 = RRDM

data RigidFlexibleDirectMatcher g1 g2 = RFDM

data RigidRigidReverseMatcher g1 g2 = RRRM

data RigidFlexibleReverseMatcher g1 g2 = RFRM

instance (Matcher m g1 g2) => Matcher (FlipMatcher (m g1 g2)) g2 g1 where
  isMatch (FlipMatcher matcher) g1 g2 = isMatch matcher g2 g1
  isDirectMatch (FlipMatcher matcher) g1 g2 = isDirectMatch matcher g2 g1
  isReverseMatch (FlipMatcher matcher) g1 g2 = isReverseMatch matcher g2 g1
  areBalanced (FlipMatcher matcher) g1 g2 = areBalanced matcher g2 g1
  allowsReverse (FlipMatcher matcher) = allowsReverse matcher

instance Matcher RigidRigidDirectMatcher GenesIRsR GenesIRsR where
  allowsReverse _ = False
  isMatch = isDirectMatch
  isDirectMatch _ g h = g == h
  isReverseMatch _ _ _ = False
  areBalanced _ (GLR (GenesIRs _ _ genesG irsG _)) (GLR (GenesIRs _ _ genesH irsH _)) =
    List.sort (fmap canonicOri (Vec.toList genesG)) == List.sort (fmap canonicOri (Vec.toList genesH))
      && sum (fmap irToInt irsG) == sum (fmap irToInt irsH)

instance Matcher RigidFlexibleDirectMatcher GenesIRsR GenesIRsF where
  allowsReverse _ = False
  isMatch = isDirectMatch
  isDirectMatch RFDM (GLR (GenesIRs chtypeG _ genesG irsG _)) (GLF (GenesIRs chtypeH _ genesH irsH _)) =
    chtypeG == chtypeH
      && Vec.length genesG == Vec.length genesH
      && case chtypeG of
        Linear -> matchAt genesH irsH
        Circular -> any (\k -> matchAt (rot k genesH) (rot k irsH)) [0 .. Vec.length genesH - 1]
    where
      matchAt gs is = genesG == gs && and (Vec.zipWith check irsG is)
      rot k v = Vec.drop k v Vec.++ Vec.take k v
      check :: IR -> IR -> Bool
      check (R ir) (F irl irr) = (irl <= ir) && (ir <= irr)
      check _ _ = error patternError
  isReverseMatch _ _ _ = False
  areBalanced _ (GLR (GenesIRs _ _ genesG irsG _)) (GLF (GenesIRs _ _ genesH irsH _)) =
    List.sort (fmap canonicOri (Vec.toList genesG)) == List.sort (fmap canonicOri (Vec.toList genesH))
      && sum (fmap irToIntLow irsH) <= sum (fmap irToInt irsG)
      && sum (fmap irToInt irsG) <= sum (fmap irToIntHigh irsH)
    where
      irToIntLow ir = case ir of
        R _ -> error patternError
        F l _ -> l
      irToIntHigh ir = case ir of
        R _ -> error patternError
        F _ r -> r

instance Matcher RigidRigidReverseMatcher GenesIRsR GenesIRsR where
  allowsReverse _ = True
  isMatch RRRM g h = isDirectMatch RRRM g h || isReverseMatch RRRM g h
  isDirectMatch _ = isDirectMatch RRDM
  isReverseMatch _ g = isDirectMatch RRDM (intergenicFullReversal g)
  areBalanced _ = areBalanced RRDM

instance Matcher RigidFlexibleReverseMatcher GenesIRsR GenesIRsF where
  allowsReverse _ = True
  isMatch RFRM g h = isDirectMatch RFRM g h || isReverseMatch RFRM g h
  isDirectMatch _ = isDirectMatch RFDM
  isReverseMatch _ g = isDirectMatch RFDM (intergenicFullReversal g)
  areBalanced _ = areBalanced RFDM

-- | check if @g1@ is compatible with a subchromosome of @g2@ and return the index
-- of the first gene from the first compatible subchromosome
isCompatibleWithSubGenome :: (Matcher m c1 c2, Chromosome c1, Chromosome c2) => m c1 c2 -> c1 -> c2 -> Maybe (GIdx c2)
isCompatibleWithSubGenome matcher g1 g2 = List.find testSubgenome (allGIdx g2)
  where
    testSubgenome beg = maybe False (isMatch matcher g1) (subGenomeBySize beg (size g1) g2)

randomGenome :: (MonadRandom mon) => Bool -> Int -> Int -> Sign -> mon GenesIRsR
randomGenome zeros n lim signed = do
  coins <- getRandoms
  ls <- case signed of
    Unsigned -> take n <$> getRandomRs (1 :: Gene, coerce lim)
    Signed -> zipWith swaps coins . take n <$> getRandomRs (1 :: Gene, coerce lim)
  li <- take (n + 1) <$> if zeros then return (repeat 0) else getRandomRs (1, 100)
  return $ mkRGenome Linear (Just (0, maxBound)) signed (coerce ls) li
  where
    swaps b v = if b then v else invOri v

randomGenomeWithReplicas :: (MonadRandom mon) => Bool -> Int -> Int -> Int -> Int -> Sign -> mon GenesIRsR
randomGenomeWithReplicas zeros n rep low high signed
  -- the repeated genes alone may need up to @rep * high@ positions
  | low > high || rep * high > n = error inputError
  | otherwise = do
      coins <- getRandoms
      occs <- getRandomRs (low, high)
      let ls =
            ( case signed of
                Unsigned -> id
                Signed -> zipWith swaps coins
            )
              . (\l -> l ++ coerce [length l + 1 .. n])
              . concat
              . zipWith replicate occs
              $ [Gene 1 .. Gene rep]
      li <- take (n + 1) <$> if zeros then return (repeat 0) else getRandomRs (1, 100)
      return $ mkRGenome Linear (Just (0, maxBound)) signed (coerce ls) li
  where
    swaps b v = if b then v else invOri v

shuffleGenome :: (MonadRandom mon) => GenesIRsR -> mon GenesIRsR
shuffleGenome g@(GLR (GenesIRs chtype sign ls_ li _)) = do
  coins <- getRandoms
  let s = sum . fmap irToInt $ li
      -- the caps of a linear chromosome stay in place
      (ext, ls) = case chtype of
        Linear -> (Just (coerce (Vec.head ls_), coerce (Vec.last ls_)), Vec.toList (Vec.slice 1 (size g - 2) ls_))
        Circular -> (Nothing, Vec.toList ls_)
  ls' <- case sign of
    Unsigned -> shuffleM ls
    Signed -> zipWith swaps coins <$> shuffleM ls
  -- same number of intergenic regions, with the same total size
  x <- List.sort . take (Vec.length li - 1) <$> getRandomRs (0, s)
  let li' = zipWith (-) (x ++ [s]) (0 : x)
  return $ mkRGenome chtype ext sign (coerce ls') li'
  where
    swaps b v = if b then v else invOri v
