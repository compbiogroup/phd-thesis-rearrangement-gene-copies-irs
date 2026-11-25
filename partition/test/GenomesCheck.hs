{-# LANGUAGE ExistentialQuantification #-}
{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE TypeFamilies #-}

module GenomesCheck (tests, genGenome, genRGenome, genMultiPair, genMultiPairWith, rearrangeGenome, rearrangeAndFlexibilizeGenome, GenomeWrapper (..)) where

import Control.Monad.Random (replicateM)
import Data.ByteString.Char8 qualified as BS
import Data.List qualified as List
import Genomes
import Hedgehog
import Hedgehog.Gen qualified as Gen
import Hedgehog.Range qualified as Range
import LocalBase

genGene :: Gen Int
genGene = Gen.int (Range.linear 1 100)

genIR :: Gen Int
genIR = Gen.int (Range.linear 0 100)

data GenomeWrapper = forall g. (Genome g, Show g, GIdx g ~ Idx) => GW g

instance Show GenomeWrapper where
  show (GW g) = show g

-- Genomes will have size at least 3 (extensions plus at least one gene)
genGenome :: Int -> Gen GenomeWrapper
genGenome size_lim = do
  flex <- Gen.bool
  if flex
    then GW <$> genFGenome size_lim
    else GW <$> genRGenome size_lim

-- | Generate an empty rigid genome
genRGenome :: Int -> Gen GenesIRsR
genRGenome = genRGenomeWithSign Signed

-- | Generate an empty flexible genome
genFGenome :: Int -> Gen GenesIRsF
genFGenome = genFGenomeWithSign Signed

genRGenomeWithSign :: Sign -> Int -> Gen GenesIRsR
genRGenomeWithSign sign size_lim = do
  n <- Gen.int (Range.linear 3 size_lim)
  coins <- Gen.list (Range.singleton n) Gen.bool
  genes <-
    ( case sign of
        Signed -> zipWith swaps coins
        Unsigned -> id
      )
      <$> Gen.list (Range.singleton $ n - 2) genGene
  irs <- Gen.list (Range.singleton $ n - 1) genIR
  return $ mkRGenome Linear (Just (0, maxBound)) sign genes irs
  where
    swaps b v = if b then v else -v

genFGenomeWithSign :: Sign -> Int -> Gen GenesIRsF
genFGenomeWithSign sign size_lim = do
  n <- Gen.int (Range.linear 3 size_lim)
  coins <- Gen.list (Range.singleton n) Gen.bool
  genes <-
    ( case sign of
        Signed -> zipWith swaps coins
        Unsigned -> id
      )
      <$> Gen.list (Range.singleton $ n - 2) genGene
  irs_low <- Gen.list (Range.singleton $ n - 1) genIR
  porc <- Gen.list (Range.singleton $ n - 1) $ Gen.int (Range.linear 1 50)
  let irs_high = zipWith (\ir p -> ir + p * ir `div` 100) irs_low porc
  let irs = zip irs_low irs_high
  return $ mkFGenome Linear (Just (0, maxBound)) sign genes irs
  where
    swaps b v = if b then v else -v

-- | Generate a circular rigid chromosome (n genes and n intergenic regions)
genCircRGenome :: Int -> Gen GenesIRsR
genCircRGenome size_lim = do
  n <- Gen.int (Range.linear 3 size_lim)
  coins <- Gen.list (Range.singleton n) Gen.bool
  genes <- zipWith swaps coins <$> Gen.list (Range.singleton n) genGene
  irs <- Gen.list (Range.singleton n) genIR
  return $ mkRGenome Circular Nothing Signed genes irs
  where
    swaps b v = if b then v else -v

-- | Adjacencies (unordered pair of unsigned genes) with the IR between them
adjacencyIRs :: GenesIRsR -> [((Gene, Gene), Int)]
adjacencyIRs g = List.sort [(sortPair (abs (getGene i g), abs (getGene j g)), irToInt (getIR i g)) | (i, j) <- adjacencies g]
  where
    sortPair (a, b) = (min a b, max a b)

rearrangeGenome :: (RigidIntergenicChromosome g) => Int -> g -> Gen g
rearrangeGenome = applyReversals

rearrangeAndFlexibilizeGenome :: Int -> GenesIRsR -> Gen GenesIRsF
rearrangeAndFlexibilizeGenome k g = do
  g' <- rearrangeGenome k g
  f <- Gen.int (Range.linear 0 100)
  return (flexibilize f g')

applyReversals :: (RigidIntergenicChromosome g) => Int -> g -> Gen g
applyReversals k g =
  if size g <= 4
    then return g
    else do
      revs <- replicateM k rev
      return $ foldr (\(i, j) g' -> intergenicReversal (mkIdx i) (mkIdx j) 0 0 g') g revs
  where
    rev = do
      i <- Gen.int (Range.linear (2 :: Int) (size g - 2))
      j <- Gen.int (Range.linear (i + 1) (size g - 1))
      return (i, j)

prop_getGeneIsGene :: Property
prop_getGeneIsGene = property $ do
  (GW g) <- forAll (genGenome 100)
  idx <- mkIdx <$> forAll (Gen.int (Range.linear 1 (size g)))
  assert $ getGene idx g `isGene` g

prop_invGeneIsGene :: Property
prop_invGeneIsGene = property $ do
  (GW g) <- forAll (genGenome 100)
  idx <- mkIdx <$> forAll (Gen.int (Range.linear 1 (size g)))
  assert $ invGene g (getGene idx g) `isGene` g

prop_invGeneIsAutoinverse :: Property
prop_invGeneIsAutoinverse = property $ do
  (GW g) <- forAll (genGenome 100)
  idx <- mkIdx <$> forAll (Gen.int (Range.linear 1 (size g)))
  let gene = getGene idx g
  invGene g (invGene g gene) === gene

prop_setGenesAreGenesInIndices :: Property
prop_setGenesAreGenesInIndices = property $ do
  (GW g) <- forAll (genGenome 100)
  genes <- forAll $ Gen.list (Range.linear 1 (size g - 2)) (mkGene <$> Gen.int (Range.linear 1 200))
  indices <- forAll $ take (length genes) <$> Gen.shuffle [2 .. mkIdx (size g - 1)]
  let g' = setGenes indices genes g
  map (`getGene` g') indices === genes

prop_makeSingletonsAddSingletonsInIndices :: Property
prop_makeSingletonsAddSingletonsInIndices = property $ do
  (GW g) <- forAll (genGenome 100)
  k <- forAll $ Gen.int (Range.linear 1 (min 5 (size g - 2)))
  indices <- forAll $ take k <$> Gen.shuffle [2 .. mkIdx (size g - 1)]
  let (g', genes) = makeSingletons g indices g
  let posMap = positionMap g'
      testSingleton gene =
        case geneMapLookup gene posMap of
          Nothing -> False
          Just pos -> length pos == 1
  map (`getGene` g') indices === genes
  assert $ all testSingleton genes

prop_alphabetAreGenes :: Property
prop_alphabetAreGenes = property $ do
  (GW g) <- forAll (genGenome 100)
  let alp = alphabet g
  assert $ all (`isGene` g) alp

prop_newGeneIsNotGeneAftersetGenes :: Property
prop_newGeneIsNotGeneAftersetGenes =
  property $ do
    (GW g) <- forAll (genGenome 100)
    genes <- forAll $ Gen.list (Range.linear 1 (size g - 2)) (mkGene <$> Gen.int (Range.linear 1 200))
    indices <- forAll $ take (length genes) <$> Gen.shuffle [2 .. mkIdx (size g - 1)]
    let g' = setGenes indices genes g
    assert $ not (getNewGene g' `isGene` g')

prop_newGeneIsNotGeneAfterMakeSingletons :: Property
prop_newGeneIsNotGeneAfterMakeSingletons = property $ do
  (GW g) <- forAll (genGenome 100)
  k <- forAll $ Gen.int (Range.linear 1 (min 5 (size g - 2)))
  indices <- forAll $ take k <$> Gen.shuffle [2 .. mkIdx (size g - 1)]
  let (g', _) = makeSingletons g indices g
  assert $ not (getNewGene g' `isGene` g')

prop_reversalsKeepBalanced :: Property
prop_reversalsKeepBalanced = property $ do
  g <- forAll (genRGenome 100)
  k <- forAll $ Gen.int (Range.linear 0 (size g))
  h <- forAll $ applyReversals k g
  assert $ areBalanced RRRM g h

prop_occurrenceWithGeneMapImplementationWorks :: Property
prop_occurrenceWithGeneMapImplementationWorks = property $ do
  g <- forAll (genRGenome 100)
  i <- forAll (Gen.int (Range.linear 1 (size g)))
  let a = canonicOri (getGene (mkIdx i) g)
  occurrence (positionMap g) a === (sum . fmap (\i' -> if abs (getGene (mkIdx i') g) == a then 1 else 0) $ [1 .. size g])

prop_joinIsInverseOfCut :: Property
prop_joinIsInverseOfCut = property $ do
  g <- toMC <$> forAll (genRGenome 100)
  chr_i <- mkIdx <$> forAll (Gen.int (Range.linear (1 :: Int) (numChromosomes g)))
  -- a linear chromosome can not be cut after its last gene
  i <- mkIdx <$> forAll (Gen.int (Range.linear (1 :: Int) (size (getChromosome chr_i g) - 1)))
  (join chr_i (chr_i + 1) False False False . cut (chr_i,i) 0 $ g) === g

prop_fullReversalKeepsAdjacencyIRs :: Property
prop_fullReversalKeepsAdjacencyIRs = property $ do
  circ <- forAll Gen.bool
  g <- forAll (if circ then genCircRGenome 100 else genRGenome 100)
  adjacencyIRs (intergenicFullReversal g) === adjacencyIRs g

prop_circularSubGenomeHasGenesBetween :: Property
prop_circularSubGenomeHasGenesBetween = property $ do
  g <- forAll (genCircRGenome 100)
  i <- mkIdx <$> forAll (Gen.int (Range.linear 1 (size g)))
  j <- mkIdx <$> forAll (Gen.int (Range.linear 1 (size g)))
  let sub = subGenome i j g
  map (`getGene` sub) (allGIdx sub) === map (`getGene` g) (betweenGIdx g i j)

prop_turnCircularIsInverseOfCut :: Property
prop_turnCircularIsInverseOfCut = property $ do
  chr <- forAll (genCircRGenome 100)
  let (bs_genes, bs_irs) = writeRGenome False chr
      g = readMCRGenome Signed (BS.pack "C " <> bs_genes) bs_irs
  i <- mkIdx <$> forAll (Gen.int (Range.linear (1 :: Int) (size g)))
  (turnCircularAtIndex 1 . cut (1, i) 0 $ g) === g

-- | Pair of multichromosome genomes (1 to 3 linear or circular chromosomes,
-- capped with 0 as in readMCGenome) over a small alphabet, so that there are
-- many repeated genes. H is G after k random DCJs, made flexible. If unbal is
-- True, some positions of G and of H receive exclusive genes before the DCJs.
genMultiPair :: Bool -> Sign -> Int -> Gen (ChromList GenesIRsR, ChromList GenesIRsF)
genMultiPair unbal = genMultiPairWith unbal False 50

-- | As genMultiPair. If shared is True, a replaced position may also receive a gene
-- of the alphabet, so that some labels occur in both genomes with distinct numbers
-- of occurrences. The flexibility of H is at most max_flex (0 gives the equality
-- relation).
genMultiPairWith :: Bool -> Bool -> Int -> Sign -> Int -> Gen (ChromList GenesIRsR, ChromList GenesIRsF)
genMultiPairWith unbal shared max_flex sign size_lim = do
  nchr <- Gen.int (Range.linear 1 3)
  alph <- Gen.int (Range.linear 1 6)
  chrs <- replicateM nchr $ do
    chtype <- Gen.element [Linear, Circular]
    n <- Gen.int (Range.linear 1 (max 1 (size_lim `div` nchr)))
    genes <- replicateM n (signGene =<< Gen.int (Range.linear 1 alph))
    irs <- replicateM (if chtype == Linear then n + 1 else n) genIR
    return (chtype, genes, irs)
  chrsG <- if unbal then mapM (exclusive alph 100) chrs else return chrs
  chrsH <- if unbal then mapM (exclusive alph 200) chrs else return chrs
  k <- Gen.int (Range.linear 0 (2 * nchr + 2))
  h <- applyDCJs k (mkMC chrsH)
  f <- Gen.int (Range.linear 0 max_flex)
  return (mkMC chrsG, flexibilizeMC f h)
  where
    mkMC = ChromList [0] . map (\(chtype, genes, irs) -> mkRGenome chtype (Just (0, 0)) sign genes irs)
    signGene a = case sign of
      Signed -> Gen.element [a, -a]
      Unsigned -> return a
    -- replace each gene by a new exclusive gene with probability 1/5 (or, if
    -- shared, by a new gene or a gene of the alphabet)
    exclusive alph base (chtype, genes, irs) = do
      genes' <- mapM (\(i, a) -> Gen.frequency [(4, return a), (1, signGene =<< newLabel alph (base + i))]) (zip [1 ..] genes)
      return (chtype, genes', irs)
    newLabel alph new = if shared then Gen.choice [return new, Gen.int (Range.linear 1 alph)] else return new

    -- same choice of positions as DB.applyDCJs
    applyDCJs :: Int -> ChromList GenesIRsR -> Gen (ChromList GenesIRsR)
    applyDCJs 0 g = return g
    applyDCJs nop g = do
      let chrs_idxs = filter (\idx -> let chr = getChromosome idx g in size chr > 1 || isCircular chr && size chr > 0) [1 .. mkIdx (numChromosomes g)]
      if null chrs_idxs
        then return g
        else do
          chri_ <- Gen.element chrs_idxs
          let chrs_idxs' = filter (\idx -> let chr = getChromosome idx g in idx /= chri_ || size chr > 2 || size chr == 2 && isCircular chr) chrs_idxs
          if null chrs_idxs'
            then return g
            else do
              chrj_ <- Gen.element chrs_idxs'
              let chri = min chri_ chrj_
                  chrj = max chri_ chrj_
                  chr_i = getChromosome chri g
                  chr_j = getChromosome chrj g
              i <- mkIdx <$> Gen.int (Range.linear 1 (size chr_i - (if isCircular chr_i then 0 else 1) - (if chri == chrj then 1 else 0)))
              j <- mkIdx <$> Gen.int (Range.linear (if chri == chrj then idxToInt i + 1 else 1) (size chr_j - (if isCircular chr_j then 0 else 1)))
              x <- Gen.int (Range.linear 0 (irToInt $ getIR i chr_i))
              y <- Gen.int (Range.linear 0 (irToInt $ getIR j chr_j))
              connectAC <- Gen.bool
              applyDCJs (nop - 1) (intergenicDCJ (chri, i) (chrj, j) connectAC x y g)

tests :: IO Bool
tests = checkSequential $$(discover)
