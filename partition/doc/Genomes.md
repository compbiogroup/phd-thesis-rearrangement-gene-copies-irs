# `src/Genomes.hs`

This module covers how genomes are represented, how they are indexed and changed, the rearrangement operations used to generate instances, reading and writing, and the *matchers* that decide when two blocks are equivalent.

## Basic types

| Type | Definition | Notes |
|------|-----------|-------|
| `Gene` | `newtype Gene = Gene Int` | The sign is the orientation. Built with `mkGene`, read with `geneToInt` |
| `IR` | `R Int \| F Int Int` | Rigid size or flexible interval `[l, u]`. `irToInt` works only for `R`. `sumIR` adds two rigid IRs |
| `Idx` | `newtype Idx = Idx Int` | 1-based gene index inside a chromosome |
| `CIdx` | `newtype CIdx = CIdx Int` | 1-based chromosome index inside a genome |
| `Sign` | `Signed \| Unsigned` | With `Unsigned`, inverting a gene keeps it the same (`invGene`) |
| `ChromType` | `Circular \| Linear` | |

**`class GenIdx`** abstracts indices: `mkIdx`, `idxToInt`, `incIdx`, `decIdx`, `incKIdx`, `decKIdx`, `idxDist`. `Idx` and `CIdx` are instances. These functions don't wrap around; for chromosome-aware stepping use `nextIdx`/`prevIdx` (below).

## Chromosome representation: `GenesIRs`

```haskell
data GenesIRs = GenesIRs ChromType Sign (Vector Gene) (Vector IR) Gene
```

- **Fields**: chromosome type, sign, genes, IRs, and a *fresh gene*. The fresh gene is larger than every gene ever stored, and `getNewGene` uses it to create singletons.
- **Linear**: *n* genes and *n − 1* IRs. IR *k* (vector position *k − 1*) sits between genes *k* and *k + 1*. Linear chromosomes read from files start and end with cap genes.
- **Circular**: *n* genes and *n* IRs. IR *n* sits between gene *n* and gene 1.

The `Show` instance prints `(g1) - ir1 - (g2) …`, with `maxBound` shown as `inf`.

Two newtypes fix whether the IRs are rigid or flexible: `GenesIRsR` (`GLR`) and `GenesIRsF` (`GLF`). Both derive every class through `GenesIRs`.

`Eq GenesIRsR` compares sign, genes and IRs after `canonicPosition`, so two circular chromosomes are equal if one is a rotation of the other. A linear chromosome is never equal to a circular one.

### `canonicPosition`

For a circular chromosome it returns the rotation with the lexicographically smallest `(genes, IRs)`, rotating genes and IRs together. Linear chromosomes are returned unchanged.

## Genome representation: `ChromList c`

```haskell
data ChromList c = ChromList {caps :: [Gene], chromosomes :: [c]}
```

A list of chromosomes plus the cap genes. `Show` separates chromosomes with `|`. `Eq` compares only the chromosomes. The constructor is exported, and [Balancer](Balancer.md) uses it to build reduced genomes.

`toMC` wraps a single linear `GenesIRsR` into a `ChromList`, taking its first and last genes as caps. It fails on circular chromosomes.

## Type classes

### `Genome g`

The interface shared by single chromosomes (`GenesIRs`, index `Idx`) and whole genomes (`ChromList c`, index `(CIdx, GIdx c)`).

| Method | Meaning |
|--------|---------|
| `allGIdx`, `firstGIdx`, `lastGIdx` | All indices in order, first, last |
| `betweenGIdx g i j` | Indices from `i` to `j`. On a circular chromosome, `i > j` wraps around the origin. For `ChromList`, it can span chromosomes |
| `isGene a g` | `a` occurs in `g` (either sign when signed) |
| `size` | Number of genes |
| `getGene`, `setGenes` | Read / overwrite genes at indices (`setGenes` updates the fresh gene) |
| `invGene g a` | `-a` if signed, `a` otherwise |
| `subGenome i j g` | The block from `i` to `j` as a **linear** chromosome (see below) |
| `alphabet` | Set of absolute gene values |
| `getNewGene` | A gene value absent from the genome |
| `makeSingletons h idxs g` | Replaces the genes at `idxs` with new genes that don't occur in `g` or `h`; returns the new genome and the genes used |
| `hasSigns` | Whether the genome is signed |

`subGenome` on `GenesIRs`:
- If `i <= j`, it slices genes `i..j` and the `j − i` IRs between them.
- If `i > j` on a circular chromosome, it concatenates `i..n` and `1..j`, together with their IRs (including IR *n*). The block crosses the origin.
- The result is always `Linear`, so a block of *k* genes has *k − 1* IRs.

`subGenome` on `ChromList` returns a `ChromList`. If both ends are in the same chromosome it holds one block. Otherwise it holds the suffix of the first chromosome, the whole chromosomes in between, and the prefix of the last one.

### `Chromosome g` (extends `Genome`)

| Method | Meaning |
|--------|---------|
| `isLinear`, `isCircular` | Chromosome type |
| `turnCircular` | Linear → circular: drops both caps and adds the first and last IRs into the closing IR |
| `turnLinear i` | Circular → linear by cutting after gene `i` (rotates genes and IRs, drops the cut IR; doesn't add caps) |
| `combine g1 ir g2` | Concatenates two linear chromosomes with `ir` between them (fails on circular ones or different signs) |
| `subGenomeBySize i k` | The `k` genes starting at `i`, wrapping on circular chromosomes; `Nothing` if they don't fit |
| `getDuo i` | `subGenomeBySize i 2` |
| `nextIdx`, `prevIdx` | Next / previous index; wraps on circular, `Nothing` past the end of a linear chromosome |
| `adjacencies` | All `(i, nextIdx i)`; includes `(n, 1)` on circular chromosomes |

The partition algorithms walk chromosomes only through `nextIdx`/`prevIdx`/`adjacencies`, which is how they support circular chromosomes.

### `MultiChromosome g`

Implemented by `ChromList`:
- `getChromosome c` takes a 1-based `c`.
- `numChromosomes`, `getChromosomes`.
- `getIdxOfChrom` / `getIdxInChrom` split a genome position into its two parts, and `setIdx` builds a position from them.

### `IntergenicChromosome g`

- `getIR i` returns the IR after gene `i`.
- `intergenicFullReversal` reverses the whole chromosome, inverting signs if signed.
  - Linear: the IR vector is reversed.
  - Circular: the first *n − 1* IRs are reversed and the closing IR stays last, which keeps every adjacency's IR.

`Orientable GenesIRs` uses `intergenicFullReversal` as `invOri`. `getOri` compares a chromosome with its reversal: lexicographically for linear ones, and after `canonicPosition` for circular ones.

### `RigidIntergenicChromosome` (rigid single chromosome, used by `DB-exe`)

Each operation splits the IRs at its cut points: `x`, `y`, `z` nucleotides stay on the left, and the rest goes to the right.

| Operation | Effect |
|-----------|--------|
| `intergenicReversal i j x y` | Reverses genes `i..j`. The IR before `i` becomes `x + y` and the IR after `j` becomes `(ir_x − x) + (ir_y − y)` |
| `intergenicTransposition i j k x y z` | Swaps the segments `[i, j−1]` and `[j, k−1]`, recombining the three cut IRs |
| `intergenicDeletion i j x` | Deletes genes `i..j−1`; the two cut IRs become one IR of size `x` |
| `intergenicInsertion i g'` | Inserts the open genome `g'` (*m* genes, *m + 1* IRs) after gene `i`, replacing IR `i` |

On circular chromosomes the asserts allow `i = 1`, where the IR before gene 1 is the closing IR. They forbid the cases that would cut the same IR twice or delete the whole chromosome:
- a reversal of genes `1..n`;
- a transposition with `i = 1` and `k = n + 1`;
- a deletion with `i = 1` and `j = n + 1`.

### `RigidIntergenicMultiChromGenome` (rigid multi-chromosome, used by `DB-exe`)

- `cut (c, i) x`: cuts chromosome `c` after gene `i`, leaving `x` nucleotides on the left.
  - A linear chromosome becomes two linear chromosomes, with new caps at the cut.
  - A circular one becomes a single linear chromosome with caps on both ends.
- `join ci cj moveFirst invi invj`: concatenates two linear chromosomes, optionally reversing either. It removes the inner caps and adds their IRs together. The result goes to the position of the first chromosome, or of the second if `moveFirst`.
- `turnCircularAtIndex c`: `turnCircular` applied to chromosome `c`.
- `intergenicDCJ (ci, i) (cj, j) connectAC x y`: a double-cut-and-join.
  1. Cuts after `i` and after `j`. If both cuts are on the same chromosome, `j` is re-indexed into the second piece after the first cut.
  2. Rejoins the four ends either as (a, c)(b, d) or as (a, d)(b, c). There are 12 cases, depending on whether the chromosomes are linear or circular and whether the cuts are on the same chromosome. Each case is drawn in an ASCII diagram in the source.

## Gene maps

`GeneMap a` is an `IntMap` keyed by the **absolute** gene value.

| Function | Meaning |
|----------|---------|
| `positionMap g` | Gene → list of positions |
| `occurrence m a` | Number of positions of `a` in `m` |
| `geneMapLookup`, `geneMapAdjust` | Lookup / update a gene's entry |
| `geneMapUnion`, `geneMapDifference`, `geneMapIntersection` | Set operations on the keys |
| `geneMapDelete` | Removes a gene's entry |
| `occurrenceMax` | Largest number of occurrences of a gene |
| `singletonOnBoth _ pmG pmH a` | `a` occurs exactly once in each genome |

## Construction and I/O

- **`mkRGenome chtype ext sign genes irs`** and **`mkFGenome …`** (IR pairs) build a chromosome.
  - If `ext = Just (b, e)` and the chromosome is linear, caps `b` and `e` are added at the ends.
  - The fresh gene is `max |gene| + 1` over all genes, caps included, except `maxBound`. That value may be used as a cap, and adding 1 to it would overflow. It is 1 for an empty chromosome.
- **`mkRGenome0` / `mkFGenome0`**: chromosomes with zero IRs and caps `0` / `maxBound`.
- **`readRGenome` / `readFGenome`** parse one chromosome. Values are separated by spaces or commas, and flexible IRs are written `l:u`.
- **`readMCRGenome` / `readMCFGenome sign genesLine irsLine`** parse a genome:
  - chromosomes are separated by `;`;
  - each one may start with `L` (linear, the default) or `C` (circular);
  - every linear chromosome gets caps `0 … 0`, and the genome's cap list is `[0]`.
- **`writeRGenome rext` / `writeFGenome rext`** write the genes and IRs of one chromosome. With `rext`, the caps of linear chromosomes are omitted.
- **`writeMultiC writer`** writes a whole genome in the format read by `readMC*`, keeping the order of the chromosomes.
- **`writeIR`**: `x` or `l:u`.
- **`flexibilize p g` / `flexibilizeMC`** turn each rigid IR `x` into `F (x − p% x) (x + p% x)`.
- **`equalizeLinear sign g h`** adds empty linear chromosomes (the two caps and one IR 0, or `0:0` in H) to the genome with fewer linear chromosomes, so that both have the same number. Used by `MCFISP-exe` on every instance.

## Matchers

```haskell
class Matcher m g1 g2 where
  isMatch, isDirectMatch, isReverseMatch, areBalanced :: m g1 g2 -> g1 -> g2 -> Bool
```

| Matcher | Types | `isMatch` |
|---------|-------|-----------|
| `RRDM` | rigid × rigid | `==` (direct only) |
| `RRRM` | rigid × rigid | direct, or direct against the reversed first argument |
| `RFDM` | rigid × flexible | same chromosome type and gene list, every rigid IR inside its interval. For circular chromosomes, every rotation of H is tried |
| `RFRM` | rigid × flexible | `RFDM` direct, or direct against the reversed first argument. **This is the matcher used by `MCFISP-exe`** |
| `FlipMatcher m` | swaps the arguments of `m` | used to run an algorithm with H as the first genome |

`areBalanced` compares the multisets of canonical genes. For the IRs it checks that the rigid sum equals the other rigid sum, or lies between the sums of the interval bounds.

`isCompatibleWithSubGenome matcher b c` returns the first index of `c` where a substring of `c` (as long as `b`, wrapping on circular chromosomes) matches `b`.

## Random genomes (used by `DB-exe`)

- **`randomGenome zeros n lim sign`**: *n* genes drawn from `1..lim`, random signs if signed, and caps `0`/`maxBound`. The IRs are 0 if `zeros`, otherwise drawn from 1..100.
- **`randomGenomeWithReplicas zeros n rep low high sign`**: genes `1..rep` repeated `low..high` times each, then singletons up to *n* (see `BUGS.md` item 16).
- **`shuffleGenome`**: shuffles the genes (and flips signs at random). It also redistributes the IR total at random over the IRs.
