# Tests (`test/`)

The tests are Hedgehog property tests. `test/Spec.hs` runs `GenomesCheck.tests` and `PartitionCheck.tests` and prints whether all of them passed (`stack test`, or `make test`).

> **Status:** the suite compiles and all 26 properties pass (14 in `GenomesCheck`, 12 in `PartitionCheck`).

## `test/GenomesCheck.hs`

### Generators

| Generator | Output |
|-----------|--------|
| `genRGenome n` / `genFGenome n` | Signed linear rigid / flexible chromosome with 3..n genes, caps included (`0`/`maxBound`). Flexible IRs are `[x, x + p% x]` |
| `genGenome n` | Either of the two, wrapped in `GenomeWrapper` |
| `genCircRGenome n` | Signed circular rigid chromosome (*n* genes, *n* IRs) |
| `rearrangeAndFlexibilizeGenome k g` | `g` after `k` random reversals, made flexible |

`adjacencyIRs g` lists, for each adjacency, the unordered pair of genes with its IR. It's used to check that operations keep adjacencies.

### Properties

- **Genes and gene maps**:
  - `getGene` and `invGene` give genes of the genome, and `invGene` is its own inverse;
  - `setGenes` writes the given genes;
  - `makeSingletons` creates singletons, and the fresh gene isn't in the genome afterwards;
  - `alphabet` contains only genes of the genome;
  - `occurrence` agrees with a direct count.
- **Operations**:
  - reversals keep the genomes balanced;
  - `join` is the inverse of `cut` (the cut position is drawn in `[1, size − 1]`, since a linear chromosome can't be cut after its last gene);
  - `turnCircularAtIndex` after `cut` gives back the original circular chromosome;
  - `intergenicFullReversal` keeps the IR of every adjacency, for linear and circular chromosomes.
- **Circular subgenomes**: `subGenome i j` on a circular chromosome contains exactly the genes of `betweenGIdx i j`, including blocks that wrap.

## `test/PartitionCheck.hs`

### Helpers

| Helper | Meaning |
|--------|---------|
| `toMCF h` | Wraps a flexible chromosome into a `ChromList`, taking its first and last genes as caps (the flexible version of `toMC`) |
| `genPair n` | A random genome G (`toMC (genRGenome n)`) and H, made from G by up to `size g` random reversals and then made flexible |
| `subChr q a b` | The block of `q` between positions `a` and `b` of the same chromosome, as a chromosome |

### Properties

- `mkPartitionFromBreakpoints` and `mkPartitionFromBlocks` build the same partition from a random set of breakpoints.
- `commonPrefix` (from a random position of G and an occurrence of the same gene in H) and `longestSubstring` return blocks that match.
- The SOAR reduction rules (`suboptimalRuleInterval`, `suboptimalRulePairs`) keep every original singleton and keep the genomes balanced.
- For identical genomes, every component of the block match graph is perfect.
- Greedy, GreedySin, SOAR and Approx produce partitions that, after `combine`, are balanced (`CGPbal`) and satisfy `checkCommonBal`. The property fails if `combine` returns `CGPunbal`, because the genomes are balanced.
- A disabled (commented) FPT property does the same.

## `test/PFptCheck.hs`

Placeholder with no properties. It isn't included in `Spec.hs`.

## Suggested additions

- Port `heuristic/test/BalancerCheck.hs`: after balancing, the reduced genomes have the same gene content, and the final partition is common.
- Add a generator for multi-chromosome genomes with circular chromosomes, and run the partition properties on it.
