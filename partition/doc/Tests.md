# Tests (`test/`)

The tests are Hedgehog property tests. `test/Spec.hs` runs `GenomesCheck.tests`, `PartitionCheck.tests` and `PFptCheck.tests` and prints whether all of them passed (`stack test`, or `make test`).

> **Status:** the suite compiles and all properties pass.

## `test/GenomesCheck.hs`

### Generators

| Generator | Output |
|-----------|--------|
| `genRGenome n` / `genFGenome n` | Signed linear rigid / flexible chromosome with 3..n genes, caps included (`0`/`maxBound`). Flexible IRs are `[x, x + p% x]` |
| `genGenome n` | Either of the two, wrapped in `GenomeWrapper` |
| `genCircRGenome n` | Signed circular rigid chromosome (*n* genes, *n* IRs) |
| `rearrangeAndFlexibilizeGenome k g` | `g` after `k` random reversals, made flexible |
| `genMultiPair unbal sign n` | Pair of genomes with 1 to 3 linear or circular chromosomes over a small alphabet; H is G after random DCJs, made flexible. With `unbal`, some genes of G and of H are replaced by exclusive genes |
| `genMultiPairWith unbal shared max_flex sign n` | The same; with `shared`, a replaced gene may also get a label of the alphabet (labels in both genomes with distinct occurrences), and the flexibility is at most `max_flex` (0 gives the equality relation) |

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

- `prop_fptIsValidAndNotWorseThanOthers` (reverse, signed) and `prop_fptDirectIsValidAndNotWorseThanOthers` (direct, unsigned): on balanced genomes (`genMultiPair`), the FPT finishes, gives a balanced common partition of the input genomes, and its cost is never larger than the costs of Greedy, SOARComb and Approx.
- `prop_fptUnbalIsValidAndNotWorseThanOthers` and `prop_fptUnbalDirectIsValidAndNotWorseThanOthers`: the same on unbalanced genomes with labels in both genomes with distinct occurrences (`genMultiPairWith True True 50`), which use branching rule 2. The partition may be unbalanced.
- `prop_fptIsExact` and `prop_fptDirectIsExact`: with the equality relation (`max_flex` 0), on small balanced and unbalanced genomes, the cost of the FPT is the minimum cost of a common partition, found by enumerating all pairs of breakpoint sets (`bruteForceCost`).
- In these properties, an instance in which the FPT reaches its time limit (20 s) is discarded: a few small instances, for example with many caps, take longer.

## Suggested additions

- Port `heuristic/test/BalancerCheck.hs`: after balancing, the reduced genomes have the same gene content, and the final partition is common.
- Add a generator for multi-chromosome genomes with circular chromosomes, and run the partition properties on it.
