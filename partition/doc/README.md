# Code documentation

These pages describe how the code in `partition/` is implemented, one page per source file. The code solves the **Signed Minimum Common Flexible Intergenic String Partition** problem and its unbalanced variant. Given an origin genome G, whose intergenic regions (IRs) have exact sizes, and a target genome H, whose IRs are intervals, it finds a partition of both genomes into blocks such that each block of G matches a block of H. The goal is the fewest blocks.

Known defects are tracked separately in [`../BUGS.md`](../BUGS.md).

## Pages

| Page | Source | Role |
|------|--------|------|
| [LocalBase](LocalBase.md) | `src/LocalBase.hs` | Small generic helpers and error messages |
| [Genomes](Genomes.md) | `src/Genomes.hs` | Genome/chromosome representation, rearrangement operations, I/O, matchers |
| [Partition](Partition.md) | `src/Partition.hs` | Partitions, common partitions, matching between blocks, costs, output |
| [PGreedy](PGreedy.md) | `src/PGreedy.hs` | Greedy algorithm (longest common substring first) |
| [PSOAR](PSOAR.md) | `src/PSOAR.hs` | SOAR-framework algorithm (independent set on a pair-match graph) |
| [PApprox](PApprox.md) | `src/PApprox.hs` | Approximation algorithm based on the block match graph |
| [PFpt](PFpt.md) | `src/PFpt.hs` | Fixed-parameter (exact) algorithm based on sample graphs |
| [Balancer](Balancer.md) | `src/Balancer.hs` | Heuristic that removes exclusive blocks so a balanced algorithm can be used |
| [Main](Main.md) | `app/Main.hs` | `MCFISP-exe`: reads instances, runs an algorithm, writes partitions |
| [DB](DB.md) | `app/DB.hs` | `DB-exe`: generates random instances |
| [Tests](Tests.md) | `test/*.hs` | Hedgehog property tests |

## Module dependencies

```
LocalBase
   └── Genomes
          └── Partition
                 ├── PGreedy
                 ├── PSOAR
                 ├── PApprox
                 ├── PFpt
                 └── Balancer
app/Main.hs ── all algorithm modules + Balancer
app/DB.hs   ── Genomes
```

## Core concepts

- **Gene**: an `Int`. Its sign is its orientation. Gene `0` is used as the *cap* that marks both ends of every linear chromosome, and `maxBound` was used as the right cap in older code.
- **Chromosome** (`GenesIRs`): a vector of genes plus a vector of IRs.
  - A linear chromosome with *n* genes, caps included, has *n − 1* IRs.
  - A circular chromosome with *n* genes has *n* IRs. The last IR sits between gene *n* and gene 1.
- **Rigid / flexible**:
  - `GenesIRsR` holds IRs as `R x`, an exact size; G is always rigid.
  - `GenesIRsF` holds IRs as `F l u`, an interval; H is always flexible.
- **Genome** (`ChromList c`): a list of chromosomes plus the list of cap genes.
- **Index**:
  - `Idx` is 1-based and indexes genes inside a chromosome.
  - `CIdx` is 1-based and indexes chromosomes.
  - A position in a genome is the pair `(CIdx, Idx)`.
- **Breakpoint** `(c, i)`: a cut between gene *i* and the next gene of chromosome *c*. On a circular chromosome, the next gene after *n* is 1. The end of a linear chromosome is never a breakpoint.
- **Block**: a maximal segment between breakpoints. On a circular chromosome a block may cross the origin. A circular chromosome without breakpoints is a single "whole" block, which stays circular.
- **Matcher**: decides whether a block of G matches a block of H.
  - Direct match: same genes, and every rigid IR lies inside the corresponding interval.
  - Reverse match: the same comparison against the reversed block.
- **Common partition**: two partitions with a one-to-one matching between their blocks.
  - **Balanced** (`CGPbal`): every block of G is matched to a block of H.
  - **Unbalanced** (`CGPunbal`): some blocks are *exclusive* and stay unmatched.

## Typical data flow (`MCFISP-exe`)

1. Read four lines per instance: genes of G, IRs of G, genes of H, IRs of H (see [Main](Main.md)).
2. `readMCRGenome` / `readMCFGenome` build `ChromList GenesIRsR` / `ChromList GenesIRsF`, adding caps to linear chromosomes.
3. An algorithm computes the breakpoint sets of both genomes. `mkCommonPartition2` turns them into a `CommonPartition`, checking it is balanced or unbalanced.
4. `writePartition` writes the blocks and the IRs at the breakpoints. The block correspondences and the running time are appended.
