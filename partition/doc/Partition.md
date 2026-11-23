# `src/Partition.hs`

This module represents a partition of a genome into blocks and the common partition of two genomes. It also checks whether a pair of partitions is common, computes costs, merges blocks (`combine`), and writes partitions.

## `Partition c`

```haskell
GP :: Chromosome c => ChromList c -> [(CIdx, c)] -> Breakpoints (ChromList c) -> Int -> Partition c
```

The fields are the underlying genome, the blocks (each tagged with its chromosome), the breakpoints, and the number of blocks. `Breakpoints g = Set (GIdx g)`; for a genome a breakpoint is `(CIdx, Idx)`, a cut after that gene.

Accessors: `underlineGenome`, `blocks`, `cindexedBlocks`, `breakpoints`, `blockDelims`, and the internal `size`.

### Block delimiters

```haskell
data BlockDelim c = BlockDelim { b_chr :: CIdx, b_beg, b_end :: GIdx c, b_whole :: Bool }
```

A block is described by its chromosome and its first and last genes. `b_whole` marks a circular chromosome without breakpoints, which is kept as a circular block. `delimToBlock` extracts the block: the chromosome itself if `b_whole`, otherwise `subGenome b_beg b_end`, which is linear.

**`bpsToBlockDelims g bps`** produces the blocks of each chromosome:
- **Linear**: breakpoints at the last gene are ignored. The blocks are `[1, k1]`, `[k1+1, k2]`, …, `[km+1, n]`.
- **Circular without breakpoints**: one whole block.
- **Circular with breakpoints `k1 < … < km`**: the blocks are `[next(km), k1]`, `[k1+1, k2]`, …, `[km−1+1, km]`.
  - The first block contains gene 1 and crosses the origin, unless `km = n`.
  - A breakpoint at *n* is a real cut between genes *n* and 1.

### Constructors

| Function | Input |
|----------|-------|
| `mkPartitionFromBreakpoints g bps` | Breakpoint set; blocks come from `bpsToBlockDelims` |
| `mkPartitionFromBlocks g bs` | Blocks in genome order. Each chromosome's blocks must start at gene 1, so they can't cross the origin |
| `mkPartitionFromDelims g ds` (internal) | Delimiters; every block end is a breakpoint, except the end of a linear chromosome |
| `trivialPartition g` | One block per chromosome, no breakpoints |

## `CommonPartition c1 c2`

- **`CGPbal pg ph`**: pg and ph have the same number of blocks, and there is a perfect matching between them.
- **`CGPunbal pg ph`**: some blocks are exclusive.

**`mkCommonPartition matcher pg ph`** tries `checkCommonBal`, then `checkCommonUnbal`, and fails with `logicError` if neither holds. **`mkCommonPartition2 matcher g bg h bh`** builds both partitions from breakpoints first. All algorithms finish through one of these two functions.

### Block correspondence

`getBlocksCorrespondence_ matcher pg ph` returns, for each block of pg (0-based), the list of blocks of ph that it matches. `getBlocksCorrespondence` does the same for a `CommonPartition`.

### Balanced check: `checkCommonBal` / `findUncommon`

`checkCommonBal` requires an equal number of blocks and `findUncommon` returning `Nothing`.

`findUncommon` builds a perfect matching with augmenting paths (Kuhn's algorithm):
- The blocks of G are processed from largest to smallest.
- `fixCorres` maps each block of H to the block of G it is assigned to.
- For each block of G, it tries every candidate. A candidate that is already taken is freed by recursively reassigning its current owner.
- If a block can't be assigned, it returns that block together with the blocks of G that hold its candidates. [PApprox](PApprox.md) uses this list to decide what to break.

### Unbalanced check: `getBlocksFit` / `checkCommonUnbal`

`getBlocksFit` computes a maximum matching, block of G → block of H (both 0-based), also with Kuhn's algorithm:
- `assigned` maps blocks of H to blocks of G.
- `augment seen b_g assigned` searches for an augmenting path. It returns the updated `seen` even when the search fails, so each block of H is visited at most once per search, which keeps it polynomial.

`checkCommonUnbal` holds when the unmatched (exclusive) blocks of G and of H share no gene. The cap genes (`caps` of either genome) are ignored, because every linear chromosome contains them.

### Costs

- `costBalanced` gives the number of blocks of G.
- `costUnbalanced matcher` gives the number of blocks of H plus the number of exclusive blocks of G (the size of G's partition minus the size of the maximum matching). For `CGPbal` it is the number of blocks of H.

## `combine`

A post-processing step (used by `SOARComb`) that merges neighbouring blocks while the partition stays common.

1. Each block is tagged with whether it contains a gene that is a singleton in its genome (`testSingletons`).
2. `combine_` runs three rounds with progressively looser rules:
   - `CombineSinSin`: both blocks of the pair contain a singleton.
   - `CombineSinRep`: at least one does.
   - `CombineRepRep`: no requirement.
3. `go1` walks the two block sequences, using `Seq` zippers, looking for a consecutive pair `(a, b)` in G and a consecutive pair `(c, d)` in H, on the same chromosome, such that `a+b` matches `c+d`. Each pair is replaced by its merge, and the round restarts.
4. When no consecutive pair merges, `wrapMerge` tries pairs that cross the origin of a circular chromosome: the last and first blocks of a circular chromosome with at least two blocks (`adjPairs`).
5. The result is rebuilt with `mkPartitionFromDelims` and `mkCommonPartition`.

## Output: `writePartition`

Returns four `ByteString` lines: genes of G, IRs of G, genes of H, IRs of H. Each chromosome is written as `L`/`C`, its blocks separated by `|`, and ` ;` at the end. In the IR line, the IRs inside each block alternate with the IRs at the breakpoints, also separated by `|`. A circular chromosome with *k* breakpoints has *k* blocks and *k* breakpoint IRs. A linear one with *k* breakpoints has *k + 1* blocks.
