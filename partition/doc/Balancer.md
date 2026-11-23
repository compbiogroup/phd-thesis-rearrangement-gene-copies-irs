# `src/Balancer.hs`

A heuristic that turns two unbalanced genomes into two balanced ones, so that an algorithm for balanced genomes (currently [PApprox](PApprox.md)) can be used. It was ported from the older project in `heuristic/src/Balancer.hs` (`reduceGenesForDeletion`).

## Idea

When a gene occurs more times in one genome than in the other, the extra copies must be in exclusive blocks. The heuristic removes those copies and prefers long runs of consecutive removable genes, so the exclusive blocks are few and large. The remaining genomes have the same gene content. They are partitioned with the balanced algorithm, and the result is mapped back to the original genomes. Each removed segment becomes an exclusive block.

## Functions

### `occBalance g h`

For each absolute gene value, the occurrences in `g` minus the occurrences in `h`, as an `IntMap Int`.

### `markExclusive sign q occ`

Selects the positions of `q` to remove. `sign` is `1` for G (genes with balance > 0 are in excess) and `-1` for H (balance < 0).

1. For each chromosome, `chromRuns` scans the positions in order and builds runs of consecutive positions whose gene is in excess and not selected yet.
   - Inside a run the excess is decreased as genes are added, so a run never takes more copies of a gene than its excess. When a run ends, the counter is reset.
   - On a circular chromosome the scan starts right after a gene that isn't in excess, so a run can cross the origin.
   - Runs never cross chromosome boundaries.
2. The largest run is selected; `sortWith` is stable, so the leftmost one wins ties. Its positions are added to the selected set and the excess of their genes is decreased.
3. Steps 1 and 2 repeat until no gene is in excess.

Each iteration is O(n log n), and there are at most *n* iterations.

`balancedPartition` runs it first on G with the initial balance, then on H with the balance left after G. At that point every gene has balance ≤ 0, so only the copies H has in excess are removed.

### `segments p chr` (internal)

Splits the indices of `chr` into maximal segments of consecutive indices with the same value of `p`. On circular chromosomes the first and last segments are merged when they have the same value.

### `reduce marks q` (internal)

Builds the reduced genome:
- A chromosome without selected positions is kept as is, circular ones included.
- A chromosome with selected positions is split into its unselected segments, each one a linear chromosome (`subGenome`, which handles segments that cross the origin). The IRs at the cuts are dropped, because those adjacencies are breakpoints in the final partition anyway.
- The caps of `q` are kept.

It also returns:
- `posMap`: position in the reduced genome → position in `q`;
- `cuts`: the breakpoints of `q` around each removed segment, namely the gene before the segment (if any) and the segment's last gene. If a chromosome is removed entirely, no cut is added, so it becomes a single exclusive block.

### `balanceGenomes g h`

Runs `markExclusive` on G and then on H, and `reduce` on both. It returns `Nothing` when nothing is in excess. Otherwise it returns a `Balanced` record with the reduced genomes (`bal_g`, `bal_h`), the position maps (`bal_posMapG/H`), and the cuts (`bal_cutsG/H`). Both functions below use it.

### `balancedPartition partAlg matcher g h`

1. If nothing is in excess, it returns `partAlg matcher g h`.
2. Otherwise it reduces both genomes and runs `partAlg` on the reduced pair. Both reduced genomes have the same genes, so they are either both empty or both non-empty.
3. It maps the breakpoints of the result back with `posMap` and adds the `cuts`.
4. It builds the final partition with `mkCommonPartition2`. That function checks that the partition is common; normally it is `CGPunbal`, with the removed segments as the exclusive blocks.

### `balancedPartitionWith def partAlg matcher g h`

The same as `balancedPartition`, for an algorithm that also returns a value (for example, `approxPartitionPhase2`, with the number of uses of the second stage of Approx). The value is the one of the algorithm on the balanced genomes, or `def` when the balanced genomes are empty and the algorithm is not run.

### `balancedLowerBound lbAlg matcher g h`

Applies a lower bound for balanced genomes, such as `approxLowerBound`, to genomes that may be unbalanced. It computes the bound on the reduced genomes (0 if they are empty) and adds the number of cuts that are real breakpoints. The end of a linear chromosome doesn't count.

It bounds the number of breakpoints, in both genomes, of partitions whose exclusive blocks are the segments chosen by `markExclusive`. That isn't a lower bound for every partition of the unbalanced genomes, because a different choice of exclusive blocks could need fewer breakpoints.

## Differences from the original heuristic

- The original replaced each removed interval with a new marker gene. Here the intervals are removed and the chromosome is split, so the balanced algorithm never sees genes that exist in only one genome.
- Multiple and circular chromosomes are supported: runs stay inside one chromosome and can cross the origin.
- The cap genes (`0`) are ordinary genes here. `MCFISP-exe` gives both genomes the same number of linear chromosomes before running any algorithm (`equalizeLinear`, see [Main](Main.md)), so the caps are always balanced there. If the Balancer is called directly with different numbers of linear chromosomes, the extra caps are selected for removal.

## Usage

`app/Main.hs` runs `-a Approx` as `balancedPartition approxPartition RFRM g h` and `-a ApproxLB` as `balancedLowerBound approxLowerBound RFRM g h`. Any function with type `m c1 c2 -> ChromList c1 -> ChromList c2 -> CommonPartition c1 c2` can be passed as `partAlg`.
