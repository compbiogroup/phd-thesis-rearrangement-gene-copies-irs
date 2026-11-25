# `src/PFpt.hs`

A fixed-parameter (exact) algorithm based on *sample graphs* (Section `sec:fpt` of the thesis). Its running time is polynomial for a fixed partition size, maximum gene occurrence and total difference between the occurrences of the labels in the two genomes (Δ, the number of genes in exclusive blocks). According to the module header, it is meant for **direct** and **signed** partition problems. It accepts balanced and unbalanced genomes.

## Data structures

- **`Vertex i = Vtx { vInG, vIdx }`**: a gene position in G or in H.
- **`BlockWitness i = BW { w_left_pair, w_black_pair, w_right_pair, w_inverted }`**: a *witness* of a common block. The black pair `(posG, posH)` is a pair of matched genes. The left and right pairs are the leftmost and rightmost pairs (looking at G) that are parallel to the black pair, i.e. the extent of the block. `w_inverted` means H is read in the opposite direction.
- **`Edge i = Edge { eOther, eVpair, eParent, eBlack }`**: an edge between a vertex of G and a vertex of H. The black edge is the witness's own pair; the green edges are its parallel pairs. `eParent` is the witness that created the edge.
- **`SampleGraph`** (`SG`): the matcher, both genomes, their position maps, the *sample* (list of witnesses), the edges (`MultiMap Vertex Edge`) and the *removed* genes (`removed`, the set D of the thesis): genes chosen to be in exclusive blocks. A removed gene is not a vertex of the graph and has no edges.

## Flow: `fptPartition time_limit matcher g h`

1. **`sample0`**: one pair for each gene that occurs once in each genome, kept only if the two occurrences match directly or in reverse (`genesMatch`).
2. **`createSampleGraph`** adds a witness for each pair (`findLimits` + `addBlackEdge`). A pair that is parallel to, or conflicts with, a witness already in the sample is skipped, because it is already covered by that witness.
3. **`removeAbsent`** removes the genes whose label does not occur in the other genome. They are in exclusive blocks in every partition (branching rule 2 would have a single branch for each of them), and removing them from the start lets the pruning count their blocks.
4. **`initial`**: the partitions of Greedy (`greedyPartition False`) and SOARComb (`soarPartition True`) are computed, and the one with the smaller cost is the best partition found so far. Its cost is used by the pruning from the start, so the search prunes from the first node, it never returns a partition worse than the heuristics, and it always returns a partition, even when the time limit is reached before any leaf. The search is still exact, since it only discards branches that can not have a smaller cost.
5. **`runFpt`** branches (see below) and keeps the cost of the best partition and its sample graph (`Nothing` while the best one is the initial partition) in an `MVar`. It runs in a `race` against `threadDelay time_limit`. The returned `Bool` says whether the search finished before the time limit.
6. The best graph goes through `finalizeCorrespondence`, then `sampleGraphToPartition`. If no graph was better than the initial partition, that partition is returned.

### `findLimits`

Starting at the black pair, it extends to the right and to the left while the duos keep matching. It moves forward or backward in H depending on whether the first genes match directly or in reverse. The total extent is limited to `min sizeG sizeH − 1` steps, so a witness never covers more than a whole circle. Stepping uses `nextIdx`/`prevIdx`, so witnesses can cross the origin.

`findLimitsDir` computes the largest extension on each side (each one limited to that number of steps). If together they pass the limit, there are several maximal windows (rotations of a circular chromosome) and all of them are returned; otherwise the window is unique. `findAllLimits` returns the windows of both directions in which the genes match, and branching rule 1 tries all of them. `findLimits`, used for the initial sample, takes the first one, and `sample0` leaves out the singletons with more than one witness, so that they are handled by branching (BUGS #28 and #29).

### `addBlackEdge`

1. Removes the edges that are no longer valid at both new vertices (`cleanEdges`).
2. Inserts the black edge.
3. Adds green edges from the black pair towards both limits (`addEdges`), stopping at a vertex that already has a black edge or is removed.

### `cleanEdges v`

For each edge at `v` from a witness `bw'`, it removes the edges of `bw'` at `v` and at every vertex further from `bw'`'s black pair. The direction is decided with `isAhead`, which uses `fwdOffset`, the offset modulo the size on circular chromosomes. The walk stops before reaching the black pair again.

### `removeVertex v`

Puts `v` in the set of removed genes after `cleanEdges v`, which deletes the green edges that have `v` between them and their parent.

### Branching (`runFpt.branch`)

- `toBig` prunes a branch when a lower bound on the cost of its partitions is at least the cost of the best one found: `2|A| + exclusiveRuns − n_c`, where n_c is the number of chromosomes of both genomes. Each witness is at least one common block on each side, and `exclusiveRuns` (computed only for unbalanced genomes) counts the exclusive blocks that every partition of the branch has (see below).
- **Rule 1**: if there is a *rare odd path* (`getRareOddPath`), the algorithm branches on every candidate witness (`candidateWitnesses`) of the vertices on the larger side of the path (`largerSide`, the genome of its ends).
- **Rule 2**: otherwise, if there is an *abundant odd path* (`getAbundantOddPath`), for each vertex u on the larger side it branches into removing u (`removeVertex`), if fewer than |occ(a, G) − occ(a, H)| genes with the label a of u were removed (`canRemove`), and into adding each candidate witness of u.
- With no odd path, the graph is a candidate (`selectSG`) and is compared with the best one using `cost`. At this point each label a has exactly |occ(a, G) − occ(a, H)| removed genes, all in the genome where a is abundant.

Each branch adds a witness or a removed gene, so the search ends.

### `candidateWitnesses u`

The witnesses (`findAllLimits`, one per direction in which the genes match and per window) of the pairs of u and a gene v of the other genome with the same label, where v is not removed. A witness is skipped if it is parallel to or in conflict with the sample (`checkParallel`, `checkConflict`; the conflict test also leaves out the genes of black edges).

### `exclusiveRuns`

A lower bound on the number of exclusive blocks of the partitions of a branch. Each gene gets a status: removed, never removable (in a black edge, or its label is not abundant in its genome or already has all its genes in excess removed) or undecided. A label whose undecided genes are exactly the ones it still needs to remove has them counted as removed, since at the end of the branch it has exactly |occ(a, G) − occ(a, H)| removed genes. The bound is the number of maximal runs of removed genes of each chromosome, ignoring the undecided genes (two runs separated only by undecided genes may become one block); on a circular chromosome the first and the last runs are joined. The number of runs of removed genes is not monotone along a branch (removing the genes between two runs joins them), so they can not be counted directly.

### `getOddPaths`

Walks the components of the graph, where each vertex has degree ≤ 2, and returns the paths with an odd number of vertices, each one starting at one of its ends. The paths are found from the non-removed vertices of G, plus the isolated vertices of H (every other path has a vertex in G). A path is *rare* when the label of its first vertex (an end, in the genome with more vertices in the path) occurs at most as often in that genome as in the other one (`isRare`), and *abundant* otherwise. With balanced genomes every odd path is rare.

### `checkParallel bw bw'`

True when both witnesses are on the same pair of chromosomes with the same orientation, and the black pair of `bw'` is at the same offset from the black pair of `bw` in G and in H (negated if inverted, modulo the size on circular chromosomes), within the extent of `bw`.

### `finalizeCorrespondence`

At the end, every vertex that is not removed has a black edge or is on an even path or cycle of green edges. This function deletes alternating green edges so that every such vertex has exactly one edge, which gives a one-to-one correspondence between the gene positions in common blocks.

### `sampleGraphToPartition`

- The removed genes form the exclusive blocks, one per maximal run of removed genes: an adjacency with exactly one removed gene is a breakpoint, and one with two removed genes is not.
- Otherwise, an adjacency `(i, i')` is a breakpoint when `i` and `i'` come from different witnesses (`parentOf`).
- It is also a breakpoint when `i'` is the first gene of the witness (`startOf`). This is the closing adjacency of a circular chromosome covered by a single witness.
- The exception is `closesCircle`: the witness covers two whole circular chromosomes of the same size, and their closing duos match in the witness's direction. In that case the block stays a whole circular block.

The cost is computed with `costUnbalanced` (breakpoints of both partitions), which also holds for balanced partitions.

## Limitations

- The pruning sees the exclusive blocks only as the genes are removed, and rule 2 is applied only when there are no rare odd paths. On unbalanced genomes with labels of many occurrences (for example the caps of several linear chromosomes) the search can be much slower than on balanced genomes of the same size.
