# `src/PFpt.hs`

A fixed-parameter (exact) algorithm based on *sample graphs*. Its running time is polynomial for a fixed partition size and maximum gene occurrence. According to the module header, it is meant for **direct** and **signed** partition problems.

## Data structures

- **`Vertex i = Vtx { vInG, vIdx }`**: a gene position in G or in H.
- **`BlockWitness i = BW { w_left_pair, w_black_pair, w_right_pair, w_inverted }`**: a *witness* of a common block. The black pair `(posG, posH)` is a pair of matched genes. The left and right pairs are the leftmost and rightmost pairs (looking at G) that are parallel to the black pair, i.e. the extent of the block. `w_inverted` means H is read in the opposite direction.
- **`Edge i = Edge { eOther, eVpair, eParent, eBlack }`**: an edge between a vertex of G and a vertex of H. The black edge is the witness's own pair; the green edges are its parallel pairs. `eParent` is the witness that created the edge.
- **`SampleGraph`** (`SG`): the matcher, both genomes, their position maps, the *sample* (list of witnesses), and the edges (`MultiMap Vertex Edge`).

## Flow: `fptPartition time_limit matcher g h`

1. **`sample0`**: one pair for each gene that occurs once in each genome, kept only if the two occurrences match directly or in reverse (`genesMatch`).
2. **`createSampleGraph`** adds a witness for each pair (`findLimits` + `addBlackEdge`). A pair that is parallel to, or conflicts with, a witness already in the sample is skipped, because it is already covered by that witness.
3. **`runFpt`** branches (see below) and keeps the best sample graph in an `MVar`. It runs in a `race` against `threadDelay time_limit`. The returned `Bool` says whether the search finished before the time limit.
4. The best graph goes through `finalizeCorrespondence`, then `sampleGraphToPartition`.

### `findLimits`

Starting at the black pair, it extends to the right and to the left while the duos keep matching. It moves forward or backward in H depending on whether the first genes match directly or in reverse. The total extent is limited to `min sizeG sizeH − 1` steps, so a witness never covers more than a whole circle. Stepping uses `nextIdx`/`prevIdx`, so witnesses can cross the origin.

`findLimitsDir` computes the largest extension on each side (each one limited to that number of steps). If together they pass the limit, there are several maximal windows (rotations of a circular chromosome) and all of them are returned; otherwise the window is unique. `findAllLimits` returns the windows of both directions in which the genes match, and branching rule 1 tries all of them. `findLimits`, used for the initial sample, takes the first one, and `sample0` leaves out the singletons with more than one witness, so that they are handled by branching (BUGS #28 and #29).

### `addBlackEdge`

1. Removes the edges that are no longer valid at both new vertices (`cleanEdges`).
2. Inserts the black edge.
3. Adds green edges from the black pair towards both limits (`addEdges`), stopping at a vertex that already has a black edge.

### `cleanEdges v`

For each edge at `v` from a witness `bw'`, it removes the edges of `bw'` at `v` and at every vertex further from `bw'`'s black pair. The direction is decided with `isAhead`, which uses `fwdOffset`, the offset modulo the size on circular chromosomes. The walk stops before reaching the black pair again.

### Branching (`runFpt.branch1`)

- `toBig` prunes branches whose sample is already at least as large as the best one found.
- **Rule 1**: if there is a *rare odd path* (`getRareOddPath`), the algorithm branches on every new witness that `edgesFromRareOddPath` derives from the vertices on the larger side. Only pairs whose genes match (`genesMatch`) are candidates. A witness is skipped if it is parallel to or in conflict with the sample (`checkParallel`, `checkConflict`).
- With no odd path, the graph is a candidate (`selectSG`) and is compared with the best one using `cost`.
- **Rule 2** (abundant odd paths, `removeVertex`) is commented out.

### `getOddPaths`

Walks the components of the graph, where each vertex has degree ≤ 2, and returns the paths with an odd number of vertices. A path is *rare* when its starting vertex's gene occurs at most as often in its own genome as in the other one (`isRare`).

### `checkParallel bw bw'`

True when both witnesses are on the same pair of chromosomes with the same orientation, and the black pair of `bw'` is at the same offset from the black pair of `bw` in G and in H (negated if inverted, modulo the size on circular chromosomes), within the extent of `bw`.

### `finalizeCorrespondence`

At the end, every vertex has a black edge or is on an even path or cycle of green edges. This function deletes alternating green edges so that every vertex has exactly one edge, which gives a one-to-one correspondence between gene positions.

### `sampleGraphToPartition`

- An adjacency `(i, i')` is a breakpoint when `i` and `i'` come from different witnesses (`parentOf`).
- It is also a breakpoint when `i'` is the first gene of the witness (`startOf`). This is the closing adjacency of a circular chromosome covered by a single witness.
- The exception is `closesCircle`: the witness covers two whole circular chromosomes of the same size, and their closing duos match in the witness's direction. In that case the block stays a whole circular block.

## Limitations

- Unbalanced genomes trigger the `checkDegrees` assertion (a TODO in the code).
