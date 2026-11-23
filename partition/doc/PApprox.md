# `src/PApprox.hs`

An approximation algorithm whose ratio depends on the maximum number of occurrences of a gene. It assumes **balanced** genomes; `MCFISP-exe` runs it through [Balancer](Balancer.md) so it never sees unbalanced input. It starts from the trivial partitions and adds breakpoints until the partitions are common.

## Vertices and the block match graph

```haskell
data Vertex i = Vtx { vSize :: Int, vInG :: Bool, vChr :: CIdx, vIdx :: i, vWhole :: Bool }
```

A vertex is a substring of a current block:
- it has `vSize` genes, starting at `vIdx` of chromosome `vChr`;
- `vInG` says whether it is in G or in H;
- it wraps around the origin on circular chromosomes;
- `vWhole` marks a whole circular chromosome without breakpoints, which also contains its closing adjacency.

The derived `Ord` compares `vSize` first, so ordered traversals visit small substrings first.

Position helpers:
- `shiftPos`, `offsetOf` and `posIn` move along a chromosome, using modular arithmetic on circular ones;
- `pieceIn` / `vSubGenome` extract the substring;
- `pairAt v o` returns the canonical gene pair at offsets `o`, `o + 1`;
- `adjOffsets v` lists the offsets of the adjacencies inside `v`, including the closing one of a whole vertex;
- `rightOverlap v u` returns the size of the overlap when `u` starts inside `v` and ends after it.

**`partVertices`** lists every substring of size ≥ 2 of every block. For a whole circular block it also adds the whole block itself, even with a single gene. Every rotation of every length is included.

**`mkBlockMatchGraph`** (`BG pg ph corres`) connects a vertex of G to a vertex of H when they match. Compatible substrings have the same genes up to orientation and rotation, so `isMatch` is only tested between vertices with the same key: size, whether they are whole circular chromosomes, and the sorted list of genes without orientation. The map `corres` stores, for each vertex, its component (filled in later) and its neighbours.

**`subgraphOf`**: the compatibility of two substrings does not depend on the partitions, so the graph of any partitions is the subgraph of the initial graph G_0 induced by the vertices of their blocks (all of them are vertices of G_0). `runApprox` builds G_0 once and uses this at the start of each iteration.

**`removeAdjacency`**: after a breakpoint is added, the vertices that contain its adjacency are removed. Only the components of these vertices change, and only their vertices are labelled again (`labelComponents`). This is how `addBreaks` updates the graph.

**`getConnectedComponents`**:
- labels the components with a DFS;
- for each component, records a representative vertex (`cVtx`) and the number of vertices of G and of H;
- `cPerfect` is true when the component has the same number of vertices on each side and a perfect matching between them, found with the same augmenting-path method as `findUncommon`.

## Main loop: `approxPartition matcher g h`

It follows Algorithm 1 of the thesis (`teo:aprox_part`), with a set of **active classes**.

```
G_0 := compatibility graph of the trivial partitions
classes, B_C := mkClasses G_0
(pg, ph) := trivial partitions;  active := {};  q := 0
repeat
   while some component is not perfect: addBreaks
   compatibilityCompletion
until the partition is common   (error if an iteration adds no breakpoint)
```

`runApprox` returns the partitions, q (the number of components whose sets B_C were activated) and the number of blocks D broken by the second stage. `approxLowerBound` returns q, and `approxPartitionPhase2` returns the partition with the number of blocks D broken (used by `MCFISP-exe` for the file `OFILE_phase2`). These components do not share adjacencies, so a common partition has a distinct breakpoint in each of them, and the cost is at most 2tk·q (Theorem `teo:aprox_part`, for symmetric and transitive relations).

### `mkClasses`

- The **class** of an adjacency is the component of G_0 that contains it, seen as a substring with two genes (`clsOf`).
- **Tmin(G, H)**: the non-perfect components that don't contain a smaller non-perfect component. A vertex contains a vertex of a non-perfect component if its prefix or suffix with one gene less is such a vertex or contains one (`containsNotPerf`, computed once for all vertices); a whole circular chromosome contains its substrings with all its genes. `clsTminComp` maps their vertices to the component.
- `interMap` records, for each pair of Tmin vertices from different components that overlap on the right (`rightOverlap`), the smallest overlap on each side and which components intersect. `updateEntry` flips the side when the vertex's orientation differs from the orientation of the component's representative. `sharedMap` has the classes of the shared adjacencies.
- The set B_C of each component (`clsB`) follows the cases of the lemma `tmin_comp`:
  - **i**, no intersections: the class of the first adjacency;
  - **ii**, `case2Cands` not empty (classes of C that are intersection edges with respect to every component that intersects C): one of them, preferring the adjacency common by position to all the intersections;
  - **iii**, otherwise: the classes of the first and the last adjacencies.

### `addBreaks`

- C is the component of Tmin(P) of the first candidate. The candidates are the vertices of non-perfect components with fewest adjacencies, sorted by preferring the genome with more vertices of the component, then higher degree. Counting adjacencies (and not genes) matters for whole circular chromosomes: one with n genes has n adjacencies and contains substrings with n genes.
- If a vertex of C has an adjacency of an active class, the breakpoint goes there.
- Otherwise, C is a component of Tmin(G, H) (Lemma `tmin_ativa`). The classes of B_C are activated, q increases, and the breakpoint goes in an adjacency of B_C of a vertex of C.
- If neither applies (it should not happen while every breakpoint is in an active class), the first adjacency of the first vertex is cut and its class is activated, so that the algorithm still terminates.

### `compatibilityCompletion`

Once every component is perfect, there may still be blocks without a partner, because the matchings of different components can conflict. Each call breaks at most one block and returns, so that the main loop makes Tmin empty again before the next one:

1. The unmatched blocks are the blocks without a pair in a maximum matching (`getBlocksFit`). B is the largest of them, in G or in H.
2. `findBreakpoint` looks for the first block `b2` of the other genome that contains a substring compatible with B (`isCompatibleWithSubGenome`); the blocks of the alternating path of `findUncommon` are used only if there is none. It also takes a block `b1` that matches `b2` as a whole.
   - For two whole circular blocks, `alignment` finds the rotation and the direction.
3. `blockCuts` breaks `b2` around the substring and breaks `b1` at the corresponding offsets. In reverse, offset `o` maps to `len − o − k`. The classes of the cuts are activated, so that every breakpoint stays in an active class. This stage only adds breakpoints with relations that are not symmetric and transitive, where the thesis does not guarantee a factor.

When no block can be broken, the partitions are returned unchanged. If an iteration adds no breakpoint and the partition is still not common, the loop stops with `logicError` (`noProgress`) instead of repeating forever. Through `MCFISP-exe`, `Approx` only receives balanced genomes (see [Balancer](Balancer.md)).
