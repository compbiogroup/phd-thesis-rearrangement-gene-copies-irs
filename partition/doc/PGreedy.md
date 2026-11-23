# `src/PGreedy.hs`

The greedy algorithm for common string partition, adapted to IRs. It repeatedly picks the longest pair of matching substrings that don't overlap any block chosen so far, and turns that pair into a block.

## `greedyPartition withSingleton matcher g h`

- **State**: the positions already in blocks (`genesOnBlocksG/H`) and the breakpoints collected so far (`breaksG/H`).
- **`getAllLS`** loops:
  1. If `withSingleton` (algorithm `GreedySin`), first go through `all_singletons`, the genes that occur once in each genome. For each one, look only for blocks that contain it. A singleton with no such block is skipped.
  2. Then look for the longest match with no restriction (`longestSubstring Nothing …`).
  3. For each match `(g_beg, g_end)` / `(h_beg, h_end)`:
     - add the breakpoint `g_end`, and the breakpoint before `g_beg` (`prevChromIdx`, which wraps on circular chromosomes and is `Nothing` at the start of a linear one);
     - do the same in H;
     - mark the genes in between as used.
  4. Stop when no match exists.
- **`removeRedundantBreaks`** drops breakpoints at the last gene of a linear chromosome.
- The result is `mkCommonPartition2 matcher g bg h bh`.

## `longestSubstring maybe_singleton matcher usedG usedH g h`

For every unused position `ig` of G and every unused position `ih` of H, it calls `commonPrefix` and keeps the longest result (`maxWith`, so the last one on ties). It returns the four end indices and the updated sets of used positions. `betweenGIdx` handles blocks that wrap around the origin.

## `commonPrefix maybe_singleton matcher usedG usedH g h ig ih`

- **Direction**: if the single genes at `ig` and `ih` match directly, it extends forward in both genomes. If they match in reverse, it extends forward in G and backward in H. Otherwise there is no match.
- **`commonPrefix' rev k jg jh`** advances one gene at a time, as long as all of these hold:
  - the next position exists (`nextIdx`, or `prevIdx` in H when reversed; wraps on circular chromosomes);
  - neither next position is used;
  - the prefix length `k` is still smaller than both chromosome sizes, so a block never goes all the way around a circle;
  - the *duo* (two genes and the IR between them) matches, directly or in reverse.

  Checking each consecutive duo is enough, because a block matches exactly when all its genes and IRs match.
- **Singleton**: with `Just a`, the match is accepted only if `a` occurs in the block of G.
- **Result**: `((g_beg, g_end), (h_beg, h_end), length)`. In H the pair is ordered left to right, so for a reverse match it is `(endCPH, ih)`.

## Notes

- Because the length is limited to *n* genes, a circular chromosome matched entirely becomes a linear block with one breakpoint. The block count is the same as a whole circular block.
- The complexity is roughly O(n²) calls to `commonPrefix` per selected block, so O(n³) or more in total.
