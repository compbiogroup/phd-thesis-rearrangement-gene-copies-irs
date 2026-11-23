# `src/PSOAR.hs`

An algorithm in the SOAR framework, adapted to IRs. It considers every pair of matching *duos* (a duo is two adjacent genes and the IR between them) and builds a conflict graph whose vertices are these pairs. It then chooses an independent set: the chosen pairs are the adjacencies that are kept, and every other adjacency becomes a breakpoint.

## `soarPartition use_combine matcher g_ h_`

1. Optional reduction rules: `suboptimalRuleInterval`, then `suboptimalRulePairs`. They rewrite `g_`, `h_` into `g`, `h` by turning some repeated genes into singletons. **At the moment they're disabled in the source** (`(g, h) = (g_, h_)`; the commented line above it enables them).
2. `makePMGraph4 matcher g h` builds the pair-match graph.
3. `independentSet` takes the complement of a 2-approximate vertex cover. `getBpsFromIS` turns the independent set into breakpoints.
4. The partition is built on the **original** genomes `g_`, `h_`. The rules only replace genes at fixed positions, and they do it the same way in both genomes, so breakpoints computed on `g`, `h` remain valid for `g_`, `h_`.
5. With `use_combine` (algorithm `SOARComb`), `Partition.combine` merges blocks afterwards.

## Pair-match graph (`PMGraph4`)

- **Vertices**: `PairMatch = (duoG, duoH, idxG, idxH, vertexId)` for every duo of G and duo of H that match (`isMatch`). Duos come from `getDuo` at every index, including the duo `(n, 1)` of circular chromosomes. Vertices are numbered from 1.
- **Edges**: `isConflict` decides whether two pairs conflict. They do if:
  - they use the same duo of G with different duos of H;
  - the two G duos don't overlap but the two H duos do;
  - the G duos are consecutive (`followsG`, via `nextIdx`) but the H duos aren't consistent. Consistent means that for a direct match the second H duo follows the first, and for a reverse match it precedes it.

### `vertexCover` / `independentSet`

The vertex set is a `BitMask` (an `Integer`). The vertex cover is built greedily from the edges: when neither endpoint of an edge is covered, both are added. This gives a maximal matching, hence a 2-approximation. The independent set is its complement.

### `getBpsFromIS`

Every adjacency of a genome (`adjacencies`, so the end of a linear chromosome isn't one) becomes a breakpoint, unless its duo starts at the G (or H) index of a pair in the independent set.

## Reduction rules

### `suboptimalRuleInterval`

For each gene `a^l` that occurs once in each genome:
- it looks for an interval `a^l … a^r` in G, where `a^r` is the next such singleton;
- the interval must correspond to `b^l … b^r` in H in one of three ways:
  - `checkDirMatch`: every duo matches directly, walking forward in both genomes;
  - `checkRevMatch`: every duo matches in reverse, walking backward in H;
  - `checkMiddleMatch`: the ends match directly and the inside matches in reverse (`a^l a^{l+1} … a^{r−1} a^r = b^l −b^{r−1} … −b^{l+1} b^r`).

The genes strictly inside the interval are replaced by fresh singletons with `makeSingletons`, the same values in G and H (inverted and in reverse order for a reverse match). `stepsFrom` walks at most *size − 1* steps with `nextIdx`/`prevIdx`, so intervals can cross the origin of a circular chromosome.

### `suboptimalRulePairs`

This rule looks for a singleton `a^h` followed by a repeated gene `a^t` (`match_ah`), or a repeated gene followed by a singleton (`match_at`). If the duo also occurs, directly or in reverse, next to the singleton's unique occurrence in H, the repeated gene is replaced by the same fresh singleton in both genomes. Each check runs on the current (already rewritten) genomes, so earlier replacements are taken into account.
