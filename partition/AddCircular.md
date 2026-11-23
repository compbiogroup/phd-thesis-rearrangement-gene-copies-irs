I read every source file and checked each place that touches a chromosome against your convention: in a circular chromosome, gene 1 comes after gene *n*, and IR *n* sits between gene *n* and gene 1. Nothing was changed.

## The main pattern

Almost all position-walking code uses `incIdx`/`decIdx` and stops at `firstGIdx`/`lastGIdx`. That's right for linear chromosomes, but in a circular one it never sees the adjacency (*n*, 1), and a block can't cross the origin. Rather than patching each place, I'd add a few helpers to the `Chromosome` class and use them everywhere:

```haskell
nextIdx :: g -> GIdx g -> Maybe (GIdx g)    -- circular: n -> 1;  linear: n -> Nothing
prevIdx :: g -> GIdx g -> Maybe (GIdx g)    -- circular: 1 -> n;  linear: 1 -> Nothing
adjacencies :: g -> [(GIdx g, GIdx g)]      -- circular also includes (n, 1)
subGenomeWrap :: GIdx g -> GIdx g -> g -> g -- allows beg > end on circular
```

## `src/Genomes.hs`

**Probably wrong for circular chromosomes:**
- **`intergenicFullReversal` (489–492):** `Vec.reverse vi` is off by one rotation. Genes `g1..gn` become `gn..g1`, so the first IR must be `r(n-1)`, not `rn`. The correct value is `reverse (init irs) ++ [last irs]`. Everything that reverses a circular chromosome depends on this: `isReverseMatch`, `getOri` and DCJ inversions.
- **`getOri` circular (495–498):** this compares `canonicPosition g` with `invOri` of that result. But `canonicPosition` has already tagged it `Linear`, so the linear reversal is used. It should compare `canonicPosition g` with `canonicPosition (invOri g)`.
- **`canonicPosition` (505–510):** it returns a chromosome tagged `Linear` that has *n* IRs, which breaks the linear invariant of *n* − 1 IRs. Every caller then sees a "linear" chromosome.
- **`isDirectMatch` RFDM (919–926):**
  - `chtypeG == chtypeH` is always true, because both sides went through `canonicPosition`. So a circular block can match a linear one.
  - The two sides choose their minimum rotation separately, one from rigid IRs and one from flexible IRs. When genes repeat, they can pick different rotations, so equal chromosomes fail to match. For circular chromosomes, try every rotation of one side instead.
- **`mkRGenome` / `mkFGenome` (529–547):** the `ext` caps are added even when `chtype == Circular`.
- **`readMCGenome` (564):** it passes `Just (0, 0)` for every chromosome, so circular chromosomes get cap genes.
- **The four rigid operations (629–740), marked TODO:**
  - The asserts allow `i = 1` for circular chromosomes, but the code then reads `vi ! (i - 2)`, which is index −1. For a circular chromosome, the IR before gene 1 is at index *n* − 1.
  - No operation can wrap around the origin.
- **`subGenome` (381):** `Vec.slice` can't express a block that crosses the origin. It also always returns `Linear`, so a whole circular chromosome turned into a block loses its circularity. `trivialPartition` keeps it circular, so the two disagree.
- **`subGenomeBySize` / `getDuo` (425–429, 268):** there's no duo (*n*, 1) for circular chromosomes.
- **`isCompatibleWithSubGenome` (954–959):** it never tests windows that wrap around.
- **`ChromList` `betweenGIdx` (438–449):** returns `[]` when `beg > end` in a circular chromosome.

**Worth checking:**
- **`turnCircular` / `turnLinear` (403–412):** I traced both and they look right. `turnLinear` doesn't add caps, which is fine because `cut` adds them.
- **`toMC` (526):** it takes the caps from the first and last genes, which is only valid for linear chromosomes. DB and the tests only use it with linear genomes today.
- **`intergenicDCJ` (750–856):** there are 12 linear/circular cases and the `j'`/`chrj'` remapping after the first `cut`. The cases look plausible, but I'd test each one, for example with a property test that checks gene content and IR sums are preserved.
- **`cut` (858–869):** the circular branch looks right. It can't handle a linear chromosome with `i = n`, but the assert already prevents that.
- **`join` (870–883):** errors if either chromosome is circular, because `combine` rejects circular chromosomes. Every caller has to make sure both are linear.

## `src/Partition.hs`

- **`bpsToBlockDelims` (82–103):**
  - In a circular chromosome with *k* ≥ 1 breakpoints, the segment from the last breakpoint through *n* → 1 up to the first breakpoint is a single block. Right now it is split into two.
  - A breakpoint at (*c*, *n*) is a real cut between genes *n* and 1. Right now it is treated like a chromosome end.
- **`mkPartitionFromBreakpoints` (106):** the block count `|bps| + numChromosomes` is one too high for each circular chromosome that has breakpoints. Such a chromosome has *k* blocks, not *k* + 1.
- **`mkPartitionFromBlocks` (112–126):** assumes the blocks start at position 1 and run in order, so it can't hold a block that wraps.
- **`combine` / `go1` (333–372):** only merges blocks that are next to each other in the list. The last and first blocks of a circular chromosome are also neighbours, but they're never merged.
- **`testSingletons` (325):** `betweenGIdx` returns nothing for a block that wraps.
- **`writePartition` (280–294):** interleaving blocks with breakpoint IRs assumes #IRs = #blocks − 1. For a circular chromosome they are equal. The output also has no chromosome separators or `L`/`C` markers.

## `src/PGreedy.hs`

- **`removeRedundantBreaks` (38–43):** drops breaks at `lastGIdx`. For circular chromosomes, that is the real adjacency (*n*, 1) and must be kept.
- **`decChromIdx` (58):** should wrap from 1 to *n* for circular chromosomes instead of producing index 0.
- **`commonPrefix'` (122–134):** stops at `lastGIdx chr_g` / `firstGIdx chr_h`. To extend past the origin, it also needs a stop after `size` steps, or it will loop forever.
- **`longestSubstring` (78–79):** `betweenGIdx g g_beg g_end` has the same problem with blocks that wrap.

## `src/PSOAR.hs`

- **`getDuos` (100–103):** misses the duo (*n*, 1).
- **`isConflict` (68–82):** adjacency is checked with `incIdx i == i'` and `[j, incIdx j]`, so it misses (*n*, 1).
- **`getBpsFromIS` (116–117):** `allGIdx \\ duos` always includes (*c*, *n*) as a breakpoint.
  - Linear: this is the cap position, so it creates an empty block (*n* + 1 .. *n*). Check this one even without circular chromosomes.
  - Circular: it's a real breakpoint, and should be one unless the duo (*n*, 1) was chosen.
- **`suboptimalRuleInterval` (167–229):** the `checkDirMatch`, `checkRevMatch` and `checkMiddleMatch` ranges stop at the chromosome ends, so they don't wrap.
- **`suboptimalRulePairs` (258, 265–266, 273, 281–282):** the checks `i /= lastGIdx`, `i > firstGIdx`, `j < lastGIdx` and `j > firstGIdx` all leave out the wrap neighbour.

## `src/PApprox.hs`

- **`vEnd` (54–55) and every `vBeg`/`vEnd` interval test:** these are `selectTmin` (284–289), `updateInterMap` (297–309), `separateInterCases` (328–334) and `addBreaks_ii`/`addBreaks_iii`. They all assume intervals that don't wrap.
- **`mkBlockMatchGraph` / `corresList` / `getCorresIn` (230–252):** positions are computed as `prev_sizes + beg` from the start of the chromosome. That breaks once a block wraps around the origin.
- **`addBreaks.testIR` (144–153), `getGenePair` (364–371) and `findBreakpoint'` (421):** `i + 1` doesn't wrap.
- **`compatibilityCompletion` (412–444):**
  - Splitting a block into `[1..i-1] | b' | [rest]` assumes a linear block. For a whole circular chromosome from `trivialPartition`, the first and last pieces are adjacent, so you should get 2 blocks, not 3.
  - The reverse-index formula on line 437 also assumes a linear block.

## `src/PFpt.hs`

- **`findLimit` (352–373):** stops at the chromosome bounds. If you make it wrap, it also needs a stop once the witness covers the whole circle.
- **`checkParallel` (265–284):** `idxDist` and the `<=`/`>=` tests in `intersected` assume a linear order. On a circle, distances and the left/right order need modular arithmetic.
- **`addEdges` (395–407) and `cleanEdges.move`/`go` (424–440):** these step with `second incIdx`/`decIdx` and decide direction with `vIdx u < idxG`. Neither works on a circle. The bounds check uses the whole-genome `firstGIdx sg_g`/`lastGIdx sg_g`, not the chromosome's own bounds. It only stops correctly because looking up a missing vertex returns `[]`.
- **`sampleGraphToPartition` (180–188):** never checks the pair (*n*, 1).

## `app/` and `test/`

- **`app/Main.hs:120–121`:** everything is read as `Linear` with `readRGenome`. To accept circular input, use `readMCRGenome`/`readMCFGenome`.
- **`app/DB.hs` `applyDCJs` (279–294):** the index ranges look right for circular chromosomes. The other operations (reversal, transposition, indels) are only ever applied to linear genomes.
- **`test/GenomesCheck.hs`:** the generators only make linear genomes. A circular generator would let property tests catch most of the problems above, for example "reversing twice gives back the original" (which the circular `intergenicFullReversal` will currently fail) and `turnLinear`/`turnCircular` round-trips.

## Bugs unrelated to circular chromosomes

While reading I also found these. They affect linear multi-chromosome genomes too:
- **`Genomes.hs:453` and `:467`:** `(!! idxToInt chr_i)` uses the 1-based `CIdx` to index a 0-based list, so it picks the wrong chromosome. `getChromosome` subtracts 1 correctly.
- **`Genomes.hs:471`:** `drop chr_i . take (chr_j + 1)` includes `chr_j` and `chr_j + 1`. To get only the chromosomes strictly between them, it should be `take (chr_j - 1)`.
- **`Genomes.hs:426`:** `fpow n incIdx` gives `beg + n`, so `subGenomeBySize _ 2` returns 3 genes. This means SOAR's duos are actually triples. It should be `n - 1`.
- **`Genomes.hs:538` and `:550`:** `mkRGenome0`/`mkFGenome0` add 2 caps but create only `length genes - 1` IRs; they need `length genes + 1`. They aren't used anywhere at the moment.
- **`PApprox.hs:143`:** `chr_h = getChromosome (vChr v) g` should use `h`.
- **`PApprox.hs:279, 284–298`:** pairs of vertices are compared by index without checking `vChr`, so vertices on different chromosomes can look like they overlap.
- **`PSOAR.hs:280`:** `checkPair_at` requires `chr_i == chr_j`, but G and H chromosome numbers don't correspond, and `checkPair_ah` has no such check.
- **`PFpt.hs:180–188`:** `skipFirst` only skips (1, 1). For the first gene of every other chromosome, `decIdx` produces (*c*, 0), the edge lookup returns `[]`, and the code hits `error patternError`.
