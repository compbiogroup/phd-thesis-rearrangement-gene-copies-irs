# Naming Conventions

These are the naming conventions found in the Haskell sources (`src/`, `app/`, `test/`).

## Case style

| Kind | Convention | Examples |
|---|---|---|
| Types, constructors, classes | `PascalCase` | `GenesIRs`, `ChromList`, `BlockWitness`, `RigidIntergenicChromosome` |
| Short constructor tags | Upper-case abbreviations | `GP`, `BW`, `SG`, `BG`, `GM`, `GLR`/`GLF`, `RRDM`, `RFRM`, `CGPbal`/`CGPunbal` |
| Top-level functions | `camelCase` | `mkPartitionFromBlocks`, `getConnectedComponents`, `suboptimalRulePairs` |
| Local variables / `where` bindings | Mostly `snake_case` | `chr_i`, `ir_x`, `blocks_fit`, `posMapG`, `maybe_oddPath` |
| Record fields | Two styles: a short prefix + `_` (`b_chr`, `w_black_pair`, `sg_posG`, `db_num_pairs`), or a camelCase prefix letter (`vInG`, `vIdx`, `eOther`, `cPerfect`) | |

## Function name prefixes

- **`mk`** builds something: `mkGene`, `mkIdx`, `mkRGenome`, `mkCommonPartition`
- **`get`/`set`** read and write: `getGene`, `getIR`, `getChromosome`, `setGenes`, `setIdx`
- **`is`/`has`/`are`** return a Boolean: `isGene`, `isLinear`, `isMatch`, `hasSigns`, `areBalanced`
- **`check`/`test`** run a validity test: `checkCommonBal`, `checkDegrees`, `testSingleton`, `testMatch`
- **`read`/`write`** handle ByteString I/O: `readRGenome`, `writeFGenome`, `writeMultiC`
- **`inc`/`dec`** step an index: `incIdx`, `decKIdx`
- **`xToY`** converts: `geneToInt`, `irToInt`, `bpsToBlockDelims`, `sampleGraphToPartition`
- A trailing **`_`** on a function marks an internal helper or variant: `getBlocksCorrespondence_`, `combine_`
- A trailing **`'`** or **`''`** marks a local helper version: `selectBlock'`, `findBreakpoint''`

## Domain-specific names

- **`g` and `h`** are always the two genomes (origin G and target H). Related names get a `G`/`H` or `_g`/`_h` suffix: `posMapG`, `posMapH`, `bsG`, `bs_h`, `breaksG`, `chr_g`.
- **`pg`, `ph`** are the partitions of g and h. **`bg`, `bh`** are their breakpoint sets.
- **`c1`, `c2`** are the chromosome type variables of the two genomes. **`m`** is the matcher type and **`matcher`** is its value.
- **`mc`** is a multi-chromosome genome (a `ChromList`).
- **`chr`** is a chromosome. **`chri`/`chrj`** are chromosome indices (`CIdx`) and **`chr_i`/`chr_j`** are the chromosomes at those indices. This is used consistently in `Genomes.hs` and `DB.hs`, but other files mix the two (`chr_i` is sometimes an index in `PSOAR.hs` and `GenomesCheck.hs`).
- **`i`, `j`, `k`** are gene positions. **`x`, `y`, `z`** are nucleotide amounts cut inside intergenic regions, matched to `i`, `j`, `k`.
- **`idx`**, **`idxs`**, **`gidx`** are indices. Variants: `idx_c`/`idx_i` (chromosome / inner index), `idxG`/`idxH`, `idxG0`, `idxGL`/`idxGR` (left/right).
- **`ir`/`irs`** are intergenic regions. **`vs`/`vi`** are the vectors of genes and IRs inside `GenesIRs`.
- **`a`, `b`** (or `gene`) are single genes. **`al`, `ar`, `bl`, `br`** stand for aˡ, aʳ, bˡ, bʳ, as in the paper's notation.
- **`bs`** is a list of blocks, **`bps`** a set of breakpoints, **`bw`** a block witness, **`sg`** a sample graph, **`bmg`** a block match graph.
- **`comps`/`comp`/`comp_id`** are components, **`corres`/`cors`** correspondences, **`neigs`** neighbours, **`cands`** candidates.
- **`v`, `u`** are vertices and **`e`, `e1`, `e2`** are edges.
- **`chtype`** is a `ChromType` and **`sign`** a `Sign`.
- **`posMap`** maps genes to positions (`positionMap`). **`occ`** is an occurrence count (`occG`, `occH`).

## Suffix conventions

- **`'`, `''`, `'''`**: successive updated versions of a value (`g'`, `edges''`, `corres'''`, `fixCorres'`).
- **`_`** suffix: the raw input before processing (`genes_`, `g_`, `h_`, `ls_`, `chri_`), or a pattern binding meant to be renamed (`pg_`, `ph_`).
- **`_0`** or **`0`** suffix: the initial value (`pg_0`, `breakSet_0`, `sg0`, `sample0`, `interMap_0`).
- **`_aux`** or **`aux_`** prefix: an accumulator inside a fold or recursion (`aux_seen`, `aux_edges`, `sg_aux`, `tminSet_aux`).
- **`_rec`** suffix: the result of a recursive call (`seen_rec`, `fixCorres_rec`).
- **Plural `s`** for lists: `genes`, `irs`, `bs`, `us`, `b2s`, `cs1`. **`bss`** is a list of lists.
- **`maybe_`** prefix for `Maybe` values: `maybe_sg`, `maybe_part`, `maybe_singleton`, `maybe_breakSet`.
- **`l_`** prefix for something converted to a list: `l_genes`, `l_irs`.
- **`bs_`** prefix for ByteStrings in I/O code: `bs_genes`, `bs_irs`. This clashes with `bs` meaning "blocks" in the partition code.
- **`new_`** for freshly made items (`new_genes`, `new_breaks`). **`old_`** for prior values (`old_b2_chr`).
- **`beg`/`end`** for interval bounds (`g_beg`, `h_end`, `idx_beg`). **`l`/`u`** (or `irl`/`irr`) for the lower and upper bounds of a flexible IR.

## Other patterns

- **Error strings** are named `xxxError` (`patternError`, `logicError`, `indexError`, `inputError`) and passed to `error`.
- **Property tests** are named `prop_camelCaseDescription`, as Hedgehog's `discover` expects.
- **Generators** in tests are `genX` (`genGene`, `genRGenome`, `genFGenomeWithSign`).
- **Rigid vs. flexible** is marked by `R` or `F` in names: `GenesIRsR`/`GenesIRsF`, `mkRGenome`/`mkFGenome`, `readRGenome`.
- **Multi-chromosome variants** use an `MC` suffix or infix: `readMCRGenome`, `flexibilizeMC`, `toMC`.

## Inconsistencies

- `snake_case` and `camelCase` are mixed in local names, sometimes in one function: `fixCorres_rec'`, `genesOnBlocksG`, `final_breaksG`, `initalVtxInG`.
- The `G`/`H` suffix is sometimes glued on (`posMapG`) and sometimes joined with an underscore (`bs_g`, `chr_g`).
- Record field prefixes come in two styles: `b_`/`w_`/`sg_`/`db_` versus `v`/`e`/`c`.
- `bs` means "blocks" in the partition code but "ByteString" in the I/O code.
- There are typos in identifiers: `geneMapUniom`, `separetInterCases`, `compatibilyCompletion`, `sucess`, `crhs_idxs`, `initalVtxInG`.
