# `app/DB.hs` (`DB-exe`)

Generates random instances in the input format of [Main](Main.md). It builds an ancestor genome and derives G and H from it with rearrangements and insertions/deletions (indels).

## Options

| Option | Meaning |
|--------|---------|
| `RepDB -l L -h H -r REP` | Genes `1..REP` each repeated between `L` and `H` times; the remaining genes are singletons |
| `RandDB -l ALPH` | Genes drawn uniformly from `1..ALPH` |
| `-k K` | Number of instances |
| `-n N` | Number of genes of the ancestor |
| `-r R` | Number of operations (`-1`: shuffle the genome instead) |
| `-p P` | Percentage of reversals among the operations (the rest are transpositions); default 100 |
| `-s` | Signed genes |
| `-d D` | Number of deletions, followed by the same number of insertions, applied to H's side |
| `-c` | Multiple chromosomes: every operation is a DCJ |
| `-z` | IRs equal to 0 |
| `-f LF`, `-F HF` | Make the IRs of H flexible by ±LF%. With both options, one file is written per value in `LF, LF+10, …, HF`, named `OFILE_<value>` |
| `-o OFILE` | Output file |

Without `-f`/`-F` the IRs of H are rigid, and `MCFISP-exe` can't read the file (it expects flexible IRs for H).

## `genPair`

1. Builds the ancestor `g` with `randomGenomeWithReplicas` or `randomGenome`. It is linear, with caps `0`/`maxBound`.
2. With `-c`, `g_final` is `g` after ⌊R/2⌋ random DCJs (`applyDCJs`); otherwise `g_final = g`. `g_final` is written as G.
3. `g_indel` is `g` after `applyIndels`:
   - D deletions of one gene each (`intergenicDeletion`), with the merged IR drawn at random;
   - then D insertions of new genes (`intergenicInsertion`), with two random IRs each.
4. H is derived from `g_indel` in one of three ways:
   - with `-c`: ⌈R/2⌉ DCJs;
   - with `-r -1`: `shuffleGenome`;
   - otherwise: R reversals/transpositions (`applyOperations`, choosing each operation at random until the counts are used up).
5. `genBS` writes G and H with `writeMultiC`, omitting the caps, and makes H flexible for each requested level.

## `applyDCJs`

1. Picks two chromosomes, possibly the same one, that are big enough.
2. Picks a cut position in each one. Linear chromosomes can't be cut after their last gene. On the same chromosome, `j > i`.
3. Picks the amounts `x`, `y` left on the left of each cut IR, and the reconnection type `connectAC`.
4. Applies `intergenicDCJ` (see [Genomes](Genomes.md)).

## Known issues

See `BUGS.md`:
- `RepDB` may produce invalid genomes (item 16).
