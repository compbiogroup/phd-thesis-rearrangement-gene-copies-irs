# Signed Minimum Common Flexible Intergenic String Partition

Algorithms for the Signed Minimum Common Flexible Intergenic String Partition (SMCFISP) problem and its unbalanced variant, used in the genome partition chapter of the thesis. Given an origin genome, whose intergenic regions have exact sizes, and a target genome, whose intergenic regions are intervals, the algorithms partition both genomes into blocks such that the blocks of the origin match the blocks of the target, minimizing the number of breakpoints.

## Algorithms

The executable `MCFISP-exe` runs one algorithm (option `-a`) on every instance of an input file:

| algorithm | description |
|---|---|
| `Approx` | approximation algorithm of the thesis (active classes), with factor 2tk for symmetric and transitive compatibility relations; unbalanced genomes first go through a balancing heuristic |
| `Greedy`, `GreedySin` | greedy algorithm (longest common substring first; `GreedySin` starts with the blocks of genes that occur once in each genome) |
| `SOAR`, `SOARComb` | adjacency matching algorithm (`SOARComb` builds the minimal partition from the matching) |
| `FPT` | exact fixed-parameter algorithm, with a limit of 10 minutes per instance |
| `ApproxLB` | writes only the lower bound used in the proof of the approximation |

The executable `DB-exe` generates random instances.

## Usage

```sh
make                    # builds MCFISP-exe and DB-exe in this folder (requires stack, with Nix)
make test               # runs the property tests
./MCFISP-exe -a SOARComb -s -i INPUT -o OUTPUT_PREFIX
```

The instances of the thesis are in `../db` (see `../db/README.md`), and `../scripts/run_part_tese.sh` runs the algorithms on them. The input and output formats and the other options are described in `doc/Main.md`.

## Documentation

- `doc/`: one page per source file, describing how the code is implemented (`doc/README.md` is the index).
- `CompatibilityWithText.md`: comparison between the code and the text of the thesis (in Portuguese).
- `BUGS.md`: bugs found and fixed (in Portuguese).
- `AddCircular.md`: notes on the review made to support circular chromosomes.
- `NameConvention.md`: naming conventions of the code.
