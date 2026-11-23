# `app/Main.hs` (`MCFISP-exe`)

The command-line program that runs a partition algorithm on every instance of an input file.

## Options

| Option | Meaning |
|--------|---------|
| `-a, --algorithm ALGO` | `SOAR`, `SOARComb`, `Greedy`, `GreedySin`, `Approx`, `FPT`, or `ApproxLB` |
| `-i, --input IFILE` | Input file |
| `-o, --outfile OFILE` | Output prefix; instance *k* is written to `OFILE_000k` |
| `--no-par` | Process the instances one after another instead of in parallel (`parallel-io`) |
| `-s, --signed` | The genes are signed |

## Input format

Blank lines and lines starting with `#` (possibly after spaces) are ignored. Every four remaining lines form one instance:

```
L 1 2 3 4 ; C 5 6          # genes of G   (chromosomes separated by ';', L/C prefix)
1 2 3 4 5 ; 7 8            # IRs of G     (rigid)
L 1 -3 -2 4 ; C 6 5        # genes of H
1:1 2:3 3:3 4:4 5:6 ; 7:9 8:8   # IRs of H (flexible l:u)
```

Linear chromosomes get caps `0` at both ends when they are read. The IR line must therefore contain *n + 1* values for a linear chromosome with *n* genes, and *n* values for a circular one.

If G and H have different numbers of linear chromosomes, the genome with fewer gets empty linear chromosomes (only the two caps, with intergenic region 0, or `0:0` in H) until both have the same number, as in the definition of the thesis (`equalizeLinear` in [Genomes](Genomes.md)). This is done for every algorithm. The empty chromosomes are part of the output and of the cost.

## Algorithms

The matcher is always `RFRM` (rigid G, flexible H, direct or reverse matches).

| ALGO | Call |
|------|------|
| `SOAR` / `SOARComb` | `soarPartition False/True` ([PSOAR](PSOAR.md)) |
| `Greedy` / `GreedySin` | `greedyPartition False/True` ([PGreedy](PGreedy.md)) |
| `Approx` | `balancedPartition approxPartition` ([Balancer](Balancer.md) + [PApprox](PApprox.md)) |
| `FPT` | `fptPartition 600000000` ([PFpt](PFpt.md)), with a 10-minute limit |
| `ApproxLB` | `balancedLowerBound approxLowerBound`: writes only the lower bound on the number of breakpoints of both genomes (see [Balancer](Balancer.md)) |

## Output

For each instance:

1. Four lines from `writePartition`: genes of G, IRs of G, genes of H, IRs of H, with blocks separated by `|` (see [Partition](Partition.md)).
2. The block correspondences: for each block of G, the indices of the blocks of H it matches, separated by `|`.
3. `# Exact solution` or `# Partial solution`. Partial only happens with FPT when the time limit is reached.
4. `# Time: <seconds>`.

With `-a Approx`, the file `OFILE_phase2` is also written, with one line per instance (in the order of the input): the number of times the second stage of the approximation broke a block D, or `NA` if the instance failed. The name does not match `OFILE_[0-9]*`, so scripts that count the outputs of the instances are not affected.

If FPT finds no partition, the four partition lines are replaced by a comment.
