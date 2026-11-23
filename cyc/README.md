# Heuristics for Cycle Packing of Adjacency Graphs and Approximation Algorithms to Rearrangement Distance on Permutations

This repository has the implementation of two heuristics for the Adjacency Graph Packing problem. One is called Random Packings and tries to solve the problem by randomly generating cycle packings of the Adjacency Graph. The other is based on the Genetic Algorithm metaheuristic.

The following algorithms to approximate rearrangement distances are also implemented:
- A 3 asymptotic approximation algorithm for Sorting by Intergenic Transposition in unsigned permutations.
- An 4-approximation algorithm for Sorting by Intergenic Reversals [[1]](#1) in unsigned permutations.
- An 4.5-approximation algorithm for Sorting by Intergenic Reversals and Transposition [[1]](#1).
- An 4-approximation algorithm for Sorting by Intergenic Reversals and Transposition in unsigned permutations (`reversal_transposition4`, a C++ version of `external/reversal_transposition4.py`). Only the lower bound of the intergenic regions is used and the genomes must have the same genes and the same total of nucleotides.
- An 2-approximation for Reversal and Indel Distance in signed strings. This algorithm uses an implementation of an exact algorithm for Sorting by Reversal [[2]](#2).
- An 2.5-approximation for Reversal and Indel Distance considering intergenic regions in signed strings [[4]](#4) (`reversal_indel25`, requires extended genomes). The intergenic regions of the target may be flexible (intervals), only the lower bound of the origin intergenic regions is used.
- An 6-approximation for Reversal, Transposition and Indel Distance considering intergenic regions [[5]](#5) adapted to signed strings (`reversal_transposition_indel6`, requires extended genomes). Signed breakpoints are used, and reversals on elements with consecutive absolute values and different signs replace the reversals on decreasing strips. The intergenic regions of the target may be flexible (intervals), only the lower bound of the origin intergenic regions is used.
- An heuristic for Reversal, Transposition and Indels Distance in signed strings. This algorithm uses an implementation of a 2-approximation algorithm for Sorting by Reversal and Transposition [[3]](#3).
- Any executable in the external folder. This executable must be called with three arguments corresponding to the permutation, the target intergenic regions, and the source intergenic regions (as comma separated lists of integers). The permutation starts with 0 and ends with the size of the target (caps), genes to be deleted are represented by 0 and genes to be inserted are the labels missing in the permutation. Only the lower bound of flexible intergenic regions is sent. The result of the executable must be a single line containing the distance.

The algorithms can be applied to genomes with multiple genes, in that case random mappings of the genomes into permutations will be used to produce the distances.

With the options `-f FIX_FILE` and `-p`, `dec` fixes the adjacencies (and, with `-p`, the correspondence of genes) of genome partitions, given in the FIX files of `scripts/run_part_fix_tese.sh` (see `scripts/PARTICAO_CYC.md`): the adjacencies that are not breakpoints of the partitions cannot become breakpoints, so the packings induce the partition (or one with fewer breakpoints) while the heuristics still choose the bijection between compatible blocks.

With the option `-F` (`--fix-flexible`), `dist` replaces the flexible intergenic regions of the target by fixed sizes inside their intervals before calculating each distance, adapting the reduction of Brito et al. [[6]](#6) to genomes with distinct gene sets: the stable cycles of the breakpoint graph that do not need to give or receive nucleotides become balanced, and the unstable cycles get the minimum (maximum) sizes when they need more (less) nucleotides, which are taken from (given to) auxiliary cycles. A solution for the fixed instance is also a solution for the flexible one, so algorithms that only accept fixed intergenic regions (such as `reversal_transposition4`) can be used with flexible ones.

## Usage

Compile the code by running `make` and see the running options with `./dec --help` (for the cycle packing) or `./dist --help` (for the rearrangement distances). To use other algorithm when calculating the distances include a executable in the `external` folder and pass its name as the algorithm parameter. In that case, the executable should receive one instance as command line arguments (the three comma separated lists described above) and produce the distance in the standard output.

## References

<a id="1">[1]</a> 
Klairton Lima Brito, Géraldine Jean, Guillaume Fertin, Andre Rodrigues Oliveira, Ulisses Dias and Zanoni Dias. Sorting by genome rearrangements on both gene order and intergenic sizes. Journal of Computational Biology 2020;27(2):156–74.

<a id="2">[2]</a> 
Anne Bergeron, Julia Mixtacki and Jens Stoye. “Reversal Distance without Hurdles and Fortresses”. Combinatorial Pattern Matching. Berlin, Heidelberg: Springer.

<a id="3">[3]</a> 
Maria Emília M.T. Walter, Zanoni Dias and João Meidanis. “Reversal and transposition distance of linear chromosomes”. Proceedings of the String Processing and Information Retrieval: A South American Symposium (SPIRE’1998). Los Alamitos, CA, USA, 1998, pp. 96–102.

<a id="4">[4]</a> 
Alexsandro Oliveira Alexandrino, Klairton Lima Brito, Andre Rodrigues Oliveira, Ulisses Dias and Zanoni Dias. “Reversal and Indel Distance With Intergenic Region Information”. IEEE/ACM Transactions on Computational Biology and Bioinformatics 2023;20(3):1628–40.

<a id="5">[5]</a> 
Alexsandro Oliveira Alexandrino, Andre Rodrigues Oliveira, Ulisses Dias and Zanoni Dias. “Incorporating intergenic regions into reversal and transposition distances with indels”. Journal of Bioinformatics and Computational Biology 2021;19(6):2140011.

<a id="6">[6]</a> 
Klairton Lima Brito, Alexsandro Oliveira Alexandrino, Andre Rodrigues Oliveira, Ulisses Dias and Zanoni Dias. “Genome Rearrangement Distance with a Flexible Intergenic Regions Aspect”. IEEE/ACM Transactions on Computational Biology and Bioinformatics 2022.
