#pragma once
#include "../cycle/cycles.hpp"
#include "../heur/ga.hpp"
#include "genome.hpp"
#include "permutation.hpp"
#include <getopt.h>
#include <iostream>

struct CommonArgs {
  string input_file;    // path of input file
  string output_folder; // path of output folder
  int iterations = 100; // number of decompositions/maps to generate
  bool extend =
      false; // whether to extend (cap) the genomes before apply the algorithm
  bool fill_zero = false; // whether to fill intergenic regions with zero
                          // instead of reading them
  bool unsign = false;    // whether to use unsigned genomes (default is signed)
  bool unbalanced =
      false; // whether to allow for genomes with distinct sets of genes
};

const struct option common_longopts[] = {{"input", 1, NULL, 'i'},
                                         {"output", 1, NULL, 'o'},
                                         {"iterations", 1, NULL, 'k'},
                                         {"extend", 0, NULL, 'e'},
                                         {"help", 0, NULL, 'h'},
                                         {"zero", 0, NULL, 'z'},
                                         {"unsign", 0, NULL, 'u'},
                                         {"unbalanced", 0, NULL, 'd'},
                                         {NULL, 0, NULL, 0}};

const string common_str = "i:o:k:ezudh";

template <class T> void print_vec(ostream &os, vector<T> vec) {
  for (auto a : vec) {
    os << a << " ";
  }
}

void common_help();
bool read_common_args(CommonArgs &args, int op, char *optarg);

vector<string> *read_lines(istream &is);
/* Read the decompositions printed by dec (one "Dec: ..." record per instance,
 * e.g. the -best files or their concatenation). */
vector<string> read_decompositions(const string &file);
/* Read the partitions of a FIX file of run_part_fix_tese.sh (3 lines per
 * instance: breakpoints of the origin and of the target, as indices of
 * intergenic regions of the genomes without caps with the chromosomes
 * separated by ';' and "-" for no breakpoint, and the position in the target
 * of each gene of the origin; three lines "NA" for an instance without a
 * partition). */
vector<PartitionFix> read_partition_fixes(const string &file);
GenomePair input(string &line1, string &line2, bool extend);
GenomePair input(string &line1, string &line2, string &line3, string &line4,
                 bool extend);
MultiGenomePair input_mg(string &line1, string &line2, bool extend);
MultiGenomePair input_mg(string &line1, string &line2, string &line3,
                         string &line4, bool extend);
void output(ostream &os, int dist, double time);
void output(ostream &os, string &dist_etc, double time);
void output(ostream &os, const Permutation &perm, bool indel_zero);
void output(ostream &os, const CycleGraph &, double time);
void output(ostream &os, GenomePair &data);
void output(ostream &os, PermsIrs);
void output(ostream &os, CommonPartition);
