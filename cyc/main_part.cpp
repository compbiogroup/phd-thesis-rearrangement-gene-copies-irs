#include <iostream>
#include <random>

#include "cycle/cycles.hpp"
#include "heur/heur.hpp"
#include "misc/io.hpp"
#include "misc/main_common.hpp"
#include "misc/reduction_rules.hpp"
#include "misc/timer.hpp"
using namespace std;

void help(char *name) {
  cout << "usage: Find a cycle packing for the generalized adjacency graph "
          "formed "
          "from the origin and target signed genomes and produce a partition "
          "from it."
       << endl
       << "\t" << name << " HEUR DEC_ALGO [OPTIONS]" << endl
       << endl
       << "positional arguments:" << endl
       << "\tHEUR                    the algorithm to use (ga|rand)" << endl
       << "\tDEC_ALGO                the decomposition algorithm to use "
          "(bfs|1c2c3c|map). So far 1c2c3c is only for unichromossomal genomes "
          "with a single copy of each gene."
       << endl
       << endl;
  common_help();
  heur_help();
  cout << endl << endl;

  exit(EXIT_SUCCESS);
}

// TODO: adapt partition for multiple chromosomes
int main(int argc, char *argv[]) {
  return run_main([&]() {
    CommonArgs args;
    HeurArgs heur_args;
    vector<string> pos_args;
    parse_args(args, pos_args, argc, argv, 2, help, heur_str, heur_longopts,
               [&](int op, char *optarg) {
                 return read_heur_args(heur_args, op, optarg);
               });
    const string &heur = pos_args[0], &dec_algo = pos_args[1];

    random_device rd;
    mt19937 rand_gen(rd());

    unique_ptr<vector<string>> input_lines = read_input(args);
    int div = instance_lines(args);

    for (size_t i = 0; i < input_lines->size(); i += div) {
      Timer timer;
      ofstream os;
      int name_idx = i / div;

      GenomePair data = read_instance(args, *input_lines, i);
      MultiGenome g(std::move(data.g)), h(std::move(data.h));

      // Application of suboptimal rules from the MSOAR framework
      // Uncomment the following lines to apply the rules
      // suboptimal_rule_interval(g, h);
      // suboptimal_rule_pairs(g, h);

      CycleGraph cg(g, h, rand_gen, args.unsign, true);

      ostream *all_os =
          heur_args.all ? &open_output(args, os, "-all", name_idx) : nullptr;
      unique_ptr<CycleGraph> cg_best =
          best_decomposition(cg, heur, dec_algo, args.iterations, heur_args,
                             all_os, timer, rand_gen);

      output(open_output(args, os, "-best", name_idx), *cg_best,
             timer.elapsed_time());

      /* write the breakpoints in the partition correspondent to the best
       * decomposition. Attention: for now this answer is only valid if the
       * genomes are balanced */
      output(open_output(args, os, "-part", name_idx), cg_best->get_part(true));
    }
  });
}
