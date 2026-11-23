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
          "from the origin and target signed genomes."
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
  cout << "\t-f, --fix FIX_FILE      fix the adjacencies of the partitions in "
          "FIX_FILE (written by run_part_fix_tese.sh, one partition per "
          "instance, in the order of the input): the adjacencies that are not "
          "breakpoints are kept (their genes are paired with genes that are "
          "adjacent in the other genome, or deleted/inserted together), so the "
          "packings induce the partition or one with fewer breakpoints. The "
          "reduction rules are not applied. The linear chromosomes must be "
          "extended (-e). Instances without a partition (NA) get the result "
          "NA."
       << endl
       << "\t-p, --fix-pairs         with -f, also fix the correspondence of "
          "genes of FIX_FILE (the packing is then determined)."
       << endl;
  cout << endl << endl;

  exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[]) {
  return run_main([&]() {
    CommonArgs args;
    HeurArgs heur_args;
    vector<string> pos_args;
    string fix_file;
    bool fix_pairs = false;
    vector<option> longopts = heur_longopts;
    longopts.push_back({"fix", 1, NULL, 'f'});
    longopts.push_back({"fix-pairs", 0, NULL, 'p'});
    parse_args(args, pos_args, argc, argv, 2, help, heur_str + "f:p", longopts,
               [&](int op, char *optarg) {
                 if (op == 'f') {
                   fix_file = optarg;
                   return true;
                 }
                 if (op == 'p') {
                   fix_pairs = true;
                   return true;
                 }
                 return read_heur_args(heur_args, op, optarg);
               });
    const string &heur = pos_args[0], &dec_algo = pos_args[1];
    if (fix_pairs && fix_file == "") {
      throw invalid_argument("--fix-pairs requires --fix");
    }

    random_device rd;
    mt19937 rand_gen(rd());

    unique_ptr<vector<string>> input_lines = read_input(args);
    int div = instance_lines(args);
    vector<PartitionFix> fixes;
    if (fix_file != "") {
      fixes = read_partition_fixes(fix_file);
      if (fixes.size() * div != input_lines->size()) {
        throw invalid_argument(to_string(fixes.size()) + " partitions in " +
                               fix_file + " for " +
                               to_string(input_lines->size() / div) +
                               " instances");
      }
    }

    for (size_t i = 0; i < input_lines->size(); i += div) {
      Timer timer;
      ofstream os;
      int name_idx = i / div;

      if (missing_instance(args, *input_lines, i) ||
          (!fixes.empty() && fixes[name_idx].missing)) {
        if (heur_args.all) {
          open_output(args, os, "-all", name_idx) << "Dec: NA" << endl;
        }
        open_output(args, os, "-best", name_idx) << "Dec: NA" << endl;
        open_output(args, os, "-dcj", name_idx) << "Dist: NA" << endl;
        continue;
      }

      MultiGenomePair data = read_instance_mg(args, *input_lines, i);

      // Application of suboptimal rules from the MSOAR framework
      // Uncomment the following lines to apply the rules
      // (not applied with a partition, whose pairs they could contradict)
      if (fixes.empty()) {
        suboptimal_rule_interval(*data.g, *data.h);
        suboptimal_rule_pairs(*data.g, *data.h);
      }

      CycleGraph cg(*data.g, *data.h, rand_gen, args.unsign);
      if (!fixes.empty()) {
        cg.fix_partition(fixes[name_idx], fix_pairs);
      }

      ostream *all_os =
          heur_args.all ? &open_output(args, os, "-all", name_idx) : nullptr;
      unique_ptr<CycleGraph> cg_best =
          best_decomposition(cg, heur, dec_algo, args.iterations, heur_args,
                             all_os, timer, rand_gen);

      if (cg_best->released_adjacencies() > 0) {
        cerr << "instance " << name_idx << ": "
             << cg_best->released_adjacencies()
             << " fixed adjacencies released (dead end of the heuristic)"
             << endl;
      }
      output(open_output(args, os, "-best", name_idx), *cg_best,
             timer.elapsed_time());
      timer.mark_time();

      // Add output for dcj distance
      // output(open_output(args, os, "-dcj", name_idx),
      //        cg_best->dcj_ir_distance(args.unbalanced),
      //        timer.since_last_mark());
      output(open_output(args, os, "-dcj", name_idx),
             cg_best->dcj_indel_ir_distance(), timer.since_last_mark());
    }
  });
}
