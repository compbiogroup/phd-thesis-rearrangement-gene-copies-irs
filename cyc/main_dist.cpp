#include "distance_algorithms/decomposition.hpp"
#include "distance_algorithms/r_or_rt_noir.hpp"
#include "distance_algorithms/reversal4.hpp"
#include "distance_algorithms/reversal_indel25.hpp"
#include "distance_algorithms/reversal_transposition4.hpp"
#include "distance_algorithms/reversal_transposition45.hpp"
#include "distance_algorithms/reversal_transposition_indel6.hpp"
#include "distance_algorithms/transposition3.hpp"
#include "external/external.hpp"
#include "misc/fix_flexible.hpp"
#include "misc/genome.hpp"
#include "misc/io.hpp"
#include "misc/main_common.hpp"
#include "misc/permutation.hpp"
#include "misc/seeds.hpp"
#include "misc/timer.hpp"
#include <iostream>
#include <sstream>
using namespace std;

// TODO: adapt all distance algorithms to consider the flexible intergenic regions, as they are they only work correctly if the bounds of the intergenic regions are the same
void help(char *name) {
  cout << "usage: Calculate distances for pairs of strings. If replicas are "
       << "presents multiple random mappings are generated." << endl
       << "\t" << name << " ALG [OPTIONS]" << endl
       << endl
       << "positional arguments:" << endl
       << "\tALG      the algorithm to use, the options are:" << endl
       << "\t\t - transposition3: Asymptotic 3-approximation for sorting by "
          "intergenic transpositions in unsigned permutations"
       << endl
       << "\t\t - reversal4: 4-approximation for sorting by intergenic "
          "reversals in unsigned permutations"
       << endl
       << "\t\t - reversal_indel25: 2.5-approximation for sorting by "
          "reversals and indels considering intergenic regions in signed "
          "strings (requires extended genomes)"
       << endl
       << "\t\t - reversal_transposition4: 4-approximation for sorting by "
          "intergenic reversals and transpositions in unsigned permutations"
       << endl
       << "\t\t - reversal_transposition45: 4.5-approximation for sorting by "
          "intergenic reversals and transpositions in unsigned permutations"
       << endl
       << "\t\t - reversal_transposition_indel6: 6-approximation for sorting "
          "by reversals, transpositions and indels considering intergenic "
          "regions in signed strings (requires extended genomes)"
       << endl
       << "\t\t - reversal_noir: 2-approximation for sorting by reversals and "
          "indels in signed permutations (ignores intergenic regions)"
       << endl
       << "\t\t - reversal_transposition_noir: heuristic for sorting by "
          "reversals, transposition and indels in signed permutations (ignores "
          "intergenic regions)"
       << endl
       << "\t\t - dcj: 2-approximation for the DCJ and indel distance "
          "considering intergenic regions (Siqueira et al., 2026), computed on "
          "the cycles of a decomposition of the adjacency graph (the "
          "decompositions of -D, or random mappings of the genes otherwise). "
          "Accepts genomes with multiple chromosomes; flexible intergenic "
          "regions are fixed cycle by cycle."
       << endl
       << "\t\t - the name of an executable in the external folder" << endl
       << endl;
  common_help();
  cout << "\t-D, --dec DEC_FILE      use the decompositions in DEC_FILE (the "
          "output of dec, e.g. the -best files, one per instance in the order "
          "of the input) to map the genes instead of random mappings. The "
          "input and the options -e, -z and -u must be the ones used to run "
          "dec. Only for genomes with a single chromosome."
       << endl
       << "\t-F, --fix-flexible      before calculating each distance, replace "
          "the flexible intergenic regions of the target by fixed sizes inside "
          "their intervals, as in Brito et al. (2022), making balanced as many "
          "cycles of the breakpoint graph as possible (requires extended "
          "genomes without replicas)."
       << endl
       << endl;

  exit(EXIT_SUCCESS);
}

/* Genomes of the instance starting at line i relabeled according to the
 * decomposition dec (printed by dec), so that each gene of the origin has the
 * label of the gene of the target assigned to it. */
GenomePair genomes_from_dec(const CommonArgs &args, vector<string> &lines,
                            size_t i, const string &dec, mt19937 &rand_gen) {
  MultiGenomePair mg = read_instance_mg(args, lines, i);
  CycleGraph cg(*mg.g, *mg.h, rand_gen, args.unsign);
  cg.read_decomposition(dec);
  cg.make_assignment();
  MultiGenomePair assigned = cg.get_genomes();
  if (assigned.g->n_chrom() != 1 || assigned.h->n_chrom() != 1) {
    throw invalid_argument(
        "--dec only supports genomes with a single chromosome");
  }
  GenomePair data;
  data.g = make_unique<Genome>((*assigned.g)[1]);
  data.h = make_unique<Genome>((*assigned.h)[1]);
  return data;
}

/* Return the algorithm with the given name and set duplicate to indicate the
 * cases where genes have to be duplicated to remove signs. */
unique_ptr<DistAlg> make_alg(const string &name, bool unsign,
                             bool &duplicate) {
  duplicate = false;
  if (name == "transposition3") {
    duplicate = !unsign;
    return make_unique<Transposition3>();
  } else if (name == "reversal4") {
    duplicate = !unsign;
    return make_unique<Reversal4>();
  } else if (name == "reversal_transposition4") {
    duplicate = !unsign;
    return make_unique<ReversalTransposition4>();
  } else if (name == "reversal_transposition45") {
    duplicate = !unsign;
    return make_unique<ReversalTransposition45>();
  } else if (name == "reversal_transposition_indel6") {
    return make_unique<ReversalTranspositionIndel6>();
  } else if (name == "reversal_indel25") {
    return make_unique<ReversalIndel25>();
  } else if (name == "decomposition") {
    return make_unique<Decomposition>();
  } else if (name == "reversal_noir") {
    return make_unique<ReversalNOIR>();
  } else if (name == "reversal_transposition_noir") {
    return make_unique<ReversalTranspositionNOIR>();
  }
  return make_unique<ExternalDistAlg>("external/" + name);
}

int main(int argc, char *argv[]) {
  return run_main([&]() {
    CommonArgs args;
    vector<string> pos_args;
    string dec_file;
    bool fix_flexible = false;
    parse_args(args, pos_args, argc, argv, 1, help, "D:F",
               {{"dec", 1, NULL, 'D'}, {"fix-flexible", 0, NULL, 'F'}},
               [&](int op, char *optarg) {
                 if (op == 'D') {
                   dec_file = optarg;
                   return true;
                 }
                 if (op == 'F') {
                   fix_flexible = true;
                   return true;
                 }
                 return false;
               });

    random_device rd;
    mt19937 rand_gen(rd());

    bool duplicate = false;
    bool dcj = pos_args[0] == "dcj";
    unique_ptr<DistAlg> alg;
    if (!dcj) {
      alg = make_alg(pos_args[0], args.unsign, duplicate);
    }

    unique_ptr<vector<string>> input_lines = read_input(args);
    int div = instance_lines(args);
    int n_inst = input_lines->size() / div;
    vector<string> decs;
    if (dec_file != "") {
      decs = read_decompositions(dec_file);
      if (int(decs.size()) != n_inst) {
        throw invalid_argument(to_string(decs.size()) +
                               " decompositions in " + dec_file + " for " +
                               to_string(n_inst) + " instances");
      }
    }

    // Each instance has its own generator, its output (when there is no output
    // folder) is buffered and printed in input order after the loop, and its
    // errors are reported after the loop (an exception must not leave an
    // OpenMP loop).
    vector<mt19937::result_type> seeds = draw_seeds(rand_gen, n_inst);
    vector<string> outputs(n_inst), errors(n_inst);

#pragma omp parallel for schedule(dynamic) // Process each instance in parallel
    for (int k = 0; k < n_inst; ++k) {
      try {
        mt19937 inst_gen(seeds[k]);
        Timer timer;
        ofstream os;
        ostringstream buf;
        auto out = [&](const string &suffix) -> ostream & {
          if (args.output_folder == "") {
            return buf;
          }
          return open_output(args, os, suffix, k);
        };
        // a missing instance or a decomposition NA (dec with a FIX file
        // without the partition of the instance)
        if (missing_instance(args, *input_lines, k * div) ||
            (dec_file != "" && decs[k].compare(0, 2, "NA") == 0)) {
          out("-all") << "Dist: NA" << endl;
          out("-best") << "NA" << endl << "Dist: NA" << endl;
          outputs[k] = buf.str();
          continue;
        }

        int dist, dist_best = std::numeric_limits<int>::max();

        if (dcj) {
          /* DCJ and indel distance on the cycles of a decomposition: the one
           * of -D or random mappings of the genes (best of the iterations) */
          MultiGenomePair mg = read_instance_mg(args, *input_lines, k * div);
          ostream &os_all = out("-all");
          string dec_best;
          int iterations = (dec_file == "") ? args.iterations : 1;
          for (int j = 1; j <= iterations; ++j) {
            timer.mark_time();
            CycleGraph cg(*mg.g, *mg.h, inst_gen, args.unsign);
            if (dec_file == "") {
              cg.decompose(true, "map");
            } else {
              cg.read_decomposition(decs[k]);
              cg.rebuild_cycles();
            }
            dist = cg.dcj_indel_ir_distance();
            output(os_all, dist, timer.since_last_mark());
            if (dist < dist_best) {
              dist_best = dist;
              ostringstream dec;
              dec << cg;
              dec_best = dec.str();
            }
          }
          ostream &os_best = out("-best");
          os_best << "Dec: " << dec_best << endl;
          output(os_best, dist_best, timer.elapsed_time());
          outputs[k] = buf.str();
          continue;
        }

        unique_ptr<Permutation> pi, pi_best;

        GenomePair data =
            (dec_file == "")
                ? read_instance(args, *input_lines, k * div)
                : genomes_from_dec(args, *input_lines, k * div, decs[k],
                                   inst_gen);

        ostream &os_all = out("-all");
        for (int j = 1; j <= args.iterations; ++j) {
          timer.mark_time();
          pi.reset(new Permutation(*data.g, *data.h, duplicate, inst_gen));
          if (fix_flexible) {
            fix_flexible_irs(*pi);
          }
          dist = alg->estimate_distance(*pi, inst_gen);
          output(os_all, dist, timer.since_last_mark());
          if (dist < dist_best) {
            dist_best = dist;
            pi_best = std::move(pi);
          }
        }

        ostream &os_best = out("-best");
        os_best << *pi_best << endl;
        output(os_best, dist_best, timer.elapsed_time());
        outputs[k] = buf.str();
      } catch (const exception &e) {
        errors[k] = e.what();
      }
    }

    for (int k = 0; k < n_inst; ++k) {
      if (errors[k] != "") {
        throw runtime_error("instance " + to_string(k) + ": " + errors[k]);
      }
      cout << outputs[k];
    }
  });
}
