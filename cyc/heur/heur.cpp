#include "heur.hpp"
#include "ga.hpp"
#include <iostream>
#include <stdexcept>

void heur_help() {
  cout << "\t-m, --mutation MUT      mutation rate in percentage for GA "
          "(default 50%)"
       << endl
       << "\t-c, --crossover CROS    crossover rate in percentage for GA "
          "(default 50%)"
       << endl
       << "\t-t, --tournament TOUR   number of elements in each tournament for "
          "GA (default 2)"
       << endl
       << "\t-a, --all               print all cycles" << endl;
}

bool read_heur_args(HeurArgs &args, int op, char *optarg) {
  switch (op) {
  case 'm':
    args.mutation_rate = atoi(optarg);
    break;
  case 'c':
    args.crossover_rate = atoi(optarg);
    break;
  case 't':
    args.tournament_size = atoi(optarg);
    break;
  case 'a':
    args.all = true;
    break;
  default:
    return false; // Unknown option
  }
  return true;
}

/* Decomposition found in iteration idx of the rand heuristic. */
struct RandCand {
  CycleGraph *cg;
  int idx;
};

/* Keep the candidate with the highest gain (ties go to the lowest iteration,
 * so the result does not depend on the order in which threads finish) and
 * delete the other. */
static RandCand best_cand(RandCand a, RandCand b) {
  if (a.cg == nullptr) {
    return b;
  }
  if (b.cg == nullptr) {
    return a;
  }
  bool a_wins = a.cg->gain() > b.cg->gain() ||
                (a.cg->gain() == b.cg->gain() && a.idx < b.idx);
  delete (a_wins ? b.cg : a.cg);
  return a_wins ? a : b;
}

unique_ptr<CycleGraph> best_decomposition(const CycleGraph &cg,
                                          const string &heur,
                                          const string &dec_algo,
                                          int iterations, const HeurArgs &args,
                                          ostream *all_os, Timer &timer,
                                          mt19937 &rand_gen) {
  if (heur == "rand") {
    RandCand best = {new CycleGraph(cg), 0};
    best.cg->decompose(false, dec_algo);

    vector<mt19937::result_type> seeds = draw_seeds(rand_gen, iterations);
#pragma omp declare reduction(select_cand : RandCand : omp_out =              \
                                  best_cand(omp_in, omp_out))                  \
    initializer(omp_priv = RandCand{nullptr, 0})
#pragma omp parallel for reduction(                                            \
        select_cand : best) // Process each decomposition in parallel
    for (int i = 1; i < iterations; ++i) {
      CycleGraph *cand = new CycleGraph(cg);
      cand->reseed(seeds[i]);
      cand->decompose(true, dec_algo);
      best = best_cand({cand, i}, best);
    }
    return unique_ptr<CycleGraph>(best.cg);
  }

  if (heur == "ga") {
    int start = max(2, iterations / 10);
    if (start % 2 == 1) {
      start += 1;
    }
    int gens = max(0, (iterations - start) / start);
    GA ga = GA(make_unique<Chromossome>(cg), dec_algo,
               args.mutation_rate / 100.0, args.crossover_rate / 100.0,
               args.tournament_size, start, start, gens, all_os, timer,
               rand_gen);
    ga.solve(timer);
    return ga.get_best_chr();
  }

  throw invalid_argument("Unknown heuristic: " + heur);
}
