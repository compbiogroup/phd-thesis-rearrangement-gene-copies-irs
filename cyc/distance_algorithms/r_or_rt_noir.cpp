#include "r_or_rt_noir.hpp"

#include <cassert>
#include <cmath>
#include <sstream>

#include "../cycle/cycles.hpp"
#include "../misc/io.hpp"
#include "../misc/permutation.hpp"
#include "auxiliar.hpp"

extern "C" {
#include "../GRIMM/grimm.h"
#include "perm/perm_rearrange.h"
}

/* Position of genome g after which run is inserted. The genomes here have a
 * single linear chromosome, so a run after the left cap goes after the first
 * gene. */
static int insertion_pos(const Genome &g, const Run &run) {
  if (run.after_cap) {
    return 1;
  }
  vector<int> pos = g.pos(run.fst_gene);
  if (pos.empty()) {
    throw runtime_error("insertion position not found: gene " +
                        to_string(run.fst_gene));
  }
  return pos[0];
}

int lower_bound_gen(GenomePair &data, unique_ptr<CycleGraph> &cg, int div) {
  unordered_map<int,int> alp;
  data.g->alphabet(alp,false);
  data.h->alphabet(alp,true);
  int inter = 0;
  for (auto el : alp) {
    if (el.second == 0) inter++;
  }
  return ceil((inter - cg->dec_size() + cg->potation() - 1) / double(div));
}

int ReversalNOIR::lower_bound(GenomePair &data, unique_ptr<CycleGraph> &cg) {
    return lower_bound_gen(data, cg, 1);
}

int ReversalTranspositionNOIR::lower_bound(GenomePair &data, unique_ptr<CycleGraph> &cg) {
    return lower_bound_gen(data, cg, 2);
}

int TranspositionNOIR::lower_bound(GenomePair &data, unique_ptr<CycleGraph> &cg) {
    return lower_bound_gen(data, cg, 2);
}


int R_OR_RT_NOIR::estimate_distance(Permutation &pi, mt19937 &rand_gen) {
  unique_ptr<CycleGraph> cg;
  GenomePair data;
  int dist = 0;

  assert(pi.occ_max() == 1);
  data = pi.split_iota(); // create a permutation iota (except the last gene is bigger than any gene in pi)
  cg = make_unique<CycleGraph>(MultiGenome(make_unique<Genome>(*data.g)), make_unique<Genome>(*data.h), rand_gen);
  cg->decompose_with_bfs(true);
  // Track if the number of cycles is increasing or the potential in decreasing
  int last_lb = lower_bound_gen(data, cg, 1);

  bool found_run;
  do {
    found_run = false;
    for (auto c : cg->cycle_list()) {
      Run run = cg->cycle_run(c.start_vtx);
      if (!run.genes_to_add.empty()) {
        if (run.genome == 'G') {
          int i = insertion_pos(*data.g, run);
          data.g->insertion(i, run.genes_to_add);
        } else {
          int i = insertion_pos(*data.h, run);
          data.h->insertion(i, run.genes_to_add);
        }
        dist++;
        found_run = true;
        break;
      }
    }
    if (found_run) {
      cg.reset(new CycleGraph(MultiGenome(make_unique<Genome>(*data.g)), make_unique<Genome>(*data.h), rand_gen));
      cg->decompose_with_bfs(true);
      int lb = lower_bound_gen(data, cg, 1);
      assert(lb < last_lb);
      last_lb = lb;
    }
  } while (found_run);

  assert(data.g->size() == data.h->size());
  assert(data.g->balanced(*data.h));
  assert(data.g->occ_max() == 1);


  pair<int *, int *> perms;
  perms.first = data.g->get_perm();
  perms.second = data.h->get_perm();
  dist += dist_aux(perms.first, perms.second, data.g->size());
  delete[] perms.first;
  delete[] perms.second;

  return dist;
}

unique_ptr<string> R_OR_RT_NOIR::estimate_distance_str(Permutation &pi, mt19937 &rand_gen) {
  unique_ptr<CycleGraph> cg;
  string out = "(";
  GenomePair data;
  int dist = 0;

  assert(pi.occ_max() == 1);
  data = pi.split_iota(); // create a permutation iota (except the last gene is bigger than any gene in pi)
  cg = make_unique<CycleGraph>(MultiGenome(make_unique<Genome>(*data.g)), make_unique<Genome>(*data.h), rand_gen);
  cg->decompose_with_bfs(true);
  // Each insertion of a run must decrease the numerator of the lower bound
  // (the bound itself may not decrease when it is divided by 2).
  int last_lb = lower_bound_gen(data, cg, 1);

  out += to_string(cg->dec_size());
  out += ", " + to_string(cg->potation());

  bool found_run;
  do {
    found_run = false;
    for (auto c : cg->cycle_list()) {
      Run run = cg->cycle_run(c.start_vtx);
      if (!run.genes_to_add.empty()) {
        if (run.genome == 'G') {
          int i = insertion_pos(*data.g, run);
          data.g->insertion(i, run.genes_to_add);
        } else {
          int i = insertion_pos(*data.h, run);
          data.h->insertion(i, run.genes_to_add);
        }
        dist++;
        found_run = true;
        break;
      }
    }
    if (found_run) {
      cg.reset(new CycleGraph(MultiGenome(make_unique<Genome>(*data.g)), make_unique<Genome>(*data.h), rand_gen));
      cg->decompose_with_bfs(true);
      int lb = lower_bound_gen(data, cg, 1);
      assert(lb < last_lb);
      last_lb = lb;
    }
  } while (found_run);

  assert(data.g->size() == data.h->size());
  assert(data.g->balanced(*data.h));
  assert(data.g->occ_max() == 1);

  out += ", " + to_string(dist);

  pair<int *, int *> perms;
  perms.first = data.g->get_perm();
  perms.second = data.h->get_perm();
  dist += dist_aux(perms.first, perms.second, data.g->size());
  delete[] perms.first;
  delete[] perms.second;

  out += ", " + to_string(dist);
  out += ")";
  return make_unique<string> (out);
}

GenomePair R_OR_RT_NOIR::make_genomes(Permutation pi, mt19937 &rand_gen) {
  unique_ptr<CycleGraph> cg;
  GenomePair data;

  assert(pi.occ_max() == 1);
  data = pi.split_iota(); // create a permutation iota (except the last gene is bigger than any gene in pi)
  cg = make_unique<CycleGraph>(MultiGenome(make_unique<Genome>(*data.g)), make_unique<Genome>(*data.h), rand_gen);
  cg->decompose_with_bfs(true);
  // Each insertion of a run must decrease the numerator of the lower bound
  // (the bound itself may not decrease when it is divided by 2).
  int last_lb = lower_bound_gen(data, cg, 1);

  bool found_run;
  do {
    found_run = false;
    for (auto c : cg->cycle_list()) {
      Run run = cg->cycle_run(c.start_vtx);
      if (!run.genes_to_add.empty()) {
        if (run.genome == 'G') {
          int i = insertion_pos(*data.g, run);
          data.g->insertion(i, run.genes_to_add);
        } else {
          int i = insertion_pos(*data.h, run);
          data.h->insertion(i, run.genes_to_add);
        }
        found_run = true;
        break;
      }
    }
    if (found_run) {
      cg.reset(new CycleGraph(MultiGenome(make_unique<Genome>(*data.g)), make_unique<Genome>(*data.h), rand_gen));
      cg->decompose_with_bfs(true);
      int lb = lower_bound_gen(data, cg, 1);
      assert(lb < last_lb);
      last_lb = lb;
    }
  } while (found_run);

  /* assert(data.g->size() == data.h->size()); */
  /* assert(data.g->balanced(*data.h)); */
  assert(data.g->occ_max() == 1);

  return data;
}


int ReversalNOIR::dist_aux(int *g1, int *g2, int size) {
  return grimm(g1, g2, size);
}

int ReversalTranspositionNOIR::dist_aux(int *g1, int *g2, int size) {
  return dist(g1, g2, size, 1, 2, 1);
}

int TranspositionNOIR::dist_aux(int *g1, int *g2, int size) {
  int dist;
  stringstream ss;

  ss << "distance_algorithms/perm/perm_rearrange.py" << " ";

  for (size_t i = 0; i < size-1; i++) {
    ss << g1[i] << ",";
  }
  ss << g1[size - 1] << " ";
  for (size_t i = 0; i < size-1; i++) {
    ss << g2[i] << ",";
  }
  ss << g2[size - 1] << " ";
  FILE *prog_pipe;
  if ((prog_pipe = popen(ss.str().c_str(), "r")) == NULL) {
    throw runtime_error("Error calling child process.");
  }
  auto ret = fscanf(prog_pipe, "%d", &dist);
  if (ret != 1 || pclose(prog_pipe) != EXIT_SUCCESS) {
    throw runtime_error("Error on child process.");
  }
  return dist;
}
