#pragma once

#include "../misc/dist.hpp"
#include "../misc/permutation.hpp"

class R_OR_RT_NOIR : public DistAlg {
public:
  int estimate_distance(Permutation &pi, mt19937 &rand_gen) override;
  unique_ptr<string> estimate_distance_str(Permutation &pi, mt19937 &rand_gen) override;
  GenomePair make_genomes(Permutation pi, mt19937 &rand_gen);
  virtual int dist_aux(int *g1, int *g2, int size) = 0;
  virtual int lower_bound(GenomePair &data, unique_ptr<CycleGraph> &cg) = 0;
};

class ReversalNOIR : public R_OR_RT_NOIR {
public:
  int dist_aux(int *g1, int *g2, int size) override;
  int lower_bound(GenomePair &data, unique_ptr<CycleGraph> &cg) override;
};

class ReversalTranspositionNOIR : public R_OR_RT_NOIR {
public:
  int dist_aux(int *g1, int *g2, int size) override;
  int lower_bound(GenomePair &data, unique_ptr<CycleGraph> &cg) override;
};

class TranspositionNOIR : public R_OR_RT_NOIR {
public:
  int dist_aux(int *g1, int *g2, int size) override;
  int lower_bound(GenomePair &data, unique_ptr<CycleGraph> &cg) override;
};
