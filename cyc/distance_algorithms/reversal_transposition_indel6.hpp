#pragma once

#include "../misc/dist.hpp"
#include "../misc/permutation.hpp"

/* 6-approximation for the reversal, transposition and indel distance
 * considering intergenic regions (Algorithm 1 of Alexandrino et al., 2021,
 * "Incorporating intergenic regions into reversal and transposition distances
 * with indels", using Lemma 9 when there are only increasing strips) adapted
 * to signed strings. */
class ReversalTranspositionIndel6 : public DistAlg {
public:
  int estimate_distance(Permutation &pi, mt19937 &rand_gen) override;
  unique_ptr<string> estimate_distance_str(Permutation &pi, mt19937 &rand_gen) override {
      return unique_ptr<string> (new string(to_string(estimate_distance(pi, rand_gen))));
  }
};
