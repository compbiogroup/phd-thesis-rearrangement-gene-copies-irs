#pragma once

#include "../misc/dist.hpp"
#include "../misc/permutation.hpp"

/* 4-approximation for sorting by intergenic reversals and transpositions in
 * unsigned permutations (C++ version of external/reversal_transposition4.py,
 * without the greedy strategy, which is disabled in the Python script). */
class ReversalTransposition4 : public DistAlg {
public:
  int estimate_distance(Permutation &pi, mt19937 &rand_gen) override;

  unique_ptr<string> estimate_distance_str(Permutation &pi, mt19937 &rand_gen) override {
      return unique_ptr<string> (new string(to_string(estimate_distance(pi, rand_gen))));
  }
};
