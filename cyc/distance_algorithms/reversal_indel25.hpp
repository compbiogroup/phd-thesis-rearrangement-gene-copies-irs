#pragma once

#include "../misc/dist.hpp"
#include "../misc/permutation.hpp"

/* 2.5-approximation for the reversal and indel distance considering
 * intergenic regions in signed strings (Alexandrino et al., 2022, "Reversal
 * and Indel Distance With Intergenic Region Information"). */
class ReversalIndel25 : public DistAlg {
public:
  int estimate_distance(Permutation &pi, mt19937 &rand_gen) override;
  unique_ptr<string> estimate_distance_str(Permutation &pi, mt19937 &rand_gen) override {
      return unique_ptr<string> (new string(to_string(estimate_distance(pi, rand_gen))));
  }
};
