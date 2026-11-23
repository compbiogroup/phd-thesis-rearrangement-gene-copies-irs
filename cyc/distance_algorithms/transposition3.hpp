#pragma once

#include "../misc/dist.hpp"
#include "../misc/permutation.hpp"

void transposition_move(Permutation &, Gene, Gene, Gene, IR, IR, IR);

class Transposition3 : public DistAlg {
public:
  int estimate_distance(Permutation &pi, mt19937 &rand_gen) override;

  unique_ptr<string> estimate_distance_str(Permutation &pi, mt19937 &rand_gen) override {
      return unique_ptr<string> (new string(to_string(estimate_distance(pi, rand_gen))));
  }
};
