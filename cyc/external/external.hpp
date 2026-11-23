#pragma once
#include "../misc/dist.hpp"

/* Call to an external program to calculate the distances. The program must
 * receive three comma separated list as command line arguments (the
 * permutation, the target intergenic regions, and the origin intergenic
 * regions). The permutation starts with 0 and ends with the size of the target
 * (caps), genes to be deleted are 0 and genes to be inserted are the labels
 * missing in the permutation. The program must print the distance to its
 * stdout.*/
class ExternalDistAlg : public DistAlg {
  string prog;

public:
  ExternalDistAlg(string prog) : prog(prog){};
  int estimate_distance(Permutation &pi, mt19937 &rand_gen) override;

  unique_ptr<string> estimate_distance_str(Permutation &pi, mt19937 &rand_gen) override {
      return make_unique<string> (to_string(estimate_distance(pi, rand_gen)));
  }
};
