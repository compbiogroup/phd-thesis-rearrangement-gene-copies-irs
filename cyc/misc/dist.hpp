#pragma once

#include "permutation.hpp"

class DistAlg {
public:
    virtual ~DistAlg() = default;
    virtual int estimate_distance(Permutation &pi, mt19937 &rand_gen) = 0;
    virtual unique_ptr<string> estimate_distance_str(Permutation &pi, mt19937 &rand_gen) = 0;
};
