#include "decomposition.hpp"
#include "../misc/permutation.hpp"
#include "auxiliar.hpp"
#include <cassert>
#include <limits>
#include <vector>

int Decomposition::estimate_distance(Permutation &pi, mt19937 &rand_gen) {
    GenomePair data = pi.split_iota();
    unique_ptr<CycleGraph> cg = unique_ptr<CycleGraph>(new CycleGraph(MultiGenome(std::move(data.g)), MultiGenome(std::move(data.h)), rand_gen));
    cg->decompose_with_bfs(true);
    return cg->gain()[0];
}
