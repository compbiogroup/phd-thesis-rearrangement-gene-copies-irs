#pragma once
#include <random>
#include <vector>
using namespace std;

/* Draw n seeds from rand_gen. Used before a parallel loop so each iteration
 * gets its own generator: there is no shared state between threads and the
 * result does not depend on the number of threads. */
inline vector<mt19937::result_type> draw_seeds(mt19937 &rand_gen, size_t n) {
  vector<mt19937::result_type> seeds(n);
  for (auto &s : seeds) {
    s = rand_gen();
  }
  return seeds;
}
