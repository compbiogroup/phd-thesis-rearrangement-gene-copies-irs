#include "reversal_transposition4.hpp"
#include <algorithm>
#include <array>
#include <climits>
#include <deque>
#include <numeric>
#include <stdexcept>
#include <vector>

/* Translation of external/reversal_transposition4.py. The instance uses the
 * same indices of the Python code: pi has the caps 0 and n+1, bpi[i] is the
 * intergenic region between pi[i-1] and pi[i] and biota[x] is the intergenic
 * region between x-1 and x in the target (bpi[0] and biota[0] are dummies).
 * Only the lower bound of the intergenic regions is used. */

namespace {

/* Position of v in vec (as list.index in Python). */
size_t index_of(const vector<int> &vec, int v) {
  auto it = find(vec.begin(), vec.end(), v);
  if (it == vec.end()) {
    throw runtime_error("reversal_transposition4: breakpoint not found");
  }
  return it - vec.begin();
}

int at(const vector<int> &vec, size_t k) {
  if (k >= vec.size()) {
    throw runtime_error("reversal_transposition4: breakpoint out of range");
  }
  return vec[k];
}

bool contains(const vector<int> &vec, int v) {
  return find(vec.begin(), vec.end(), v) != vec.end();
}

typedef array<int, 3> Params;

/* Parameters of two transpositions moving a, b, and c nucleotides into x, y,
 * and z (maximum flow with Edmonds-Karp, as in external/MF.py). */
pair<Params, Params> redistribution(int a, int b, int c, int x, int y, int z) {
  if (a + b + c != x + y + z) {
    throw runtime_error("reversal_transposition4: invalid redistribution");
  }
  const int n = 11, source = 0, sink = 10;
  int t = a + b + c;
  int g[n][n] = {{0, a, b, c, 0, 0, 0, 0, 0, 0, 0},
                 {0, 0, 0, 0, t, t, 0, 0, 0, 0, 0},
                 {0, 0, 0, 0, t, 0, t, 0, 0, 0, 0},
                 {0, 0, 0, 0, 0, t, t, 0, 0, 0, 0},
                 {0, 0, 0, 0, 0, 0, 0, t, t, 0, 0},
                 {0, 0, 0, 0, 0, 0, 0, t, 0, t, 0},
                 {0, 0, 0, 0, 0, 0, 0, 0, t, t, 0},
                 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x},
                 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, y},
                 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, z},
                 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  vector<int> parent(n, -1);
  auto bfs = [&]() {
    vector<bool> visited(n, false);
    deque<int> queue = {source};
    visited[source] = true;
    while (!queue.empty()) {
      int u = queue.front();
      queue.pop_front();
      for (int v = 0; v < n; ++v) {
        if (!visited[v] && g[u][v] > 0) {
          queue.push_back(v);
          visited[v] = true;
          parent[v] = u;
        }
      }
    }
    return bool(visited[sink]);
  };
  while (bfs()) {
    int flow = INT_MAX;
    for (int v = sink; v != source; v = parent[v]) {
      flow = min(flow, g[parent[v]][v]);
    }
    for (int v = sink; v != source; v = parent[v]) {
      g[parent[v]][v] -= flow;
      g[v][parent[v]] += flow;
    }
  }
  return make_pair(Params{{g[4][1], g[6][2], g[5][3]}},
                   Params{{g[7][4], g[9][5], g[8][6]}});
}

class Instance {
public:
  Instance(vector<int> pi, vector<int> bpi, vector<int> biota);
  int run();

private:
  int size;
  int ops = 0;
  vector<int> pi, bpi, biota;
  vector<int> breakpoints, soft, overcharged, undercharged;

  void reversal(int i, int j, int x, int y);
  void transposition(int i, int j, int k, int x, int y, int z);
  void apply_reversal(int i, int j, int x, int y) {
    reversal(i, j, x, y);
    ops++;
  }
  void apply_transposition(int i, int j, int k, int x, int y, int z) {
    transposition(i, j, k, x, y, z);
    ops++;
  }
  /* Target intergenic region of the adjacency (p, p + 1) of pi. */
  int target(int p) const { return biota[max(pi[p], pi[p + 1])]; }
  bool is_block(int a, int b) const;
  void update();
  void two_transpositions(int i, int j, int k, int x, int y, int z);
  template <class F> bool soft_pairs(F f);

  bool case_i();
  bool soft_case_i(int i, int j);
  bool soft_case_ii(int i, int j);
  bool soft_case_iii(int i, int j);
  bool case_ii();
  bool case_iii();
  bool classic_reversal(int i, int j);
  bool classic_transposition(int i, int j);
  bool create_undercharged();
  bool case_iv();
};

Instance::Instance(vector<int> pi_, vector<int> bpi_, vector<int> biota_)
    : pi(pi_) {
  if (bpi_.size() != biota_.size()) {
    throw invalid_argument("The total of intergenic regions in source and "
                           "target genomes is different.");
  }
  if (accumulate(bpi_.begin(), bpi_.end(), 0) !=
      accumulate(biota_.begin(), biota_.end(), 0)) {
    throw invalid_argument("The total of nucleotides in source and target "
                           "genomes is different.");
  }
  if (*min_element(bpi_.begin(), bpi_.end()) < 0 ||
      *min_element(biota_.begin(), biota_.end()) < 0) {
    throw invalid_argument("There is a negative intergenic region size.");
  }
  if (pi.size() - 1 != bpi_.size()) {
    throw invalid_argument("There is an incompatible number of genes and "
                           "intergenic regions.");
  }
  vector<int> sorted_pi(pi);
  sort(sorted_pi.begin(), sorted_pi.end());
  for (size_t i = 0; i < sorted_pi.size(); ++i) {
    if (sorted_pi[i] != int(i)) {
      throw invalid_argument(
          "The order of the genes is not an unsigned permutation.");
    }
  }
  size = pi.size() - 2;
  bpi.push_back(0);
  bpi.insert(bpi.end(), bpi_.begin(), bpi_.end());
  biota.push_back(0);
  biota.insert(biota.end(), biota_.begin(), biota_.end());
}

void Instance::reversal(int i, int j, int x, int y) {
  while (i <= j) {
    if (x < 0 || x > bpi[i]) {
      throw runtime_error("Intergenic Reversal: Invalid Parameter x.");
    }
    if (y < 0 || y > bpi[j + 1]) {
      throw runtime_error("Intergenic Reversal: Invalid Parameter y.");
    }
    swap(pi[i], pi[j]);
    int x_prime = bpi[i] - x, y_prime = bpi[j + 1] - y;
    bpi[j + 1] = x_prime + y_prime;
    bpi[i] = x + y;
    i++;
    j--;
    x = 0;
    y = bpi[j + 1];
  }
}

void Instance::transposition(int i, int j, int k, int x, int y, int z) {
  if (x < 0 || x > bpi[i]) {
    throw runtime_error("Intergenic Transposition: Invalid Parameter x.");
  }
  if (y < 0 || y > bpi[j]) {
    throw runtime_error("Intergenic Transposition: Invalid Parameter y.");
  }
  if (z < 0 || z > bpi[k]) {
    throw runtime_error("Intergenic Transposition: Invalid Parameter z.");
  }
  int x_prime = bpi[i] - x, y_prime = bpi[j] - y;
  // Transposition simulated using three reversals
  reversal(i, k - 1, x, z);
  reversal(i, i + k - j - 1, x, y_prime);
  reversal(i + k - j, k - 1, z, x_prime);
}

bool Instance::is_block(int a, int b) const {
  return abs(a - b) == 1 && abs(pi[a] - pi[b]) == 1 &&
         bpi[max(a, b)] == biota[max(pi[a], pi[b])];
}

void Instance::update() {
  breakpoints.clear();
  soft.clear();
  overcharged.clear();
  undercharged.clear();
  for (int i = 0; i < int(pi.size()) - 1; ++i) {
    if (!is_block(i, i + 1)) {
      breakpoints.push_back(i);
      if (abs(pi[i + 1] - pi[i]) != 1) {
        soft.push_back(i);
      } else if (bpi[i + 1] > target(i)) {
        overcharged.push_back(i);
      } else {
        undercharged.push_back(i);
      }
    }
  }
}

/* Two transpositions on the breakpoints i < j < k leaving x, y, and z
 * nucleotides in them. */
void Instance::two_transpositions(int i, int j, int k, int x, int y, int z) {
  pair<Params, Params> t =
      redistribution(bpi[i + 1], bpi[j + 1], bpi[k + 1], x, y, z);
  apply_transposition(i + 1, j + 1, k + 1, t.first[0], t.first[1],
                      t.first[2]);
  apply_transposition(i + 1, (i + 1) + (k + 1) - (j + 1), k + 1, t.second[0],
                      t.second[1], t.second[2]);
}

/* Call f(i, j) for the pairs of soft breakpoints (i, i + 1) and (j, j + 1)
 * that may be connected, stopping when f returns true. */
template <class F> bool Instance::soft_pairs(F f) {
  // Inverse permutation indicating the position of each element in pi
  vector<int> inverse(pi.size());
  for (size_t i = 0; i < pi.size(); ++i) {
    inverse[pi[i]] = i;
  }
  vector<bool> is_soft(size + 1, false);
  for (int i : soft) {
    is_soft[i] = true;
  }
  for (int i : soft) {
    // It is necessary to check the adjacencies of at most four elements
    for (int e : {pi[i] - 1, pi[i] + 1, pi[i + 1] - 1, pi[i + 1] + 1}) {
      if (e < 0 || e > size + 1) {
        continue;
      }
      // Check (j, j + 1) and (j - 1, j)
      int j = inverse[e];
      if (j <= size && is_soft[j] && f(i, j)) {
        return true;
      }
      j--;
      if (j >= 0 && is_soft[j] && f(i, j)) {
        return true;
      }
    }
  }
  return false;
}

// Lemma 3.4 - There are two or more overcharged breakpoints
bool Instance::case_i() {
  if (overcharged.size() < 2) {
    return false;
  }
  int i = overcharged[0], j = overcharged[1], k;
  // If there are more than two overcharged breakpoints, the third breakpoint
  // receives the excess
  if (overcharged.size() > 2) {
    k = overcharged[2];
    int sum = bpi[i + 1] + bpi[j + 1] + bpi[k + 1];
    int x = target(i), y = target(j);
    two_transpositions(i, j, k, x, y, sum - (x + y));
    return true;
  }
  // Otherwise, a soft (prioritized) or an undercharged breakpoint receives
  // the excess
  k = soft.empty() ? at(undercharged, 0) : soft[0];
  int x, y, z;
  if (k < i) {
    int old_i = i, old_j = j;
    i = k;
    j = old_i;
    k = old_j;
    int sum = bpi[i + 1] + bpi[j + 1] + bpi[k + 1];
    y = target(j);
    z = target(k);
    x = sum - (y + z);
  } else if (i < k && k < j) {
    swap(j, k);
    int sum = bpi[i + 1] + bpi[j + 1] + bpi[k + 1];
    x = target(i);
    z = target(k);
    y = sum - (x + z);
  } else {
    int sum = bpi[i + 1] + bpi[j + 1] + bpi[k + 1];
    x = target(i);
    y = target(j);
    z = sum - (x + y);
  }
  two_transpositions(i, j, k, x, y, z);
  return true;
}

// Lemma 3.5 - Case 1 (i, j) or (i + 1, j + 1)
bool Instance::soft_case_i(int i, int j) {
  if (i > j) {
    swap(i, j);
  }
  if (abs(pi[i] - pi[j]) == 1) {
    int goal = biota[max(pi[i], pi[j])];
    if (bpi[i + 1] + bpi[j + 1] >= goal) {
      int x, y;
      if (bpi[i + 1] >= goal) {
        x = goal;
        y = 0;
      } else {
        x = bpi[i + 1];
        y = goal - x;
      }
      apply_reversal(i + 1, j, x, y);
      return true;
    }
  }
  if (abs(pi[i + 1] - pi[j + 1]) == 1) {
    int goal = biota[max(pi[i + 1], pi[j + 1])];
    if (bpi[i + 1] + bpi[j + 1] >= goal) {
      int x, y;
      if (bpi[j + 1] >= goal) {
        x = bpi[i + 1];
        y = bpi[j + 1] - goal;
      } else {
        x = bpi[i + 1] - (goal - bpi[j + 1]);
        y = 0;
      }
      apply_reversal(i + 1, j, x, y);
      return true;
    }
  }
  return false;
}

// Lemma 3.5 - Case 2 (i + 1, j)
bool Instance::soft_case_ii(int i, int j) {
  if (i > j) {
    swap(i, j);
  }
  if (is_block(i + 1, j)) {
    return false;
  }
  if (abs(pi[i + 1] - pi[j]) == 1) {
    int goal = biota[max(pi[i + 1], pi[j])];
    if (bpi[i + 1] + bpi[j + 1] >= goal) {
      int k = j;
      j = at(breakpoints, index_of(breakpoints, i) + 1);
      int x, y = 0, z;
      if (bpi[i + 1] >= goal) {
        x = bpi[i + 1] - goal;
        z = 0;
      } else if (bpi[k + 1] >= goal) {
        x = bpi[i + 1];
        z = goal;
      } else {
        x = 0;
        z = goal - bpi[i + 1];
      }
      apply_transposition(i + 1, j + 1, k + 1, x, y, z);
      return true;
    }
  }
  return false;
}

// Lemma 3.5 - Case 3 (i, j + 1)
bool Instance::soft_case_iii(int i, int j) {
  if (i > j) {
    swap(i, j);
  }
  if (abs(pi[i] - pi[j + 1]) == 1) {
    int goal = biota[max(pi[i], pi[j + 1])];
    if (bpi[i + 1] + bpi[j + 1] >= goal) {
      int k, x, y, z;
      if (index_of(soft, i) != 0) {
        k = j;
        j = i;
        i = soft[0];
        if (bpi[j + 1] >= goal) {
          x = 0;
          y = goal;
          z = bpi[k + 1];
        } else if (bpi[k + 1] >= goal) {
          x = 0;
          y = 0;
          z = bpi[k + 1] - goal;
        } else {
          x = 0;
          y = goal - bpi[k + 1];
          z = 0;
        }
      } else {
        k = at(soft, index_of(soft, j) + 1);
        if (bpi[i + 1] >= goal) {
          x = goal;
          y = bpi[j + 1];
          z = 0;
        } else if (bpi[j + 1] >= goal) {
          x = 0;
          y = bpi[j + 1] - goal;
          z = 0;
        } else {
          x = goal - bpi[j + 1];
          y = 0;
          z = 0;
        }
      }
      apply_transposition(i + 1, j + 1, k + 1, x, y, z);
      return true;
    }
  }
  return false;
}

// Lemma 3.5 - There is a pair of soft connected breakpoints
bool Instance::case_ii() {
  return soft_pairs([this](int i, int j) {
    return soft_case_i(i, j) || soft_case_ii(i, j) || soft_case_iii(i, j);
  });
}

// Lemma 3.8 - There is one overcharged breakpoint and at least one
// undercharged breakpoint
bool Instance::case_iii() {
  if (overcharged.size() != 1 || undercharged.empty()) {
    return false;
  }
  // If the instance has only two breakpoints, a sequence of two reversals is
  // applied
  if (undercharged.size() + soft.size() == 1) {
    int i = min(overcharged[0], undercharged[0]);
    int j = max(overcharged[0], undercharged[0]);
    if (abs(i - j) != 1) {
      apply_reversal(i + 1, j, 0, 0);
      return soft_case_i(i, j);
    }
    if (contains(overcharged, i)) {
      apply_reversal(i + 1, j, target(i), 0);
    } else {
      apply_reversal(i + 1, j, bpi[i + 1], bpi[j + 1] - target(j));
    }
    return true;
  }
  // Otherwise, a sequence of two transpositions is applied, the soft
  // breakpoint is prioritized as the third breakpoint
  array<int, 3> b = {{overcharged[0], undercharged[0],
                      soft.empty() ? at(undercharged, 1) : soft[0]}};
  sort(b.begin(), b.end());
  int i = b[0], j = b[1], k = b[2], x, y, z;
  int sum = bpi[i + 1] + bpi[j + 1] + bpi[k + 1];
  if (!soft.empty()) {
    if (contains(overcharged, i)) {
      if (contains(undercharged, j)) {
        x = target(i);
        y = target(j);
        z = sum - (x + y);
      } else {
        x = target(i);
        z = target(k);
        y = sum - (x + z);
      }
    } else if (contains(overcharged, j)) {
      if (contains(undercharged, i)) {
        y = target(j);
        x = target(i);
        z = sum - (x + y);
      } else {
        y = target(j);
        z = target(k);
        x = sum - (y + z);
      }
    } else {
      if (contains(undercharged, i)) {
        z = target(k);
        x = target(i);
        y = sum - (x + z);
      } else {
        z = target(k);
        y = target(j);
        x = sum - (y + z);
      }
    }
  } else {
    if (contains(overcharged, i)) {
      x = target(i);
      y = target(j);
      z = sum - (x + y);
    } else if (contains(overcharged, j)) {
      y = target(j);
      x = target(i);
      z = sum - (x + y);
    } else {
      z = target(k);
      y = target(j);
      x = sum - (y + z);
    }
  }
  two_transpositions(i, j, k, x, y, z);
  return true;
}

// If this method applies a reversal, then the instance remains with no soft
// connected breakpoints.
bool Instance::classic_reversal(int i, int j) {
  if (i > j) {
    swap(i, j);
  }
  // Moving all the nucleotides from bpi[j + 1] to bpi[i + 1]
  if (abs(pi[i] - pi[j]) == 1) {
    apply_reversal(i + 1, j, bpi[i + 1], bpi[j + 1]);
    return true;
  }
  // Moving all the nucleotides from bpi[i + 1] to bpi[j + 1]
  if (abs(pi[i + 1] - pi[j + 1]) == 1) {
    apply_reversal(i + 1, j, 0, 0);
    return true;
  }
  return false;
}

// If this method applies a transposition, then a super hard breakpoint is
// created.
bool Instance::classic_transposition(int i, int j) {
  if (i > j) {
    swap(i, j);
  }
  // (i + 1, j) are adjacent in iota and not in pi
  if (abs((i + 1) - j) != 1 && abs(pi[i + 1] - pi[j]) == 1) {
    int k = j;
    j = at(soft, index_of(soft, i) + 1);
    apply_transposition(i + 1, j + 1, k + 1, 0, 0, 0);
    return true;
  }
  // (i, j + 1) are adjacent in iota and not in pi
  if (abs(pi[i] - pi[j + 1]) == 1) {
    int k;
    size_t idx = index_of(soft, i);
    if (idx != 0) {
      k = j;
      j = i;
      i = soft[idx - 1];
    } else {
      k = at(soft, index_of(soft, j) + 1);
    }
    apply_transposition(i + 1, j + 1, k + 1, 0, 0, 0);
    return true;
  }
  return false;
}

bool Instance::create_undercharged() {
  // If exists a decreasing strip, then a reversal is applied, otherwise a
  // transposition is applied
  return soft_pairs([this](int i, int j) { return classic_reversal(i, j); }) ||
         soft_pairs(
             [this](int i, int j) { return classic_transposition(i, j); });
}

// Lemma 3.10 - There is one overcharged breakpoint and no undercharged
// breakpoint
bool Instance::case_iv() {
  if (!create_undercharged()) {
    return false;
  }
  update();
  if (case_ii()) {
    update();
    return case_i() || case_ii() || case_iii();
  }
  return case_iii();
}

int Instance::run() {
  update();
  // The Python code has no limit, it is here only to avoid an infinite loop
  size_t max_iter = 100 * (pi.size() + 1);
  for (size_t iter = 0; !breakpoints.empty(); ++iter) {
    if (iter > max_iter) {
      throw runtime_error("reversal_transposition4: no progress");
    }
    if (!(case_i() || case_ii() || case_iii() || case_iv())) {
      throw runtime_error("Algorithm: None of the cases I, II, III, or IV was "
                          "performed.");
    }
    update();
  }
  return ops;
}

} // namespace

int ReversalTransposition4::estimate_distance(Permutation &pi, mt19937 &) {
  vector<int> genes, irs, irs_target;
  for (size_t i = 1; i <= pi.size(); ++i) {
    genes.push_back(abs(pi[i]));
  }
  for (size_t i = 1; i < pi.size(); ++i) {
    irs.push_back(pi.get_ir(i).first);
  }
  for (int i = 1; i < pi.target_size(); ++i) {
    irs_target.push_back(pi.get_ir_target(i).first);
  }
  Instance inst(genes, irs, irs_target);
  return inst.run();
}
