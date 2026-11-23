#include "reversal_transposition_indel6.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <climits>
#include <stdexcept>
#include <vector>

/* Adaptation to signed strings: a pair (A_i, A_i+1) is a breakpoint if
 * A_i+1 - A_i != 1 (so (+x, +(x+1)) and (-(x+1), -x) are adjacencies), which
 * keeps the lower bound of Theorem 1. Decreasing strips are replaced by
 * negative elements: if there is one, a reversal on two elements with
 * consecutive absolute values and different signs creates an adjacency
 * (Lemmas 6 and 7). Otherwise all the elements are positive and Lemma 9 is
 * used with transpositions, which do not change signs.
 *
 * Elements to be deleted (alphas) and runs of genes to be inserted are seen as
 * single elements. Only the lower bound of the origin intergenic regions is
 * used, while the target intergenic regions may be intervals. */

namespace {

const int INF = INT_MAX / 4;

bool inside(int v, IR w) { return w.first <= v && v <= w.second; }

class Instance {
public:
  explicit Instance(const Permutation &pi);
  int sort();

private:
  int n;       // label of the last cap
  int ops = 0; // number of operations applied
  vector<int> el;
  vector<bool> alpha;
  vector<int> ir;   // ir[p] is the region between el[p-1] and el[p]
  vector<IR> t;     // t[x] is the region between x-1 and x in the target
  vector<int> pos;  // positions of the labels (-1 if missing)

  static bool consecutive(int a, int b) { return b - a == 1; }
  /* Target region of two consecutive elements. */
  IR target(int a, int b) const { return t[max(abs(a), abs(b))]; }
  bool breakpoint(int p) const {
    return alpha[p - 1] || alpha[p] || !consecutive(el[p - 1], el[p]);
  }
  bool intergenic_breakpoint(int p) const {
    return breakpoint(p) || !inside(ir[p], target(el[p - 1], el[p]));
  }
  int intergenic_breakpoints() const;
  void update_pos();

  void reversal(int s, int e, int x, int y);
  void transposition(int i, int j, int k, int x, int y, int z);
  void set_region(int p, int w);
  void delete_alpha(int p);
  void insert_run(int x, int y);

  void remove_indels();
  bool lemma_4();
  bool oriented_reversal();
  void reversal_regions(int s, int e, bool adj_left, bool adj_right);
  void transposition_regions(int i, int j, int k, array<int, 3> order,
                             int r1_min = 0);
  void lemma_9();
};

Instance::Instance(const Permutation &pi) {
  n = pi.target_size() - 1;
  if (pi[1] != 0 || pi[pi.size()] != n || pi.occ_max() != 1) {
    throw invalid_argument("reversal_transposition_indel6 expects extended "
                           "genomes without replicas");
  }
  t.resize(n + 1);
  for (int x = 1; x <= n; ++x) {
    t[x] = pi.get_ir_target(x);
  }
  el.push_back(0);
  alpha.push_back(false);
  ir.push_back(0); // ir[0] is not used
  for (size_t i = 2; i <= pi.size(); ++i) {
    // Consecutive alphas are a single element (the region between them is
    // deleted with them)
    if (pi.check_alpha(i) && alpha.back()) {
      continue;
    }
    el.push_back(pi[i]);
    alpha.push_back(pi.check_alpha(i));
    ir.push_back(pi.get_ir(i - 1).first);
  }
}

int Instance::intergenic_breakpoints() const {
  int bi = 0;
  for (size_t p = 1; p < el.size(); ++p) {
    if (intergenic_breakpoint(p)) {
      bi++;
    }
  }
  return bi;
}

void Instance::update_pos() {
  pos.assign(n + 1, -1);
  for (size_t p = 0; p < el.size(); ++p) {
    if (!alpha[p]) {
      pos[abs(el[p])] = p;
    }
  }
}

/* Reversal of el[s..e] leaving x nucleotides of ir[s] and y of ir[e+1]
 * outside of the reversed segment. */
void Instance::reversal(int s, int e, int x, int y) {
  assert(1 <= s && s <= e && e + 1 < int(el.size()));
  if (x < 0 || x > ir[s] || y < 0 || y > ir[e + 1]) {
    throw runtime_error("reversal_transposition_indel6: invalid reversal");
  }
  int x_prime = ir[s] - x, y_prime = ir[e + 1] - y;
  reverse(el.begin() + s, el.begin() + e + 1);
  for (int k = s; k <= e; ++k) {
    el[k] = -el[k];
  }
  reverse(ir.begin() + s + 1, ir.begin() + e + 1);
  ir[s] = x + y;
  ir[e + 1] = x_prime + y_prime;
  ops++;
}

/* Transposition exchanging el[i..j-1] and el[j..k-1] (definition of the
 * paper). */
void Instance::transposition(int i, int j, int k, int x, int y, int z) {
  assert(1 <= i && i < j && j < k && k < int(el.size()));
  if (x < 0 || x > ir[i] || y < 0 || y > ir[j] || z < 0 || z > ir[k]) {
    throw runtime_error("reversal_transposition_indel6: invalid transposition");
  }
  vector<int> new_el(el.begin(), el.begin() + i);
  new_el.insert(new_el.end(), el.begin() + j, el.begin() + k);
  new_el.insert(new_el.end(), el.begin() + i, el.begin() + j);
  new_el.insert(new_el.end(), el.begin() + k, el.end());
  vector<int> new_ir(ir.begin(), ir.begin() + i);
  new_ir.push_back(x + ir[j] - y);
  new_ir.insert(new_ir.end(), ir.begin() + j + 1, ir.begin() + k);
  new_ir.push_back(z + ir[i] - x);
  new_ir.insert(new_ir.end(), ir.begin() + i + 1, ir.begin() + j);
  new_ir.push_back(y + ir[k] - z);
  new_ir.insert(new_ir.end(), ir.begin() + k + 1, ir.end());
  el.swap(new_el);
  ir.swap(new_ir);
  ops++;
}

/* Insertion or deletion of nucleotides in the region ir[p]. */
void Instance::set_region(int p, int w) {
  assert(w >= 0 && w != ir[p]);
  ir[p] = w;
  ops++;
}

/* Deletion of the alpha el[p] (Lemma 5), the remaining nucleotides make the
 * new pair balanced if it is an adjacency. */
void Instance::delete_alpha(int p) {
  int w = ir[p] + ir[p + 1];
  if (consecutive(el[p - 1], el[p + 1])) {
    w = min(w, target(el[p - 1], el[p + 1]).second);
  }
  el.erase(el.begin() + p);
  alpha.erase(alpha.begin() + p);
  ir[p] = w;
  ir.erase(ir.begin() + p + 1);
  ops++;
}

/* Insertion of the missing genes x..y (Lemma 5): after +(x-1) or before
 * -(x-1), with the target region between x-1 and x. If the other neighbor
 * is y+1 (with the same sign), its region is also fixed when possible. */
void Instance::insert_run(int x, int y) {
  int p = pos[x - 1];
  bool positive = el[p] >= 0;
  vector<int> run, inner;
  for (int k = x; k <= y; ++k) {
    int label = positive ? k : -(x + y - k);
    run.push_back(label);
    if (k > x) {
      inner.push_back(t[positive ? k : x + y - k + 1].first);
    }
  }
  int at = positive ? p + 1 : p; // insertion index in el
  int w = ir[at];
  IR near = t[x];
  bool far_adj = positive ? consecutive(y, el[p + 1])
                          : consecutive(el[p - 1], -y);
  int keep_near = min(w, far_adj ? near.second : near.first);
  int near_w = max(keep_near, near.first);
  int far_w = w - keep_near;
  if (far_adj) {
    far_w = max(far_w, t[y + 1].first);
  }

  vector<int> regions;
  regions.push_back(positive ? near_w : far_w);
  regions.insert(regions.end(), inner.begin(), inner.end());
  regions.push_back(positive ? far_w : near_w);
  el.insert(el.begin() + at, run.begin(), run.end());
  alpha.insert(alpha.begin() + at, run.size(), false);
  ir[at] = regions[0];
  ir.insert(ir.begin() + at + 1, regions.begin() + 1, regions.end());
  ops++;
}

/* Lines 2-4 of Algorithm 1: indels until the origin has the same genes as
 * the target. */
void Instance::remove_indels() {
  for (size_t p = 1; p < el.size(); ++p) {
    if (alpha[p]) {
      delete_alpha(p);
      p--;
    }
  }
  update_pos();
  vector<pair<int, int>> runs;
  for (int x = 1; x < n; ++x) {
    if (pos[x] < 0) {
      int y = x;
      while (pos[y + 1] < 0) {
        y++;
      }
      runs.push_back(make_pair(x, y));
      x = y;
    }
  }
  for (auto r : runs) {
    update_pos();
    insert_run(r.first, r.second);
  }
  update_pos();
}

/* Lemma 4: fix an overcharged or undercharged breakpoint with one indel. */
bool Instance::lemma_4() {
  for (size_t p = 1; p < el.size(); ++p) {
    if (breakpoint(p)) {
      continue;
    }
    IR w = target(el[p - 1], el[p]);
    if (!inside(ir[p], w)) {
      set_region(p, ir[p] > w.second ? w.second : w.first);
      return true;
    }
  }
  return false;
}

/* Reversal of el[s..e] creating the adjacencies (el[s-1], -el[e]) if adj_left
 * and (-el[s], el[e+1]) if adj_right, with their target regions when
 * possible. When both are created and there are not enough nucleotides, an
 * insertion is applied first (Lemma 7). */
void Instance::reversal_regions(int s, int e, bool adj_left, bool adj_right) {
  IR tl = adj_left ? target(el[s - 1], -el[e]) : IR(0, INF);
  IR tr = adj_right ? target(-el[s], el[e + 1]) : IR(0, INF);
  if (adj_left && adj_right && ir[s] + ir[e + 1] < tl.first + tr.first) {
    set_region(s, tl.first + tr.first - ir[e + 1]);
  }
  int l = ir[s] + ir[e + 1], left;
  if (adj_left && adj_right) {
    left = min(max(tl.first, l - tr.second), min(tl.second, l));
  } else if (adj_left) {
    left = min(l, tl.first);
  } else {
    left = l - min(l, tr.first);
  }
  int x = min(ir[s], left);
  reversal(s, e, x, left - x);
}

/* Lemmas 6 and 7 for signed strings: if there is a negative element, apply a
 * reversal on two elements with consecutive absolute values and different
 * signs that turns them adjacent. Reversals that create two adjacencies are
 * preferred, then the ones that leave a negative element. */
bool Instance::oriented_reversal() {
  vector<int> neg(el.size() + 1, 0); // prefix count of negative elements
  for (size_t p = 0; p < el.size(); ++p) {
    neg[p + 1] = neg[p] + (el[p] < 0 ? 1 : 0);
  }
  if (neg.back() == 0) {
    return false;
  }
  int best_score = -1, best_s = 0, best_e = 0;
  bool best_l = false, best_r = false;
  for (int v = 0; v < n; ++v) {
    int a = pos[v], b = pos[v + 1];
    if ((el[a] >= 0) == (el[b] >= 0)) {
      continue;
    }
    int i = min(a, b), j = max(a, b), s, e;
    if (consecutive(el[i], -el[j])) {
      s = i + 1;
      e = j;
    } else {
      s = i;
      e = j - 1;
    }
    bool adj_l = consecutive(el[s - 1], -el[e]);
    bool adj_r = consecutive(-el[s], el[e + 1]);
    assert(adj_l || adj_r);
    int neg_seg = neg[e + 1] - neg[s];
    int neg_after = neg.back() - neg_seg + (e - s + 1 - neg_seg);
    int score = (adj_l && adj_r) ? 2 : (neg_after > 0 ? 1 : 0);
    if (score > best_score) {
      best_score = score;
      best_s = s;
      best_e = e;
      best_l = adj_l;
      best_r = adj_r;
      if (score == 2) {
        break;
      }
    }
  }
  assert(best_score >= 0);
  reversal_regions(best_s, best_e, best_l, best_r);
  return true;
}

/* Transposition (i, j, k) choosing the regions of the new pairs r1 =
 * (el[i-1], el[j]), r2 = (el[k-1], el[i]) and r3 = (el[j-1], el[k]) so that
 * the ones that are adjacencies have their target regions, satisfying as many
 * as possible in the given order. r1 must also have at least r1_min. */
void Instance::transposition_regions(int i, int j, int k,
                                     array<int, 3> order, int r1_min) {
  array<IR, 3> w;
  array<pair<int, int>, 3> pairs = {{make_pair(el[i - 1], el[j]),
                                     make_pair(el[k - 1], el[i]),
                                     make_pair(el[j - 1], el[k])}};
  for (int m = 0; m < 3; ++m) {
    w[m] = consecutive(pairs[m].first, pairs[m].second)
               ? target(pairs[m].first, pairs[m].second)
               : IR(0, INF);
  }
  w[0].first = max(w[0].first, r1_min);
  int a = ir[i], b = ir[j], c = ir[k];

  for (int drop = 0; drop <= 3; ++drop) {
    array<bool, 3> act = {{false, false, false}};
    for (int m = 0; m < 3 - drop; ++m) {
      act[order[m]] = true;
    }
    for (int y = 0; y <= b; ++y) {
      // r1 = x + (b - y), r2 = z + (a - x), r3 = y + (c - z)
      int x_lo = 0, x_hi = a, z_lo = 0, z_hi = c;
      if (act[0]) {
        x_lo = max(x_lo, w[0].first - b + y);
        x_hi = min(x_hi, w[0].second - b + y);
      }
      if (act[2]) {
        z_lo = max(z_lo, y + c - w[2].second);
        z_hi = min(z_hi, y + c - w[2].first);
      }
      if (x_lo > x_hi || z_lo > z_hi) {
        continue;
      }
      int d_lo = z_lo - x_hi, d_hi = z_hi - x_lo; // d = z - x
      if (act[1]) {
        d_lo = max(d_lo, w[1].first - a);
        d_hi = min(d_hi, w[1].second - a);
      }
      if (d_lo > d_hi) {
        continue;
      }
      int x = max(x_lo, z_lo - d_lo);
      transposition(i, j, k, x, y, x + d_lo);
      return;
    }
  }
  throw runtime_error("reversal_transposition_indel6: no transposition found");
}

/* Lemma 9: with only positive elements and no overcharged or undercharged
 * breakpoints, remove one intergenic breakpoint with a transposition or two
 * with three operations. */
void Instance::lemma_9() {
  int i = 1;
  while (!intergenic_breakpoint(i)) {
    i++;
  }
  int j = i;
  while (!intergenic_breakpoint(j + 1)) {
    j++;
  }
  // Strip el[i..j], el[x] = el[j] + 1 and el[y] = el[j+1] - 1
  int x = pos[el[j] + 1];
  assert(x > j + 1);
  IR tx = t[el[x]];
  if (ir[j + 1] + ir[x] >= tx.first) {
    transposition_regions(i, j + 1, x, {{2, 0, 1}});
    return;
  }

  int u = el[j + 1], y = pos[u - 1];
  int need = max(0, t[u].first - ir[y + 1]);
  set_region(j + 1, tx.first - ir[x] + need);
  transposition_regions(i, j + 1, x, {{2, 0, 1}}, need);

  // Now el[i] = u, move its strip to the right of u - 1
  update_pos();
  if (consecutive(el[i - 1], u)) {
    return;
  }
  int q = i;
  while (q + 1 < int(el.size()) && consecutive(el[q], el[q + 1])) {
    q++;
  }
  int p = pos[u - 1];
  assert(p > q);
  transposition_regions(i, q + 1, p + 1, {{1, 0, 2}});
}

int Instance::sort() {
  update_pos();
  int missing = 0; // runs of genes to insert
  for (int x = 1; x < n; ++x) {
    if (pos[x] < 0 && pos[x - 1] >= 0) {
      missing++;
    }
  }
  int bound = intergenic_breakpoints() + missing; // Lemma 8

  remove_indels();
  assert(int(el.size()) == n + 1);

  size_t max_iter = 10 * el.size() + 100;
  for (size_t iter = 0; intergenic_breakpoints() > 0; ++iter) {
    if (iter > max_iter) {
      throw runtime_error("reversal_transposition_indel6: no progress");
    }
    if (!lemma_4() && !oriented_reversal()) {
      lemma_9();
    }
    update_pos();
  }

  for (int x = 0; x <= n; ++x) {
    assert(el[x] == x);
  }
  (void)bound;
  assert(ops <= 2 * bound);
  return ops;
}

} // namespace

int ReversalTranspositionIndel6::estimate_distance(Permutation &pi,
                                                   mt19937 &) {
  Instance inst(pi);
  return inst.sort();
}
