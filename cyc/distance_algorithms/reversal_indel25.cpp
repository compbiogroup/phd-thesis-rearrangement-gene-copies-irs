#include "reversal_indel25.hpp"
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <vector>

/* The instance is represented by the labeled intergenic breakpoint graph of
 * the paper. Elements of the origin that are not in the target (alphas) are
 * kept inside the origin edges and elements of the target that are not in the
 * origin are the labels of the target edges. Each iteration recomputes the
 * cycles of the graph and applies the operations of one of the lemmas 5-10.
 * Only the lower bound of the origin intergenic regions is used, while the
 * target intergenic regions may be intervals (a cycle is balanced if the
 * weight of its origin edges is inside the sum of the target intervals). */

namespace {

/* Origin edge (intergenic region between two consecutive elements of the
 * reduced origin string). If alpha is set, the region has elements to be
 * deleted, with pre nucleotides before them and post nucleotides after them
 * (the nucleotides between two alphas are deleted with them). Otherwise all
 * the nucleotides are in pre. */
struct Gap {
  int pre = 0;
  bool alpha = false;
  int post = 0;
  int weight() const { return pre + post; }
};

Gap concat(const Gap &g1, const Gap &g2) {
  Gap g;
  if (!g1.alpha && !g2.alpha) {
    g.pre = g1.pre + g2.pre;
  } else if (!g2.alpha) {
    g = g1;
    g.post += g2.pre;
  } else if (!g1.alpha) {
    g = g2;
    g.pre += g1.pre;
  } else {
    g.pre = g1.pre;
    g.alpha = true;
    g.post = g2.post;
  }
  return g;
}

Gap rev(const Gap &g) {
  Gap r = g;
  if (g.alpha) {
    swap(r.pre, r.post);
  }
  return r;
}

/* Cut g leaving left nucleotides in the left part. The alphas go to the left
 * part if alpha_left (requires left >= g.pre) or to the right part otherwise
 * (requires left <= g.pre). */
pair<Gap, Gap> split(const Gap &g, int left, bool alpha_left) {
  Gap l, r;
  assert(left >= 0 && left <= g.weight());
  if (!g.alpha) {
    l.pre = left;
    r.pre = g.pre - left;
  } else if (alpha_left) {
    assert(left >= g.pre);
    l.pre = g.pre;
    l.alpha = true;
    l.post = left - g.pre;
    r.pre = g.weight() - left;
  } else {
    assert(left <= g.pre);
    l.pre = left;
    r.pre = g.pre - left;
    r.alpha = true;
    r.post = g.post;
  }
  return make_pair(l, r);
}

/* Nucleotides of g that can go to its left (or right) part without alphas. */
int alpha_free(const Gap &g, bool left) {
  if (!g.alpha) {
    return g.weight();
  }
  return left ? g.pre : g.post;
}

IR intersect(IR i1, IR i2) {
  return IR(max(i1.first, i2.first), min(i1.second, i2.second));
}

bool empty(IR i) { return i.first > i.second; }

struct Cycle {
  vector<int> o;    // origin edges in traversal order (o[0] is the rightmost)
  vector<bool> conv; // whether each origin edge is traversed right to left
  vector<int> t;     // target edge (by its smaller label) after each origin edge
  vector<pair<int, int>> t_vtx; // vertices of each target edge (exit, entry)
  int weight = 0;               // sum of the weights of the origin edges
  IR target = ZERO_IR;          // sum of the weights of the target edges
  bool labeled = false;
  bool trivial() const { return o.size() == 1; }
  bool balanced() const {
    return target.first <= weight && weight <= target.second;
  }
  bool divergent() const {
    return find(conv.begin(), conv.end(), false) != conv.end();
  }
};

class Instance {
public:
  explicit Instance(const Permutation &pi);
  int sort();

private:
  int n;                 // label of the last cap
  int ops = 0;           // number of operations applied
  vector<int> el;        // reduced origin string (with caps 0 and n)
  vector<Gap> gaps;      // gaps[i] is the origin edge between el[i-1] and el[i]
  vector<IR> target;     // target[x] is the region between x-1 and x in the target
  vector<int> pos, nxt, prv;
  vector<Cycle> cycles;

  void update();
  /* Vertex ids are 2p for the left and 2p+1 for the right vertex of el[p]. */
  int vertex(int label, bool head) const;
  /* Vertex connected to v by a target edge and the target edge. */
  pair<int, int> partner(int v) const;
  IR target_weight(int x) const;
  bool target_labeled(int x) const { return nxt[x] != x + 1; }
  Cycle traverse(int o0) const;
  bool sorted() const;
  bool has_divergent() const;

  void reversal(int a, int b, int left_a, bool alpha_left_a, int left_b,
                bool alpha_left_b);
  void deletion(int i, int w);
  void insert_nucleotides(int i, int d);
  void insert_run(int i, int x);
  void split_reversal(const Cycle &c, size_t x, int s);
  int split_capacity(const Cycle &c, size_t x) const;
  bool try_reversal(int a, int b);

  bool lemma_5_6();
  bool lemma_7();
  bool lemma_8();
  bool lemma_9();
  bool lemma_10();
};

Instance::Instance(const Permutation &pi) {
  n = pi.target_size() - 1;
  if (pi[1] != 0 || pi[pi.size()] != n || pi.occ_max() != 1) {
    throw invalid_argument("reversal_indel25 expects extended genomes without "
                           "replicas");
  }
  target.resize(n + 1);
  for (int x = 1; x <= n; ++x) {
    target[x] = pi.get_ir_target(x);
  }
  el.push_back(0);
  gaps.push_back(Gap()); // gaps[0] is not used
  Gap cur;
  for (size_t i = 2; i <= pi.size(); ++i) {
    int ir = pi.get_ir(i - 1).first;
    if (cur.alpha) {
      cur.post += ir;
    } else {
      cur.pre += ir;
    }
    if (pi.check_alpha(i)) {
      cur.alpha = true;
      cur.post = 0;
    } else {
      el.push_back(pi[i]);
      gaps.push_back(cur);
      cur = Gap();
    }
  }
}

void Instance::update() {
  pos.assign(n + 1, -1);
  for (size_t p = 0; p < el.size(); ++p) {
    pos[abs(el[p])] = p;
  }
  nxt.assign(n + 1, -1);
  prv.assign(n + 1, -1);
  int last = 0;
  for (int x = 1; x <= n; ++x) {
    if (pos[x] >= 0) {
      nxt[last] = x;
      prv[x] = last;
      last = x;
    }
  }

  cycles.clear();
  vector<bool> visited(el.size(), false);
  for (size_t o = 1; o < el.size(); ++o) {
    if (visited[o]) {
      continue;
    }
    Cycle c = traverse(o);
    int o_max = *max_element(c.o.begin(), c.o.end());
    if (o_max != int(o)) {
      c = traverse(o_max);
    }
    for (size_t k = 0; k < c.o.size(); ++k) {
      visited[c.o[k]] = true;
      c.weight += gaps[c.o[k]].weight();
      c.target = ir_sum(c.target, target_weight(c.t[k]));
      c.labeled = c.labeled || gaps[c.o[k]].alpha || target_labeled(c.t[k]);
    }
    cycles.push_back(c);
  }
}

int Instance::vertex(int label, bool head) const {
  int p = pos[label];
  bool positive = el[p] >= 0;
  return 2 * p + ((head == positive) ? 1 : 0);
}

pair<int, int> Instance::partner(int v) const {
  int p = v / 2;
  int x = abs(el[p]);
  bool head = ((v % 2) == 1) == (el[p] >= 0);
  if (head) {
    return make_pair(vertex(nxt[x], false), x);
  }
  return make_pair(vertex(prv[x], true), prv[x]);
}

/* For a labeled target edge only the regions adjacent to the genes of the
 * origin count, the run of genes to insert is seen as a single element. */
IR Instance::target_weight(int x) const {
  IR w = target[x + 1];
  if (target_labeled(x)) {
    w = ir_sum(w, target[nxt[x]]);
  }
  return w;
}

Cycle Instance::traverse(int o0) const {
  Cycle c;
  int o = o0;
  bool conv = true;
  do {
    c.o.push_back(o);
    c.conv.push_back(conv);
    int exit = conv ? 2 * (o - 1) + 1 : 2 * o;
    pair<int, int> ve = partner(exit);
    c.t.push_back(ve.second);
    c.t_vtx.push_back(make_pair(exit, ve.first));
    int p = ve.first / 2;
    if (ve.first % 2 == 0) { // right extremity of origin edge p
      o = p;
      conv = true;
    } else { // left extremity of origin edge p + 1
      o = p + 1;
      conv = false;
    }
  } while (o != o0);
  assert(conv);
  return c;
}

bool Instance::sorted() const {
  for (auto &c : cycles) {
    if (!c.trivial() || c.labeled || !c.balanced()) {
      return false;
    }
  }
  return true;
}

bool Instance::has_divergent() const {
  for (auto &c : cycles) {
    if (c.divergent()) {
      return true;
    }
  }
  return false;
}

/* Reverse el[a..b-1], leaving left_a nucleotides of gap a and left_b of gap b
 * out of the reversed segment. */
void Instance::reversal(int a, int b, int left_a, bool alpha_left_a,
                        int left_b, bool alpha_left_b) {
  assert(0 < a && a < b && b < int(el.size()));
  pair<Gap, Gap> sa = split(gaps[a], left_a, alpha_left_a);
  pair<Gap, Gap> sb = split(gaps[b], left_b, alpha_left_b);
  reverse(el.begin() + a, el.begin() + b);
  for (int k = a; k < b; ++k) {
    el[k] = -el[k];
  }
  reverse(gaps.begin() + a + 1, gaps.begin() + b);
  for (int k = a + 1; k < b; ++k) {
    gaps[k] = rev(gaps[k]);
  }
  gaps[a] = concat(sa.first, rev(sb.first));
  gaps[b] = concat(rev(sa.second), sb.second);
  ops++;
}

/* Delete the alphas of gap i (if any) and nucleotides, leaving w. */
void Instance::deletion(int i, int w) {
  assert(w <= gaps[i].weight() && (gaps[i].alpha || w < gaps[i].weight()));
  gaps[i] = Gap();
  gaps[i].pre = w;
  ops++;
}

void Instance::insert_nucleotides(int i, int d) {
  assert(d > 0);
  gaps[i].pre += d;
  ops++;
}

/* Insert the genes of the label of target edge x in the clean origin edge i,
 * which forms a non-negative trivial cycle with x, making the new trivial
 * cycles clean and balanced (Lemma 5). */
void Instance::insert_run(int i, int x) {
  int nx = nxt[x];
  bool positive = el[i - 1] >= 0;
  IR first = positive ? target[x + 1] : target[nx];
  IR last = positive ? target[nx] : target[x + 1];
  int w = gaps[i].weight();
  assert(!gaps[i].alpha && w <= first.second + last.second);

  int p = min(w, first.second);
  vector<int> run;
  vector<Gap> new_gaps(1);
  new_gaps[0].pre = max(p, first.first);
  for (int k = 1; k < nx - x; ++k) {
    int label = positive ? x + k : nx - k;
    run.push_back(positive ? label : -label);
    if (k > 1) {
      Gap g;
      g.pre = target[positive ? label : label + 1].first;
      new_gaps.push_back(g);
    }
  }
  Gap g;
  g.pre = max(w - p, last.first);
  new_gaps.push_back(g);

  el.insert(el.begin() + i, run.begin(), run.end());
  gaps[i] = new_gaps[0];
  gaps.insert(gaps.begin() + i + 1, new_gaps.begin() + 1, new_gaps.end());
  ops++;
}

/* Maximum weight that the reversal on the divergent pair (o[x], o[x+1]) of c
 * can leave in the trivial cycle it creates without alphas. */
int Instance::split_capacity(const Cycle &c, size_t x) const {
  int a = min(c.o[x], c.o[(x + 1) % c.o.size()]);
  int b = max(c.o[x], c.o[(x + 1) % c.o.size()]);
  bool trivial_a = c.t_vtx[x].first % 2 == 1;
  return alpha_free(gaps[a], trivial_a) + alpha_free(gaps[b], trivial_a);
}

/* Apply the reversal on the divergent pair (o[x], o[x+1]) of c, which breaks
 * c into a trivial cycle (with the target edge t[x]) of weight s and without
 * alphas, and another cycle with the remaining weight and alphas. */
void Instance::split_reversal(const Cycle &c, size_t x, int s) {
  int a = min(c.o[x], c.o[(x + 1) % c.o.size()]);
  int b = max(c.o[x], c.o[(x + 1) % c.o.size()]);
  assert(c.t_vtx[x].first % 2 == c.t_vtx[x].second % 2);
  // The target edge connects the right vertices of the origin edges (and the
  // trivial cycle uses the new edge a) or their left vertices (new edge b).
  bool trivial_a = c.t_vtx[x].first % 2 == 1;
  int cap_a = alpha_free(gaps[a], trivial_a);
  int w_a = gaps[a].weight(), w_b = gaps[b].weight();
  assert(s <= split_capacity(c, x));
  int s_a = min(cap_a, s);
  if (trivial_a) {
    reversal(a, b, s_a, false, s - s_a, false);
  } else {
    reversal(a, b, w_a - s_a, true, w_b - (s - s_a), true);
  }
}

/* Reversal on origin edges a and b that keeps their weights, it is applied
 * only if it creates a divergent cycle. */
bool Instance::try_reversal(int a, int b) {
  if (a == b) {
    return false;
  }
  if (a > b) {
    swap(a, b);
  }
  Instance copy(*this);
  copy.reversal(a, b, copy.gaps[a].weight(), true, 0, false);
  copy.update();
  if (copy.has_divergent()) {
    *this = std::move(copy);
    return true;
  }
  return false;
}

/* Turn an unbalanced or labeled trivial cycle into balanced and clean cycles
 * (one operation by Lemma 5 or two by Lemma 6). */
bool Instance::lemma_5_6() {
  for (auto &c : cycles) {
    if (!c.trivial() || (!c.labeled && c.balanced())) {
      continue;
    }
    int i = c.o[0], x = c.t[0];
    if (gaps[i].alpha || c.weight > c.target.second) {
      // Lemma 6: the deletion turns c into a non-negative cycle without
      // alphas (balanced if possible), Lemma 5 is applied next if needed
      deletion(i, min(c.weight, c.target.second));
    } else if (!target_labeled(x)) {
      insert_nucleotides(i, c.target.first - c.weight);
    } else {
      insert_run(i, x);
    }
    return true;
  }
  return false;
}

/* Break a labeled divergent cycle creating a non-negative trivial cycle with
 * a clean origin edge (Lemma 5 is applied next if needed). */
bool Instance::lemma_7() {
  for (auto &c : cycles) {
    if (!c.labeled || !c.divergent()) {
      continue;
    }
    size_t x = 0;
    while (c.conv[x] == c.conv[x + 1]) {
      x++;
    }
    int s = min(target_weight(c.t[x]).second, split_capacity(c, x));
    split_reversal(c, x, s);
    return true;
  }
  return false;
}

/* Turn a clean divergent cycle into a balanced cycle (if it is positive) or
 * break it increasing the number of balanced cycles. */
bool Instance::lemma_8() {
  for (auto &c : cycles) {
    if (c.labeled || !c.divergent()) {
      continue;
    }
    if (c.weight < c.target.first) {
      insert_nucleotides(c.o[0], c.target.first - c.weight);
      return true;
    }
    // Weight s of the trivial cycle for the best divergent pair, two balanced
    // cycles are preferred over one.
    size_t k = c.o.size(), best_x = k;
    int best_score = -1, best_s = 0;
    for (size_t x = 0; x < k; ++x) {
      if (c.conv[x] == c.conv[(x + 1) % k]) {
        continue;
      }
      IR tw = target_weight(c.t[x]);
      IR range(0, split_capacity(c, x));
      IR triv = intersect(tw, range);
      IR rest = intersect(IR(c.weight - c.target.second + tw.second,
                             c.weight - c.target.first + tw.first),
                          range);
      IR both = intersect(triv, rest);
      int score = 0, s = min(tw.second, range.second);
      if (!empty(both)) {
        score = 2;
        s = both.first;
      } else if (!empty(triv)) {
        score = 1;
        s = triv.first;
      } else if (!empty(rest)) {
        score = 1;
        s = rest.first;
      }
      if (score > best_score) {
        best_score = score;
        best_x = x;
        best_s = s;
      }
    }
    assert(best_x < k);
    split_reversal(c, best_x, best_s);
    return true;
  }
  return false;
}

/* Turn an oriented cycle into a divergent cycle. */
bool Instance::lemma_9() {
  for (auto c : cycles) {
    size_t k = c.o.size();
    if (k < 3) {
      continue;
    }
    for (size_t i = 0; i < k; ++i) {
      int o0 = c.o[i], o1 = c.o[(i + 1) % k], o2 = c.o[(i + 2) % k];
      if (min(o0, o1) <= o2 && o2 <= max(o0, o1) && try_reversal(o0, o2)) {
        return true;
      }
    }
  }
  return false;
}

/* With only non-oriented cycles, apply a reversal on the origin edges
 * adjacent to a target edge that closes an open gate of a cycle C, turning C
 * into a divergent cycle. */
bool Instance::lemma_10() {
  // Target edges of each cycle by the indices of their adjacent origin edges
  vector<vector<IR>> te(cycles.size());
  for (size_t ci = 0; ci < cycles.size(); ++ci) {
    const Cycle &c = cycles[ci];
    for (size_t k = 0; k < c.o.size(); ++k) {
      int x = c.o[k], y = c.o[(k + 1) % c.o.size()];
      te[ci].push_back(IR(min(x, y), max(x, y)));
    }
  }
  auto cross = [](IR e, IR f) {
    return (e.first < f.first && f.first <= e.second && e.second < f.second) ||
           (f.first < e.first && e.first <= f.second && f.second < e.second);
  };

  for (size_t ci = 0; ci < cycles.size(); ++ci) {
    if (cycles[ci].trivial()) {
      continue;
    }
    for (auto gate : te[ci]) {
      bool open = true;
      for (auto e : te[ci]) {
        open = open && !cross(gate, e);
      }
      if (!open) {
        continue;
      }
      for (size_t di = 0; di < cycles.size(); ++di) {
        if (di == ci) {
          continue;
        }
        for (auto e : te[di]) {
          if (cross(gate, e) && try_reversal(e.first, e.second)) {
            return true;
          }
        }
      }
    }
  }

  // Not expected, but try every target edge before giving up
  for (size_t ci = 0; ci < cycles.size(); ++ci) {
    for (auto e : te[ci]) {
      if (try_reversal(e.first, e.second)) {
        return true;
      }
    }
  }
  return false;
}

int Instance::sort() {
  update();
  int clean_balanced = 0;
  for (auto &c : cycles) {
    if (!c.labeled && c.balanced()) {
      clean_balanced++;
    }
  }
  int lower_bound = int(el.size()) - 1 - clean_balanced; // Lemma 4
  size_t max_iter = 10 * (el.size() + n) + 100;

  for (size_t iter = 0; !sorted(); ++iter) {
    if (iter > max_iter) {
      throw runtime_error("reversal_indel25: no progress");
    }
    if (!(lemma_5_6() || lemma_7() || lemma_8() || lemma_9() || lemma_10())) {
      throw runtime_error("reversal_indel25: no case found");
    }
    update();
  }

  for (int x = 0; x <= n; ++x) {
    assert(el[x] == x);
  }
  (void)lower_bound;
  assert(ops <= 2.5 * lower_bound);
  return ops;
}

} // namespace

int ReversalIndel25::estimate_distance(Permutation &pi, mt19937 &) {
  Instance inst(pi);
  return inst.sort();
}
