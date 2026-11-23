#include "fix_flexible.hpp"
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <vector>

/* The cycles are the ones of the flexible weighted cycle graph (breakpoint
 * graph of the signed permutation). Elements that are only in the origin
 * (alphas) are kept inside the black edges (the nucleotides between two
 * alphas are deleted with them) and the gray edges skip the elements that are
 * only in the target, so a gray edge (x, next(x)) has the regions between x and
 * x+1 and between next(x)-1 and next(x) (the regions between the inserted
 * genes do not depend on the origin and are fixed at their minimum).
 *
 * Free cycles are made balanced. If the unstable cycles need more (case i) or
 * less (case ii) nucleotides than they have, they get their minimum (maximum)
 * sizes, as the auxiliary cycles (the stable cycles with largest gap_min
 * (gap_max) that can provide (receive) the nucleotides), except for the last
 * auxiliary cycle that keeps the remaining nucleotides. Otherwise, the weight of
 * the unstable cycles is distributed among them. */

namespace {

struct FlexCycle {
  int wb = 0, wmin = 0, wmax = 0;
  vector<int> gray; // gray edges by their smaller label
  bool stable() const { return wmin <= wb && wb <= wmax; }
  int gap_min() const { return wb - wmin; }
  int gap_max() const { return wmax - wb; }
};

/* Sizes inside the intervals w with sum total (or the closest possible sum):
 * the minimums and then the remaining nucleotides up to the maximums. */
vector<int> distribute(const vector<IR> &w, int total) {
  vector<int> vals;
  for (auto i : w) {
    vals.push_back(i.first);
    total -= i.first;
  }
  for (size_t k = 0; k < w.size() && total > 0; ++k) {
    int add = min(total, w[k].second - w[k].first);
    vals[k] += add;
    total -= add;
  }
  return vals;
}

inline bool inside(int v, IR w) { return w.first <= v && v <= w.second; }

} // namespace

void fix_flexible_irs(Permutation &pi) {
  int n = pi.target_size() - 1;
  if (pi[1] != 0 || pi[pi.size()] != n || pi.occ_max() != 1) {
    throw invalid_argument("fixing flexible intergenic regions requires "
                           "extended genomes without replicas");
  }
  vector<IR> t(n + 1);
  for (int x = 1; x <= n; ++x) {
    t[x] = pi.get_ir_target(x);
  }

  // Reduced origin (without alphas) and weights of the black edges, wo[p] is
  // the weight between el[p-1] and el[p]
  vector<int> el = {0}, wo = {0};
  int pre = 0, post = 0;
  bool alpha = false;
  for (size_t i = 2; i <= pi.size(); ++i) {
    (alpha ? post : pre) += pi.get_ir(i - 1).first;
    if (pi.check_alpha(i)) {
      alpha = true;
      post = 0;
    } else {
      el.push_back(pi[i]);
      wo.push_back(pre + post);
      pre = post = 0;
      alpha = false;
    }
  }
  vector<int> pos(n + 1, -1), nxt(n + 1, -1), prv(n + 1, -1);
  for (size_t p = 0; p < el.size(); ++p) {
    pos[abs(el[p])] = p;
  }
  for (int x = 1, last = 0; x <= n; ++x) {
    if (pos[x] >= 0) {
      nxt[last] = x;
      prv[x] = last;
      last = x;
    }
  }
  auto gray_weight = [&](int x) {
    return (nxt[x] == x + 1) ? t[x + 1] : ir_sum(t[x + 1], t[nxt[x]]);
  };

  // Vertices are 2p (left) and 2p+1 (right) of el[p]
  auto vertex = [&](int label, bool head) {
    int p = pos[label];
    return 2 * p + ((head == (el[p] >= 0)) ? 1 : 0);
  };
  auto gray = [&](int v) {
    int x = abs(el[v / 2]);
    bool head = ((v % 2) == 1) == (el[v / 2] >= 0);
    return head ? make_pair(vertex(nxt[x], false), x)
                : make_pair(vertex(prv[x], true), prv[x]);
  };

  vector<FlexCycle> cycles;
  vector<bool> visited(el.size(), false);
  for (size_t o = 1; o < el.size(); ++o) {
    if (visited[o]) {
      continue;
    }
    FlexCycle c;
    int start = 2 * (o - 1) + 1, v = start;
    do {
      int e = (v % 2 == 1) ? v / 2 + 1 : v / 2; // black edge from v
      visited[e] = true;
      c.wb += wo[e];
      pair<int, int> g = gray((v % 2 == 1) ? v + 1 : v - 1);
      c.gray.push_back(g.second);
      c.wmin += gray_weight(g.second).first;
      c.wmax += gray_weight(g.second).second;
      v = g.first;
    } while (v != start);
    cycles.push_back(c);
  }

  // Size of each gray edge
  vector<int> val(n + 1, -1);
  auto assign = [&](const vector<int> &edges, int total) {
    vector<IR> w;
    for (int x : edges) {
      w.push_back(gray_weight(x));
    }
    vector<int> vals = distribute(w, total);
    for (size_t k = 0; k < edges.size(); ++k) {
      val[edges[k]] = vals[k];
    }
  };

  vector<FlexCycle *> unstable, stable;
  int wb_u = 0, wmin_u = 0, wmax_u = 0;
  for (auto &c : cycles) {
    if (c.stable()) {
      stable.push_back(&c);
    } else {
      unstable.push_back(&c);
      wb_u += c.wb;
      wmin_u += c.wmin;
      wmax_u += c.wmax;
    }
  }
  bool case_i = wb_u < wmin_u, case_ii = wb_u > wmax_u;
  size_t n_aux = 0;
  if (case_i || case_ii) {
    // Auxiliary cycles: the first ones in the order of the gaps
    sort(stable.begin(), stable.end(), [&](FlexCycle *a, FlexCycle *b) {
      return case_i ? a->gap_min() > b->gap_min() : a->gap_max() > b->gap_max();
    });
    int need = case_i ? wmin_u - wb_u : wb_u - wmax_u, gaps = 0;
    while (n_aux < stable.size() && gaps < need) {
      gaps += case_i ? stable[n_aux]->gap_min() : stable[n_aux]->gap_max();
      n_aux++;
    }
    for (auto c : unstable) {
      assign(c->gray, case_i ? c->wmin : c->wmax);
    }
    for (size_t k = 0; k < n_aux; ++k) {
      FlexCycle *c = stable[k];
      int excess = (k + 1 == n_aux) ? max(0, gaps - need) : 0;
      assign(c->gray, case_i ? c->wmin + excess : c->wmax - excess);
    }
  } else {
    vector<int> edges;
    for (auto c : unstable) {
      edges.insert(edges.end(), c->gray.begin(), c->gray.end());
    }
    assign(edges, wb_u);
  }
  for (size_t k = n_aux; k < stable.size(); ++k) {
    assign(stable[k]->gray, stable[k]->wb);
  }

  // Split the size of each gray edge into the regions of the target
  for (int x = 0; x < n; ++x) {
    if (pos[x] < 0) {
      continue;
    }
    int nx = nxt[x], s = val[x];
    assert(inside(s, gray_weight(x)));
    if (nx == x + 1) {
      pi.set_ir_target(x + 1, IR(s, s));
      continue;
    }
    int s1 = max(t[x + 1].first, min(t[x + 1].second, s - t[nx].first));
    pi.set_ir_target(x + 1, IR(s1, s1));
    pi.set_ir_target(nx, IR(s - s1, s - s1));
    for (int k = x + 2; k < nx; ++k) {
      pi.set_ir_target(k, IR(t[k].first, t[k].first));
    }
  }
  for (int x = 1; x <= n; ++x) {
    assert(inside(pi.get_ir_target(x).first, t[x]) &&
           pi.get_ir_target(x).first == pi.get_ir_target(x).second);
  }
}
