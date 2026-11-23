// #define DEBUG_CYCLEGRAPH
#include "cycles.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <regex>
#include <sstream>
#include <vector>

#include "../misc/genome.hpp"
#include "../misc/permutation.hpp"

CycleGraph::CycleGraph(const MultiGenome &origin, const MultiGenome &target,
                       mt19937 &rand_gen_init, bool unsign,
                       bool only_trivial_init)
    : only_trivial(only_trivial_init), vertices(), cycles(), indel_count(),
      rand_gen(rand_gen_init()) {
  int a = max(origin.full_size() + 2 * origin.n_chrom(),
              target.full_size() + 2 * target.n_chrom());
  int b = max(origin.get_op_max(), target.get_op_max());
  op_max = max(a, b) + 1;
  int alt_cap, cap = op_max;
  op_max++;
  int v_id = 0;
  size_t beg, end;

  /* If both genomes have only one chromosome we use distinct values for the
   * caps */
  if (origin.n_chrom() == 1 && target.n_chrom() == 1) {
    alt_cap = op_max;
    op_max++;
  } else {
    alt_cap = cap;
  }

  /* Count the number of linear chromosomes in the origin and target genomes */
  int n_lin_origin = origin.num_lin();
  int n_lin_target = target.num_lin();

  /* Using gene values to label the vertices, also indicate how to update indels
   * of each vertice. Furthermore, include the black and indel edges. */
  // vertices and edges of origin and target chromosomes
  int new_chroms = n_lin_target - n_lin_origin;
  for (bool positive : {false, true}) {
    const MultiGenome &genome = (!positive) ? origin : target;
    for (size_t i = 1; i <= genome.n_chrom(); ++i) {
      if (genome[i].is_linear()) {
        vertices.push_back(
            Vertex(cap, true, Cap, v_id + 1, genome[i].get_ir(1), positive));
        vertices[v_id].indel_update = (positive) ? 1 : -1;
        indel_count[cap] -= vertices[v_id].indel_update;
        v_id++;
        // The first gene was replaced by the cap, star from 2 and end in size -
        // 1
        beg = 2;
        end = genome[i].size() - 1;
      } else {
        beg = 1;
        end = genome[i].size();
      }
      int fst_v_id = v_id;
      for (size_t j = beg; j <= end; ++j) {
        if (genome[i][j] >= 0) {
          vertices.push_back(Vertex(abs(genome[i][j]), true, Tail, v_id - 1,
                                    genome[i].safe_get_ir(j - 1), positive));
          vertices.push_back(Vertex(abs(genome[i][j]), true, Head, v_id + 2,
                                    genome[i].get_ir(j), positive));
        } else {
          vertices.push_back(Vertex(abs(genome[i][j]), false, Head, v_id - 1,
                                    genome[i].safe_get_ir(j - 1), positive));
          vertices.push_back(Vertex(abs(genome[i][j]), false, Tail, v_id + 2,
                                    genome[i].get_ir(j), positive));
        }
        vertices[v_id].indel_update = (positive) ? 1 : -1;
        vertices[v_id].red = v_id + 1;
        vertices[v_id + 1].indel_update = vertices[v_id].indel_update;
        vertices[v_id + 1].red = v_id;
        indel_count[abs(genome[i][j])] -= vertices[v_id].indel_update;
        v_id += 2;
      }
      if (genome[i].is_linear()) {
        vertices.push_back(Vertex(alt_cap, true, Cap, v_id - 1,
                                  genome[i].get_ir(genome[i].size() - 1),
                                  positive));
        vertices[v_id].indel_update = (positive) ? 1 : -1;
        indel_count[alt_cap] -= vertices[v_id].indel_update;
        v_id++;
      } else {
        /* In circular chromosomes the first and last vertices are connected by
         * a black edge */
        vertices[v_id - 1].black = fst_v_id;
        vertices[fst_v_id].black = v_id - 1;
        vertices[v_id - 1].weight =
            add_sign_to_ir(genome[i].get_ir(genome[i].size()), positive);
        vertices[fst_v_id].weight = vertices[v_id - 1].weight;
      }
    }

    // Add empty chromosomes
    for (int i = 0; i < new_chroms; ++i) {
      vertices.push_back(Vertex(cap, true, Cap, v_id + 1, ZERO_IR, positive));
      vertices[v_id].indel_update = (positive) ? 1 : -1;
      indel_count[cap] -= vertices[v_id].indel_update;
      v_id++;
      vertices.push_back(
          Vertex(alt_cap, true, Cap, v_id - 1, ZERO_IR, positive));
      vertices[v_id].indel_update = (positive) ? 1 : -1;
      indel_count[alt_cap] -= vertices[v_id].indel_update;
      v_id++;
    }

    /* Save index of last vertex of origin and target genomes */
    if (!positive) {
      last_origin_vertex = v_id - 1;
    } else {
      last_target_vertex = v_id - 1;
    }

    new_chroms = n_lin_origin - n_lin_target;
  }

  /* Gray Edges */
  for (int i = 0; i <= last_origin_vertex; i++) {
    for (int j = last_origin_vertex + 1; j <= last_target_vertex; j++) {
      if (vertices[i].gene_val == vertices[j].gene_val &&
          (unsign || vertices[i].vtx_type == vertices[j].vtx_type)) {
        vertices[i].grays.push_back(j);
        vertices[j].grays.push_back(i);
      }
    }
  }
}

Vtx_id CycleGraph::partner(Vtx_id x, const map<Vtx_id, Vtx_id> &fix,
                           const map<Vtx_id, Vtx_id> *base) const {
  auto it = fix.find(x);
  if (it != fix.end()) {
    return it->second;
  }
  if (base != nullptr && (it = base->find(x)) != base->end()) {
    return it->second;
  }
  return vertices[x].is_indel ? vertices[x].red : vertices[x].fix_gray;
}

bool CycleGraph::gray_closure(Vtx_id v, Vtx_id u, map<Vtx_id, Vtx_id> &fix,
                              const map<Vtx_id, Vtx_id> *base) const {
  vector<pair<Vtx_id, Vtx_id>> work = {make_pair(v, u)};
  while (!work.empty()) {
    Vtx_id a = work.back().first, b = work.back().second;
    work.pop_back();
    Vtx_id pa = partner(a, fix, base), pb = partner(b, fix, base);
    if (pa == b && pb == a) {
      continue;
    }
    const vector<Vtx_id> &grays = vertices[a].grays;
    if (pa != NO_EDGE || pb != NO_EDGE ||
        find(grays.begin(), grays.end(), b) == grays.end()) {
      return false;
    }
    fix[a] = b;
    fix[b] = a;
    if (vertices[a].vtx_type != Cap) {
      work.push_back(make_pair(vertices[a].red, vertices[b].red));
    }
    if (vertices[a].fix_adj || vertices[b].fix_adj) {
      if (!ir_intersect(ir_sum(vertices[a].weight, vertices[b].weight),
                        ZERO_IR)) {
        return false;
      }
      work.push_back(make_pair(vertices[a].black, vertices[b].black));
    }
  }
  return true;
}

bool CycleGraph::indel_closure(Vtx_id v, const map<Vtx_id, Vtx_id> &fix,
                               vector<Vtx_id> &out) const {
  set<Vtx_id> seen;
  vector<Vtx_id> work = {v};
  while (!work.empty()) {
    Vtx_id x = work.back();
    work.pop_back();
    if (!seen.insert(x).second) {
      continue;
    }
    Vtx_id p = partner(x, fix);
    if (vertices[x].vtx_type == Cap ||
        (p != NO_EDGE && p != vertices[x].red)) {
      return false;
    }
    out.push_back(x);
    work.push_back(vertices[x].red);
    if (vertices[x].fix_adj) {
      work.push_back(vertices[x].black);
    }
  }
  return true;
}

void CycleGraph::commit_gray(Vtx_id v, Vtx_id u) {
  map<Vtx_id, Vtx_id> fix;
  if (!gray_closure(v, u, fix)) {
    throw runtime_error("gray edge (" + to_string(v) + ", " + to_string(u) +
                        ") conflicts with the fixed adjacencies");
  }
  for (auto &e : fix) {
    vertices[e.first].fix_gray = e.second;
  }
}

void CycleGraph::commit_indel(Vtx_id v, bool fixed) {
  vector<Vtx_id> out;
  if (!indel_closure(v, map<Vtx_id, Vtx_id>(), out)) {
    throw runtime_error("indel of vertex " + to_string(v) +
                        " conflicts with the fixed adjacencies");
  }
  for (Vtx_id x : out) {
    if (!vertices[x].is_indel) {
      vertices[x].is_indel = true;
      vertices[vertices[x].red].is_indel = true;
      indel_count[vertices[x].gene_val] += vertices[x].indel_update;
    }
    vertices[x].is_indel_fix = vertices[x].is_indel_fix || fixed;
  }
}

void CycleGraph::fix_partition(const PartitionFix &fix, bool fix_chars) {
  // First vertex of the chromosomes of each genome (the left cap of a linear
  // chromosome, the first vertex of gene 1 of a circular one), number of genes
  // in each and whether it is circular (its last vertex is the black neighbor
  // of the first one)
  auto chromosomes = [&](int first, int last, vector<Vtx_id> &bases,
                         vector<int> &sizes, vector<bool> &circular) {
    for (int v = first; v <= last;) {
      int w;
      if (vertices[v].vtx_type == Cap) {
        w = v + 1;
        while (vertices[w].vtx_type != Cap) {
          w += 2;
        }
        sizes.push_back((w - v - 1) / 2);
      } else {
        w = vertices[v].black;
        sizes.push_back((w - v + 1) / 2);
      }
      bases.push_back(v);
      circular.push_back(vertices[v].vtx_type != Cap);
      v = w + 1;
    }
  };
  vector<Vtx_id> base[2];
  vector<int> n_genes[2];
  vector<bool> circular[2];
  chromosomes(0, last_origin_vertex, base[0], n_genes[0], circular[0]);
  chromosomes(last_origin_vertex + 1, last_target_vertex, base[1], n_genes[1],
              circular[1]);

  // Fixed adjacencies (the regions that are not breakpoints) and the first
  // vertex of each gene (genes counted over all chromosomes). Region r is the
  // one after gene r: 0..n in a linear chromosome (0 and n are the regions of
  // the caps), 1..n in a circular one (n is between gene n and gene 1).
  vector<Vtx_id> gene_vtx[2];
  for (int side = 0; side < 2; ++side) {
    const vector<vector<int>> &breaks = (side == 0) ? fix.breaks_g : fix.breaks_h;
    for (size_t k = 0; k < base[side].size(); ++k) {
      if (k >= breaks.size()) {
        if (n_genes[side][k] > 0) {
          throw invalid_argument("FIX: fewer chromosomes than the genomes");
        }
        continue; // empty chromosome added by the graph
      }
      int n = n_genes[side][k];
      bool circ = circular[side][k];
      Vtx_id b0 = base[side][k];
      set<int> bps(breaks[k].begin(), breaks[k].end());
      for (int r : bps) {
        if (r < (circ ? 1 : 0) || r > n) {
          throw invalid_argument(
              "FIX: breakpoint " + to_string(r) + " out of the chromosome" +
              (circ ? "" : " (use -e)"));
        }
      }
      for (int r = circ ? 1 : 0; r <= n; ++r) {
        // vertex after gene r (or the left cap) and its black neighbor
        Vtx_id a = circ ? b0 + 2 * r - 1 : b0 + 2 * r;
        Vtx_id b = (circ && r == n) ? b0 : a + 1;
        assert(vertices[a].black == b);
        if (bps.count(r) == 0) {
          vertices[a].fix_adj = vertices[b].fix_adj = true;
        }
      }
      for (int j = 1; j <= n; ++j) {
        gene_vtx[side].push_back(circ ? b0 + 2 * (j - 1) : b0 + 2 * j - 1);
      }
    }
    if (breaks.size() > base[side].size()) {
      throw invalid_argument("FIX: more chromosomes than the genomes");
    }
  }
  // With a single linear chromosome in each genome the caps are distinct (the
  // left cap of one genome is only paired with the left cap of the other), but
  // the partition has no orientation and may pair the block of the left cap of
  // the origin with the block of the right cap of the target in reverse (or the
  // block of the right cap of the origin with the one of the left cap of the
  // target). These cap adjacencies cannot be kept in the graph, so they become
  // breakpoints (otherwise the blocks are left without a pair).
  bool single_linear = base[0].size() == 1 && base[1].size() == 1 &&
                       !circular[0][0] && !circular[1][0];
  int ng = single_linear ? n_genes[0][0] : 0;
  int nh = single_linear ? n_genes[1][0] : 0;
  if (ng > 0 && int(fix.corr.size()) == ng) {
    auto unfix = [&](Vtx_id cap) {
      vertices[cap].fix_adj = vertices[vertices[cap].black].fix_adj = false;
    };
    if (fix.corr[0] == -nh) {
      unfix(base[0][0]);
      unfix(last_target_vertex);
    }
    if (fix.corr[ng - 1] == -1) {
      unfix(last_origin_vertex);
      unfix(base[1][0]);
    }
  }
  use_fix = true;
  if (!fix_chars) {
    return;
  }

  // Fixed characters
  if (fix.corr.size() != gene_vtx[0].size()) {
    throw invalid_argument("FIX: " + to_string(fix.corr.size()) +
                           " genes in the correspondence, the origin has " +
                           to_string(gene_vtx[0].size()) + " (use -e)");
  }
  auto fix_edge = [&](Vtx_id x, Vtx_id y) {
    const vector<Vtx_id> &grays = vertices[x].grays;
    if (find(grays.begin(), grays.end(), y) == grays.end()) {
      throw invalid_argument("FIX: vertices " + to_string(x) + " and " +
                             to_string(y) + " cannot be paired");
    }
    vertices[x].fix_gray = y;
    vertices[y].fix_gray = x;
    vertices[x].fix_char = vertices[y].fix_char = true;
  };
  auto fix_indel = [&](Vtx_id x) {
    for (Vtx_id y : {x, vertices[x].red}) {
      vertices[y].is_indel = vertices[y].is_indel_fix = true;
      vertices[y].fix_char = true;
    }
    indel_count[vertices[x].gene_val] += vertices[x].indel_update;
  };
  vector<bool> paired(gene_vtx[1].size(), false);
  for (size_t j = 0; j < fix.corr.size(); ++j) {
    Vtx_id g = gene_vtx[0][j];
    int p = abs(fix.corr[j]);
    if (p == 0) {
      fix_indel(g);
      continue;
    }
    if (p > int(gene_vtx[1].size()) || paired[p - 1]) {
      throw invalid_argument("FIX: invalid position " + to_string(p) +
                             " in the correspondence");
    }
    paired[p - 1] = true;
    Vtx_id h = gene_vtx[1][p - 1];
    for (Vtx_id x : {g, g + 1}) {
      fix_edge(x, (vertices[h].vtx_type == vertices[x].vtx_type) ? h : h + 1);
    }
  }
  for (size_t p = 0; p < paired.size(); ++p) {
    if (!paired[p]) {
      fix_indel(gene_vtx[1][p]);
    }
  }
  // With a single linear chromosome in each genome the caps are paired
  if (single_linear) {
    fix_edge(base[0][0], base[1][0]);
    fix_edge(last_origin_vertex, last_target_vertex);
  }

  // The correspondence must keep the fixed adjacencies
  for (Vtx_id x = 0; x <= last_target_vertex; ++x) {
    Vtx_id y = vertices[x].black;
    if (!vertices[x].fix_adj) {
      continue;
    }
    bool ok = vertices[x].is_indel == vertices[y].is_indel;
    Vtx_id px = vertices[x].fix_gray, py = vertices[y].fix_gray;
    if (ok && !vertices[x].is_indel && px != NO_EDGE && py != NO_EDGE) {
      ok = vertices[px].black == py &&
           ir_intersect(ir_sum(vertices[x].weight, vertices[px].weight),
                        ZERO_IR);
    }
    if (!ok) {
      throw invalid_argument("FIX: the correspondence breaks the fixed "
                             "adjacency of vertices " +
                             to_string(x) + " and " + to_string(y));
    }
  }
}

vector<Vtx_id> CycleGraph::produce_random_idxs(bool is_random) {
  vector<Vtx_id> idxs(size());
  iota(idxs.begin(), idxs.end(), 0);
  if (is_random) {
    shuffle(idxs.begin(), idxs.end(), rand_gen);
  }
  return idxs;
}

/* Whether the gene of vertex v can still become an indel: neither extremity is
 * an indel, is in a cycle or has a fixed gray edge (the graph may be partially
 * decomposed, as in the crossover and the mutation of ga) */
bool CycleGraph::gene_free(Vtx_id v) const {
  Vtx_id r = vertices[v].red;
  return !vertices[v].is_indel && !vertices[v].in_cycle &&
         vertices[v].fix_gray == NO_EDGE && !vertices[r].is_indel &&
         !vertices[r].in_cycle && vertices[r].fix_gray == NO_EDGE;
}

// Greed heuristic to mark indel edges, always select the largest sequence of
// indels
void CycleGraph::mark_indel_edges(bool is_random) {
  vector<Vtx_id> idxs = produce_random_idxs(is_random);
  /* genes whose block cannot be an indel (with fixed adjacencies); this does
   * not change while indels are marked, since it depends only on the caps and
   * on the fixed gray edges */
  vector<bool> blocked(size(), false);

  int indel_size, indel_size_max;
  do {
    vector<Vtx_id> indel_idxs;
    vector<Vtx_id> indel_idxs_max;
    indel_size_max = 0;
    indel_size = 0;
    for (Vtx_id i : idxs) {
      bool is_indel =
          vertices[i].vtx_type != Cap && !blocked[i] && gene_free(i) &&
          indel_count[vertices[i].gene_val] * vertices[i].indel_update < 0;
      if (is_indel) {
        indel_size++;
        indel_idxs.push_back(i);
      }
      if (i == last_origin_vertex || i == last_target_vertex || !is_indel) {
        if (indel_size > indel_size_max) {
          indel_size_max = indel_size;
          indel_idxs_max = indel_idxs;
        }
        indel_size = 0;
        indel_idxs.clear();
      }
    }
    if (indel_size_max == 0) {
      indel_size_max = indel_size;
      indel_idxs_max = indel_idxs;
    }

    for (Vtx_id i : indel_idxs_max) {
      if (!gene_free(i) ||
          indel_count[vertices[i].gene_val] * vertices[i].indel_update >= 0) {
        continue;
      }
      if (use_fix) {
        // the whole exclusive block or nothing
        vector<Vtx_id> out;
        if (indel_closure(i, map<Vtx_id, Vtx_id>(), out)) {
          commit_indel(i, true);
        } else {
          blocked[i] = blocked[vertices[i].red] = true;
        }
        continue;
      }
      vertices[i].is_indel = true;
      vertices[i].is_indel_fix = true;
      vertices[vertices[i].red].is_indel = true;
      vertices[vertices[i].red].is_indel_fix = true;
      indel_count[vertices[i].gene_val] += vertices[i].indel_update;
    }
  } while (indel_size_max > 0);
}

void CycleGraph::find_1_cycles(bool is_random) {
  vector<Vtx_id> idxs = produce_random_idxs(is_random);
  for (Vtx_id i : idxs) {
    get_trivial_bal_cycle(i, is_random);
  }
}

bool check_cycles_intersect(vector<Vtx_id> c1, vector<Vtx_id> c2) {
  for (int element : c1) {
    if (std::find(c2.begin(), c2.end(), element) != c2.end()) {
      return true;
    }
  }
  return false;
}

/* Two candidate cycles conflict if they share a black edge or if the gray
 * edges they fix (including the implied ones) are incompatible. */
bool CycleGraph::cycles_conflict(
    const pair<vector<Vtx_id>, vector<Vtx_id>> &c1,
    const pair<vector<Vtx_id>, vector<Vtx_id>> &c2) {
  if (check_cycles_intersect(c1.first, c2.first)) {
    return true;
  }
  map<Vtx_id, Vtx_id> fix;
  return !(cycle_fixes(c1.second, fix) && cycle_fixes(c2.second, fix));
}

/* Select a set of pairwise compatible cycles from the candidates. This is the
 * maximum independent set problem on the conflict graph, solved with the
 * greedy heuristic: take the candidate with fewest conflicts and discard
 * every candidate that conflicts with it. Ties keep the (random) order of the
 * candidates. */
vector<vector<Vtx_id>> CycleGraph::select_cycles(
    const vector<pair<vector<Vtx_id>, vector<Vtx_id>>> &cands) {
  size_t n = cands.size();
  vector<set<size_t>> conflicts(n);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      if (cycles_conflict(cands[i], cands[j])) {
        conflicts[i].insert(j);
        conflicts[j].insert(i);
      }
    }
  }

  vector<vector<Vtx_id>> selected;
  vector<bool> removed(n, false);
  while (true) {
    size_t best = n;
    for (size_t i = 0; i < n; i++) {
      if (!removed[i] &&
          (best == n || conflicts[i].size() < conflicts[best].size())) {
        best = i;
      }
    }
    if (best == n) {
      break;
    }
    selected.push_back(cands[best].second);
    removed[best] = true;
    for (size_t j : conflicts[best]) {
      removed[j] = true;
      for (size_t k : conflicts[j]) {
        if (k != best) {
          conflicts[k].erase(j);
        }
      }
    }
  }
  return selected;
}

void CycleGraph::find_2_cycles(bool is_random, vector<vector<Vtx_id>> &cycles) {
  cycles = select_cycles(get_all_tuptok_cycles(4, 4, is_random));
}

/* Whether the black edge of x is free to be in a clean cycle */
bool CycleGraph::clean_edge_free(Vtx_id x) const {
  Vtx_id y = vertices[x].black;
  return y != NO_EDGE && y != x && !vertices[x].in_cycle &&
         !vertices[y].in_cycle && !vertices[x].is_indel &&
         !vertices[y].is_indel;
}

/* Paths u1 c1 u2 c2 ... ur cr with r gray edges starting at c0 = c, where ui
 * is a gray neighbor of c(i-1) and ci is the black neighbor of ui, indexed by
 * cr. The black edges of a path are free and distinct. */
void CycleGraph::back_half_paths(
    Vtx_id c, int r, bool is_random, vector<Vtx_id> &path,
    unordered_map<Vtx_id, vector<vector<Vtx_id>>> &paths) {
  if (r == 0) {
    paths[c].push_back(path);
    return;
  }
  for (Vtx_id u : get_neighbors(c, is_random)) {
    if (!clean_edge_free(u) || find(path.begin(), path.end(), u) != path.end())
      continue;
    path.push_back(u);
    path.push_back(vertices[u].black);
    back_half_paths(vertices[u].black, r - 1, is_random, path, paths);
    path.pop_back();
    path.pop_back();
  }
}

/* Extend path = v b1 g1 b2 ... bk (b = bk) with f more gray edges and join it
 * with a path of back (the last gray neighbor g is the cr of the back path).
 * Adds the first clean, balanced and compatible cycle found. */
bool CycleGraph::forward_half_paths(
    Vtx_id b, int f, bool is_random, vector<Vtx_id> &path,
    const unordered_map<Vtx_id, vector<vector<Vtx_id>>> &back) {
  for (Vtx_id g : get_neighbors(b, is_random)) {
    if (!clean_edge_free(g) || find(path.begin(), path.end(), g) != path.end())
      continue;
    if (f > 1) {
      path.push_back(g);
      path.push_back(vertices[g].black);
      if (forward_half_paths(vertices[g].black, f - 1, is_random, path, back))
        return true;
      path.pop_back();
      path.pop_back();
      continue;
    }
    auto it = back.find(g);
    if (it == back.end())
      continue;
    for (const auto &bp : it->second) {
      /* cycle v b1 g1 ... bk g ur c(r-1) u(r-1) ... c1 u1 */
      vector<Vtx_id> cycle = path;
      cycle.push_back(g);
      for (int i = int(bp.size()) - 2; i >= 0; i -= 2) {
        cycle.push_back(bp[i]);
        if (i > 0)
          cycle.push_back(bp[i - 1]);
      }
      bool simple = true;
      for (size_t i = 0; i < cycle.size() && simple; i++)
        for (size_t j = i + 1; j < cycle.size() && simple; j++)
          simple = cycle[i] != cycle[j];
      if (!simple)
        continue;
      IR weight = ZERO_IR;
      for (size_t i = 0; i < cycle.size(); i += 2) {
        weight.first += vertices[cycle[i]].weight.first;
        weight.second += vertices[cycle[i]].weight.second;
      }
      if (ir_intersect(weight, ZERO_IR) && check_cycle(cycle)) {
        add_cycle(cycle);
        return true;
      }
    }
  }
  return false;
}

void CycleGraph::add_small_cycles(int s, bool is_random) {
  for (Vtx_id v : produce_random_idxs(is_random)) {
    if (vertices[v].black < v || !clean_edge_free(v))
      continue;
    unordered_map<Vtx_id, vector<vector<Vtx_id>>> back;
    vector<Vtx_id> path;
    back_half_paths(v, s, is_random, path, back);
    if (back.empty())
      continue;
    path = {v, vertices[v].black};
    forward_half_paths(vertices[v].black, s, is_random, path, back);
  }
}

vector<Vtx_id> CycleGraph::get_neighbors(Vtx_id v, bool is_random) {
  vector<Vtx_id> neigs;
  if (vertices[v].fix_gray != NO_EDGE) {
    neigs.push_back(vertices[v].fix_gray);
  } else if (!vertices[v].is_indel) {
    for (Vtx_id u : vertices[v].grays) {
      if (vertices[u].fix_gray == NO_EDGE && !vertices[u].is_indel) {
        neigs.push_back(u);
      }
    }
    if (is_random) {
      shuffle(neigs.begin(), neigs.end(), rand_gen);
    }
  }
  return neigs;
}

bool CycleGraph::is_only_indel(Cycle c) const {
  bool ok = false;
  Vtx_id v = c.start_vtx;
  if (vertices[v].is_indel) {
    v = vertices[vertices[v].red].black;
    while (vertices[v].is_indel) {
      v = vertices[vertices[v].red].black;
      if (v == c.start_vtx) {
        ok = true;
        break;
      }
    }
  }
  return ok;
}

bool CycleGraph::is_trivial_cycle(Cycle c) const {
  Vtx_id v = c.start_vtx;
  vector<Vtx_id> vs;

  if (is_only_indel(c)) {
    return false;
  }
  v = next_noindel(v, vs);
  v = vertices[v].fix_gray;
  v = next_noindel(v, vs);
  v = vertices[v].fix_gray;
  if (v == c.start_vtx) {
    return true;
  }
  v = vertices[v].black;
  if (v == c.start_vtx) {
    return true;
  }
  while (vertices[v].is_indel) {
    v = vertices[v].red;
    if (v == c.start_vtx) {
      return true;
    }
    v = vertices[v].black;
    if (v == c.start_vtx) {
      return true;
    }
  }

  return false;
}

void CycleGraph::decompose(bool is_random, string dec_algo) {
  if (dec_algo == "bfs") {
    this->decompose_with_bfs(is_random);
  } else if (dec_algo == "1c2c3c") {
    this->decompose_1c_2c_3c(is_random);
  } else if (dec_algo == "map") {
    this->decompose_random_map(is_random);
  } else {
    throw runtime_error("Unknown decomposition algorithm: " + dec_algo);
  }
}

void CycleGraph::decompose_1c_2c_3c(bool is_random) {
  this->mark_indel_edges(is_random);
  /* clean and balanced cycles with 2, 4 and 6 black edges, smallest first */
  for (int s = 1; s <= 3; s++) {
    this->add_small_cycles(s, is_random);
  }
  /* Finish the decomposition with bfs */
  this->decompose_with_bfs(is_random);
}

/* Decompose the graph by using a random assignment of the genes */
void CycleGraph::decompose_random_map(bool is_random) {
  Vtx_id headtailcorresp_v, headtailcorresp_u;

  for (int v = 0; v <= last_origin_vertex; v++) {
    /* Skip vertices marked as indels */
    if (vertices[v].is_indel) {
      continue;
    }
    if (use_fix) {
      /* take the first gray edge whose closure has no conflict */
      if (vertices[v].fix_gray != NO_EDGE) {
        continue;
      }
      for (Vtx_id u : get_neighbors(v, is_random)) {
        map<Vtx_id, Vtx_id> fix;
        if (gray_closure(v, u, fix)) {
          for (auto &e : fix) {
            vertices[e.first].fix_gray = e.second;
          }
          break;
        }
      }
      continue;
    }
    vector<Vtx_id> neigs = get_neighbors(v, is_random);
    /* If there are no neighbors than this character will have no correspondent,
     * so we live it to be eliminated by an indel. */
    if (neigs.empty()) {
      continue;
    }
    Vtx_id u = neigs.front();
    vertices[v].fix_gray = u;
    vertices[u].fix_gray = v;
    if (vertices[v].vtx_type != Cap) {
      headtailcorresp_v = vertices[v].red;
      headtailcorresp_u = vertices[u].red;
      assert(vertices[headtailcorresp_u].fix_gray == NO_EDGE ||
             vertices[headtailcorresp_u].fix_gray == headtailcorresp_v);
      assert(vertices[headtailcorresp_v].fix_gray == NO_EDGE ||
             vertices[headtailcorresp_v].fix_gray == headtailcorresp_u);
      vertices[headtailcorresp_v].fix_gray = headtailcorresp_u;
      vertices[headtailcorresp_u].fix_gray = headtailcorresp_v;
      v++;
    }
  }
  this->decompose_with_bfs(true);
}

/* Decompose the remaning graph */
/* TODO: Think about how the weights affect indels in the bfs */
void CycleGraph::decompose_with_bfs(bool is_random) {
  vector<Vtx_id> idxs = produce_random_idxs(is_random);
  if (only_trivial) {
    for (Vtx_id i : idxs) {
      get_trivial_bal_cycle(i, is_random);
    }
  }
  for (Vtx_id i : idxs) {
    bfs(i, is_random);
  }
  if (use_fix) {
    /* The fixed adjacencies may leave genes without an available pair (the
     * pairs and the indels of whole blocks do not follow the counts of each
     * label). These genes are deleted (inserted) with their blocks: first the
     * ones without a gray edge that can be fixed, then any one left. */
    for (int pass = 0; pass < 2; ++pass) {
      bool marked = false;
      for (Vtx_id v = 0; v < Vtx_id(size()); ++v) {
        if (vertices[v].in_cycle || vertices[v].is_indel ||
            vertices[v].fix_gray != NO_EDGE || vertices[v].vtx_type == Cap) {
          continue;
        }
        bool has_gray = false;
        for (Vtx_id u : (pass == 0) ? get_neighbors(v, false)
                                    : vector<Vtx_id>()) {
          map<Vtx_id, Vtx_id> fix;
          if (gray_closure(v, u, fix)) {
            has_gray = true;
            break;
          }
        }
        vector<Vtx_id> out;
        if (!has_gray && indel_closure(v, map<Vtx_id, Vtx_id>(), out)) {
          commit_indel(v);
          marked = true;
        }
      }
      if (marked) {
        for (Vtx_id i : idxs) {
          bfs(i, is_random);
        }
      }
    }
    /* Dead end of the greedy choices (e.g. a block of a genome that ends in
     * a cap whose only possible pair lost its neighbor to an indel): release
     * the fixed adjacencies of the vertices left so that the decomposition is
     * complete (the partition induced may then have more breakpoints) */
    bool released = false;
    for (Vtx_id v = 0; v < Vtx_id(size()); ++v) {
      if (!vertices[v].in_cycle && vertices[v].fix_adj) {
        vertices[v].fix_adj = vertices[vertices[v].black].fix_adj = false;
        released = true;
        released_adj++;
      }
    }
    if (released) {
      for (Vtx_id i : idxs) {
        bfs(i, is_random);
      }
    }
  }
}

void CycleGraph::get_trivial_bal_cycle(Vtx_id start, bool is_random) {
  vector<Vtx_id> cycle;
  Vtx_id b1 = vertices[start].black;
  Vtx_id b1_ = b1;
  map<Gene, int> indel_count = this->indel_count;

  if (vertices[start].in_cycle)
    return;

  vector<Vtx_id> neigs = get_neighbors(b1_, is_random);

  for (Vtx_id g1 : neigs) {
    if (vertices[g1].in_cycle) {
      continue;
    }
    Vtx_id b2 = vertices[g1].black;
    Vtx_id b2_ = b2;

    vector<Vtx_id> neigs2 = get_neighbors(b2_, is_random);
    for (Vtx_id g2 : neigs2) {
      if (vertices[g2].in_cycle) {
        continue;
      }
      if (g2 == start &&
          ir_intersect(ir_sum(vertices[b1].weight, vertices[b2].weight),
                       ZERO_IR) &&
          (!use_fix || check_cycle({start, b1, g1, b2}))) {
        cycle.resize(0);
        cycle.push_back(start);
        cycle.push_back(b1);
        cycle.push_back(g1);
        cycle.push_back(b2);
      }
    }
  }

  if (cycle.size() > 0) {
    add_cycle(cycle);
  }
}

/* Follow indel edges from the last entry of queue and add all vertices reached
 * to the queue */
void CycleGraph::follow_indel_edge(set<Vtx_id> &vizited_in_level,
                                   vector<unique_ptr<QEntry>> &q,
                                   unique_ptr<QEntry> &entry) {
  Vtx_id v = vertices[entry->vtx].black;
  /* For debug purposes */
  // cout << v << " " << (vertices[v].fix_gray == NO_EDGE) << " "
  //      << (vertices[v].vtx_type != Cap) << " "
  //      << (!vertices[vertices[v].red].in_cycle) << endl;
  /* Stop at a fixed gray edge at the end of the chromosome (linear or
   * circular), or if the red neighbor is in a cycle */
  /* TODO: See how to consider the fact that more red edges is better (maybe
   * starting by changing this if to a while and adding the commented
   * follow_indel_edge in the next function. Also I may need to add some
   * condition to stop in the end of a circular chromosome) If I change this I
   * have to remember to describe on the dissertation. */
  if (vertices[v].fix_gray == NO_EDGE && vertices[v].vtx_type != Cap &&
      !vertices[vertices[v].red].in_cycle) {
    Vtx_id u = vertices[v].red;
    /* Check if we have an indel available and if the vertex has not been
     * visited (there is a special excepion for circular genomes of size one,
     * where the vertex is already visited). */
    if ((entry->indel_count[vertices[u].gene_val] * vertices[u].indel_update <
             0 ||
         vertices[v].is_indel) &&
        (entry->not_yet_visited(vizited_in_level, u) ||
         vertices[v].red == vertices[v].black)) {
      if (use_fix) {
        /* A fixed adjacency has both extremities indels or none, and the
         * whole exclusive block must be deleted (inserted) */
        Vtx_id p = entry->vtx, pp = partner(p, entry->fixed);
        vector<Vtx_id> out;
        if ((vertices[v].fix_adj && pp != NO_EDGE && pp != vertices[p].red) ||
            !indel_closure(v, entry->fixed, out)) {
          return;
        }
      }
      q.push_back(make_unique<QEntry>(v, u, entry->fixed, entry->vizited,
                                      entry->indel_count));
      if (!vertices[v].is_indel) {
        q.back()->indel_count[vertices[u].gene_val] += vertices[u].indel_update;
      }
      vizited_in_level.insert(u);
    }
  }
}

void print_queue(const vector<unique_ptr<QEntry>> &q, int pos) {
  cout << "Queue: ";
  for (size_t i = pos; i < q.size(); ++i) {
    cout << q[i]->vtx << " ";
  }
  cout << endl;
}

void CycleGraph::bfs(Vtx_id start, bool is_random) {
  unique_ptr<QEntry> entry;
  vector<unique_ptr<QEntry>> q1;
  vector<unique_ptr<QEntry>> q2;
  set<Vtx_id> vizited_in_level;
  Vtx_id headtailcorresp_v, headtailcorresp_u;
  size_t pos = 0;

  if (vertices[start].in_cycle)
    return;

  if (use_fix && vertices[start].fix_gray == NO_EDGE &&
      !vertices[start].is_indel) {
    /* skip the search if start cannot be closed: no gray edge with a valid
     * closure and no indel allowed (it would explore the whole graph) */
    bool can_close = false;
    for (Vtx_id u : get_neighbors(start, false)) {
      map<Vtx_id, Vtx_id> fix;
      if (gray_closure(start, u, fix)) {
        can_close = true;
        break;
      }
    }
    vector<Vtx_id> out;
    if (!can_close && vertices[start].vtx_type != Cap &&
        indel_count[vertices[start].gene_val] * vertices[start].indel_update <
            0 &&
        indel_closure(start, map<Vtx_id, Vtx_id>(), out)) {
      can_close = true;
    }
    if (!can_close) {
      return;
    }
  }

  q1.push_back(make_unique<QEntry>(start, this->indel_count));
  vizited_in_level.insert(start);
  /* follow_indel_edge(vizited_in_level, q1, q1.back()); */
#ifdef DEBUG_CYCLEGRAPH
  cout << *this << endl;
#endif

  do {
    entry = std::move(q1[pos]);
    pos++;
#ifdef DEBUG_CYCLEGRAPH
    cout << "current entry vtx: " << entry->vtx << endl;
    cout << "q1: ";
    print_queue(q1, pos);
#endif

    Vtx_id v = vertices[entry->vtx].black;

    /* Follow each gray edge */
    vector<Vtx_id> neigs;
    auto vu = entry->fixed.find(v);
    if (vu != entry->fixed.end()) {
      neigs.push_back(vu->second);
    } else {
      neigs = get_neighbors(v, is_random);
    }

    vector<Vtx_id> notDone;
    for (Vtx_id u : neigs) {
      auto uv = entry->fixed.find(u);
      if ((uv == entry->fixed.end() || uv->second == v) &&
          !vertices[u].in_cycle &&
          entry->not_yet_visited(vizited_in_level, u)) {
        notDone.push_back(u);
      }
    }

    if (use_fix) {
      /* keep only the gray edges whose closure has no conflict (the cycle
       * found is checked with the whole closures in the end) */
      for (Vtx_id u : notDone) {
        map<Vtx_id, Vtx_id> add;
        if (!gray_closure(v, u, add, &entry->fixed)) {
          continue;
        }
        if (vertices[v].fix_gray == NO_EDGE && vertices[v].vtx_type != Cap) {
          q2.push_back(make_unique<QEntry>(
              v, u, vertices[v].red, vertices[u].red, entry->fixed,
              entry->vizited, entry->indel_count));
        } else {
          q2.push_back(make_unique<QEntry>(v, u, entry->fixed, entry->vizited,
                                           entry->indel_count));
        }
        vizited_in_level.insert(u);
      }
    } else if (vertices[v].fix_gray == NO_EDGE &&
               vertices[v].vtx_type != Cap) {
      // v doest have a fixed gray edge and is not a cap vertex
      for (Vtx_id u : notDone) {
        headtailcorresp_v = vertices[v].red;
        headtailcorresp_u = vertices[u].red;
        q2.push_back(make_unique<QEntry>(v, u, headtailcorresp_v,
                                         headtailcorresp_u, entry->fixed,
                                         entry->vizited, entry->indel_count));
        vizited_in_level.insert(u);
        /* follow_indel_edge(vizited_in_level, q2, q2.back()); */
      }
    } else {
      for (Vtx_id u : notDone) {
        q2.push_back(make_unique<QEntry>(v, u, entry->fixed, entry->vizited,
                                         entry->indel_count));
        vizited_in_level.insert(u);
        /* follow_indel_edge(vizited_in_level, q2, q2.back()); */
      }
    }

    follow_indel_edge(vizited_in_level, q2, entry);

#ifdef DEBUG_CYCLEGRAPH
    cout << "q2: ";
    print_queue(q2, 0);
#endif
    if (pos == q1.size()) {
      q1.clear();
      for (unique_ptr<QEntry> &e : q2) {
        q1.push_back(std::move(e));
      }
      if (is_random) {
        shuffle(q1.begin(), q1.end(), rand_gen);
      }
      q2.clear();
      vizited_in_level.clear();
      pos = 0;
    }

  } while (pos < q1.size() && q1[pos]->vtx != start);

  if (pos >= q1.size()) {
    return;
  }

  /* Find balanced cycle in level if it exists */
  vector<Vtx_id> cycle, cand;
  for (; pos < q1.size(); pos++) {
    if (q1[pos]->vtx == start) {
      entry = std::move(q1[pos]);
      cand.clear();
      Vtx_id u, v;
      IR weight = ZERO_IR;
      v = start;
      do {
        u = vertices[v].black;
        cand.push_back(v);
        cand.push_back(u);
        v = entry->fixed[u];
        if (v != vertices[u].red) {
          weight.first += vertices[v].weight.first;
          weight.second += vertices[v].weight.second;
        }
      } while (v != start);
      /* with fixed adjacencies, the cycle must also respect them where it
       * closes (at start) */
      if (use_fix && !check_cycle(cand)) {
        continue;
      }
      cycle = cand;
      if (weight.first <= 0 && 0 <= weight.second) {
        break;
      }
    }
  }

  /* Add the balanced cycle, if there is none take any cycle in level */
  /* TODO: test if there is a way of looking for balanced cycles of other
   * sizes as well (for instance with a parameter to keep running the search
   * until a certain number of cycles is reached) */
  if (!cycle.empty()) {
    add_cycle(cycle);
  }
}

string CycleGraph::show_cycles() const {
  ostringstream ss;

  ss << '[';
  for (size_t i = 0; i < cycles.size(); ++i) {
    Vtx_id v = cycles[i].start_vtx;
    ss << '[' << v;
    v = vertices[v].black;
    ss << "," << v;
    v = vertices[v].is_indel ? vertices[v].red : vertices[v].fix_gray;
    while (v != cycles[i].start_vtx) {
      ss << "," << v;
      v = vertices[v].black;
      ss << "," << v;
      v = vertices[v].is_indel ? vertices[v].red : vertices[v].fix_gray;
    }
    ss << ']';
    if (i != cycles.size() - 1)
      ss << ',';
  }
  ss << ']';

  return ss.str();
}

void CycleGraph::read_decomposition(const string &str) {
  // a vertex as written by serialize: id(sign gene,black,ir,in_cycle,type|[edges])
  static const regex vtx_re(
      R"((\d+)\([+-]\d+,(-?\d+),-?\d+:-?\d+,[o_],([HTC])\|\[([^\]]*)\]\))");
  size_t n = 0;
  for (sregex_iterator it(str.begin(), str.end(), vtx_re), end; it != end;
       ++it, ++n) {
    Vtx_id v = stoi((*it)[1]);
    VtxType type = (*it)[3] == "H" ? Head : ((*it)[3] == "T" ? Tail : Cap);
    if (v != Vtx_id(n) || v > last_target_vertex ||
        stoi((*it)[2]) != vertices[v].black || type != vertices[v].vtx_type) {
      throw runtime_error("decomposition does not match the genomes (vertex " +
                          to_string(n) + ")");
    }
    string edge = (*it)[4];
    if (edge.compare(0, 2, "!!") == 0) {
      vertices[v].is_indel = true;
    } else if (edge.compare(0, 1, "*") == 0) {
      vertices[v].fix_gray = stoi(edge.substr(1));
    } else {
      throw runtime_error("vertex " + to_string(v) +
                          " is not in the decomposition");
    }
    vertices[v].in_cycle = true;
  }
  if (n != size()) {
    throw runtime_error("decomposition has " + to_string(n) +
                        " vertices, the genomes have " + to_string(size()));
  }
  for (size_t v = 0; v < size(); ++v) {
    Vtx_id u = vertices[v].fix_gray;
    if (!vertices[v].is_indel &&
        (u < 0 || u >= Vtx_id(size()) || vertices[u].fix_gray != Vtx_id(v))) {
      throw runtime_error("invalid gray edge in the decomposition (vertex " +
                          to_string(v) + ")");
    }
  }
}

void CycleGraph::rebuild_cycles() {
  vector<bool> seen(size(), false);
  for (size_t v = 0; v < size(); ++v) {
    if (seen[v]) {
      continue;
    }
    vector<Vtx_id> cy = get_cycle(v);
    for (Vtx_id u : cy) {
      seen[u] = true;
      vertices[u].in_cycle = false;
    }
    add_cycle(cy);
  }
}

void CycleGraph::read_cycles(string str) {
  string str_cycle, str_el;
  vector<Vtx_id> cycle;

  /* Remove '[[' and ']]'. */
  str = str.substr(2, str.size() - 4);

  /* Read each cycle. */
  string delimiter = "],[";
  while (not str.empty()) {
    size_t idx = str.find(delimiter);
    if (idx != str.npos) {
      str_cycle = str.substr(0, idx);
      str.erase(0, idx + delimiter.length());
    } else {
      str_cycle = str;
      str = "";
    }

    /* Read each element. */
    cycle.clear();
    stringstream ss(str_cycle);
    while (getline(ss, str_el, ',')) {
      cycle.push_back(stoi(str_el));
    }

    add_cycle(cycle);
  }
}

void CycleGraph::make_assignment() {
  int val = 1;
  for (size_t i = last_origin_vertex + 1; i <= last_target_vertex; ++i) {
    if (vertices[i].is_indel) {
      vertices[i].gene_val = 0;
    } else {
      vertices[i].gene_val = val;
      if (vertices[i].vtx_type == Cap || vertices[i].red < i) {
        val++;
      }
    }
  }
  for (size_t i = 0; i <= last_origin_vertex; ++i) {
    if (vertices[i].is_indel) {
      vertices[i].gene_val = 0;
      vertices[i].sign_positive = true;
    } else {
      vertices[i].gene_val = vertices[vertices[i].fix_gray].gene_val;
      vertices[i].sign_positive = (vertices[vertices[i].fix_gray].sign_positive)
                                      ? vertices[i].sign_positive
                                      : !vertices[i].sign_positive;
    }
  }
  for (size_t i = last_origin_vertex + 1; i <= last_target_vertex; ++i) {
    vertices[i].sign_positive = true;
  }
}

Genea vertex_to_gene(Vertex v) {
  return Genea(((v.sign_positive) ? 1 : -1) * v.gene_val, v.is_indel);
}

/* TODO: Verify if genes to be deleted or added should be included in the
 * genome. Currently we delete those genes.
 */
MultiGenomePair CycleGraph::get_genomes() const {
  MultiGenomePair data;
  vector<pair<bool, vector<Genea>>> gss;
  vector<vector<IR>> irss;

  for (size_t i = 0; i <= last_target_vertex;) {
    if (vertices[i].vtx_type == Cap) {
      gss.push_back(pair<bool, vector<Genea>>(false, vector<Genea>()));
      irss.push_back(vector<IR>());
      assert(!vertices[i].is_indel);
      gss.back().second.push_back(vertex_to_gene(vertices[i]));
      i++;
      while (vertices[i].vtx_type != Cap) {
        gss.back().second.push_back(vertex_to_gene(vertices[i]));
        irss.back().push_back(ir_abs(vertices[i].weight));
        if (vertices[i].is_indel) {
          while (vertices[i].is_indel) {
            i += 2;
          }
        } else {
          i += 2;
        }
      }
      assert(!vertices[i].is_indel);
      gss.back().second.push_back(vertex_to_gene(vertices[i]));
      irss.back().push_back(ir_abs(vertices[i].weight));
      i++;
    } else {
      gss.push_back(pair<bool, vector<Genea>>(true, vector<Genea>()));
      irss.push_back(vector<IR>());
      while (i <= last_target_vertex) {
        gss.back().second.push_back(vertex_to_gene(vertices[i]));
        irss.back().push_back(ir_abs(vertices[i].weight));
        if (vertices[i + 1].black <= i) {
          i += 2;
          break;
        }
        if (vertices[i].is_indel) {
          while (i <= last_target_vertex && vertices[i].is_indel &&
                 vertices[i + 1].black > i) {
            i += 2;
          }
        } else {
          i += 2;
        }
      }
    }

    if (i == last_origin_vertex + 1) {
      data.g = make_unique<MultiGenome>(gss, irss);
      gss.clear();
      irss.clear();
    }
    if (i == last_target_vertex + 1) {
      data.h = make_unique<MultiGenome>(gss, irss);
    }
  }

  return data;
}

PermsIrs CycleGraph::get_perms() const {
  PermsIrs perms_irs;

  MultiGenomePair data = get_genomes();
  if (data.g->n_chrom() != 1 || data.h->n_chrom() != 1) {
    throw runtime_error("get_perms can only be used if the cycle was produce "
                        "from single chromosomes.");
  }

  perms_irs.s_n = (*data.g)[1].size();
  perms_irs.p_n = (*data.h)[1].size();
  perms_irs.s = (int *)malloc(perms_irs.s_n * sizeof(int));
  perms_irs.s_ir_low = (int *)malloc((perms_irs.s_n - 1) * sizeof(int));
  perms_irs.s_ir_high = (int *)malloc((perms_irs.s_n - 1) * sizeof(int));
  perms_irs.p = (int *)malloc(perms_irs.p_n * sizeof(int));
  perms_irs.p_ir_low = (int *)malloc((perms_irs.p_n - 1) * sizeof(int));
  perms_irs.p_ir_high = (int *)malloc((perms_irs.p_n - 1) * sizeof(int));

  for (int i = 1; i <= perms_irs.s_n; i++) {
    perms_irs.s[i - 1] = (*data.g)[1][i];
  }
  for (int i = 1; i < perms_irs.s_n; i++) {
    perms_irs.s_ir_low[i - 1] = (*data.g)[1].get_ir(i).first;
    perms_irs.s_ir_high[i - 1] = (*data.g)[1].get_ir(i).second;
  }
  for (int i = 1; i <= perms_irs.p_n; i++) {
    perms_irs.p[i - 1] = (*data.h)[1][i];
  }
  for (int i = 1; i < perms_irs.p_n; i++) {
    perms_irs.p_ir_low[i - 1] = (*data.h)[1].get_ir(i).first;
    perms_irs.p_ir_high[i - 1] = (*data.h)[1].get_ir(i).second;
  }

  return perms_irs;
}

PermsIrs flip_perms_irs(PermsIrs perms_irs) {
  PermsIrs perms_irs_flip;
  perms_irs_flip.s = perms_irs.p;
  perms_irs_flip.p = perms_irs.s;
  perms_irs_flip.s_ir_low = perms_irs.p_ir_low;
  perms_irs_flip.s_ir_high = perms_irs.p_ir_high;
  perms_irs_flip.p_ir_low = perms_irs.s_ir_low;
  perms_irs_flip.p_ir_high = perms_irs.s_ir_high;
  perms_irs_flip.s_n = perms_irs.p_n;
  perms_irs_flip.p_n = perms_irs.s_n;
  return perms_irs_flip;
}

/* TODO: change function so it does not need to create a Permutation
 * (which needs non trivial cycles), after that remove the bfs from the
 * decomposition when only_trivial is true. */
CommonPartition CycleGraph::get_part(bool allow_rev) const {
  PermsIrs perms_irs = get_perms();
  Permutation pi_s(perms_irs, rand_gen);
  CommonPartition cp;

  for (size_t i = 1; i < pi_s.size(); ++i) {
    Gene a = pi_s[i];
    if (pi_s.breakpoint(a, allow_rev)) {
      cp.breaks_s.push_back(i);
    }
  }

  Permutation pi_p(flip_perms_irs(perms_irs), rand_gen);

  for (size_t i = 1; i < pi_p.size(); ++i) {
    Gene a = pi_p[i];
    if (pi_p.breakpoint(a, allow_rev)) {
      cp.breaks_p.push_back(i);
    }
  }

  free(perms_irs.s);
  free(perms_irs.s_ir_low);
  free(perms_irs.s_ir_high);
  free(perms_irs.p);
  free(perms_irs.p_ir_low);
  free(perms_irs.p_ir_high);
  return cp;
}

void CycleGraph::serialize(ostream &os) const {
  for (size_t i = 0; i <= last_target_vertex; ++i) {
    os << i << "(" << (vertices[i].sign_positive ? "+" : "-")
       << vertices[i].gene_val << "," << vertices[i].black << ","
       << vertices[i].weight.first << ":" << vertices[i].weight.second << ","
       << (vertices[i].in_cycle ? "o" : "_") << ",";
    switch (vertices[i].vtx_type) {
    case Head:
      os << "H";
      break;
    case Tail:
      os << "T";
      break;
    case Cap:
      os << "C";
      break;
    }
    os << "|[";
    if (vertices[i].is_indel) {
      os << "!!" << vertices[i].red;
    } else if (vertices[i].fix_gray == NO_EDGE) {
      for (auto gray : vertices[i].grays) {
        os << gray << ",";
      }
      os << "!" << vertices[i].red;
    } else {
      os << "*" << vertices[i].fix_gray;
    }
    os << "])"
       << " ";
    if (int(i) == last_origin_vertex)
      os << endl;
  }
}

ostream &operator<<(std::ostream &os, const CycleGraph &cg) {
  cg.serialize(os);
  return os;
}

ostream &operator<<(ostream &os, const Run &r) {
  os << r.genome << r.fst_gene << ":";
  for (Gene a : r.genes_to_add) {
    os << " " << a;
  }
  os << endl;
  return os;
}

void CycleGraph::rem_cycle(Cycle c) {
  Vtx_id headtailcorresp_v, headtailcorresp_u;
  Vtx_id i = c.start_vtx;

  indel_potation -= cycle_potation(i);
  vector<Vtx_id> cycle = get_cycle(i);

  if (ir_intersect(cycle_weight(i), ZERO_IR)) {
    balanced_cycles--;
    if (c.pot == 0) {
      clean_balanced_cycles--;
    }
    if (cycle.size() == 4) {
      trivial_balanced_cycles--;
    }
  }

  if (is_divergent(i)) {
    divergent_cycles--;
    if (cycle_weight(i).first <= 0) {
      divergent_balanced_or_neg_cycles--;
    }
  }

  assert(cycle.size() % 2 == 0);
  for (size_t i = 0; i < cycle.size(); ++i) {
    Vtx_id v = cycle[i];
    assert(vertices[v].in_cycle);
    vertices[v].in_cycle = false;
    if (i % 2 == 1) {
      Vtx_id u = (i == cycle.size() - 1) ? cycle[0] : cycle[i + 1];
      if (vertices[v].is_indel) {
        if (!vertices[v].is_indel_fix) {
          vertices[v].is_indel = false;
          vertices[u].is_indel = false;
          this->indel_count[vertices[v].gene_val] -= vertices[v].indel_update;
        }
      } else if (vertices[v].fix_char) {
        continue; // fixed by fix_partition
      } else if (vertices[v].vtx_type == Cap) {
        assert(vertices[v].fix_gray != NO_EDGE);
        assert(vertices[u].fix_gray != NO_EDGE);
        vertices[v].fix_gray = NO_EDGE;
        vertices[u].fix_gray = NO_EDGE;
      } else {
        headtailcorresp_v = vertices[v].red;
        headtailcorresp_u = vertices[u].red;
        if (not vertices[headtailcorresp_v].in_cycle &&
            not vertices[headtailcorresp_u].in_cycle) {
          assert(vertices[v].fix_gray != NO_EDGE);
          assert(vertices[u].fix_gray != NO_EDGE);
          assert(vertices[headtailcorresp_v].fix_gray != NO_EDGE);
          assert(vertices[headtailcorresp_u].fix_gray != NO_EDGE);
          vertices[v].fix_gray = NO_EDGE;
          vertices[u].fix_gray = NO_EDGE;
          vertices[headtailcorresp_v].fix_gray = NO_EDGE;
          vertices[headtailcorresp_u].fix_gray = NO_EDGE;
        }
      }
    }
  }
  cycles.erase(find(cycles.begin(), cycles.end(), c));
}

bool CycleGraph::cycle_fixes(const vector<Vtx_id> &cycle,
                             map<Vtx_id, Vtx_id> &fix) {
  auto get = [&](Vtx_id x) {
    auto it = fix.find(x);
    return it != fix.end() ? it->second : vertices[x].fix_gray;
  };
  for (Vtx_id v : cycle)
    if (vertices[v].in_cycle)
      return false;
  if (use_fix) {
    size_t n = cycle.size();
    for (size_t i = 1; i < n; i += 2) {
      Vtx_id v = cycle[i], u = cycle[(i + 1) % n];
      if (u == vertices[v].red) {
        vector<Vtx_id> out;
        if (!indel_closure(v, fix, out))
          return false;
        for (Vtx_id x : out)
          fix[x] = vertices[x].red;
      } else if (!gray_closure(v, u, fix)) {
        return false;
      }
    }
    // a fixed adjacency has both extremities indels or none
    for (size_t i = 0; i < n; i += 2) {
      Vtx_id x = cycle[i], y = cycle[i + 1];
      if (vertices[x].fix_adj &&
          (cycle[(i + n - 1) % n] == vertices[x].red) !=
              (cycle[(i + 2) % n] == vertices[y].red))
        return false;
    }
    return true;
  }
  for (size_t i = 1; i < cycle.size(); i += 2) {
    Vtx_id v = cycle[i];
    Vtx_id u = (i == cycle.size() - 1) ? cycle[0] : cycle[i + 1];
    if (vertices[v].is_indel) {
      if (vertices[v].red != u)
        return false;
      continue;
    }
    if (get(v) != NO_EDGE || get(u) != NO_EDGE) {
      if (get(v) != u || get(u) != v)
        return false;
      continue;
    }
    if (u == vertices[v].red) // add_cycle marks this edge as an indel
      continue;
    fix[v] = u;
    fix[u] = v;
    // add_cycle also fixes the edge between the other extremities of the genes
    if (vertices[v].vtx_type != Cap) {
      Vtx_id rv = vertices[v].red, ru = vertices[u].red;
      if (get(rv) != NO_EDGE || get(ru) != NO_EDGE)
        return false;
      fix[rv] = ru;
      fix[ru] = rv;
    }
  }
  return true;
}

bool CycleGraph::check_cycle(vector<Vtx_id> cycle) {
  map<Vtx_id, Vtx_id> fix;
  return cycle_fixes(cycle, fix);
}

void CycleGraph::add_cycle(vector<Vtx_id> cycle) {
  Vtx_id headtailcorresp_v, headtailcorresp_u;

  assert(cycle.size() % 2 == 0);
  if (use_fix) {
    /* fix the closures of the edges of the cycle (the whole blocks) */
    for (size_t i = 1; i < cycle.size(); i += 2) {
      Vtx_id v = cycle[i], u = cycle[(i + 1) % cycle.size()];
      if (u == vertices[v].red) {
        commit_indel(v);
      } else {
        commit_gray(v, u);
      }
    }
  }
  for (size_t i = 0; i < cycle.size(); ++i) {
    Vtx_id v = cycle[i];
    assert(not vertices[v].in_cycle);
    vertices[v].in_cycle = true;
    Vtx_id u = (i == cycle.size() - 1) ? cycle[0] : cycle[i + 1];
    if (i % 2 == 1 && vertices[v].fix_gray != NO_EDGE &&
        vertices[v].fix_gray != u && u != vertices[v].red) {
      throw runtime_error("cycle does not follows the fixed edge.");
    }
    if (i % 2 == 1 && vertices[v].fix_gray == NO_EDGE) {
      if (u == vertices[v].red) {
        if (!vertices[v].is_indel) {
          vertices[v].is_indel = true;
          vertices[u].is_indel = true;
          this->indel_count[vertices[v].gene_val] += vertices[v].indel_update;
        }
      } else {
        if (vertices[u].fix_gray == NO_EDGE) {
          vertices[v].fix_gray = u;
          vertices[u].fix_gray = v;
          if (vertices[v].vtx_type != Cap) {
            headtailcorresp_v = vertices[v].red;
            headtailcorresp_u = vertices[u].red;
            assert(vertices[headtailcorresp_u].fix_gray == NO_EDGE);
            assert(vertices[headtailcorresp_v].fix_gray == NO_EDGE);
            vertices[headtailcorresp_v].fix_gray = headtailcorresp_u;
            vertices[headtailcorresp_u].fix_gray = headtailcorresp_v;
          }
        } else {
          if (vertices[u].fix_gray != v || vertices[v].fix_gray != u) {
            throw runtime_error(
                "Atempting to add a cycle with gray edge (" + to_string(v) +
                " <-> " + to_string(u) + ") distinct from the fixed one (" +
                to_string(v) + " <-> " + to_string(vertices[v].fix_gray) +
                " or " + to_string(u) + " <-> " +
                to_string(vertices[u].fix_gray) + ").");
          }
        }
      }
    }
  }
  int pot = cycle_potation(cycle[0]);
  IR weight = cycle_weight(cycle[0]);
  cycles.push_back(Cycle(cycle, weight, pot));

  indel_potation += pot;
  if (ir_intersect(weight, ZERO_IR)) {
    balanced_cycles++;
    if (pot == 0) {
      clean_balanced_cycles++;
    }
    if (cycle.size() == 4) {
      trivial_balanced_cycles++;
    }
  }

  if (is_divergent(cycle[0])) {
    divergent_cycles++;
    if (weight.first <= 0) {
      divergent_balanced_or_neg_cycles++;
    }
  }
}

void CycleGraph::check_and_add_cycle(vector<Vtx_id> cycle) {
  if (check_cycle(cycle)) {
    add_cycle(cycle);
  }
}

vector<Vtx_id> CycleGraph::get_cycle(Vtx_id i) const {
  vector<Vtx_id> cycle;

  Vtx_id v = i;
  cycle.push_back(v);
  v = vertices[v].black;
  cycle.push_back(v);
  v = vertices[v].is_indel ? vertices[v].red : vertices[v].fix_gray;
  while (v != cycle[0]) {
    cycle.push_back(v);
    v = vertices[v].black;
    cycle.push_back(v);
    v = vertices[v].is_indel ? vertices[v].red : vertices[v].fix_gray;
  }

  return cycle;
}

int CycleGraph::is_divergent(Vtx_id i) const {
  bool ans = false;
  int flips = 0;
  int count = 0;
  bool dir = i < vertices[i].black;

  Vtx_id v = i;
  int end = v;
  do {
    if (dir != (v < vertices[v].black)) {
      flips += 1;
      if (flips == 1) {
        end = v;
      }
      if (flips >= 3) {
        ans = true;
      }
    }
    if (flips == 1) {
      count += 1;
    }
    if (flips == 2) {
      count -= 1;
    }

    v = vertices[v].black;
    if (vertices[v].is_indel) {
      v = vertices[v].red;
    } else {
      v = vertices[v].fix_gray;
    }
  } while (v != end && !ans);

  if (ans || count != 0)
    return true;
  else
    return false;
}

IR CycleGraph::cycle_weight(Vtx_id i) const {
  IR weight = ZERO_IR;

  Vtx_id v = i;
  do {
    v = vertices[v].black;
    if (vertices[v].is_indel) {
      v = vertices[v].red;
    } else {
      v = vertices[v].fix_gray;
      weight.first += vertices[v].weight.first;
      weight.second += vertices[v].weight.second;
    }
  } while (v != i);

  return weight;
}

int CycleGraph::cycle_potation(int i) const {
  int runs = 0;
  int first_state = 0, last_state = 0, state = 0;

  Vtx_id v = i;
  do {
    v = vertices[v].black;
    if (vertices[v].is_indel) {
      state = vertices[v].indel_update;
      if (last_state != state) {
        if (last_state != 0) {
          runs++;
        } else {
          first_state = state;
        }
        last_state = state;
      }
      v = vertices[v].red;
    } else {
      v = vertices[v].fix_gray;
    }
  } while (v != i);
  if (last_state != 0 && (last_state != first_state || runs == 0)) {
    runs++;
  }

  if (runs == 0) {
    return 0;
  } else {
    return ceil((runs + 1) / 2.0);
  }
}

Run CycleGraph::cycle_run(int i) const {
  Run run;
  int run_state = 0,
      state = 0; // States indicate if a run is a insertion or a deletion run
  Vtx_id v = i;

  // Find one end of the run
  int first_end = -1;
  do {
    v = vertices[v].black;
    if (vertices[v].is_indel) {
      state = vertices[v].indel_update;
      if (run_state != state) {
        if (run_state == 0) {
          run_state = state;
        } else {
          first_end = v;
          break;
        }
      }
      v = vertices[v].red;
    } else {
      v = vertices[v].fix_gray;
    }
  } while (v != i);

  if (run_state == 0)
    return run; // No run on the cycle

  bool revert = false;
  if (first_end == -1) { // single run
    while (vertices[v].indel_update == run_state) {
      v = vertices[v].black;
      if (vertices[v].is_indel) {
        v = vertices[v].red;
      } else {
        v = vertices[v].fix_gray;
      }
    }
    first_end = v;
  }

  // Get insertion position and check if the list of insertions must be
  // inverted
  assert(vertices[first_end].indel_update != run_state);
  Vtx_id left; // left vertex of the black edge where the run is inserted
  if (first_end > vertices[first_end].black) {
    left = vertices[v].black;
  } else {
    left = v;
    revert = true;
  }
  run.fst_gene = vertices[left].gene_val;
  run.after_cap = vertices[left].vtx_type == Cap;

  // Get genome where the insertion will occur
  if (run_state == -1) {
    run.genome = 'H';
  } else {
    run.genome = 'G';
  }

  // Get genes to insert
  v = first_end;
  do {
    v = vertices[v].black;
    if (vertices[v].is_indel) {
      state = vertices[v].indel_update;
      if (run_state != state) {
        break;
      }

      if (vertices[v].sign_positive == (v < vertices[v].red)) {
        run.genes_to_add.push_back(vertices[v].gene_val);
      } else {
        run.genes_to_add.push_back(-vertices[v].gene_val);
      }

      v = vertices[v].red;
    } else {
      v = vertices[v].fix_gray;
    }
  } while (v != first_end);

  if (revert) {
    reverse(run.genes_to_add.begin(), run.genes_to_add.end());
    for (int i = 0; i < int(run.genes_to_add.size()); i++) {
      run.genes_to_add[i] = -run.genes_to_add[i];
    }
  }

  return run;
}

/* Get cycles with between t and k black edges, counting black edges connected
 * with indel edges as only one. The return is a vector of pair representing
 * each cycle. The first element of the pair is a representation with one vertex
 * per edge, and the second element is the full representation with all vertices
 * of the cycle. */
vector<pair<vector<Vtx_id>, vector<Vtx_id>>>
CycleGraph::get_all_tuptok_cycles(int t, int k, bool is_random) {
  vector<pair<vector<Vtx_id>, vector<Vtx_id>>> tuptok_cycles;
  vector<Vtx_id> cycle_compact;
  vector<vector<Vtx_id>> edges;
  vector<Vtx_id> idxs = produce_random_idxs(is_random);

  for (Vtx_id i : idxs) {
    if (vertices[i].black < i && !vertices[i].is_indel) {
      get_tuptok_cycle(i, NO_EDGE, t, k, cycle_compact, edges, tuptok_cycles,
                       is_random);
    }
  }

  return tuptok_cycles;
}

/* Find next vertex in the cycle that is not an indel. If all vertices are
 * indels returns v. It follows the first and last black edge of the sequence.
 * All visited vertices are recorded in vs. */
Vtx_id CycleGraph::next_noindel(Vtx_id v, vector<Vtx_id> &vs) const {
  Vtx_id next = vertices[v].black;
  vs.push_back(next);
  while (vertices[next].is_indel && vertices[next].red != v) {
    Vtx_id r = vertices[next].red;
    vs.push_back(r);
    next = vertices[r].black;
    vs.push_back(next);
  }
  return next;
}

/* Recursive function to find all cycles with 2 up to k remaining black edges
 * (combined with indel edges) starting at vertex start and currently at vertex
 * v. */
void CycleGraph::get_tuptok_cycle(
    Vtx_id start, Vtx_id v, int t, int k, vector<Vtx_id> &cycle_compact,
    vector<vector<Vtx_id>> &edges,
    vector<pair<vector<Vtx_id>, vector<Vtx_id>>> &tuptok_cycles,
    bool is_random) {

  if (v == NO_EDGE) {
    v = start;
  } else {

    if (vertices[v].in_cycle) {
      return; // Already in a cycle
    }

    if (v > start) {
      return; // Avoid cycles with vertices after the start vertex
    }
    if (v == start) {
      if (cycle_compact.size() >= t) {
        vector<Vtx_id> cycles;
        for (const auto &e : edges) {
          cycles.insert(cycles.end(), e.begin(), e.end());
        }
        tuptok_cycles.push_back(
            pair<vector<Vtx_id>, vector<Vtx_id>>(cycle_compact, cycles));
      }
      return;
    }
  }
  if (k == 0) {
    return;
  }

  vector<Vtx_id> edge_vs;
  edge_vs.push_back(v);
  Vtx_id b = next_noindel(v, edge_vs);
  if (v == b)
    return;

  cycle_compact.push_back(max(b, v));
  edges.push_back(edge_vs);
  vector<Vtx_id> neigs = get_neighbors(b, is_random);
  for (Vtx_id g : neigs) {
    get_tuptok_cycle(start, g, t, k - 1, cycle_compact, edges, tuptok_cycles,
                     is_random);
  }
  cycle_compact.pop_back();
  edges.pop_back();
}

int CycleGraph::get_gray_degree(Vtx_id i) const {
  if (vertices[i].fix_gray != NO_EDGE) {
    return 1;
  } else {
    return vertices[i].grays.size();
  }
}

int CycleGraph::get_black_degree(Vtx_id i) const {
  return (vertices[i].black == NO_EDGE) ? 0 : 1;
}

int CycleGraph::get_red_degree(Vtx_id i) const {
  return (vertices[i].red == NO_EDGE) ? 0 : 1;
}

void CycleGraph::dcj(Vtx_id a, Vtx_id b, int x_noir, int y_noir, bool cross) {
  IR x = IR(x_noir, x_noir);
  IR y = IR(y_noir, y_noir);
  Vtx_id c = vertices[a].black;
  Vtx_id d = vertices[b].black;
  // magnitudes with sign (a previous dcj may leave a negative remainder, which
  // ir_abs would turn positive, changing the total weight of the cycle)
  IR x_ = ir_sub(add_sign_to_ir(vertices[a].weight, vertices[a].indel_update == 1), x);
  IR y_ = ir_sub(add_sign_to_ir(vertices[b].weight, vertices[b].indel_update == 1), y);

  if (cross) {
    vertices[a].black = b;
    vertices[a].weight =
        add_sign_to_ir(ir_sum(x, y), vertices[a].indel_update == 1);
    vertices[b].black = a;
    vertices[b].weight = vertices[a].weight;
    vertices[c].black = d;
    vertices[c].weight =
        add_sign_to_ir(ir_sum(x_, y_), vertices[a].indel_update == 1);
    vertices[d].black = c;
    vertices[d].weight = vertices[c].weight;
  } else {
    vertices[a].black = d;
    vertices[a].weight =
        add_sign_to_ir(ir_sum(x, y_), vertices[a].indel_update == 1);
    vertices[d].black = a;
    vertices[d].weight = vertices[a].weight;
    vertices[c].black = b;
    vertices[c].weight =
        add_sign_to_ir(ir_sum(x_, y), vertices[a].indel_update == 1);
    vertices[b].black = c;
    vertices[b].weight = vertices[c].weight;
  }
}

bool CycleGraph::split_trivial_cycle(Cycle cy) {
  if (is_trivial_cycle(cy)) {
    return false;
  }

  int weight_one = 0;
  int weight_two = 0;
  bool ok_one = true;  // first edge found
  bool ok_two = false; // second edge found
  bool ok = false;     // cycle can be split

  vector<Vtx_id> vs;
  Vtx_id i = next_noindel(cy.start_vtx, vs);
  Vtx_id v = i;
  Vtx_id a = 0;
  Vtx_id b = 0;
  do {
    v = vertices[v].black;
    if (ok_two)
      weight_two =
          add_sign_to_ir(vertices[v].weight, vertices[v].indel_update == 1).first;
    else if (ok_one)
      a = v;
    weight_one =
        add_sign_to_ir(vertices[v].weight, vertices[v].indel_update == 1).first;
    if (vertices[v].is_indel) {
      v = vertices[v].red;
      ok_one = false;
      ok_two = false;
    } else {
      v = vertices[v].fix_gray;
      if (ok_two) {
        if (weight_one >= weight_two) {
          ok = true;
          b = v;
          break;
        } else {
          weight_one = weight_two;
        }
      }
      if (ok_one) {
        ok_two = true;
      }
      ok_one = true;
    }
  } while (v != i);

  if (ok) {
    // c is in the other cycle after the dcj (black of a before the dcj)
    Vtx_id c = vertices[a].black;
    dcj(a, b, weight_two, 0, true);
    cycles.erase(find(cycles.begin(), cycles.end(), cy));
    cycles.push_back(Cycle(get_cycle(a), cycle_weight(a), cycle_potation(a)));
    cycles.push_back(Cycle(get_cycle(c), cycle_weight(c), cycle_potation(c)));
  }

  return ok;
}

void CycleGraph::fix_flexible_weights() {
  for (Cycle &cy : cycles) {
    // black edges of the cycle whose weights are summed by cycle_weight (the
    // ones reached by a gray edge), as in cycle_weight
    vector<Vtx_id> target;
    int origin = 0; // nucleotides of the origin edges of the cycle
    Vtx_id v = cy.start_vtx;
    do {
      v = vertices[v].black;
      if (vertices[v].is_indel) {
        v = vertices[v].red;
      } else {
        v = vertices[v].fix_gray;
        if (vertices[v].indel_update == 1) {
          target.push_back(v);
        } else {
          origin -= vertices[v].weight.first;
        }
      }
    } while (v != cy.start_vtx);
    // values inside the intervals with sum closest to origin
    int rest = origin;
    for (Vtx_id u : target) {
      rest -= vertices[u].weight.first;
    }
    for (Vtx_id u : target) {
      IR w = vertices[u].weight;
      int val = w.first + max(0, min(rest, w.second - w.first));
      rest -= val - w.first;
      vertices[u].weight = vertices[vertices[u].black].weight = IR(val, val);
    }
  }
  // fix also the target edges not counted in any cycle weight (after an indel)
  for (size_t u = 0; u < size(); ++u) {
    if (vertices[u].indel_update == 1) {
      IR w = vertices[u].weight;
      vertices[u].weight = IR(w.first, w.first);
    }
  }
  for (Cycle &cy : cycles) {
    cy.weight = cycle_weight(cy.start_vtx);
  }
}

int CycleGraph::dcj_indel_ir_distance() const {
  CycleGraph cg = *this;
  cg.fix_flexible_weights();
  int dist = 0;
  bool ok = true;

  set<Vtx_id> balanced;
  vector<Vtx_id> vs;

  do {
    ok = true;
    for (auto cy : cg.cycles) {
      if (!cg.is_only_indel(cy)) {
        if (!cg.is_trivial_cycle(cy)) {
          dist += 1;

          if (!cg.split_trivial_cycle(cy)) {
            Vtx_id a = cg.next_noindel(cy.start_vtx, vs);
            Vtx_id b = cg.vertices[a].fix_gray;
            b = cg.next_noindel(b, vs);
            b = cg.vertices[b].fix_gray;
            // c is in the other cycle after the dcj (black of a before it)
            Vtx_id c = cg.vertices[a].black;
            cg.dcj(a, b, 0, 0, true);
            cg.cycles.erase(find(cg.cycles.begin(), cg.cycles.end(), cy));
            cg.cycles.push_back(Cycle(cg.get_cycle(a), cg.cycle_weight(a),
                                      cg.cycle_potation(a)));
            cg.cycles.push_back(Cycle(cg.get_cycle(c), cg.cycle_weight(c),
                                      cg.cycle_potation(c)));
          }
          ok = false;
          break;
        } else if (balanced.find(cy.start_vtx) == balanced.end() &&
                   (!ir_intersect(cy.weight, ZERO_IR) || cy.pot != 0)) {
          bool has_ins =
              cy.weight.first > 0 || cy.pot > 1 ||
              (cy.pot == 1 && cg.cycle_run(cy.start_vtx).genome == 'G');
          bool has_del =
              cy.weight.second < 0 || cy.pot > 1 ||
              (cy.pot == 1 && cg.cycle_run(cy.start_vtx).genome == 'H');
          dist += 1;
          if (has_ins && has_del) {
            dist += 1;
          }
          ok = false;
          balanced.insert(cy.start_vtx);
          break;
        }
      }
    }
  } while (!ok);

  return dist;
}

int CycleGraph::dcj_ir_distance(bool use_ir_indel) const {
  int num_bal_cycles = this->balanced_cycles;
  int merges = 0;
  int splits = 0;

  if (!use_ir_indel) {
    vector<IR> unbalanced_weights;
    for (auto c : cycles) {
      if (!ir_intersect(c.weight, ZERO_IR)) {
        unbalanced_weights.push_back(c.weight);
      }
    }
    int num_unbal_cycles = unbalanced_weights.size();

    // Merge possible pairs of unbalanced cycles
    for (size_t i = 0; i < unbalanced_weights.size(); i++) {
      for (size_t j = i + 1; j < unbalanced_weights.size(); j++) {
        if (unbalanced_weights[i] != ZERO_IR &&
            unbalanced_weights[j] != ZERO_IR &&
            ir_intersect(ir_sum(unbalanced_weights[i], unbalanced_weights[j]),
                         ZERO_IR)) {
          unbalanced_weights[i] = ZERO_IR;
          unbalanced_weights[j] = ZERO_IR;
          merges += 1;
          num_bal_cycles += 1;
          num_unbal_cycles -= 2;
        }
      }
    }

    // 1/3-approximation to merge possible triples of unbalanced cycles
    for (int i = 0; i < unbalanced_weights.size(); i++) {
      for (int j = 0; j < unbalanced_weights.size(); j++) {
        if (j == i)
          continue;
        for (int k = j + 1; k < unbalanced_weights.size(); k++) {
          if (k == i)
            continue;
          if (unbalanced_weights[i] != ZERO_IR &&
              unbalanced_weights[j] != ZERO_IR &&
              unbalanced_weights[k] != ZERO_IR &&
              ir_intersect(
                  ir_sum(unbalanced_weights[i],
                         ir_sum(unbalanced_weights[j], unbalanced_weights[k])),
                  ZERO_IR)) {
            unbalanced_weights[i] = ZERO_IR;
            unbalanced_weights[j] = ZERO_IR;
            unbalanced_weights[k] = ZERO_IR;
            merges += 2;
            num_bal_cycles += 1;
            num_unbal_cycles -= 3;
          }
        }
      }
    }

    // All remaining unbalanced cycles are merged into one
    if (num_unbal_cycles > 0) {
      merges += num_unbal_cycles - 1;
      num_bal_cycles += 1;
    }
  }

  splits = (last_origin_vertex + 1) / 2 - num_bal_cycles;
  return merges + splits;
}
