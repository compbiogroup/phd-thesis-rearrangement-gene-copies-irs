#pragma once

#include "../misc/genome.hpp"
#include <bitset>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <unordered_map>

typedef int Vtx_id;
#define NO_EDGE -1

typedef enum { Head, Tail, Cap } VtxType;

inline IR add_sign_to_ir(IR ir, bool positive) {
  return (positive) ? ir : IR(-ir.second, -ir.first);
}

struct Vertex {
  Gene gene_val;              // value of the gene correspondent to this vertex
  bool sign_positive = false; // record the sign of the correspondent gene
  VtxType vtx_type; // if it is the head or tail vertex of the correspondent
                    // gene, or if it is a gene added during capping
  Vtx_id black = NO_EDGE;
  IR weight = ZERO_IR; // weight range of the incident black edge, negative
                       // values for origin genome and positive for target
  Vtx_id fix_gray = NO_EDGE;
  Vtx_id red = NO_EDGE;
  vector<Vtx_id> grays;
  bool in_cycle = false; // indicate if the vertex is in a cycle
  bool is_indel = false;
  bool is_indel_fix = false; // indicate if the vertex was marked indel by mark_indel_edges
  bool fix_adj = false;  // the incident black edge is a fixed adjacency
  bool fix_char = false; // the gray edge (or indel) was fixed by fix_partition
  int indel_update = 0; // how to update the indel count if adding the red edge
                        // incident to this vertex. Can be used to indicate if
                        // the correspondent gene is from G (-1) or H (1).
  Vertex(Gene gene_val, bool sign_positive, VtxType vtx_type, Vtx_id black,
         IR weight, bool ir_positive)
      : gene_val(gene_val), sign_positive(sign_positive), vtx_type(vtx_type),
        black(black), weight(add_sign_to_ir(weight, ir_positive)), grays() {}
  IR get_ir() { return ir_abs(weight); }
};

/* Raw representation of permutations to interact with other programs */
struct PermsIrs {
  size_t s_n, p_n;
  int *s;
  int *s_ir_low;
  int *s_ir_high;
  int *p;
  int *p_ir_low;
  int *p_ir_high;
};

struct CommonPartition {
  vector<int> breaks_s;
  vector<int> breaks_p;
};

/* Partition given by a FIX file (see read_partition_fixes), in positions of the
 * genomes without caps. */
struct PartitionFix {
  bool missing = false;           // instance without a partition (NA)
  vector<vector<int>> breaks_g;   // breakpoints (regions) of each chromosome
  vector<vector<int>> breaks_h;
  vector<int> corr; // position in the target of each gene of the origin
                    // (negative if reversed, 0 if exclusive)
};

struct QEntry {
  Vtx_id vtx;
  map<Vtx_id, Vtx_id> fixed;
  set<Vtx_id> vizited;
  map<Gene, int> indel_count;

  // Constructor for the first entry
  QEntry(Vtx_id start, map<Gene, int> indel_count)
      : fixed(), vizited(), indel_count(indel_count) {
    vtx = start;
  }

  // Constructor for entry with head-tail correspondence
  QEntry(Vtx_id v, Vtx_id u, Vtx_id alter_v, Vtx_id alter_u,
         map<Vtx_id, Vtx_id> fixed_old, set<Vtx_id> vizited_old,
         map<Gene, int> indel_count_old)
      : fixed(fixed_old), vizited(vizited_old), indel_count(indel_count_old) {
    vtx = u;
    fixed[v] = u;
    fixed[u] = v;
    fixed[alter_v] = alter_u;
    fixed[alter_u] = alter_v;
    vizited.insert(u);
    vizited.insert(v);
  }

  // Constructor for entry without head-tail correspondence
  QEntry(Vtx_id v, Vtx_id u, map<Vtx_id, Vtx_id> fixed_old,
         set<Vtx_id> vizited_old, map<Gene, int> indel_count_old)
      : fixed(fixed_old), vizited(vizited_old), indel_count(indel_count_old) {
    vtx = u;
    fixed[v] = u;
    fixed[u] = v;
    vizited.insert(u);
    vizited.insert(v);
  }

  bool not_yet_visited(const set<Vtx_id> &vizited_in_level, Vtx_id v) {
    return vizited_in_level.find(v) == vizited_in_level.end() &&
           vizited.find(v) == vizited.end();
  }
};

struct Run {
  Gene fst_gene; // gene after which the run is inserted
  // the run is inserted right after the left cap of the chromosome; fst_gene
  // is then the value of the cap, which is not a gene of the genome
  bool after_cap = false;
  char genome;   // G or H
  vector<Gene> genes_to_add;
};

struct Cycle {
  Gene start_vtx; // one of the vertices of the cycle
  size_t size;    // size of the cycle (number of black edges)
  IR weight;      // weight of the cycle (sum of weights of black edges)
  int pot;        // indel potential of the cycle

  Cycle(vector<Vtx_id> cycle, IR weight, int pot)
      : start_vtx(cycle[0]), size(cycle.size()), weight(weight), pot(pot) {}

  bool operator<(const Cycle &other) const {
    return this->size < other.size;
    // return this->size + this->pot < other.size + other.pot;
  }

  bool operator==(const Cycle &other) const {
    return this->start_vtx == other.start_vtx;
  }
};

class CycleGraph {

private:
  // Counters for types of cycles
  int balanced_cycles = 0;
  int clean_balanced_cycles = 0;
  int divergent_cycles = 0;
  int divergent_balanced_or_neg_cycles = 0;
  int trivial_balanced_cycles = 0;

  int indel_potation = 0;
  /* id of the last vertex correspondent to a gene in
   * the origin/target genome */
  int last_origin_vertex;
  int last_target_vertex;

  /* Whether to consider only trivial cycles in the gain and decomposition
   * functions */
  bool only_trivial = false;
  /* Whether there are fixed adjacencies or characters (fix_partition), the
   * heuristics then use the closures below when fixing gray edges or indels */
  bool use_fix = false;
  /* Number of fixed adjacencies released to complete the decomposition */
  int released_adj = 0;
  /* Value bigger then all the labels */
  int op_max;
  /* Map from position on the genomes to the vertex id */
  vector<Vertex> vertices;
  vector<Cycle> cycles; // List of cycles of a decomposition.
  /* Count the number of genes with positive sign for origin genome and negative
   * sign for target. The final map will indicate the number of each gene that
   * must be inserted or deleted. */
  map<Gene, int> indel_count;
  /* Each graph has its own generator, so copies can be decomposed in parallel
   * (mutable because const methods such as get_part use it). */
  mutable mt19937 rand_gen;

  vector<Vtx_id> produce_random_idxs(bool is_random);
  void follow_indel_edge(set<Vtx_id> &vizited_in_level,
                         vector<unique_ptr<QEntry>> &q,
                         unique_ptr<QEntry> &entry);
  /* Select a cycle with a bfs
   * Arguments:
   *     start - initial vertex
   *     is_random - whether to use random approach
   */
  /* Get string with the cycles from the decomposition */
  void bfs(Vtx_id start, bool is_random);
  // Get trivial balanced cycles
  void get_trivial_bal_cycle(Vtx_id start, bool is_random);
  void
  get_tuptok_cycle(Vtx_id start, Vtx_id v, int t, int k,
                   vector<Vtx_id> &cycle_compact, vector<vector<Vtx_id>> &edges,
                   vector<pair<vector<Vtx_id>, vector<Vtx_id>>> &tuptok_cycles,
                   bool is_random);
  vector<vector<int>> get_2uptok_cycle(Vtx_id start, int k) const;
  Vtx_id next_noindel(Vtx_id v, vector<Vtx_id> &vs) const;
  void mark_indel_edges(bool is_random);
  void find_1_cycles(bool is_random);
  void find_2_cycles(bool is_random, vector<vector<Vtx_id>> &cycles);
  /* Add greedily the clean (without indel edges) and balanced cycles with 2s
   * black edges, taking the black edges in random order. Each cycle is
   * found by joining paths with s gray edges from each side of the black
   * edge (meet in the middle). */
  void add_small_cycles(int s, bool is_random);
  bool clean_edge_free(Vtx_id x) const;
  bool gene_free(Vtx_id v) const;
  void back_half_paths(Vtx_id c, int r, bool is_random, vector<Vtx_id> &path,
                       unordered_map<Vtx_id, vector<vector<Vtx_id>>> &paths);
  bool forward_half_paths(
      Vtx_id b, int f, bool is_random, vector<Vtx_id> &path,
      const unordered_map<Vtx_id, vector<vector<Vtx_id>>> &back);
  /* Add to fix the gray edges that adding cycle would fix (including the
   * implied ones). Returns false if they conflict with fix or with the
   * current state of the graph. */
  bool cycle_fixes(const vector<Vtx_id> &cycle, map<Vtx_id, Vtx_id> &fix);
  bool cycles_conflict(const pair<vector<Vtx_id>, vector<Vtx_id>> &c1,
                       const pair<vector<Vtx_id>, vector<Vtx_id>> &c2);
  vector<vector<Vtx_id>>
  select_cycles(const vector<pair<vector<Vtx_id>, vector<Vtx_id>>> &cands);
  vector<Vtx_id> get_neighbors(Vtx_id v, bool is_random);
  bool is_only_indel(Cycle c) const;
  // Check trivial cycles, must have fixed gray edges
  bool is_trivial_cycle(Cycle c) const;
  // If possible uses a dcj to separate a trivial cycle from the given cycle
  // This operation messes up the signs of the genomes and the counters for
  // specific types of cycles, if we need it we have to fix them later.
  // TODO: Update the signs and counters during the operations
  bool split_trivial_cycle(Cycle c);
  /* Gray edge or red edge (indel) of x in the overlay fix (as in the QEntry of
   * bfs) or in the graph (NO_EDGE if x is free). */
  Vtx_id partner(Vtx_id x, const map<Vtx_id, Vtx_id> &fix,
                 const map<Vtx_id, Vtx_id> *base = nullptr) const;
  /* Add to fix the gray edge (v, u) and the ones it implies: the edge between
   * the other extremities of the genes and, for each fixed adjacency of v or u,
   * the edge between their black neighbors (the two black edges and the two
   * gray edges form a balanced trivial cycle, so the adjacency is not a
   * breakpoint). Returns false on a conflict with fix, base (read only,
   * consulted after fix) or the graph. */
  bool gray_closure(Vtx_id v, Vtx_id u, map<Vtx_id, Vtx_id> &fix,
                    const map<Vtx_id, Vtx_id> *base = nullptr) const;
  /* Vertices that must be indels if v is (the extremities of the genes of the
   * same exclusive block, following fixed adjacencies). Returns false if one of
   * them has a gray edge in fix or in the graph, or is a cap. */
  bool indel_closure(Vtx_id v, const map<Vtx_id, Vtx_id> &fix,
                     vector<Vtx_id> &out) const;
  /* Fix in the graph the closure of a gray edge or of an indel. */
  void commit_gray(Vtx_id v, Vtx_id u);
  void commit_indel(Vtx_id v, bool fixed = false);

public:
  /* Initial Constructor. */
  CycleGraph(const MultiGenome &origin, const MultiGenome &target,
             mt19937 &rand_gen, bool unsign = false, bool only_trivial = false);
  /* Copy Constructor. */
  CycleGraph(const CycleGraph &that)
      : balanced_cycles(that.balanced_cycles),
        clean_balanced_cycles(that.clean_balanced_cycles),
        divergent_cycles(that.divergent_cycles),
        divergent_balanced_or_neg_cycles(that.divergent_balanced_or_neg_cycles),
        trivial_balanced_cycles(that.trivial_balanced_cycles),
        indel_potation(that.indel_potation),
        last_origin_vertex(that.last_origin_vertex),
        last_target_vertex(that.last_target_vertex),
        only_trivial(that.only_trivial), use_fix(that.use_fix),
        released_adj(that.released_adj),
        op_max(that.op_max),
        vertices(that.vertices), cycles(that.cycles),
        indel_count(that.indel_count), rand_gen(that.rand_gen) {}
  /* A copy starts with the same generator state as the original, so reseed
   * copies that must be decomposed independently. */
  void reseed(mt19937::result_type seed) { rand_gen.seed(seed); }
  size_t size() const { return last_target_vertex + 1; };
  void decompose_1c_2c_3c(bool is_random);
  void decompose_random_map(bool is_random);
  void decompose_with_bfs(bool is_random);
  void decompose(bool is_random, string dec_algo = "bfs");
  /* Fix the adjacencies of the partition (the ones that are not breakpoints)
   * and, if fix_chars, the correspondence of genes. The graph must be built
   * from the genomes of the FIX file, with the linear chromosomes extended
   * (-e). */
  void fix_partition(const PartitionFix &fix, bool fix_chars);
  /* Fixed adjacencies released because the heuristic reached a dead end */
  int released_adjacencies() const { return released_adj; }
  string show_cycles() const;
  /* Recover cycles from string */
  void read_cycles(string);
  /* Recover the fixed gray edges and indels of a decomposition from the output
   * of serialize (as printed by dec). The graph must be built from the same
   * genomes. Only the assignment is recovered (enough for make_assignment and
   * get_genomes), not the list of cycles. */
  void read_decomposition(const string &);
  /* Build the list of cycles (and the counters of the gain) of a decomposition
   * recovered by read_decomposition. */
  void rebuild_cycles();
  /* Make assigment correspondent to the current decomposition by changing the
   * values of each vertex */
  /* indels are represented with a zero. */
  void make_assignment();
  /* Get genomes correspondent to the graph */
  MultiGenomePair get_genomes() const;
  /* Recover permutations from decomposition (only if we are using unichromosome
   * genomes) */
  PermsIrs get_perms() const;
  /* Recover common partition from decomposition */
  CommonPartition get_part(bool allow_rev) const;
  /* Get number of cycles from the decomposition */
  int dec_size() const { return cycles.size(); }
  void rem_cycle(Cycle c);
  /* Get indel potation */
  int potation() const { return indel_potation; }
  /* Gain function, represents the value of the cycle packing */
  vector<int> gain() const {
    vector<int> val(1);
    /* val[1] = this->dec_size() - this->potation(); */
    /* val[0] = this->dec_size() - this->potation(); */
    // val[0] = this->dec_size() + this->balanced_cycles;
    /* val[0] = this->dec_size() + this->divergent_balanced_or_neg_cycles; */
    if (only_trivial) {
      val[0] = this->trivial_balanced_cycles;
    } else {
      val[0] = this->clean_balanced_cycles;
    }
    return val;
  }

  /* Verify if a cycle can be added */
  bool check_cycle(vector<Vtx_id> cycle);
  void add_cycle(vector<Vtx_id> cycle);
  void check_and_add_cycle(vector<Vtx_id> cycle);

  /* Deal with cycles from a decomposition */
  vector<Cycle> cycle_list() const { return cycles; }
  vector<Vtx_id> get_cycle(Vtx_id i) const;
  int is_divergent(Vtx_id) const;
  IR cycle_weight(Vtx_id i) const;
  int cycle_potation(Vtx_id i) const;
  Run cycle_run(int i) const;

  /* Deal with cycles before a decomposition */
  vector<pair<vector<Vtx_id>, vector<Vtx_id>>>
  get_all_tuptok_cycles(int t, int k, bool is_random);

  void serialize(ostream &) const;
  // Get vertex info by index (starting at zero)
  int get_gray_degree(Vtx_id i) const;
  int get_black_degree(Vtx_id i) const;
  int get_red_degree(Vtx_id i) const;
  int is_indel(Vtx_id i) const { return vertices[i].is_indel; }
  VtxType get_type(Vtx_id i) const { return vertices[i].vtx_type; }
  int get_value(Vtx_id i) const { return vertices[i].gene_val; }
  Vtx_id get_correspondent(Vtx_id i) const { return vertices[i].fix_gray; }
  bool in_G(Vtx_id i) const { return vertices[i].indel_update == -1; }

  // Apply a DCJ operation on the cycle graph
  // a and b are vertices adjacent to the black edges to be cut,
  // x and y indicate the point of the cuts in the intergenic regions
  // (x nucleotides after a, following from a, and y nucleotides after b,
  // following from b), and cross indicates how the new black edges will be
  // formed (true: ab and cd, false: ad and cb).
  // Note that x and y are always positive.
  // This operation messes up the signs of the genomes and the counters for
  // specific types of cycles, if we need it we have to fix them later.
  // TODO: Update the signs and counters during the operations
  void dcj(Vtx_id a, Vtx_id b, int x, int y, bool cross);

  // Distance algorithms based on the cycle graph with gene indels
  // Requires that all gray edges are fixed
  /* Replace the flexible weights (intervals) of the target black edges by
   * fixed values inside their intervals, cycle by cycle, so that every
   * balanced cycle of the decomposition stays balanced (each target black edge
   * belongs to exactly one cycle). The unbalanced cycles get the values closest
   * to balance them. Requires a complete decomposition. */
  void fix_flexible_weights();
  int dcj_indel_ir_distance() const;

  // Distance algorithms based on the cycle graph
  // 2017-fertin-etal algorithm
  // if ir indel are considered, it balances cycles with indels instead of using
  // merges
  int dcj_ir_distance(bool use_ir_indel) const;
};

ostream &operator<<(std::ostream &os, const CycleGraph &cg);
ostream &operator<<(std::ostream &os, const Run &r);

inline CycleGraph *get_best_cg(CycleGraph *cg1, CycleGraph *cg2) {
  if (cg1 == nullptr) {
    return cg2;
  }
  if (cg2 == nullptr) {
    return cg1;
  }
  if (cg1->gain() > cg2->gain()) {
    delete cg2;
    return cg1;
  } else {
    delete cg1;
    return cg2;
  }
}
