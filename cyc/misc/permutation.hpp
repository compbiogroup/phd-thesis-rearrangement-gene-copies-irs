#pragma once
#include "../cycle/cycles.hpp"
#include "genome.hpp"
#include <random>
using namespace std;

class Permutation : public Genome {
  unique_ptr<vector<IR>> intergenic_regions_target;

private:
  IR target_ir(int i, bool use_abs) const;
  mt19937 &rand_gen;

public:
  Permutation(mt19937 &rand_gen) : rand_gen(rand_gen){};
  /* Copy constructor */
  Permutation(const Permutation &pi);
  /* Construct a permutation replicas are mapped randomly */
  Permutation(const Genome &g, const Genome &h, bool duplicate,
              mt19937 &rand_gen);
  /* Construct a permutation from the raw representation */
  Permutation(PermsIrs perms_irs, mt19937 &rand_gen)
      : Permutation(Genome(perms_irs.s, perms_irs.s_ir_low, perms_irs.s_ir_high, perms_irs.s_n),
                    Genome(perms_irs.p, perms_irs.p_ir_low, perms_irs.p_ir_high, perms_irs.p_n), false,
                    rand_gen) {}
  /* Separate iota from the permutation */
  GenomePair split_iota() const;
  /* Verify if the permutation is sorted and with the correct intergenic region
   * sizes */
  bool is_iota() const;
  /* Only read access to intergenic regions of target. */
  IR get_ir_target(int i) const { return (*intergenic_regions_target)[i - 1]; }
  /* Change an intergenic region of target. */
  void set_ir_target(int i, IR ir) { (*intergenic_regions_target)[i - 1] = ir; }
  /* Check breakpoints starting with the given gene */
  bool breakpoint(Gene, bool use_abs) const;
  bool hard_breakpoint(Gene, bool use_abs) const;
  bool soft_breakpoint(Gene, bool use_abs) const;
  bool overcharged_breakpoint(Gene, bool use_abs) const;
  bool undercharged_breakpoint(Gene, bool use_abs) const;
  bool is_block(Gene, Gene, bool use_abs) const;
  int target_size() const { return intergenic_regions_target->size() + 1; }
  virtual void serialize(ostream &) const override;
};
