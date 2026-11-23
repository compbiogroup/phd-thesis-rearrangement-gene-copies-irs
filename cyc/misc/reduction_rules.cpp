#include "reduction_rules.hpp"
#include "genome.hpp"
#include <map>

template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

// Match two correspondent intervals a^l...a^r from g and b^l...b^r from h.
void suboptimal_rule_interval(MultiGenome &g, MultiGenome &h) {
  map<pair<int, int>, bool> original_singletons_pos;
  for (int c = 1; c <= int(g.n_chrom()); c++) {
    for (int i = 1; i <= g[c].size(); i++) {
      // alphas (genes to be deleted or inserted) are never anchors
      if (!g[c].check_alpha(i) && g.occ(g[c][i]) == 1 &&
          h.occ(g[c][i]) == 1) {
        original_singletons_pos[pair<int, int>(c, i)] = true;
      } else {
        original_singletons_pos[pair<int, int>(c, i)] = false;
      }
    }
  }

  for (int cg = 1; cg <= int(g.n_chrom()); cg++) {
    for (int i = 1; i <= g[cg].size(); i++) {
      int idx_al = -1, idx_ar = -1; // index from the left and the right of the
                                    // interval a^l...a^r in g
      int idx_bl = -1, idx_br = -1; // index from the left and the right of the
                                    // interval b^l...b^r in h

      // a^l is a singleton in both genomes
      if (original_singletons_pos[pair<int, int>(cg, i)]) {
        pair<int, int> chj =
            h.pos(g[cg][i]).back(); // get correspondent of g[i] in h
        int ch = chj.first;
        int j = chj.second;
        bool dir = false, rev = false;
        // pos() matches on the absolute value, so check the orientation of the
        // anchor: a^l = b^l in the direct and middle cases, a^l = -b^l in the
        // reverse case.
        bool same_sign = g[cg][i] == h[ch][j];

        // check if the interval is a direct match a^l...a^r = b^l...b^r
        idx_al = i;
        idx_bl = j;
        for (int k = 1; same_sign && i + k <= int(g[cg].size()) &&
                        j + k <= int(h[ch].size());
             k++) {
          // I am comparing the intergenic region before the match
          if (g[cg][i + k] != h[ch][j + k] ||
              g[cg].get_ir(i + k - 1) != h[ch].get_ir(j + k - 1)) {
            break;
          }
          if (original_singletons_pos[pair<int, int>(cg, i + k)]) {
            if (k > 1) {
              idx_ar = i + k;
              idx_br = j + k;
              dir = true;
            }
            break;
          }
        }

        // check if the interval is a reverse match a^l...a^r = -b^r...-b^l
        if (!dir && !same_sign) {
          idx_al = i;
          idx_br = j;
          for (int k = 1; i + k <= int(g[cg].size()) && j - k >= 1; k++) {
            // I am comparing the intergenic region before the match
            if (g[cg][i + k] != -h[ch][j - k] ||
                g[cg].get_ir(i + k - 1) != h[ch].get_ir(j - k)) {
              break;
            }
            if (original_singletons_pos[pair<int, int>(cg, i + k)]) {
              if (k > 1) {
                idx_ar = i + k;
                idx_bl = j - k;
                rev = true;
              }
              break;
            }
          }
        }

        // check if the middle is a reverse match a^la^{l+1}...a^{r-1}a^r =
        // b^l-b^{r-1}...-b^{l+1}b^r
        if (!dir && !rev && same_sign) {
          idx_al = i;
          idx_ar = -1;
          for (int k = i + 1; k <= int(g[cg].size()); k++) {
            if (original_singletons_pos[pair<int, int>(cg, k)]) {
              idx_ar = k;
              break;
            }
          }

          if (idx_ar != -1) {
            idx_bl = j;
            for (int k = 1; idx_ar - k >= 1 && j + k <= int(h[ch].size());
                 k++) {
              if (idx_ar - k == idx_al) {
                // check if j+k is idx_br (ie. g[idx_ar] equals to h[j+k])
                // and check last ir in the middle.
                if (g[cg][idx_ar] == h[ch][j + k] &&
                    g[cg].get_ir(idx_al) == h[ch].get_ir(j + k - 1) && k > 1) {
                  idx_br = j + k;
                  rev = true;
                }
                break;
              }
              // I am comparing the intergenic region before the match
              if (g[cg][idx_ar - k] != -h[ch][j + k] ||
                  g[cg].get_ir(idx_ar - k) != h[ch].get_ir(j + k - 1)) {
                break;
              }
            }
          }
        }

        if (dir || rev) {
          for (int k = 1; idx_al + k < idx_ar; k++) {
            // Apply new assignment
            Gene new_gene = max(g.get_op_max(), h.get_op_max());
            g.replace_label(cg, idx_al + k, new_gene);
            if (dir) {
              h.replace_label(ch, idx_bl + k, new_gene);
            } else {
              h.replace_label(ch, idx_br - k, new_gene);
            }
          }
          i = idx_ar - 1;
        }
      }
    }
  }
}

// Match pair of singleton and replica
void suboptimal_rule_pairs(MultiGenome &g, MultiGenome &h) {
  map<pair<int, int>, bool> original_singletons_pos;
  for (int c = 1; c <= int(g.n_chrom()); c++) {
    for (int i = 1; i <= g[c].size(); i++) {
      // alphas (genes to be deleted or inserted) are never anchors
      if (!g[c].check_alpha(i) && g.occ(g[c][i]) == 1 &&
          h.occ(g[c][i]) == 1) {
        original_singletons_pos[pair<int, int>(c, i)] = true;
      } else {
        original_singletons_pos[pair<int, int>(c, i)] = false;
      }
    }
  }

  for (int cg = 1; cg <= int(g.n_chrom()); cg++) {
    for (int i = 1; i < int(g[cg].size()); i++) {
      int idx1 = -1, idx2 = -1;
      Gene ah = g[cg][i], at = g[cg][i + 1]; // duo a^h a^t
      int ch;

      // a^h is a singleton in both genomes and a^t is a replica
      if (original_singletons_pos[pair<int,int>(cg,i)] && g.occ(at) > 1) {
        pair<int,int> chj = h.pos(ah).back();
        ch = chj.first;
        int b_idx = chj.second;
        Gene b = h[ch][b_idx];
        idx1 = i + 1;
        // a^t matches pred of b
        if (b_idx > 1 && sgn(ah) != sgn(b) && abs(at) == abs(h[ch][b_idx - 1]) &&
            sgn(at) != sgn(h[ch][b_idx - 1]) &&
            g[cg].get_ir(i) == h[ch].get_ir(b_idx - 1)) {
          idx2 = b_idx - 1;
          // a^t matches succ of b
        } else if (b_idx < int(h[ch].size()) && sgn(ah) == sgn(b) &&
                   abs(at) == abs(h[ch][b_idx + 1]) &&
                   sgn(at) == sgn(h[ch][b_idx + 1]) &&
                   g[cg].get_ir(i) == h[ch].get_ir(b_idx)) {
          idx2 = b_idx + 1;
        }
        // a^t is a singleton in both genomes and a^h is a replica
      } else if (original_singletons_pos[pair<int,int>(cg,i + 1)] && g.occ(ah) > 1) {
        pair<int,int> chj = h.pos(at).back();
        ch = chj.first;
        int b_idx = chj.second;
        Gene b = h[ch][b_idx];
        idx1 = i;
        // a^h matches pred of b
        if (b_idx > 1 && sgn(at) == sgn(b) && abs(ah) == abs(h[ch][b_idx - 1]) &&
            sgn(ah) == sgn(h[ch][b_idx - 1]) &&
            g[cg].get_ir(i) == h[ch].get_ir(b_idx - 1)) {
          idx2 = b_idx - 1;
          // a^h matches succ of b
        } else if (b_idx < int(h[ch].size()) && sgn(at) != sgn(b) &&
                   abs(ah) == abs(h[ch][b_idx + 1]) &&
                   sgn(ah) != sgn(h[ch][b_idx + 1]) &&
                   g[cg].get_ir(i) == h[ch].get_ir(b_idx)) {
          idx2 = b_idx + 1;
        }
      }

      if (idx1 != -1 && idx2 != -1) {
        // Skip the next pair so that the new assignment is not taken into
        // account
        if (original_singletons_pos[pair<int,int>(cg,i)]) {
          i++;
        }
        // Apply new assignment
        Gene new_gene = max(g.get_op_max(), h.get_op_max());
        g.replace_label(cg, idx1, new_gene);
        h.replace_label(ch, idx2, new_gene);
      }
    }
  }
}
