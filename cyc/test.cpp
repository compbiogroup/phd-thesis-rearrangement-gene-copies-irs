#include <algorithm>
#include <functional>
#include <memory>
#include <vector>

#include "cycle/cycles.hpp"
#include "distance_algorithms/reversal4.hpp"
#include "distance_algorithms/reversal_transposition45.hpp"
#include "distance_algorithms/transposition3.hpp"
#include "external/external.hpp"
#include "heur/ga.hpp"
#include "misc/genome.hpp"
#include "misc/io.hpp"
#include "misc/permutation.hpp"
#include "misc/timer.hpp"
#include "quickcheck/quickcheck/quickcheck.hh"
using namespace quickcheck;

mt19937 rand_gen;

/* Property that reports an exception as a failure of the current test case
 * instead of aborting the whole test run. Subclasses implement safeHoldsFor. */
template <class T> class SafeProperty : public Property<T> {
  bool holdsFor(const T &x) override {
    try {
      return safeHoldsFor(x);
    } catch (const exception &e) {
      cout << endl << "exception: " << e.what() << endl;
      return false;
    }
  }
  virtual bool safeHoldsFor(const T &x) = 0;
};

struct TestPerm {
  unique_ptr<Permutation> perm;
  TestPerm() {}
  TestPerm(const TestPerm &perm_t) {
    perm = make_unique<Permutation>(*perm_t.perm);
  }
};

ostream &operator<<(std::ostream &os, const TestPerm &t_perm) {
  os << *t_perm.perm;
  return os;
}

struct TestCycleGraph {
  unique_ptr<CycleGraph> cg;
  unique_ptr<MultiGenome> g;
  unique_ptr<MultiGenome> h;

  TestCycleGraph() {}
  TestCycleGraph(const TestCycleGraph &cg_t) {
    cg = make_unique<CycleGraph>(*cg_t.cg);
    g = make_unique<MultiGenome>(*cg_t.g);
    h = make_unique<MultiGenome>(*cg_t.h);
  }
};

struct TestBreakGraph : public TestCycleGraph {
  TestBreakGraph() {}
  TestBreakGraph(const TestBreakGraph &bg_t) : TestCycleGraph(bg_t) {}
};

ostream &operator<<(std::ostream &os, const TestCycleGraph &t_cg) {
  os << *t_cg.cg;
  return os;
}

GenomePair generate_genomes(size_t n, bool add_sign, bool add_indels) {
  string sg, ig, sh, ih;
  int ir;
  int sum = 0;
  bool positive;
  GenomePair gs;

  ir = generateInRange(50, 100);
  ig.append(to_string(ir));
  ig.append(" ");
  sum += ir;
  for (size_t i = 0; i < n; i++) {
    ir = generateInRange(0, 100);
    generate(n, positive);
    if (!add_sign) {
      sg.append("1");
    } else if (positive) {
      sg.append("+1");
    } else {
      sg.append("-1");
    }
    sg.append(" ");
    ig.append(to_string(ir));
    ig.append(" ");
    generate(n, positive);
    if (!add_sign) {
      sh.append("1");
    } else if (positive) {
      sh.append("+1");
    } else {
      sh.append("-1");
    }
    sh.append(" ");
    sum += ir;
  }

  vector<int> ir_list;
  ir_list.push_back(0);
  for (size_t i = 0; i < n; i++) {
    ir_list.push_back(generateInRange(0, sum));
  }
  ir_list.push_back(sum);
  sort(ir_list.begin(), ir_list.end());
  for (size_t i = 1; i <= n + 1; i++) {
    ir = ir_list[i] - ir_list[i - 1];
    ih.append(to_string(ir));
    ih.append(" ");
  }

  gs.g = make_unique<Genome>(sg, ig, true);
  gs.h = make_unique<Genome>(sh, ih, true);

  if (add_indels) {
    for (int i = 2; i < int(gs.g->size()); i++) {
      if (rand_gen() % 100 < 10) {
        IR s = ir_sum(gs.g->get_ir(i - 1), gs.g->get_ir(i));
        IR x(rand_gen() % (s.first + 1), s.second);
        gs.g->deletion(i, i + 1, x);
        i--;
      }
    }

    for (int i = 2; i < int(gs.h->size()); i++) {
      if (rand_gen() % 100 < 10) {
        IR s = ir_sum(gs.h->get_ir(i - 1), gs.h->get_ir(i));
        IR x(rand_gen() % (s.first + 1), s.second);
        gs.h->deletion(i, i + 1, x);
        i--;
      }
    }
  }

  return gs;
}

void generate(size_t n, TestPerm &t_perm) {
  bool add_sign;
  generate(n, add_sign);
  GenomePair gs = generate_genomes(n, add_sign, false);
  t_perm.perm = make_unique<Permutation>(*gs.g, *gs.h, add_sign, rand_gen);
}

void generate_graphs(size_t n, TestCycleGraph &t_cg, bool add_sign,
                     bool only_perm, bool add_multichrom) {
  bool circular;

  GenomePair gs = generate_genomes(n, add_sign, true);
  if (only_perm) {
    Permutation perm(*gs.g, *gs.h, false, rand_gen);
    gs = perm.split_iota();
  }
  t_cg.g = make_unique<MultiGenome>(std::move(gs.g));
  t_cg.h = make_unique<MultiGenome>(std::move(gs.h));

  if (add_multichrom) {
    int full_n = t_cg.g->full_size();
    int cut;
    do {
      cut =
          rand_gen() %
          full_n; // 0 does not cut the genome and the last cut is in full_n - 1
      t_cg.g->break_genome(t_cg.g->n_chrom(), cut);
      full_n = full_n - cut + 1;
    } while (cut != 0 && full_n >= 0);
    for (size_t i = 1; i <= t_cg.g->n_chrom(); i++) {
      generate(n, circular);
      if (circular && t_cg.g->get_chrom(i).size() > 2) {
        t_cg.g->get_chrom(i).make_circular();
      }
    }

    full_n = t_cg.h->full_size();
    do {
      cut = rand_gen() % full_n;
      t_cg.h->break_genome(t_cg.h->n_chrom(), cut);
      full_n = full_n - cut + 1;
    } while (cut != 0 && full_n >= 0);
    for (size_t i = 1; i <= t_cg.h->n_chrom(); i++) {
      generate(n, circular);
      if (circular && t_cg.h->get_chrom(i).size() > 2) {
        t_cg.h->get_chrom(i).make_circular();
      }
    }
  }

  t_cg.cg = make_unique<CycleGraph>(*t_cg.g, *t_cg.h, rand_gen, !add_sign);
}

void generate(size_t n, TestCycleGraph &t_cg) {
  generate_graphs(n, t_cg, true, false, true);
}

void generate(size_t n, TestBreakGraph &t_cg) {
  generate_graphs(n, t_cg, false, true, false);
}

class PBounds3T : public SafeProperty<TestPerm> {
  bool safeHoldsFor(const TestPerm &t_perm) {
    Transposition3 alg;

    int breaks_count = 0;
    for (size_t b = 1; b <= t_perm.perm->size() - 1; b++) {
      if (t_perm.perm->breakpoint((*t_perm.perm)[b], false))
        breaks_count++;
    }

    int dist = alg.estimate_distance(*t_perm.perm, rand_gen);

    bool sucess = true;
    sucess = sucess && dist >= (breaks_count + 2) / 3;
    sucess = sucess && dist <= (4 * breaks_count / 3);
    if (!sucess) {
      cout << "pi: " << *t_perm.perm << endl;
      cout << "dist: " << dist << endl;
      cout << "b(pi):" << breaks_count << endl;
    } else {
      cout << ".";
      cout.flush();
    }
    return sucess;
  }
};

class PBounds4R : public SafeProperty<TestPerm> {
  bool safeHoldsFor(const TestPerm &t_perm) {
    Reversal4 alg;

    int breaks_count = 0;
    for (size_t b = 1; b <= t_perm.perm->size() - 1; b++) {
      if (t_perm.perm->breakpoint((*t_perm.perm)[b], true))
        breaks_count++;
    }

    int dist = alg.estimate_distance(*t_perm.perm, rand_gen);

    bool sucess = true;
    sucess = sucess && dist >= (breaks_count + 1) / 2;
    sucess = sucess && dist <= (4 * breaks_count / 2);
    if (!sucess) {
      cout << "pi: " << *t_perm.perm << endl;
      cout << "dist: " << dist << endl;
      cout << "b(pi):" << breaks_count << endl;
    } else {
      cout << ".";
      cout.flush();
    }
    return sucess;
  }
};

class PBounds45RT : public SafeProperty<TestPerm> {
  bool safeHoldsFor(const TestPerm &t_perm) {
    ReversalTransposition45 alg;

    int breaks_count = 0;
    for (size_t b = 1; b <= t_perm.perm->size() - 1; b++) {
      if (t_perm.perm->breakpoint((*t_perm.perm)[b], true))
        breaks_count++;
    }

    int dist = alg.estimate_distance(*t_perm.perm, rand_gen);

    bool sucess = true;
    sucess = sucess && dist >= (breaks_count + 2) / 3;
    sucess = sucess && dist <= (4.5 * breaks_count / 3);
    if (!sucess) {
      cout << "pi: " << *t_perm.perm << endl;
      cout << "dist: " << dist << endl;
      cout << "b(pi):" << breaks_count << endl;
    } else {
      cout << ".";
      cout.flush();
    }
    return sucess;
  }
};

bool checkGraphDegrees_generic(const TestCycleGraph &t_cg, string dec_algo,
                               bool with_sign) {
  bool sucess = true;
  for (size_t i = 0; i < t_cg.cg->size(); i++) {
    if (t_cg.cg->get_type(i) == Cap) {
      if (t_cg.g->n_chrom() == 1 && t_cg.h->n_chrom() == 1) {
        if (t_cg.cg->get_gray_degree(i) != 1) {
          sucess = false;
          cout << i << " gray_degree: " << t_cg.cg->get_gray_degree(i)
               << ", expected: 1" << endl;
        }
      } else {
        if (t_cg.cg->get_gray_degree(i) !=
            2 * (max(t_cg.g->num_lin(), t_cg.h->num_lin()))) {
          sucess = false;
          cout << i << " gray_degree: " << t_cg.cg->get_gray_degree(i)
               << ", num_lin_max: "
               << 2 * (max(t_cg.g->num_lin(), t_cg.h->num_lin())) << endl;
        }
      }
      if (t_cg.cg->get_red_degree(i) != 0) {
        sucess = false;
        cout << i << " red_degree: " << t_cg.cg->get_red_degree(i)
             << ", expected: 0" << endl;
      }
    } else {
      int mult = (with_sign) ? 1 : 2;
      if (t_cg.cg->in_G(i)) {
        if (t_cg.cg->get_gray_degree(i) !=
            mult * t_cg.h->occ(t_cg.cg->get_value(i))) {
          sucess = false;
          cout << i << " gray_degree: " << t_cg.cg->get_gray_degree(i)
               << ", occ: " << t_cg.h->occ(t_cg.cg->get_value(i)) << endl;
        }
      } else {
        if (t_cg.cg->get_gray_degree(i) !=
            mult * t_cg.g->occ(t_cg.cg->get_value(i))) {
          sucess = false;
          cout << i << " gray_degree: " << t_cg.cg->get_gray_degree(i)
               << ", occ: " << t_cg.g->occ(t_cg.cg->get_value(i)) << endl;
        }
      }
      if (t_cg.cg->get_red_degree(i) != 1) {
        sucess = false;
        cout << i << " red_degree: " << t_cg.cg->get_red_degree(i)
             << ", expected: 1" << endl;
      }
    }
    if (t_cg.cg->get_black_degree(i) != 1) {
      sucess = false;
      cout << i << " black_degree: " << t_cg.cg->get_black_degree(i) << endl;
    }
  }

  if (!sucess) {
    cout << "g: " << *t_cg.g << endl;
    cout << "h: " << *t_cg.h << endl;
    cout << "cg: " << *t_cg.cg << endl;
  } else {
    cout << ".";
    cout.flush();
  }

  // decompose a copy, so the input printed on failure is the original graph
  CycleGraph cg(*t_cg.cg);
  cg.decompose(true, dec_algo);

  bool sucess_d = true;
  for (size_t i = 0; i < cg.size(); i++) {
    if (cg.is_indel(i)) {
      if (cg.get_type(i) == Cap) {
        sucess_d = false;
        cout << i << " cap marked for deletion" << endl;
      } else {
        if (cg.get_red_degree(i) != 1 || cg.get_gray_degree(i) == 1) {
          sucess_d = false;
          cout << i << " red_degree: " << cg.get_red_degree(i)
               << ", expected 1" << endl;
          cout << i << " gray_degree: " << cg.get_gray_degree(i)
               << ", expected != 1" << endl;
        }
      }
    } else {
      if (cg.get_gray_degree(i) != 1) {
        sucess_d = false;
        cout << i << " gray_degree: " << cg.get_gray_degree(i)
             << ", expected 1" << endl;
      }
    }
    if (cg.get_black_degree(i) != 1) {
      sucess_d = false;
      cout << "black_degree: " << cg.get_black_degree(i) << ", expected 1"
           << endl;
    }
  }

  if (!sucess_d) {
    cout << "g: " << *t_cg.g << endl;
    cout << "h: " << *t_cg.h << endl;
    cout << "cg_dec: " << cg << endl;
  } else {
    cout << ".";
    cout.flush();
  }
  sucess = sucess && sucess_d;

  return sucess;
}

class CheckGraphDegreesMap : public SafeProperty<TestCycleGraph> {
  bool safeHoldsFor(const TestCycleGraph &t_cg) {
    return checkGraphDegrees_generic(t_cg, "map", true);
  }
};

class CheckGraphDegreesBFS : public SafeProperty<TestCycleGraph> {
  bool safeHoldsFor(const TestCycleGraph &t_cg) {
    return checkGraphDegrees_generic(t_cg, "bfs", true);
  }
};

class CheckGraphDegrees1C2C3C : public SafeProperty<TestBreakGraph> {
  bool safeHoldsFor(const TestBreakGraph &t_cg) {
    return checkGraphDegrees_generic(t_cg, "1c2c3c", false);
  }
};

bool verify_correspondence(const CycleGraph &cg, bool with_sign) {
  bool sucess = true;
  for (size_t i = 0; i < cg.size(); i++) {
    if (cg.is_indel(i)) {
      continue;
    }

    int j = cg.get_correspondent(i);
    if (j == NO_EDGE) {
      sucess = false;
      cout << "no correspondent for " << i << endl;
      continue;
    }
    int i_2 = cg.get_correspondent(j);
    if (i_2 == NO_EDGE) {
      sucess = false;
      cout << "no correspondent for " << j << " looking from " << i << endl;
    }
    if (i_2 != int(i)) {
      sucess = false;
      cout << "correspondence is not symetric for " << i << endl;
    }

    if (cg.get_value(i) != cg.get_value(j)) {
      sucess = false;
      cout << "different values for " << i << " and " << j << endl;
    }
    if (with_sign && cg.get_type(i) != cg.get_type(j)) {
      sucess = false;
      cout << "different types for " << i << " and " << j << endl;
    }
  }

  /* TODO: add tests for balancing of indels */

  if (!sucess) {
    cout << "cg_dec: " << cg << endl;
  } else {
    cout << ".";
    cout.flush();
  }

  return sucess;
}

class CheckGraphDecCorrespondenceMap : public SafeProperty<TestCycleGraph> {
  bool safeHoldsFor(const TestCycleGraph &t_cg) {
    CycleGraph cg(*t_cg.cg);
    cg.decompose_random_map(true);
    return verify_correspondence(cg, true);
  }
};

class CheckGraphDecCorrespondenceBFS : public SafeProperty<TestCycleGraph> {
  bool safeHoldsFor(const TestCycleGraph &t_cg) {
    CycleGraph cg(*t_cg.cg);
    cg.decompose_with_bfs(true);
    return verify_correspondence(cg, true);
  }
};

class CheckGraphDecCorrespondence1C2C3C : public SafeProperty<TestBreakGraph> {
  bool safeHoldsFor(const TestBreakGraph &t_cg) {
    CycleGraph cg(*t_cg.cg);
    cg.decompose_1c_2c_3c(true);
    return verify_correspondence(cg, false);
  }
};

class CheckGraphDecCorrespondenceGAMap : public SafeProperty<TestCycleGraph> {
  bool safeHoldsFor(const TestCycleGraph &t_cg) {
    Timer timer;

    GA ga = GA(make_unique<Chromossome>(*t_cg.cg), "map", 0.5, 0.5, 2, 2, 2, 3,
               nullptr, timer, rand_gen);
    ga.solve(timer);

    return verify_correspondence(*ga.get_best_chr(), true);
  }
};

class CheckGraphDecCorrespondenceGABFS : public SafeProperty<TestCycleGraph> {
  bool safeHoldsFor(const TestCycleGraph &t_cg) {
    Timer timer;

    GA ga = GA(make_unique<Chromossome>(*t_cg.cg), "bfs", 0.5, 0.5, 2, 2, 2, 3,
               nullptr, timer, rand_gen);
    ga.solve(timer);

    return verify_correspondence(*ga.get_best_chr(), true);
  }
};

class CheckGraphDecCorrespondenceGA1C2C3C : public SafeProperty<TestBreakGraph> {
  bool safeHoldsFor(const TestBreakGraph &t_cg) {
    Timer timer;

    GA ga = GA(make_unique<Chromossome>(*t_cg.cg), "1c2c3c", 0.5, 0.5, 2, 2, 2,
               3, nullptr, timer, rand_gen);
    ga.solve(timer);

    return verify_correspondence(*ga.get_best_chr(), false);
  }
};

bool ga_bounds_generic(const TestCycleGraph &t_cg, string dec_algo) {
  Timer timer;

  GA ga = GA(make_unique<Chromossome>(*t_cg.cg), dec_algo, 0.5, 0.5, 2, 2, 2, 3,
             nullptr, timer, rand_gen);
  ga.solve(timer);
  vector<int> obj = ga.get_best_obj();
  unique_ptr<Chromossome> best = ga.get_best_chr();

  // obj[0] counts (balanced) cycles of the decomposition. Cycles may have only
  // 2 vertices (indels, circular chromosomes with one gene), so the number of
  // cycles is at most size / 2.
  bool sucess = true;
  sucess = sucess && obj[0] >= 0;
  sucess = sucess && obj[0] <= best->dec_size();
  sucess = sucess && best->dec_size() <= int(t_cg.cg->size()) / 2;
  if (!sucess) {
    cout << "cg: " << *t_cg.cg << endl;
    cout << "cg_size: " << t_cg.cg->size() << endl;
    cout << "dec_size: " << best->dec_size() << endl;
    cout << "obj: ";
    print_vec(cout, obj);
  } else {
    cout << ".";
    cout.flush();
  }
  return sucess;
}

class GABoundsMap : public SafeProperty<TestCycleGraph> {
  bool safeHoldsFor(const TestCycleGraph &t_cg) {
    return ga_bounds_generic(t_cg, "map");
  }
};

class GABoundsBFS : public SafeProperty<TestCycleGraph> {
  bool safeHoldsFor(const TestCycleGraph &t_cg) {
    return ga_bounds_generic(t_cg, "bfs");
  }
};

class GABounds1C2C3C : public SafeProperty<TestBreakGraph> {
  bool safeHoldsFor(const TestBreakGraph &t_cg) {
    return ga_bounds_generic(t_cg, "1c2c3c");
  }
};

/* Checks that can be selected from the command line. Keep the names in sync
 * with the TEST_CHECKS list in CMakeLists.txt. */
const vector<pair<string, function<bool()>>> checks = {
    {"PBounds3T",
     [] {
       return check<PBounds3T>("upper and lower bounds hold for factor 3 "
                               "algorithm for transposition");
     }},
    {"PBounds4R",
     [] {
       return check<PBounds4R>(
           "upper and lower bounds hold for factor 4 algorithm for reversal");
     }},
    {"PBounds45RT",
     [] {
       return check<PBounds45RT>("upper and lower bounds hold for factor 4.5 "
                                 "algorithm for reversal and transposition");
     }},
    {"CheckGraphDegreesMap",
     [] {
       return check<CheckGraphDegreesMap>(
           "check degree of cycle graphs before and after decomposition with "
           "random map");
     }},
    {"CheckGraphDecCorrespondenceMap",
     [] {
       return check<CheckGraphDecCorrespondenceMap>(
           "check if decomposition produces correct correspondences with "
           "random map");
     }},
    {"CheckGraphDegreesBFS",
     [] {
       return check<CheckGraphDegreesBFS>("check degree of cycle graphs "
                                          "before and after decomposition "
                                          "with BFS");
     }},
    {"CheckGraphDecCorrespondenceBFS",
     [] {
       return check<CheckGraphDecCorrespondenceBFS>(
           "check if decomposition produces correct correspondences with BFS");
     }},
    {"CheckGraphDegrees1C2C3C",
     [] {
       return check<CheckGraphDegrees1C2C3C>(
           "check degree of cycle graphs before and after decomposition with "
           "1C2C3C");
     }},
    {"CheckGraphDecCorrespondence1C2C3C",
     [] {
       return check<CheckGraphDecCorrespondence1C2C3C>(
           "check if decomposition produces correct correspondences with "
           "1C2C3C");
     }},
    {"CheckGraphDecCorrespondenceGAMap",
     [] {
       return check<CheckGraphDecCorrespondenceGAMap>(
           "check if decomposition from GA produces correct correspondences "
           "with random map");
     }},
    {"CheckGraphDecCorrespondenceGABFS",
     [] {
       return check<CheckGraphDecCorrespondenceGABFS>(
           "check if decomposition from GA produces correct correspondences "
           "with BFS");
     }},
    {"CheckGraphDecCorrespondenceGA1C2C3C",
     [] {
       return check<CheckGraphDecCorrespondenceGA1C2C3C>(
           "check if decomposition from GA produces correct correspondences "
           "with 1C2C3C");
     }},
    {"GABoundsMap",
     [] {
       return check<GABoundsMap>("upper and lower bounds hold for cycle "
                                 "decomposition with ga and random map");
     }},
    {"GABoundsBFS",
     [] {
       return check<GABoundsBFS>("upper and lower bounds hold for cycle "
                                 "decomposition with ga and BFS");
     }},
    {"GABounds1C2C3C",
     [] {
       return check<GABounds1C2C3C>("upper and lower bounds hold for cycle "
                                    "decomposition with ga and 1C2C3C");
     }},
};

/* usage: dist_test [CHECK|all [SEED]] */
int main(int argc, char *argv[]) {
  string which = (argc > 1) ? argv[1] : "all";
  // fixed default seed so that failures can be reproduced
  unsigned long seed = (argc > 2) ? stoul(argv[2]) : 1031129599;
  cout << "seed: " << seed << endl;
  rand_gen.seed(seed);

  bool found = false, ok = true;
  for (const auto &c : checks) {
    if (which == "all" || which == c.first) {
      found = true;
      ok = c.second() && ok;
    }
  }
  if (!found) {
    cerr << "Unknown check: " << which << endl << "Available checks:";
    for (const auto &c : checks) {
      cerr << " " << c.first;
    }
    cerr << endl;
    return EXIT_FAILURE;
  }
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
