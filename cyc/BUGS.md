# Known bugs — cyc

Review of the working tree as of 2026-10-08 (uncommitted changes on top of `18e7957`).
`GRIMM/` (vendored) was not reviewed. Line numbers are approximate; search for the
quoted code if they have drifted.

The build is `Release` by default, so `NDEBUG` is set and **every `assert` is compiled
out** — several items below are "guarded" only by an assert.

Priority: items marked ⚠️ affect correctness of reported distances or crash on normal input.

---

## 1. Crashes and hangs

- [X] ⚠️ **`dist_test` does not compile** — `test.cpp:111`, `:119` use `IR + IR` and `rand_gen() % IR`.
  ```cpp
  IR s = ir_sum(gs.g->get_ir(i - 1), gs.g->get_ir(i));
  IR x(rand_gen() % (s.first + 1), s.second);
  ```
  (same at 119 with `gs.h`). Must still satisfy the asserts in `Genome::deletion`.

- [X] ⚠️ **`TranspositionNOIR::dist_aux` always throws** — `r_or_rt_noir.cpp:225` compares
  `fscanf`'s return with `EXIT_SUCCESS` (0); `fscanf` returns 1 on success. Also `pclose`
  is skipped by `||` short-circuit when the read fails (pipe leak/zombie) — same in
  `external/external.cpp:27`. Fix both:
  ```cpp
  int ret = fscanf(prog_pipe, "%d", &dist);
  int st = pclose(prog_pipe);
  if (ret != 1 || st != 0) throw runtime_error("Error on child process.");
  ```

- [X] **`reversal4.cpp`: infinite loop when case 2/3 finds no `k`** (~lines 50-71). The
  `for (auto c : breaks)` loops may finish without a move; `pi` is unchanged and
  `while (true)` repeats forever. Track whether a move happened and
  `throw runtime_error("Reversal4: no progress")` otherwise.

- [X] **`reversal_transposition45.cpp` `case_ii`/`case_iii`: dereference of `end()`**
  **Reopened 2026-10-19:** `swap_stage = false` is still set when `case_ii`/`case_iii` return 0 (~line 309), so the same breakpoints repeat forever in Release (the only guard is an `assert`).
  **Reopened 2026-10-14:** `case_ii` was fixed but `case_iii` (~line 61) still loops `breaks.begin()` → `overcharged_breaks.end()`.
  (lines 13-16, 55-58). The loop compares an iterator over `breaks` with
  `overcharged_breaks.end()` (UB: different containers); when `overcharged_breaks` is
  empty, `k_it = overcharged_breaks.begin()` (== `end()`) is dereferenced. Fix:
  ```cpp
  vector<int> cand(breaks.begin(), breaks.end());
  cand.insert(cand.end(), overcharged_breaks.begin(), overcharged_breaks.end());
  for (int g : cand) { k = pi.pos(g)[0]; ... }
  ```
  Also: if these return 0 (no `k` found), stage 2 sets `swap_stage = false` with no change
  → possible infinite loop (lines 290-301). Only set `swap_stage = false` when the case
  returned > 0, otherwise throw.

- [ ] **`transposition3.cpp`: `d` used uninitialized**
  **Reopened 2026-10-19:** only the site at ~line 92 was fixed; lines ~136 and ~180 still declare `int d, l = -1;` and rely on `assert(l != -1)`, and `softs[1]` (~line 289) is still indexed without checking `softs.size() >= 2`. (lines 91/131/174 `int d, l = -1;`
  then `assert(l != -1)`, compiler warns at 113/151/193). In Release, `l == -1` reaches
  `pi.transposition(...)` with garbage `d`. Initialize `d = -1` and handle `l == -1`
  (fall back to the plain transposition, or throw).
  Also `over[0]`/`over[1]` (46-50, 226-287) and `softs[1]` (279-291) are indexed without
  size checks → `if (over.size() < 2) throw ...`, `if (softs.size() < 2) throw ...`.

- [X] **`Genome()` leaves `circular` uninitialized** — `misc/genome.hpp:48`. Every
  `Permutation` constructor (except the copy ctor) goes through it; `serialize`,
  `is_linear()` then read garbage. Fix: `Genome() : circular(false), op_max(1) {}`.

- [ ] **Empty / one-gene genomes not rejected** — `genome.cpp:57-58` writes `(*genes)[0]`
  on an empty vector; several `size() - 1` expressions on `size_t` underflow
  (`serialize`, `permutation.cpp:30`, `io.cpp` output). After parsing:
  `if (!circular && genes->size() < 2) throw invalid_argument("linear genome needs >= 2 genes");`

- [ ] **Mains only catch `invalid_argument`**
  **Reopened 2026-10-14:** `dec`/`part` are compiled with OpenMP now, so an exception inside the parallel loops (`ga.hpp:62/88`, `heur/heur.cpp`) calls `std::terminate` (e.g. `dec ga foo`). Validate DEC_ALGO before the loop and/or catch inside the loop body and rethrow after it. (fixed: `run_main` in `misc/main_common.cpp` catches `exception`; the in-loop catch for OpenMP is still needed if `dist` is built with OpenMP) (`main_dist.cpp:203`, `main_dec.cpp:293`,
  `main_part.cpp:287`, `main_toperm.cpp:122`). `runtime_error` from `decompose` (unknown
  algorithm), `external.cpp`, `reversal4`, `rt45` aborts with no message. Catch
  `const exception &`. Once OpenMP is enabled for `dist` (see 4.1), exceptions must also be
  caught *inside* the parallel loop body, or the program calls `std::terminate`.

## 2. Wrong results — interval IRs

- [X] ⚠️ **`ir_abs` inverts origin intervals** — `misc/genome.hpp:27`. Origin weights are
  now stored as `(-hi, -lo)`; `ir_abs` returns `(hi, lo)` (first > second). Affects
  `Vertex::get_ir()` (`cycles.hpp:40`), `dcj()` remainders (`cycles.cpp:1363-1364`) and
  `split_trivial_cycle` (`cycles.cpp:1410-1421`, which reads `.first` = `hi` for origin
  edges, `lo` for target edges). Fix — undo the sign the same way it was applied
  (negation is its own inverse):
  ```cpp
  IR get_ir() { return add_sign_to_ir(weight, indel_update != -1); }
  ```
  and use the same expression in `dcj`. Then `split_trivial_cycle` should compare/cut on
  the lower bound (`.first`) on both sides.

- [ ] ⚠️ **Lexicographic `pair` comparisons used as interval comparisons.** `std::pair`'s
  `<`, `<=`, `>=`, `min`, `!=` compare lexicographically, not per component or as
  intervals. Results only match the old `int` code when `lo == hi` everywhere.
  - `distance_algorithms/auxiliar.cpp`: 13-34 (`<=` in `is_connected`), 45-69 (`min` in
    `transposition_move`), 81-97 (`>=` in `case_i`)
  - `distance_algorithms/reversal_transposition45.cpp`: 18, 32, 60 (`>=`)
  - `distance_algorithms/transposition3.cpp`: 73 (`>=`), 240
  - `misc/permutation.cpp`: `is_iota` (`!=`), `breakpoint`/`hard_breakpoint` (`!=`),
    `overcharged_breakpoint` (`>`), `undercharged_breakpoint` (`<`), `is_block` (`==`)

  Add helpers to `misc/genome.hpp`:
  ```cpp
  inline bool ir_leq(IR a, IR b) { return a.first <= b.first && a.second <= b.second; }
  inline IR   ir_min(IR a, IR b) { return IR(min(a.first, b.first), min(a.second, b.second)); }
  inline IR   ir_max(IR a, IR b) { return IR(max(a.first, b.first), max(a.second, b.second)); }
  ```
  - `transposition_move`: `min` → `ir_min`; `is_connected`: `ir_leq(target, nucleotides)`
  - `case_i`: `x = ir_min(ir_i, goal); y = ir_sub(goal, x);`
  - breakpoint predicates: breakpoint ⇔ `!ir_intersect(src, tgt)`;
    overcharged ⇔ `src.first > tgt.second`; undercharged ⇔ `src.second < tgt.first`
  - decide what a "mixed" case means for the algorithms (e.g. one component over, one
    under) and document it.

- [X] **`get_genomes` writes origin IRs negated** — `cycles.cpp:728`, `739`, `746` push the
  raw `vertices[i].weight`. Fix: push
  `add_sign_to_ir(vertices[i].weight, vertices[i].indel_update != -1)`.

- [X] **`find_2_cycles` greedy is wrong** — `cycles.cpp:~204-240` (TODO already there).
  Cycles that conflict with nothing never enter `graph`, so they are never selected; the
  greedy picks the *maximum*-degree node and does not remove its neighbours, so
  conflicting cycles reach `add_cycle`. Fix: insert every cycle as a node first
  (`graph[c.first];`), then repeatedly take the *minimum*-degree node and erase it **and
  its neighbours** (as `find_upto3_cycles` does).

- [X] **`add_cycle` silently accepts a contradicting cycle** — `cycles.cpp:994`. The check
  only runs when `vertices[v].fix_gray == NO_EDGE`; if `v` is already fixed to some
  `w != u`, the cycle is accepted. Add before that branch:
  ```cpp
  if (i % 2 == 1 && vertices[v].fix_gray != NO_EDGE && vertices[v].fix_gray != u &&
      u != vertices[v].red)
    throw runtime_error("cycle contradicts fixed gray edge");
  ```

- [X] **Wrong cap label for added empty chromosomes** — `cycles.cpp:~107`. The second cap
  uses `cap`; linear chromosomes end with `alt_cap`, so `indel_count[cap]`/`[alt_cap]`
  become unbalanced when the numbers of linear chromosomes differ. Use `alt_cap` (and
  `indel_count[alt_cap]`) for the second vertex.

- [X] **`mark_indel_edges` (random mode) can leave surplus genes unmarked** —
  `cycles.cpp:150-184`. With shuffled `idxs`, a run is only flushed at
  `last_origin_vertex`/`last_target_vertex`/a non-candidate; a trailing run is dropped.
  Flush once more after the `for` loop.

- [X] **No validation of interval IRs** — IR parsing (`genome.cpp:74-81`) accepts
  `lo > hi`; reversal/transposition splits (`genome.cpp:234-300`) check each component
  but not `x.first <= x.second`, so remainders like `[2,1]` can be produced. Throw on
  `r.first > r.second` after parsing; add the same check (as a real `if`, not an
  `assert`) after splits and before `transposition_move`.

- [ ] **External programs receive only the lower bound** — `external.cpp:15-21` sends
  `.first` of every IR. With interval targets this computes the distance to one specific
  target, not the flexible distance. At least throw when `lo != hi`, or pass `lo:hi` if
  the external program supports it.

## 3. Output

- [X] **`output(PermsIrs)` garbled**
  **Reopened 2026-10-14:** target loop (`io.cpp:~230`) still prints `p_ir_high[p_n - 2]` for every element; S/P loops also print a double space. — `misc/io.cpp:186-204`: source IRs print only
  `s_ir_low`; target loop prints `p_ir_high[p_n - 2]` for every element; line 202 prints
  an extra high value; line 203 prints the last pair again after the newline. Fix:
  ```cpp
  for (size_t i = 0; i + 1 < perms_irs.s_n; i++)
    os << perms_irs.s_ir_low[i] << ":" << perms_irs.s_ir_high[i]
       << (i + 2 < perms_irs.s_n ? " " : "\n");
  // same for p; delete line 203
  ```
- [X] **`output(Permutation)` prints only `.first` of source IRs** — `io.cpp:159`. Print
  `first << ":" << second` like the target line.
- [X] **Target-IR loops bounded by the source size**
  **Reopened 2026-10-19:** the *source* IR loop (`io.cpp:159`) is now bounded by `target_size()`; AddressSanitizer reports a heap-buffer-overflow in `Genome::get_ir` when the target is longer (`toperm -i <first 50 instances of in_full.txt>`). Bound it by `perm.size()` and end the line on `i == perm.size() - 1`.
  **Reopened 2026-10-14:** `io.cpp:~172` now loops to `target_size() - 1` but ends the line on `i == perm.size() - 1`; when sizes differ the newline is never written. — `Permutation::serialize` and
  `io.cpp:162`. With unbalanced input (`-d`) this reads past the target vector. Loop up
  to `intergenic_regions_target->size()`.
- [X] **`MultiGenome::serialize` writes to `cout`** — `genome.cpp:~533-538`. Use `os`.
- [X] **`operator<<(ostream&, const Run&)` writes `endl` to `cout`** — `cycles.cpp:902`. Use `os`.
- [X] **`print_queue` prints `q[pos]` instead of `q[i]`** — `cycles.cpp:492` (debug only).

## 4. Randomness, OpenMP, build

- [X] ⚠️ **RNG copied by value everywhere** — `Permutation` (`permutation.hpp:12-22`),
  `CycleGraph` (`cycles.hpp:154`, copy ctor), `DistAlg::estimate_distance`
  (`misc/dist.hpp:8-9`), `GA` (`heur/ga.hpp:39`). Every `-k` iteration, every
  `new CycleGraph(*cg)` in the `rand` heuristic, and the whole GA initial population start
  from the **same** generator state → identical "random" results; "best of k" is one
  sample. Pass `mt19937 &` everywhere; in the `CycleGraph` copy constructor reseed
  (`rand_gen(that.rand_gen())`) or pass the generator into `decompose`.

- [X] **GA data race** (fixed 2026-10-19: each `CycleGraph` owns its generator, seeds are drawn serially before every parallel loop via `draw_seeds`, and the `rand` reduction breaks ties by iteration; results no longer depend on the number of threads)
  **Reopened 2026-10-14:** `CycleGraph::rand_gen` is now an `mt19937 &`, so every copy shares `main`'s generator, and `dec`/`part` are now compiled with OpenMP. `decompose(true, ...)` → `shuffle(..., rand_gen)` runs concurrently in `ga.hpp:62`, `ga.hpp:88` and the `rand` heuristic (`heur/heur.cpp`); `select_parent`/`crossover`/`mutation` also call `rand_gen()` from threads. Give each thread its own generator (seeded once from the shared one under `omp critical`) and pass it into `decompose`. — `heur/ga.cpp:99-100` `shuffle(..., rand_gen)` on the shared
  member inside the `omp parallel for` (`ga.hpp:88`); `rand()` called from threads
  (`ga.cpp:14, 16, 22, 105, 133`). Give each thread its own
  `mt19937 local(seed + omp_get_thread_num())` and use
  `uniform_real_distribution<>(0, 1)(local)` instead of `rand()`.

- [X] **`dec` uses a fixed seed** — `main_dec.cpp:131` `mt19937 rand_gen(0);` (`rd` unused).
  Use `rd()`; add a `--seed` option if reproducibility is wanted.

- [X] **`dist` compiled without OpenMP** (fixed 2026-10-20: compiled with OpenMP; each instance gets its own generator from `draw_seeds` (`misc/seeds.hpp`), stdout output is buffered and printed in input order, and errors are caught inside the loop and reported after it)
  **Reopened 2026-10-19:** `CMakeLists.txt:30` still has `PRIVATE "${}"`. Before enabling it, note that every iteration shares `main`'s `rand_gen` through `Permutation`/`estimate_distance`, so it needs per-iteration seeds like the GA (`draw_seeds`). — `CMakeLists.txt:30` `PRIVATE "${}"`; linked with
  `-fopenmp` but not compiled with it, so `#pragma omp parallel for` in `main_dist.cpp`
  is ignored. Use `target_link_libraries(dist ... PRIVATE OpenMP::OpenMP_CXX)` (and the
  same for the other targets). Before enabling: buffer each instance's output and print
  it under `#pragma omp critical`, and catch exceptions inside the loop (see §1).
  `dec` is compiled without OpenMP on purpose (TODO at line 50) — add
  `-Wno-unknown-pragmas` to silence the `ga.hpp:62/88` warnings.

- [X] **`char op` holds `getopt_long`'s `int` return** — `main_dist.cpp:63`,
  `main_dec.cpp:82`, `main_part.cpp:82`, `main_toperm.cpp:35`, and
  `read_common_args(CommonArgs&, char op, ...)` (`io.hpp:41`, `io.cpp:8`). On platforms
  where `char` is unsigned (ARM), `-1` becomes 255 and the option loop never ends. Use
  `int op`.

## 5. Tests and minor CLI issues

- [X] **`dist_test` can never fail** — every `check<>` in `test.cpp:main()` (~560-606) is
  commented out. Uncomment and `return ok ? 0 : 1;` with `ok &= check<...>(...)`.
- [X] **Wrong algorithm in test** — `test.cpp:466`
  `CheckGraphDecCorrespondence1C2C3C` calls `decompose_with_bfs(true)`; should be
  `decompose_1c_2c_3c(true)`.
- [X] **Leaked buffer-less streams in GA tests** — `test.cpp:474, 487, 500, 512`
  `new ostream(0)`. Pass `nullptr` (GA already handles it).
- [X] **`toperm -o` writes nothing** — `main_toperm.cpp:85` opens `os`, then line 103 calls
  `os.open(...)` again on the open stream (sets failbit). Remove the second open (its
  `-all` name is copied from `dist`). The close at line 118 tests `input_file`; should
  test `output_folder`.
- [X] **`main_part` CLI** (fixed: now takes `HEUR DEC_ALGO`, empty `-perm` file removed) — help (line 39) advertises `HEUR DEC_ALGO` but `N_POS_ARGS` is 1
  and the algorithm is hard-coded to `"bfs"`; line ~270 creates an empty `-perm` file.

## 6. Latent (dead or unreachable code today)

- [X] `misc/reduction_rules.cpp` (now called in `main_dec`/`main_part`). Rechecked
  2026-10-20 with hand-made inputs:
  - [X] `pair<int,int>(ch, i + k)` → `(cg, i + k)`; `cg <= n_chrom()`; header signature;
    `i = idx_ar - 1` (all fixed in the current code).
  - [X] ⚠️ **every relabel uses the same label** (fixed: `MultiGenome::replace_label` updates `MultiGenome::op_max`; verified each new label occurs once per genome) (lines ~114, ~195). `new_gene` is
    `max(g.get_op_max(), h.get_op_max())`, but `Genome::replace_label` only updates the
    chromosome's `op_max`, not the `MultiGenome` one, so the value never changes. Input
    `1 2 2 3` / `1 2 2 3` turns both 2s into label 6 in both genomes: they stay
    replicas, so the rule does nothing useful and can merge unrelated matches. Fix: a
    counter `next++` per matched pair, and update `MultiGenome::op_max` too (e.g. a
    `MultiGenome::replace_label(c, i, label)`); otherwise the `CycleGraph` cap
    (`max(size + 2n, op_max) + 1`) can collide with a new label.
  - [X] **the left anchor's orientation is never checked** (fixed: `same_sign` check before the three cases) in the interval rule. `pos()`
    matches on `abs`, so `1 2 3 2` / `-1 2 3 2` is matched as a direct interval and
    `1 2 3 2` / `-3 -2 1 2` as a reverse one. Require `h[j] == g[i]` for the direct and
    middle-reversal cases and `h[j] == -g[i]` for the reverse case.
  - [X] **alpha (indel) genes can be anchors** (fixed: `check_alpha` positions are excluded from the singletons): `original_singletons_pos` only checks
    `occ == 1` in both genomes. With input zeros in both genomes, the two deleted genes
    get the same label (see the `0` parsing item), so they look like a shared
    singleton. Exclude `check_alpha` positions from the singletons.
  - [X] **`main_part.cpp:60-61` does not compile** (fixed: the rules are applied to `MultiGenome`s built before the `CycleGraph`; the calls are still commented out): it passes `Genome` (`GenomePair`)
    to functions taking `MultiGenome&`. Build the two `MultiGenome`s first, apply the
    rules to them, then construct the `CycleGraph` from them.
  - Not wrong, but missed matches: circular chromosomes do not wrap, and the
    middle-reversal case only looks to the right of `b^l` in the same orientation.
- [X] ⚠️ **input `0`s in both genomes become the same gene** (fixed 2026-10-20: `input`/`input_mg` label G's alphas above every label of both genomes and H's above G's; `MultiGenome` string constructors take `op_max_init`) — `Genome(string, ...)`
  (`genome.cpp:~43`) labels each `0` with `++op_max` of *its own* genome, so G and H
  zeros get equal labels, and `CycleGraph` adds gray edges between them (it compares
  `gene_val` only). `dec rand bfs -z -e` gives distance 0 for `1 0 2` / `1 0 2`
  (two different indel genes). Label H's zeros after G's (pass the other genome's
  `op_max` in `input`/`input_mg`) or skip alphas when adding gray edges.
- [X] `distance_algorithms/perm/misc/list.c:70` — `pop_back` sets `n->prev = NULL`;
  should be `n->next = NULL` (unused today).
- [X] `distance_algorithms/perm/perm_rearrange.c:18-23` — signed reversal returns 0
  without computing anything (unused today; throw instead).
- [X] `distance_algorithms/perm/HartmanShamir/hartman_shamir.c:243` — reads `es[idx+1]`
  past the end if the number of 2-cycles is odd → `idx + 1 < num_two_cicles`.

---

## Fixed (for reference)

- `genome.cpp:65` ternary precedence (infinite loop with `-z`)
- `MultiGenome` copy constructor `chromosomes[i + 1]`
- `record_positions` buffer size with `duplicate = true`
- `target_ir` negative index (now uses `abs`); `soft_breakpoint` treats alpha genes as breakpoints
- `make_circular` now calls `record_positions()`
- reversal assert bound for linear genomes
- `balanced()` uses `ir_intersect`
- `add_sign_to_ir` returns `(-hi, -lo)`
- `get_trivial_bal_cycle` balance test uses `ir_sum`
- `get_tuptok_cycle` `vertices[-1]` read
- `get_perms`: `p_n` set, 0-based indexing, buffers freed
- `mark_indel_edges`: double count per gene, loop condition
- `bfs`: empty-queue guard; per-path `indel_count` decremented on `q.back()`
- `add_cycle`/`rem_cycle`: pre-marked indels (`is_indel_fix`) no longer counted twice
- `dcj_ir_distance`: `merges` guard and `splits = (last_origin_vertex + 1) / 2 - ...`
- `r_or_rt_noir.cpp` use-after-move (`std::move(data.g)`)
- `external.cpp` `fscanf` check (`ret != 1`)
- `getopt_long` option tables terminated
- blank input lines (`io.cpp:83`)
- input line-count check `% div != 0`; input file open checked; `-k < 1` rejected
- `rand` heuristic OpenMP reduction (initializer, accumulate with `get_best_cg`, shared `get_best_cg` in `cycles.hpp`)
- GA elitism (`ga.cpp:45` uses `>`); GA `start` divide-by-zero
- `rt45` stage 2 `case_x == 4` now throws instead of looping; `reversal4` throws when no case found

## 7. Found in the 2026-10-14 review

- [X] ⚠️ **`reversal4` throws on every case-1 step** — `reversal4.cpp:~49`. `rev_applied`
  is never set in the `case_x == 1` branch, so any successful case-1 move ends with
  `"Error: no reversal applied."` (reproduced: `1 2 3 4` → `1 3 2 4` with `-u`). Set
  `rev_applied = true;` there.
- [X] **Indel counting in `follow_indel_edge` may be inverted** — `cycles.cpp:~486`. The
  local `indel_count` is updated only when `vertices[v].is_indel` is already true (an
  edge already counted by `mark_indel_edges`/`add_cycle`), while a newly used indel edge
  is never counted during the search. A path could then use two copies of the same extra
  gene as indels, and `add_cycle` would push `indel_count[g]` past zero. Needs checking
  against the intended semantics.
- [ ] **`read_lines` drops empty lines** — `io.cpp:~84`. An empty intergenic-region line
  is valid input (single-gene genome without `-e`; `Genome` fills zeros for `""`).
  Dropping it shifts every later instance or trips the "multiple of 4" check. Skip only
  comment lines, or skip empty lines only in `-z` mode.
- [X] ⚠️ **`1c2c3c` throws "cycle does not follows the fixed edge"** — `decompose_1c_2c_3c`
  adds the 2-/3-cycles from `find_2_cycles`/`find_upto3_cycles` without checking them.
  Fixing a gray edge v–u in `add_cycle` also fixes red(v)–red(u), and a later gray edge
  of the same (or a later) cycle can contradict that implied edge. `check_cycle` does not
  simulate the implied edges either. Verified fix (all 1C2C3C tests pass on 6 seeds):
  `check_cycle` keeps a local map of edges the cycle would fix, including
  red(v)–red(u), and rejects conflicts; call `if (check_cycle(cycle)) add_cycle(cycle);`.
- [X] **`find_2_cycles` greedy independent set is wrong** — `cycles.cpp:~230`. It keeps the
  max-degree vertex and never removes its neighbours (so overlapping 2-cycles are all
  returned), and 2-cycles with no conflicts are never put in `graph` (so they are
  dropped). Pick the min-degree vertex, remove it *and its neighbours*, and seed `graph`
  with every cycle.
- [ ] **`rt45` reaches case 4** — `reversal_transposition45.cpp:~315` throws "case iv
  should not be reached" on random inputs (test seed 1031129599). Verified fix:
  `dist += case_iv(pi, i, j); swap_stage = false;` (PBounds45RT passes on 3 seeds).

## 8. Found in the 2026-10-19 review (ASan/UBSan runs on the binaries)

- [X] ⚠️ **`dist reversal_noir` crashes** (fixed 2026-10-20: when the run goes right after the left cap, `cycle_run` returned the cap's artificial `gene_val`, which is not a gene of the `Genome`. `Run::after_cap` now flags this and `insertion_pos` in `r_or_rt_noir.cpp` inserts after the first gene; it throws instead of indexing an empty `pos()`) — `r_or_rt_noir.cpp:59`
  `data.g->pos(run.fst_gene)[0]` indexes an empty vector: `cycle_run` returns a
  `fst_gene` that is not in the genome where the run is inserted (null dereference,
  SEGV). Reproduce: `dist reversal_noir -i <first 200 lines of in_full.txt> -k 2`. Needs
  a look at how `cycle_run` picks `fst_gene` (vertex `v` vs `first_end`, and which genome
  it belongs to).
- [X] **`dist reversal_transposition_noir` fails `assert(lb < last_lb)`** (fixed 2026-10-20: the assert was too strict. With `div = 2`, ⌈x/2⌉ need not decrease when x decreases by 1; the asserts now check the numerator `lower_bound_gen(data, cg, 1)`) —
  `r_or_rt_noir.cpp:74`, same input as above. Either an insertion run does not always
  lower the bound (then the assert is wrong) or the same `cycle_run` problem.
- [X] **`cycle_fixes` accepted a red edge at a vertex with a different fixed gray edge** —
  regression from the 1C2C3C fix: the `u == red(v)` shortcut ran before the fixed-edge
  check. Fixed by moving it after (old `check_cycle` behaviour). Note `add_cycle` still
  skips its "does not follows the fixed edge" throw when `u == red(v)` (~line 1020).
- [X] **an alpha could get the label of a real gene** — `Genome(string, ...)` labeled
  alphas from the raw maximum but shifts the other genes by 1 afterwards, so without
  `-e` the first alpha got the label of the largest inner gene (`1 5 0 2 3` →
  `(0) (6) (6*) (3) (1)`). Fixed: `op_max++` before labeling the alphas.
- [X] **`dist`/`part`/`toperm` cannot read the thesis files** (fixed 2026-10-22: `input()` reads through `input_mg()` and takes the single linear chromosome; other inputs give a clear error) — they use `read_instance`
  (`Genome(string, ...)`), which does not accept the `L`/`C` chromosome markers or the
  `;` separator (`stoi` error); only `read_instance_mg` (`dec`, `dist --dec`) does.
  Workaround: `sed -E 's/^L //; s/ *; *$//'`. Fix: parse with `MultiGenome` and
  take the single chromosome (or throw if there are more).

## 9. Found on 2026-11-04 (while testing `dist -F`)

- [ ] ⚠️ **`dist reversal_noir` fails `assert(data.g->balanced(*data.h))`** —
  `r_or_rt_noir.cpp:96`, with fixed intergenic regions and indels (any thesis `F0` file
  with `I>0`, e.g. `dist reversal_noir -e -k 3 -i N200/SREV/N200_SREV100_I50_L10_F0.txt`
  in a Debug build; also with `-F` on flexible files). Not caused by `-F` (not
  investigated): after the loop of insertion runs the genomes still have different
  gene sets, probably because `cycle_run` (after `decompose_with_bfs`) does not
  produce a run for every missing gene; in Release GRIMM then receives unbalanced
  genomes without any error. Maybe related to the `fst_gene` item of section 8.

## 10. Found on 2026-11-16 (DCJ distances of the thesis)

- [X] ⚠️ **`dcj_indel_ir_distance` counts each non-trivial cycle only once** —
  `split_trivial_cycle` took `c = vertices[a].black` *after* the `dcj`, so `c` was in
  the same cycle as `a` and the other part of the split cycle was dropped from
  `cycles` (never split again). The `-dcj` outputs of `dec` were the number of
  non-trivial/unbalanced cycles (e.g. 29 instead of 105). Same pattern in the
  `else` branch of `dcj_indel_ir_distance` (pushed `cy.start_vtx`). Fixed by taking
  `c` before the `dcj` in both places. The `-dcj` files produced before the fix
  (e.g. `ANS_TESE_DEC_WITHOUTRULES/DCJ`) are wrong and must be recomputed.
- [X] ⚠️ **Signs lost in the splits of `dcj_indel_ir_distance`** — `dcj()` and
  `split_trivial_cycle` used `ir_abs` of the weights; when a split left a negative
  remainder, its sign was lost and the total weight of the cycle changed (extra
  unbalanced trivial cycles even with fixed intergenic regions). Fixed by using the
  weights with sign (`add_sign_to_ir`) (commit a0beb27).
- [X] ⚠️ **Flexible regions in `dcj_indel_ir_distance`** — each split took the lower
  bound of the interval of the target edge, consuming the flexibility and leaving the
  last trivial cycle unbalanced (the distance grew with the flexibility). Fixed by
  `fix_flexible_weights`, which fixes the target weights cycle by cycle before the
  splits, so every balanced cycle stays balanced (commit a0beb27). The choice is
  greedy (from the lower bounds): exact for balanced cycles, not necessarily the
  closest sum for unbalanced ones.
- [ ] **Dead test in `split_trivial_cycle`** — `weight_one` is overwritten with the
  weight of the second edge before `weight_one >= weight_two`, so the test is always
  true and the fallback `dcj(a, b, 0, 0, true)` of `dcj_indel_ir_distance` never runs.
  The result is correct (every split conserves the weight), but the code can be
  simplified.
- [X] **No DCJ distance in `dist`** — `run_dist_tese.sh` used `decomposition` for SDCJ,
  whose `estimate_distance` returns the gain, not a distance, and `-D` accepted only
  single-chromosome genomes. `dist dcj` (2026-11-21) computes `dcj_indel_ir_distance`
  on the decompositions of `-D` (multiple chromosomes allowed, via
  `CycleGraph::rebuild_cycles`) or on random mappings (`decompose(true, "map")`),
  and `run_dist_tese.sh` now uses it for SDCJ.

## 11. Found on 2026-11-23 (`ga 1c2c3c` on the thesis data set)

- [X] ⚠️ **`ga 1c2c3c` crashes or hangs** — e.g.
  `dec ga 1c2c3c -e -i N200/SDCJ/N200_SDCJ100_I50_L100_F0.txt --mut 20 -c 50`
  (segmentation fault; with ASan, a heap-buffer-overflow in `cycle_potation` called by
  `rem_cycle` in `GA::mutation`, and in other instances an endless walk in
  `cycle_potation`). The crossover and the mutation of `ga` decompose a graph that
  already has cycles, and `decompose_1c_2c_3c` starts with `mark_indel_edges`, which
  marked as indels genes that were in a cycle or had a fixed gray edge. The walks of
  those cycles then followed the red edge instead of the gray one and did not close.
  `bfs` and `map` do not call `mark_indel_edges`, and `rand` always starts from an
  empty graph, so only `ga 1c2c3c` was affected. Fixed by marking only the genes whose
  extremities are free (`gene_free`). With `rand -k 1` (deterministic), the outputs
  of `1c2c3c` and `bfs` did not change, with and without `-f`/`-f -p`.
- [X] ⚠️ **`1c2c3c` with `-f` hangs** — e.g. instance 5 of
  `N200/SDCJ/N200_SDCJ100_I50_L100_F0.txt` with the SOARComb `FIX` file (with any
  heuristic; `bfs` and `map` are not affected). With fixed adjacencies,
  `mark_indel_edges` marks an indel only if the whole block can be an indel
  (`indel_closure`). When the closure failed for every gene of the largest sequence
  of candidates, nothing changed and the same sequence was chosen forever. The
  closure depends only on the caps and on the fixed gray edges, which do not change
  while indels are marked, so these genes are now left out of the candidates
  (`blocked`). The outputs of the instances that finished before did not change.
