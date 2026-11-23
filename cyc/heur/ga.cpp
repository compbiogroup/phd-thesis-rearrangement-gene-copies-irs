#include "ga.hpp"

#include <queue>
#include <set>

#include "../misc/io.hpp"

bool cmp_chr(const unique_ptr<Chromossome> &c1,
             const unique_ptr<Chromossome> &c2) {
  return c1->fitness() > c2->fitness();
}

int GA::select_parent(mt19937 &local) {
  int best_index =
      uniform_int_distribution<>(0, population->size() - 1)(local) %
      population->size();
  for (int i = 0; i < tournament_size; i++) {
    int index = uniform_int_distribution<>(0, population->size() - 1)(local);
    if ((*population)[index]->fitness() >
        (*population)[best_index]->fitness()) {
      best_index = index;
    } else if ((*population)[index]->fitness() ==
               (*population)[best_index]->fitness()) {
      if (uniform_int_distribution<>(0, 1)(local) == 0) {
        best_index = index;
      }
    }
  }
  return best_index;
}

void GA::select_population(unique_ptr<Population> &mutants) {
  // Find worse mutant.
  vector<int> worse_obj = (*mutants)[0]->fitness();
  int worse_idx = 0;
  for (size_t i = 1; i < mutants->size(); ++i) {
    if ((*mutants)[i]->fitness() < worse_obj) {
      worse_obj = (*mutants)[i]->fitness();
      worse_idx = i;
    }
  }

  // Find best chromosome from the original population.
  vector<int> best_obj = (*population)[0]->fitness();
  int best_idx = 0;
  for (size_t i = 1; i < population->size(); ++i) {
    if ((*population)[i]->fitness() > best_obj) {
      best_obj = (*population)[i]->fitness();
      best_idx = i;
    }
  }

  // Replace worse mutant with best chromosome from original population
  // if the value is better
  if ((*mutants)[worse_idx]->fitness() < (*population)[best_idx]->fitness()) {
    mutants->erase(mutants->begin() + worse_idx);
    mutants->push_back(std::move((*population)[best_idx]));
  }

  /* unique_ptr<Population> offsprings = nullptr; */
  /* offsprings.reset(new Population(mutants->size())); */

  /* sort(mutants->begin(), mutants->end(), cmp_chr); */
  /* for (size_t i = 0, j = 0, k = 0; k < mutants->size(); k++) { */
  /*     if ((*mutants)[i]->fitness() > (*population)[j]->fitness()) { */
  /*         (*offsprings)[k] = std::move((*mutants)[i]); */
  /*         i++; */
  /*     } else { */
  /*         (*offsprings)[k] = std::move((*population)[j]); */
  /*         j++; */
  /*     } */
  /* } */

  /* population = std::move(offsprings); */
  population = std::move(mutants);
}

bool GA::eval_population(unique_ptr<Population> &new_population, Timer timer) {
  bool improved = false;
  for (auto &chr : *new_population) {
    if (chr->fitness() > best_obj) {
      best_obj = chr->fitness();
      best_chr.reset(new Chromossome(*chr));
      improved = true;
    }
    if (os != nullptr) {
      output(*os, *chr, timer.since_last_mark());
    }
  }
  return improved;
}

/* We include cycles from the original decompositions ignoring conflicts.
 * Afterwards we use dec_algo to find new cycles. */
Chromossome *GA::crossover(const Chromossome &chr1, const Chromossome &chr2,
                           mt19937 &local) {
  auto c_list1 = vector<Cycle>(chr1.cycle_list());
  auto c_list2 = vector<Cycle>(chr2.cycle_list());
  Chromossome *chr = new Chromossome(*original);
  chr->reseed(local());

  shuffle(c_list1.begin(), c_list1.end(), local);
  shuffle(c_list2.begin(), c_list2.end(), local);
  /* sort(c_list1.begin(), c_list1.end()); */
  /* sort(c_list2.begin(), c_list2.end()); */

  while (c_list1.size() > 0 && c_list2.size() > 0) {
    double p = (double)uniform_real_distribution<>(0, 1)(local);
    if (p < crossover_rate) {
      chr->check_and_add_cycle(chr1.get_cycle(c_list1.back().start_vtx));
      c_list1.pop_back();
    } else {
      chr->check_and_add_cycle(chr2.get_cycle(c_list2.back().start_vtx));
      c_list2.pop_back();
    }
  }

  while (c_list1.size() > 0) {
    chr->check_and_add_cycle(chr1.get_cycle(c_list1.back().start_vtx));
    c_list1.pop_back();
  }

  while (c_list2.size() > 0) {
    chr->check_and_add_cycle(chr2.get_cycle(c_list2.back().start_vtx));
    c_list2.pop_back();
  }

  chr->decompose(true, dec_algo);
  return chr;
}

/* For each cycle we have a chance equals to mutation_rate to removed.
 * Afterwards we use dec_algo to find new cycles. */
void GA::mutation(unique_ptr<Chromossome> &chr, mt19937 &local) const {
  for (auto c : chr->cycle_list()) {
    double p = (double)uniform_real_distribution<>(0, 1)(local);
    if (p < mutation_rate) {
      chr->rem_cycle(c);
    }
  }
  chr->decompose(true, dec_algo);
}
