#include "misc/genome.hpp"
#include "misc/io.hpp"
#include "misc/main_common.hpp"
#include "misc/permutation.hpp"
#include <iostream>
using namespace std;

void help(char *name) {
  cout << "usage: Read a pair of strings (maybe with intergenic regions)"
       << "and prints a single permutation correspondent to it"
       << "(uses a random map for replicas)." << endl
       << "\t" << name << " [OPTIONS]" << endl
       << endl;
  common_help();
  cout << endl;

  exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[]) {
  return run_main([&]() {
    CommonArgs args;
    vector<string> pos_args;
    parse_args(args, pos_args, argc, argv, 0, help);

    random_device rd;
    mt19937 rand_gen(rd());

    unique_ptr<vector<string>> input_lines = read_input(args);
    int div = instance_lines(args);

    ofstream os;
    ostream &out = open_output(args, os, "-perm");
    for (size_t i = 0; i < input_lines->size(); i += div) {
      GenomePair data = read_instance(args, *input_lines, i);
      Permutation pi(*data.g, *data.h, false, rand_gen);
      output(out, pi, true);
    }
  });
}
