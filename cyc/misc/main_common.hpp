#pragma once
#include "genome.hpp"
#include "io.hpp"
#include <fstream>
#include <functional>
#include <getopt.h>
#include <memory>
#include <string>
#include <vector>
using namespace std;

/* Code shared by all main programs. */

/* Parse the command line. Besides the common options, accepts the options in
 * extra_str/extra_longopts, which are handled by read_extra (it must return
 * false for unknown options). The positional arguments are stored in pos_args
 * and there must be exactly n_pos of them. Calls help(argv[0]) on -h or on
 * error (help must exit). */
void parse_args(CommonArgs &args, vector<string> &pos_args, int argc,
                char *argv[], size_t n_pos, void (*help)(char *),
                const string &extra_str = "",
                const vector<option> &extra_longopts = {},
                function<bool(int, char *)> read_extra = nullptr);

/* Number of input lines that form a single instance. */
int instance_lines(const CommonArgs &args);

/* Read and validate the input lines (from the input file or stdin). */
unique_ptr<vector<string>> read_input(const CommonArgs &args);

/* Whether the instance that starts at line i is missing (all its lines are
 * "NA", e.g. an instance of run_part_tese.sh without a partition). The mains
 * write "NA" as its result, so the results keep the order of the input. */
bool missing_instance(const CommonArgs &args, const vector<string> &lines,
                      size_t i);

/* Read the instance that starts at line i. */
GenomePair read_instance(const CommonArgs &args, vector<string> &lines,
                         size_t i);
MultiGenomePair read_instance_mg(const CommonArgs &args, vector<string> &lines,
                                 size_t i);

/* Path of an output file for the instance with index idx (or for the whole
 * input if idx < 0). */
string output_path(const CommonArgs &args, const string &suffix,
                   int idx = -1);

/* Reopen os on the output file (if an output folder was given) and return the
 * stream to write to (os or cout). */
ostream &open_output(const CommonArgs &args, ofstream &os,
                     const string &suffix, int idx = -1);

/* Run body reporting any error. Returns the program exit status. */
int run_main(const function<void()> &body);
