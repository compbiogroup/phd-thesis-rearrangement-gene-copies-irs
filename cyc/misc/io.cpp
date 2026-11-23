#include "io.hpp"
#include "genome.hpp"
#include "permutation.hpp"
#include "../cycle/cycles.hpp"
#include <fstream>
#include <iostream>
#include <sstream>

// TODO: Add single function to argument validation, also discard genomes of size 0 and 1
bool read_common_args(CommonArgs &args, int op, char *optarg) {
  switch (op) {
  case 'i':
    args.input_file = optarg;
    break;
  case 'o':
    args.output_folder = optarg;
    break;
  case 'k':
    args.iterations = atoi(optarg);
    break;
  case 'e':
    args.extend = true;
    break;
  case 'z':
    args.fill_zero = true;
    break;
  case 'u':
    args.unsign = true;
    break;
  case 'd':
    args.unbalanced = true;
    break;
  default:
    return false; // Unknown option
  }
  return true;
}

void common_help() {
  cout
      << "optional arguments:" << endl
      << "\t-h, --help              show this help message and exit" << endl
      << "\t-i, --input INPUT       "
      << "input file (if not provided stdin is "
         "used). Each 4 lines of the input file correspond to a instance, "
         "each line has a list of space separated values, and represent in "
         "order the origin strings, the origin intergenic regions lists, the "
         "target strings, and the target intergenic regions lists. Mutliple "
         "chromossomes must be separeted by the charactere ; both in the "
         "strings "
         "and intergenic regions lists. An linear chromossome may be proceeded "
         "by "
         "the charactere L and an circular chromossome must be proceeded by an "
         "charactere C (if there is no letter it is assumed linear). The"
         "intergenic region between the first and last gene of a "
         "circular chromossome must be the last on the list. If the -z"
         "option is used, each 2 lines correspond to a instance. There is two"
         "options to represent genes to be inserted or deleted, either with the"
         "number zero or with a number that appear in only one of the genomes."
      << endl
      << "\t-o, --output OUTPUT     output folder (if not provided stdout is "
         "used)"
      << endl
      << "\t-k, --iterations ITER   number of iterations (default 100)" << endl
      << endl
      << "\t-e, --extend            whether to extend (cap) the genomes before "
         "apply the algorithm. Should always be used with multiple "
         "chromossomes and is ignore in circular chromossomes. In "
         "unichromossomal linear genomes the genomes may "
         "already be extended."
      << endl
      << "\t-z, --zero              fill intergenic regions with zero instead "
         "of reading them."
      << endl
      << "\t-u, --unsign          use unsigned genomes (default is signed)."
      << endl
      << "\t-d, --unbalanced      allows the genomes to have distinct sets of genes."
      << endl;
}

vector<string> *read_lines(istream &is) {
  auto input_lines = new vector<string>();

  for (string line; getline(is, line);) {
    if (line.empty()) {
      continue;
    }
    if (line.at(0) != '#') {
      input_lines->push_back(line);
    }
  }

  return input_lines;
}

vector<string> read_decompositions(const string &file) {
  ifstream is(file);
  if (!is.is_open()) {
    throw invalid_argument("Could not open decomposition file: " + file);
  }
  stringstream ss;
  ss << is.rdbuf();
  string text = ss.str();

  const string mark = "Dec: ";
  vector<string> decs;
  size_t pos = text.find(mark);
  while (pos != string::npos) {
    size_t next = text.find(mark, pos + mark.size());
    decs.push_back(text.substr(pos + mark.size(), next - pos - mark.size()));
    pos = next;
  }
  return decs;
}

static vector<int> read_ints(const string &str) {
  vector<int> vals;
  stringstream ss(str);
  for (string tok; ss >> tok;) {
    if (tok != "-") {
      vals.push_back(stoi(tok));
    }
  }
  return vals;
}

vector<PartitionFix> read_partition_fixes(const string &file) {
  ifstream is(file);
  if (!is.is_open()) {
    throw invalid_argument("Could not open FIX file: " + file);
  }
  vector<string> lines;
  for (string line; getline(is, line);) {
    lines.push_back(line);
  }
  while (!lines.empty() && lines.back().empty()) {
    lines.pop_back();
  }
  if (lines.size() % 3 != 0) {
    throw invalid_argument("Number of lines of the FIX file " + file +
                           " is not multiple of 3.");
  }
  vector<PartitionFix> fixes;
  for (size_t i = 0; i < lines.size(); i += 3) {
    PartitionFix f;
    if (lines[i] == "NA" && lines[i + 1] == "NA" && lines[i + 2] == "NA") {
      f.missing = true;
    } else {
      for (int side = 0; side < 2; ++side) {
        vector<vector<int>> &breaks = (side == 0) ? f.breaks_g : f.breaks_h;
        stringstream ss(lines[i + side]);
        for (string chrom; getline(ss, chrom, ';');) {
          breaks.push_back(read_ints(chrom));
        }
      }
      f.corr = read_ints(lines[i + 2]);
    }
    fixes.push_back(f);
  }
  return fixes;
}

/* Genomes with a single linear chromosome, read as multichromosomal genomes so
 * the same format is accepted (L before the chromosome, ; after it). */
static GenomePair single_chromosome(MultiGenomePair mg) {
  if (mg.g->n_chrom() != 1 || mg.h->n_chrom() != 1) {
    throw invalid_argument("Expected genomes with a single chromosome.");
  }
  if ((*mg.g)[1].is_circular() || (*mg.h)[1].is_circular()) {
    throw invalid_argument("Expected linear genomes.");
  }
  GenomePair data;
  data.g = make_unique<Genome>((*mg.g)[1]);
  data.h = make_unique<Genome>((*mg.h)[1]);
  return data;
}

GenomePair input(string &line1, string &line2, bool extend) {
  return single_chromosome(input_mg(line1, line2, extend));
}

GenomePair input(string &line1, string &line2, string &line3, string &line4,
                 bool extend) {
  return single_chromosome(input_mg(line1, line2, line3, line4, extend));
}

/* Alphas (zeros in the input) get new labels, which must differ from every
 * label of both genomes. So the genomes are read once to find a bound for their
 * labels, the alphas of G are labeled above it and those of H above G's. */
MultiGenomePair input_mg(string &line1, string &line2, bool extend) {
  string str;
  MultiGenomePair data;

  Gene bound = max(MultiGenome(line1, extend).get_op_max(),
                   MultiGenome(line2, extend).get_op_max());
  data.g = make_unique<MultiGenome>(line1, extend, bound);
  data.h = make_unique<MultiGenome>(line2, extend, data.g->get_op_max());
  mark_alphas(data);

  return data;
}

MultiGenomePair input_mg(string &line1, string &line2, string &line3, string &line4, bool extend) {
  string str;
  MultiGenomePair data;

  Gene bound = max(MultiGenome(line1, line2, extend).get_op_max(),
                   MultiGenome(line3, line4, extend).get_op_max());
  data.g = make_unique<MultiGenome>(line1, line2, extend, bound);
  data.h = make_unique<MultiGenome>(line3, line4, extend, data.g->get_op_max());
  mark_alphas(data);

  return data;
}

void output(ostream &os, int dist, double time) {
  os << "Dist: " << dist;
  os.precision(5);
  os << fixed;
  os << ", Wall Time: " << time << "s" << endl;
}

void output(ostream &os, string &dist_etc, double time) {
  os << dist_etc << endl;
}

void output(ostream &os, const Permutation &perm, bool indel_zero) {
  for (size_t i = 1; i <= perm.size(); i++) {
    if (indel_zero && perm.check_alpha(i)) {
      os << "0 ";
    } else {
      os << perm[i] << " ";
    }
  }
  os << endl;
  for (size_t i = 1; i <= perm.size() - 1; i++) {
    if (perm.get_ir(i).first == perm.get_ir(i).second) {
      os << perm.get_ir(i).first;
    } else {
      os << perm.get_ir(i).first << ":" << perm.get_ir(i).second;
    }
    if (i == perm.size() - 1) {
      os << endl;
    } else {
      os << " ";
    }
  }
  for (size_t i = 1; i <= perm.target_size() - 1; i++) {
    os << perm.get_ir_target(i).first << ":" << perm.get_ir_target(i).second;
    if (i == perm.target_size() - 1) {
      os << endl;
    } else {
      os << " ";
    }
  }
}

void output(ostream &os, const CycleGraph &cyc_dec, double time) {
  os << "Dec: " << cyc_dec;
  os.precision(5);
  os << fixed;
  os << ", Wall Time: " << time << "s" << endl;
}

void output(ostream &os, GenomePair &data) {
    for (int i = 1; i <= int(data.g->size()); i++) {
            os << (*data.g)[i];
            if (i == int(data.g->size())) {
                os << endl;
            } else {
                os << " ";
            }
    }
    for (int i = 1; i <= int(data.h->size()); i++) {
            os << (*data.h)[i] << " ";
            if (i == int(data.h->size())) {
                os << endl;
            } else {
                os << " ";
            }
    }
}

void output(ostream &os, PermsIrs perms_irs) {
    for (size_t i = 0; i < perms_irs.s_n; i++) {
        os << perms_irs.s[i];
        if (i == perms_irs.s_n - 1) {
          os << endl;
        } else {
          os << " ";
        }
    }
    for (size_t i = 0; i < perms_irs.s_n - 1; i++) {
        if (perms_irs.s_ir_low[i] == perms_irs.s_ir_high[i]) {
            os << perms_irs.s_ir_low[i];
        } else {
            os << perms_irs.s_ir_low[i] << ":" << perms_irs.s_ir_high[i];
        }
        if (i == perms_irs.s_n - 2) {
          os << endl;
        } else {
          os << " ";
        }
    }
    for (size_t i = 0; i < perms_irs.p_n; i++) {
        os << perms_irs.p[i];
        if (i == perms_irs.p_n - 1) {
          os << endl;
        } else {
          os << " ";
        }
    }
    for (size_t i = 0; i < perms_irs.p_n - 1; i++) {
        os << perms_irs.p_ir_low[i] << ":" << perms_irs.p_ir_high[i];
        if (i == perms_irs.p_n - 2) {
          os << endl;
        } else {
          os << " ";
        }
    }
}

void output(ostream &os, CommonPartition cp) {
    for (auto b : cp.breaks_s) {
        os << b << " ";
    }
    os << endl;
    for (auto b : cp.breaks_p) {
        os << b << " ";
    }
    os << endl;
}
