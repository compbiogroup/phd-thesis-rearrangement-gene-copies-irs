#include "main_common.hpp"
#include <experimental/filesystem>
#include <iostream>
#include <stdexcept>
namespace fs = experimental::filesystem;

void parse_args(CommonArgs &args, vector<string> &pos_args, int argc,
                char *argv[], size_t n_pos, void (*help)(char *),
                const string &extra_str, const vector<option> &extra_longopts,
                function<bool(int, char *)> read_extra) {
  vector<option> longopts;
  for (const option *o = common_longopts; o->name != NULL; ++o) {
    longopts.push_back(*o);
  }
  longopts.insert(longopts.end(), extra_longopts.begin(),
                  extra_longopts.end());
  longopts.push_back({NULL, 0, NULL, 0});

  int op;
  while ((op = getopt_long(argc, argv, (common_str + extra_str).c_str(),
                           longopts.data(), NULL)) != -1) {
    if (!read_common_args(args, op, optarg) &&
        !(read_extra && read_extra(op, optarg))) {
      if (op != 'h') {
        cout << "Unknown option: " << op << endl;
      }
      help(argv[0]);
    }
  }

  pos_args.assign(argv + optind, argv + argc);
  if (pos_args.size() != n_pos) {
    cout << "Wrong number of positional arguments: " << pos_args.size()
         << endl;
    help(argv[0]);
  }
}

int instance_lines(const CommonArgs &args) { return args.fill_zero ? 2 : 4; }

bool missing_instance(const CommonArgs &args, const vector<string> &lines,
                      size_t i) {
  for (int k = 0; k < instance_lines(args); ++k) {
    if (lines[i + k] != "NA") {
      return false;
    }
  }
  return true;
}

unique_ptr<vector<string>> read_input(const CommonArgs &args) {
  unique_ptr<vector<string>> lines;
  if (args.input_file != "") {
    ifstream is(args.input_file);
    if (!is.is_open()) {
      throw invalid_argument("Could not open input file: " + args.input_file);
    }
    lines.reset(read_lines(is));
  } else {
    lines.reset(read_lines(cin));
  }

  if (lines->size() % instance_lines(args) != 0) {
    throw invalid_argument(
        "Number of lines is not multiple of 4 (or 2 if fill_zero).");
  }
  if (args.iterations < 1) {
    throw invalid_argument("Number of iterations must be at least 1.");
  }
  return lines;
}

GenomePair read_instance(const CommonArgs &args, vector<string> &lines,
                         size_t i) {
  if (args.fill_zero) {
    return input(lines[i], lines[i + 1], args.extend);
  }
  return input(lines[i], lines[i + 1], lines[i + 2], lines[i + 3],
               args.extend);
}

MultiGenomePair read_instance_mg(const CommonArgs &args, vector<string> &lines,
                                 size_t i) {
  if (args.fill_zero) {
    return input_mg(lines[i], lines[i + 1], args.extend);
  }
  return input_mg(lines[i], lines[i + 1], lines[i + 2], lines[i + 3],
                  args.extend);
}

string output_path(const CommonArgs &args, const string &suffix, int idx) {
  string path = (args.output_folder / fs::path(args.input_file).filename())
                    .string();
  if (idx >= 0) {
    string num = to_string(idx);
    path += string(max(0, 5 - int(num.size())), '0') + num;
  }
  return path + suffix;
}

ostream &open_output(const CommonArgs &args, ofstream &os,
                     const string &suffix, int idx) {
  if (args.output_folder == "") {
    return cout;
  }
  os.close();
  os.open(output_path(args, suffix, idx));
  return os;
}

int run_main(const function<void()> &body) {
  try {
    body();
  } catch (const exception &e) {
    cerr << "Something went wrong!!!" << endl;
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
