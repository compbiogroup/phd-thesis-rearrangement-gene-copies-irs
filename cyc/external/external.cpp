#include "external.hpp"
#include <cstdio>
#include <sstream>

// For now it only sends the lower bound of the intergenic regions. Genes to
// be deleted are sent as 0 and the target may have a different size (indels).
int ExternalDistAlg::estimate_distance(Permutation &pi, mt19937 &rand_gen) {
  int dist;
  stringstream ss;
  ss << prog << " ";
  for (size_t i = 1; i <= pi.size(); i++) {
    ss << (pi.check_alpha(i) ? 0 : pi[i]) << ((i < pi.size()) ? "," : " ");
  }
  for (int i = 1; i < pi.target_size(); i++) {
    ss << pi.get_ir_target(i).first << ((i < pi.target_size() - 1) ? "," : " ");
  }
  for (size_t i = 1; i < pi.size() - 1; i++) {
    ss << pi.get_ir(i).first << ",";
  }
  ss << pi.get_ir(pi.size() - 1).first;
  FILE *prog_pipe;
  if ((prog_pipe = popen(ss.str().c_str(), "r")) == NULL) {
    throw runtime_error("Error calling child process.");
  }
  auto ret = fscanf(prog_pipe, "%d", &dist);
  if (ret != 1 || pclose(prog_pipe) != EXIT_SUCCESS) {
    throw runtime_error("Error on child process.");
  }
  return dist;
}
