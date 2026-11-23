#pragma once
#include "permutation.hpp"

/* Replace the flexible intergenic regions of the target of pi by fixed sizes
 * inside their intervals, adapting the function F of Brito et al. (2022),
 * "Genome Rearrangement Distance with a Flexible Intergenic Regions Aspect",
 * to genomes with distinct gene sets. A solution for the fixed instance is
 * also a solution for the flexible one. Requires extended genomes without
 * replicas. */
void fix_flexible_irs(Permutation &pi);
