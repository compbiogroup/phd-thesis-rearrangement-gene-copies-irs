#include "uniinvdist.h"
#include "mcstructs.h"
#include "mcrdist.h"
#include <stdio.h>

int grimm(int *g1, int *g2, int size) {
  int dist;
  distmem_t distmem_;

  struct genome_struct *g1_str = malloc(sizeof(struct genome_struct));
  struct genome_struct *g2_str = malloc(sizeof(struct genome_struct));
  g1_str->genes = g1;
  g2_str->genes = g2;
  mcdist_allocmem(size, 0, &distmem_);

  dist = invdist_noncircular(g1_str, g2_str, 0, size, &distmem_);

  mcdist_freemem(&distmem_);
  free(g1_str);
  free(g2_str);
  return dist;
}
