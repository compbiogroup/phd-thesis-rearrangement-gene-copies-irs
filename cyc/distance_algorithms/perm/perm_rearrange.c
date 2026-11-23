#include <stdlib.h>
#include <stdio.h>
#include "perm_rearrange.h"
#include "misc/util.h"
#include "misc/perm.h"
#include "misc/list.h"
#include "HartmanShamir/hartman_shamir.h"
#include "KaceciogluSankoff/kacecioglu_sankoff.h"
#include "WalterBP/walter_bp.h"

typedef enum Alg {GreedyBP=1,Cicle} Alg;

int dist(int *g1, int *g2, int size, int with_sign, int mod_, int alg_) {
  Alg alg = alg_;
  Model mod = mod_;
  PermType type = (with_sign) ? PSign : PUnsign;
  int dist = -1;
  perm *pi = (type == PSign && mod == Rev)
    ? create_perm(size, type, mod)
    : build_and_rename_perm(g1,g2,size,type, mod);
   
  if(is_signed(pi) && mod == Rev) {
    error("Signed permutations with reversals are not supported.");
  } else {
    switch(alg) {
    case GreedyBP:
      if(mod == Rev) {
        dist = kacecioglu_sankoff(pi);
      } else {
        dist = walter_bp(pi);
      }
      break;
    case Cicle:
      if(mod == Trans) {
        dist = hartman_shamir(pi);
      } else {
        error("Algorithm not implemented for that case.");
      }
      break;
    default:
        error("Invalid Permutation Algorithm.");
    }
  }

  print_ops(pi);
  clear_perm(pi);
  return(dist);
}

// int main() {
  // int g1[5] = {5,4,3,2,1};
  // int g2[5] = {1,2,3,4,5};
//
  // printf("%d\n",dist(g1, g2, 5, 0, 2));
//
  // return 0;
// }
