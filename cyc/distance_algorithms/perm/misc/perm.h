#ifndef H_PERM
#define H_PERM

#include "../misc/util.h"

typedef struct perm perm;
typedef enum PermType {PSign, PUnsign} PermType;
typedef enum Model {Rev, Trans, TransRev} Model;

perm *build_perm(int *vet, int n, PermType type, Model mod);
perm *build_and_rename_perm(int *vet1, int *vet2, int n, PermType type, Model mod);
perm *create_perm(int size, PermType type, Model mod);
void print_perm(perm *pi);
void clear_perm(perm *pi_free);

int length_perm(perm *pi);
Model get_model(perm *pi);
int get_perm(perm *pi, int i);
int get_perm_idx(perm *pi, int val);
void set_perm(perm *pi_mut, int i, int val);

Bool is_sorted(perm *pi);
Bool is_signed(perm *pi);
Bool singleton_strip(perm *pi, int i);

Bool break_point(perm *pi, int x, int y);
int end_strip(perm *pi, int i);

void rev(perm *pi_mut, int i, int j);
void trans(perm *pi_mut, int i, int j, int k);
void print_ops(perm *pi);

#endif
