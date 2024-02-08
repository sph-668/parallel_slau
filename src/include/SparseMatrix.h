#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct SparseMatrix {
  int n;
  int *ai, *aj;
  double *di, *al, *au;
};

int spm_create(int n, unsigned int *ai, unsigned int *aj, unsigned int *al,
               unsigned int *au, struct SparseMatrix *m);

// TODO: ADD ERROR HANDLING AND RESOURCE FREEING
int spm_read(FILE *fin, struct SparseMatrix *m);

int spm_lu_factorize(struct SparseMatrix *m);

int spm_get(struct SparseMatrix *m, int i, int j, double *res);

int spm_set(struct SparseMatrix *m, int i, int j, double dbl);

// TODO: OPTIMIZE
int spm_print(struct SparseMatrix *m, FILE *fout);

void spm_destroy(struct SparseMatrix *m);

#endif  // ! SPARSE_MATRIX_H
