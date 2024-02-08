#include "SparseMatrix.h"
#include <mpi.h>

int find_y(struct SparseMatrix *m, double *b, int rank, int ranks, int root);
int find_x(struct SparseMatrix *m, double *y, int rank, int ranks, int root);