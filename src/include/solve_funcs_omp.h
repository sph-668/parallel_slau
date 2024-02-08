#include "SparseMatrix.h"
#include <omp.h>

int find_y(struct SparseMatrix *m, double *b);
int find_x(struct SparseMatrix *m, double *y);