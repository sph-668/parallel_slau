#include "SparseMatrix.h"
#include <pthread.h>

#define NTHREADS 4

struct idxs_t {
  int row;
  int sidx;
  int eidx;
  pthread_mutex_t *mutex; // shared
  struct SparseMatrix *m; // shared
  double *b; // shared
};

int find_y(struct SparseMatrix *m, double *b);
int find_x(struct SparseMatrix *m, double *y);
void *pthread_f(void* targs);