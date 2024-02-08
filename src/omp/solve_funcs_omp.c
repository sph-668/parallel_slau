#include "solve_funcs_omp.h"
// Ax = b
// LUx = b
// Ux = y; Ly = b;

int find_y(struct SparseMatrix *m, double *b) {
  // We are finding zi with bi, so previous b can be reused
  double *y = b;
  int nthreads, tid;

  // ai1*y1 + ai2*y2 + ...  + 1*yi = bi
  // yi = bi - (ai1*y1 + ai2*y2 + ...)
  for (int i = 0; i < m->n; i++) {
    double sum = 0.;
    int aih = m->ai[i+1]-m->ai[i];
    #pragma omp parallel private(nthreads, tid)
    {
      tid = omp_get_thread_num();
      nthreads = omp_get_num_threads();

      int threadh = aih / nthreads;
      int tadd = tid == nthreads-1?aih % nthreads:0;
      int startidx = m->ai[i]+threadh*tid;
      int endidx = m->ai[i]+threadh*(tid+1)+tadd;
      for (int idx = startidx; idx < endidx; idx++) {
        sum += m->al[idx] * y[m->aj[idx]];
      }
    }
    y[i] = (b[i] - sum) / m->di[i];
  }

  return 1;
}

int find_x(struct SparseMatrix *m, double *y) {
  // don't need previous y, so reusing it
  double *x = y;

  // 1*xi + aii+1*xi+1 + aii+2*xi+2 + ...  + ain*xn = yi
  // xi = yi - (aii+1*xi+1 + xii+2*xi+2 + ...)
  for (int i = m->n - 1; i >= 0; i--) {
    double sum = 0.;
    x[i] = y[i];
    for (int idx = m->ai[i]; idx < m->ai[i + 1]; idx++) {
      y[m->aj[idx]] -= m->au[idx] * x[i];
    }
  }
  return 1;
}