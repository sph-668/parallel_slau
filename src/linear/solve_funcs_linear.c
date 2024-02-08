#include "solve_funcs_linear.h"
// Ax = b
// LUx = b
// Ux = y; Ly = b;

int find_y(struct SparseMatrix *m, double *b) {
  // We are finding zi with bi, so previous b can be reused
  double *y = b;

  // ai1*y1 + ai2*y2 + ...  + 1*yi = bi
  // yi = bi - (ai1*y1 + ai2*y2 + ...)
  for (int i = 0; i < m->n; i++) {
    double sum = 0.;
    for (int idx = m->ai[i]; idx < m->ai[i + 1]; idx++) {
      sum += m->al[idx] * y[m->aj[idx]];
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