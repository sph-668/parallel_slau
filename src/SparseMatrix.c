#include "SparseMatrix.h"

int spm_create(int n, unsigned int *ai, unsigned int *aj, unsigned int *al,
               unsigned int *au, struct SparseMatrix *m) {
  return 0;
}

// TODO: ADD ERROR HANDLING
int spm_read(FILE *fin, struct SparseMatrix *m) {
  m->n = 0;
  m->ai = m->aj = NULL;
  m->al = m->au = m->di = NULL;

  if (!fin) return 0;

  if (!fscanf(fin, "%d", &m->n) || !m->n) return 0;

  m->ai = (int *)malloc(sizeof(int) * m->n);

  for (int i = 0; i <= m->n; i++) {  // i = [0; n]  because we need n+1 values
    if (!fscanf(fin, "%d", &m->ai[i])) return 0;
  }

  if (!m->ai[m->n])
  {
    m->di = (double *)malloc(sizeof(double) * m->n);
    for (int i = 0; i < m->n; i++) {
      if (!fscanf(fin, "%lf", &m->di[i])) return 0;
    }
    return 1;
  }

  m->aj = (int *)malloc(sizeof(int) * m->ai[m->n]);
  for (int i = 0; i < m->ai[m->n]; i++) {
    if (!fscanf(fin, "%d", &m->aj[i])) return 0;
  }

  m->di = (double *)malloc(sizeof(double) * m->n);
  m->al = (double *)malloc(sizeof(double) * m->ai[m->n]);
  m->au = (double *)malloc(sizeof(double) * m->ai[m->n]);

  for (int i = 0; i < m->n; i++) {
    if (!fscanf(fin, "%lf", &m->di[i])) return 0;
  }
  for (int i = 0; i < m->ai[m->n]; i++) {
    if (!fscanf(fin, "%lf", &m->al[i])) return 0;
  }
  for (int i = 0; i < m->ai[m->n]; i++) {
    if (!fscanf(fin, "%lf", &m->au[i])) return 0;
  }

  return 1;
}

int spm_lu_factorize(struct SparseMatrix *m) {
  double *L = m->al, *D = m->di, *U = m->au;
  for (int i = 0; i < m->n; i++) {
    double sd = 0;
    int i0 = m->ai[i], i1 = m->ai[i + 1];
    for (int j = i0; j < i1; j++) {
      double sl = 0, su = 0;
      int j0 = m->ai[m->aj[j]], j1 = m->ai[m->aj[j] + 1];
      int ki = i0, kj = j0;
      while (ki < j && kj < j1) {
        int jl = m->aj[ki];
        int ju = m->aj[kj];
        if (jl == ju) {
          sl += U[kj] * L[ki];
          su += L[kj] * U[ki];
          ki++;
          kj++;
        } else if (jl < ju)
          ki++;
        else
          kj++;
      }
      U[j] = (m->au[j] - su) / D[m->aj[j]];
      L[j] = (m->al[j] - sl);
      sd += U[j] * L[j];
    }
    D[i] = m->di[i] - sd;
    if (fabs(D[i]) <= 1e-15) return 0;
  }
  return 1;
}

int spm_get(struct SparseMatrix *m, int i, int j, double *res) {
  if (i == j) return m->di[i];
  double *v = i > j ? m->al : m->au;
  if (i < j) {
    int tmp = i;
    i = j;
    j = tmp;
  }
  size_t _i;
  for (_i = m->ai[i]; _i < m->ai[i + 1] && m->aj[_i] < j; _i++) {
  }
  if (_i < m->ai[i + 1] && j == m->aj[_i]) {
    *res = v[_i];
    return 1;
  }
  return 0;
}

int spm_set(struct SparseMatrix *m, int i, int j, double dbl) {
  if (i == j){
    m->di[i] = dbl;
    return 1;
  }
  double *v = i > j ? m->al : m->au;
  if (i < j) {
    int tmp = i;
    i = j;
    j = tmp;
  }
  size_t _i;
  for (_i = m->ai[i]; _i < m->ai[i + 1] && m->aj[_i] < j; _i++) {
  }
  if (_i < m->ai[i + 1] && j == m->aj[_i]) {
    v[_i] = dbl;
    return 1;
  }
  return 0;
}

// TODO: OPTIMIZE?
int spm_print(struct SparseMatrix *m, FILE *fout) {
  int res = 0;
  for (int i = 0; i < m->n; i++) {
    for (int j = 0; j < m->n; j++) {
      double val;
      if (!spm_get(m, i, j, &val)) val = 0;
      fprintf(fout, "%f ", val);
    }
    fprintf(fout, "\n");
  }
}

void spm_destroy(struct SparseMatrix *m) {
  if (!m->n) return;

  if (m->ai) {
    free(m->ai);
    m->ai = NULL;
  }
  if (m->di) {
    free(m->di);
    m->di = NULL;
  }
  // if (!m->ai[m->n]) return;
  if (m->aj) {
    free(m->aj);
    m->aj = NULL;
  }
  if (m->al) {
    free(m->al);
    m->al = NULL;
  }
  if (m->au) {
    free(m->au);
    m->au = NULL;
  }
}