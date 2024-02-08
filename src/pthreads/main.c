#include <stdio.h>
#include <pthread.h>
#include <time.h>

#include "SparseMatrix.h"
#include "solve_funcs_pthread.h"

#define OUTPUT_FNAME "output_pth.txt"

// TODO: ADD ERROR HANDLING
int read_b(FILE *fin, double **b, int n) {
  double dbg = 0;
  int ret = 0;
  *b = (double *)malloc(sizeof(double) * n);
  for (int i = 0; i < n; i++) {
    // fscanf(*fin, "%f", (*b) + i);
    ret = fscanf(fin, "%lf", &dbg);
    if (!ret) {
      return ret;
    }
    (*b)[i] = dbg;
  }
  return ret;
}

int print_result(FILE *fout, double *res, int n) {
  for (int i = 0; i < n; i++) {
    fprintf(fout, "%lf\n", res[i]);
  }
  return 1;
}

int main(void) {
  struct SparseMatrix m;
  double *b;
  FILE *fin = fopen("input.txt", "r");
  time_t before, after;
  int msecs;

  if (!fin) {
    perror("ERROR OPENING FILE \"input.txt\"\n");
    return -1;
  }

  if (!spm_read(fin, &m)) {
    perror("ERROR READING MATRIX FROM FILE\n");
    fclose(fin);
    return -1;
  }

  // ERROR READING B. WRONG VALUES FROM FILE
  if (!read_b(fin, &b, m.n)) {
    perror("ERROR READING b FROM FILE\n");
    fclose(fin);
    spm_destroy(&m);
    return -1;
  }

  fclose(fin);

  if (!spm_lu_factorize(&m)) {
    printf("SLAU CANNOT BE SOLVED WITH THIS METHOD\n");
    spm_destroy(&m);
    free(b);
    return 0;
  }

  before = clock();

  find_y(&m, b);
  find_x(&m, b);

  after = clock();
  msecs = (after-before)*1000/CLOCKS_PER_SEC;

  FILE *fout = fopen(OUTPUT_FNAME, "w");
  if (!fout) {
    perror("ERROR OPENING FILE \"input.txt\"\n");
    spm_destroy(&m);
    free(b);
    return -1;
  }

  print_result(fout, b, m.n);
  printf("[PTH]: Solved in %d ms\n", msecs);

  spm_destroy(&m);

  free(b);

  return 0;

}