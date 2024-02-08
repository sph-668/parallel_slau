#include "solve_funcs_mpi.h"
#include <mpi.h>
// Ax = b
// LUx = b
// Ux = y; Ly = b;
int find_y(struct SparseMatrix *m, double *b, int rank, int ranks, int root) {
  // We are finding yi with bi, so previous b can be reused
  double *y = b;


  // ai1*y1 + ai2*y2 + ...  + 1*yi = bi
  // yi = bi - (ai1*y1 + ai2*y2 + ...)
  for (int i = 0; i < m->n; i++) {
    double sum = 0.;
    int aih = m->ai[i+1]-m->ai[i];
    int rankh = aih / ranks;
    int radd = rank == ranks-1?aih % ranks:0;
    int startidx = m->ai[i]+rankh*rank;
    int endidx = m->ai[i]+rankh*(rank+1)+radd;
    #ifdef DEBUG_MPI
    printf("[P %d/%d]: i %d; sidx %d; eidx %d;\n",rank, ranks, i, startidx, endidx);
    #endif
    for (int idx = startidx; idx < endidx; idx++) {
      sum += m->al[idx] * y[m->aj[idx]];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == root)
    {
      b[i] -= sum;
      for (int pi = 0; pi < ranks; pi++)
      {
        if (pi != root)
        {
          MPI_Recv(&sum, 1, MPI_DOUBLE, pi, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          #ifdef DEBUG_MPI
          printf("[P %d/%d]: sum %d %lf;\n",rank, ranks, pi, sum);
          #endif
          b[i] -= sum;
        }
      }
      y[i] = b[i] / m->di[i];
    } else {
      MPI_Send(&sum, 1, MPI_DOUBLE, root, 1, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&y[i], 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
    #ifdef DEBUG_MPI
    printf("[P %d/%d]: y[%d] %lf;\n",rank, ranks, i, y[i]);
    #endif
    MPI_Barrier(MPI_COMM_WORLD);
  }

  return 1;
}

int find_x(struct SparseMatrix *m, double *y, int rank, int ranks, int root) {
  // don't need previous y, so reusing it
  double *x = y;

  // 1*xi + aii+1*xi+1 + aii+2*xi+2 + ...  + ain*xn = yi
  // xi = yi - (aii+1*xi+1 + xii+2*xi+2 + ...)
  // MPI Here will be ineffective so don't touch it ubless you root
  if (rank == root)
  {
    for (int i = m->n - 1; i >= 0; i--) {
      double sum = 0.;
      x[i] = y[i];
      for (int idx = m->ai[i]; idx < m->ai[i + 1]; idx++) {
        y[m->aj[idx]] -= m->au[idx] * x[i];
      }
    }
  }
  return 1;
}