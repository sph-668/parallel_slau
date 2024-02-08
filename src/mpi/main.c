#include <stdio.h>
#include <mpi.h>
#include <time.h>

#include "SparseMatrix.h"
#include "solve_funcs_mpi.h"


#define OUTPUT_FNAME "output_mpi.txt"

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

int main(int argc, char** argv) {
  struct SparseMatrix m;
  m.n = 0;
  double *b;
  int process_Rank, size_Of_Cluster;
  int err;
  int root = 0;
  time_t before, after;
  int msecs;
  FILE *fin;

  if ((err = MPI_Init(&argc, &argv)) != MPI_SUCCESS)
  {
    #ifdef DEBUG_MPI
      printf("[P ...]: MPI INIT ERROR;\n");
    #endif
    goto m_exit;
  }

  MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);

  #ifdef DEBUG_MPI
    printf("[P %d/%d]: Starting;\n",process_Rank, size_Of_Cluster);
  #endif

  if (process_Rank == root)
  {
    fin = fopen("input.txt", "r");
    if (!fin) {
      perror("ERROR OPENING FILE \"input.txt\"\n");
      err = 1;
      goto m_mpi_finish;
    }

    if (!spm_read(fin, &m)) {
      perror("ERROR READING MATRIX FROM FILE\n");
      fclose(fin);
      err = 1;
      goto m_mpi_finish;
    }

    // ERROR READING B. WRONG VALUES FROM FILE
    if (!read_b(fin, &b, m.n)) {
      perror("ERROR READING b FROM FILE\n");
      fclose(fin);
      err = 1;
      goto m_free_m;
    }

    fclose(fin);

    if (!spm_lu_factorize(&m)) {
      printf("SLAU CANNOT BE SOLVED WITH THIS METHOD\n");
      goto m_free_b;
    }
    
    before = clock();
  }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&m.n, 1, MPI_INT, root, MPI_COMM_WORLD);
    #ifdef DEBUG_MPI
    printf("[P %d/%d]: m.n %d;\n",process_Rank, size_Of_Cluster, m.n);
    #endif
    if (m.n)
    {
      if (process_Rank != root)
      {
        m.ai = (int*)malloc(sizeof(int)*(m.n+1));
        m.di = (double*)malloc(sizeof(double)*m.n);
        b = (double*)malloc(sizeof(double)*m.n);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(m.ai, m.n+1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(m.di, m.n, MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(b, m.n, MPI_DOUBLE, root, MPI_COMM_WORLD);
      #ifdef DEBUG_MPI
      printf("[P %d/%d]: ai %p; di %p; b %p;\n",process_Rank, size_Of_Cluster, m.ai, m.di, b);
      #endif
      if (m.ai[m.n])
      {
        #ifdef DEBUG_MPI
        printf("[P %d/%d]: has m.ai[n];\n",process_Rank, size_Of_Cluster);
        #endif
        if (process_Rank != root)
        {
          m.aj = (int*)malloc(sizeof(int)*(m.ai[m.n]));
          m.al = (double*)malloc(sizeof(double)*(m.ai[m.n]));
          m.au = (double*)malloc(sizeof(double)*(m.ai[m.n]));
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(m.aj, m.ai[m.n], MPI_INT, root, MPI_COMM_WORLD);
        MPI_Bcast(m.al, m.ai[m.n], MPI_DOUBLE, root, MPI_COMM_WORLD);
        MPI_Bcast(m.au, m.ai[m.n], MPI_DOUBLE, root, MPI_COMM_WORLD);
        #ifdef DEBUG_MPI
        printf("[P %d/%d]: aj %p; al %p; au %p;\n",process_Rank, size_Of_Cluster, m.aj, m.al, m.au);
        #endif

        MPI_Barrier(MPI_COMM_WORLD);
        find_y(&m, b, process_Rank, size_Of_Cluster, root);
        find_x(&m, b, process_Rank, size_Of_Cluster, root);
      }
    }

  if (process_Rank == root)
  {
    after = clock();
    msecs = (after-before)*1000/CLOCKS_PER_SEC;

    FILE *fout = fopen(OUTPUT_FNAME, "w");
    if (!fout) {
      perror("ERROR OPENING FILE \"input.txt\"\n");
      err = 1;
      goto m_free_b;
    }
    print_result(fout, b, m.n);
    printf("[MPI]: Solved in %d ms\n", msecs);
  }
m_broadcast_exit:
  if (err && process_Rank == root)
  {
    MPI_Bcast(&m.n, 1, MPI_INT, root, MPI_COMM_WORLD);
  }
m_free_b:
    free(b);
m_free_m:
  #ifdef DEBUG_MPI
    printf("[P %d/%d]: Destroying m;\n",process_Rank, size_Of_Cluster);
  #endif
    spm_destroy(&m);
m_mpi_finish:
  #ifdef DEBUG_MPI
    printf("[P %d/%d]: Finishing;\n",process_Rank, size_Of_Cluster);
  #endif
  MPI_Finalize();
m_exit:
  return -err;

}