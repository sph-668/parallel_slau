#include "solve_funcs_pthread.h"
// Ax = b
// LUx = b
// Ux = y; Ly = b;


int find_y(struct SparseMatrix *m, double *b) {
  // We are finding zi with bi, so previous b can be reused
  double *y = b;
  pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
  #ifdef DEBUG_PTH
    printf("[BOSS]: mutex: %p; m: %p\n",&mut, m);
  #endif

  // ai1*y1 + ai2*y2 + ...  + 1*yi = bi
  // yi = bi - (ai1*y1 + ai2*y2 + ...)
  for (int i = 0; i < m->n; i++) {
    pthread_t threads[NTHREADS];
    struct idxs_t targs[NTHREADS];
    int aih = m->ai[i+1]-m->ai[i];
    int threadh = aih / NTHREADS;
    for (int ti = 0; ti < NTHREADS; ti++) {
      int tadd = ti == NTHREADS-1?aih % NTHREADS:0;
      int startidx = m->ai[i]+threadh*ti;
      int endidx = m->ai[i]+threadh*(ti+1)+tadd;
      targs[ti].row = i;
      targs[ti].sidx = startidx;
      targs[ti].eidx =endidx;
      targs[ti].mutex = &mut;
      targs[ti].m=m;
      targs[ti].b=b;
      #ifdef DEBUG_PTH
        printf("[BOSS]: &targs[%d]: %p;\n",ti, &targs[ti]);
      #endif
      // TODO: Add error handling
      pthread_create(&threads[ti], NULL, pthread_f, &targs[ti]);
      #ifdef DEBUG_PTH
        printf("[BOSS]: threads[%d]: %lu;\n",ti, threads[ti]);
      #endif
    }
    for (int ti = 0; ti < NTHREADS; ti++) {
      pthread_join(threads[ti], NULL);
    }
    y[i] = (b[i]) / m->di[i];
  }

  return 1;
}

void *pthread_f(void* targs) {
  struct idxs_t* targ=(struct idxs_t*)targs;
  double sum = 0.;
  #ifdef DEBUG_PTH
    printf("[P %lu]: targs: %p; i %d; sidx %d; eidx %d; mutex: %p; m: %p;\n",pthread_self(), targ, targ->row, targ->sidx, targ->eidx,targ->mutex, targ->m);
  #endif
  // READONLY ACCESS?
  for (int idx = targ->sidx; idx < targ->eidx; idx++) {
      sum += targ->m->al[idx] * targ->b[targ->m->aj[idx]];
  }
  pthread_mutex_lock(targ->mutex);
  targ->b[targ->row] -= sum;
  pthread_mutex_unlock(targ->mutex);
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