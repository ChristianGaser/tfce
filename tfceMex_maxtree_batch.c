/*
 * Exact TFCE for many maps at once, one map per thread.
 *
 *   tfce = tfceMex_maxtree_batch(T, E, H, calc_neg, geom [, n_threads])
 *
 * T         - N x B matrix, one statistical map per column
 * E, H      - TFCE parameters for extent and height
 * calc_neg  - also calc neg. TFCE values
 * geom      - [nx ny nz] for volume data, or an F x 3 face list for surfaces
 * n_threads - number of worker threads (default: one per map, capped)
 *
 * tfce      - N x B matrix of TFCE maps, one per column
 *
 * Intended for permutation testing: parallelising across permutations rather
 * than within a single map keeps every thread on independent data, so there is
 * no locking on the hot path (workers only take a mutex to claim the next map
 * index). Each worker allocates its scratch arrays once and reuses them for
 * every map it processes.
 *
 * See tfce_maxtree.h for the algorithm.
 *
 * Christian Gaser
 */

#include "mex.h"
#include <pthread.h>
#include "tfce_maxtree.h"

#define MAX_THREADS 256

typedef struct {
  const double *T;
  double *out;
  int N, B;
  const Neigh *nb;
  double E, H;
  int calc_neg;
  int *next;                 /* next map to claim, shared */
  pthread_mutex_t *mtx;
  int ok;
} BatchArgs;

static void *batch_thread(void *pa)
{
  BatchArgs *A = (BatchArgs *)pa;
  Workspace ws;

  if (!ws_alloc(&ws, A->N)) { A->ok = 0; return NULL; }

  for (;;) {
    int b;

    pthread_mutex_lock(A->mtx);
    b = (*A->next)++;
    pthread_mutex_unlock(A->mtx);

    if (b >= A->B) break;

    maxtree_map(&ws, A->T + (size_t)b * A->N, A->out + (size_t)b * A->N,
                A->nb, A->E, A->H, A->calc_neg);
  }

  ws_free(&ws);
  return NULL;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double E, H;
  int calc_neg = 1, N, B, i, nthreads, next = 0;
  Neigh nb;
  BatchArgs args[MAX_THREADS];
  pthread_t th[MAX_THREADS];
  pthread_mutex_t mtx;
  int *aptr = NULL, *aidx = NULL;
  int started = 0;

  if (nrhs < 5)
    mexErrMsgTxt("Usage: tfce = tfceMex_maxtree_batch(T, E, H, calc_neg, geom [, n_threads])");
  if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
    mexErrMsgTxt("T must be real double.");
  if (mxGetNumberOfDimensions(prhs[0]) != 2)
    mexErrMsgTxt("T must be an N x B matrix, one map per column.");

  N = (int) mxGetM(prhs[0]);
  B = (int) mxGetN(prhs[0]);
  E = mxGetScalar(prhs[1]);
  H = mxGetScalar(prhs[2]);
  calc_neg = (int) mxGetScalar(prhs[3]);

  memset(&nb, 0, sizeof(nb));

  /* geom is either [nx ny nz] or an F x 3 face list */
  if (mxGetNumberOfElements(prhs[4]) == 3 && mxGetM(prhs[4]) == 1) {
    const double *g = mxGetPr(prhs[4]);
    nb.is_mesh = 0;
    nb.nx = (int)g[0]; nb.ny = (int)g[1]; nb.nz = (int)g[2];
    if (nb.nx * nb.ny * nb.nz != N)
      mexErrMsgTxt("prod(geom) does not match the number of rows of T.");
  } else if (mxGetN(prhs[4]) == 3) {
    build_adjacency(mxGetPr(prhs[4]), (int) mxGetM(prhs[4]), N, &aptr, &aidx);
    nb.is_mesh = 1; nb.ptr = aptr; nb.idx = aidx;
  } else {
    mexErrMsgTxt("geom must be [nx ny nz] or an F x 3 face list.");
  }

  nthreads = B;
  if (nrhs > 5 && !mxIsEmpty(prhs[5])) nthreads = (int) mxGetScalar(prhs[5]);
  if (nthreads < 1) nthreads = 1;
  if (nthreads > B) nthreads = B;
  if (nthreads > MAX_THREADS) nthreads = MAX_THREADS;

  plhs[0] = mxCreateDoubleMatrix((mwSize)N, (mwSize)B, mxREAL);

  if (pthread_mutex_init(&mtx, NULL) != 0)
    mexErrMsgTxt("Mutex init failed.");

  for (i = 0; i < nthreads; i++) {
    args[i].T        = mxGetPr(prhs[0]);
    args[i].out      = mxGetPr(plhs[0]);
    args[i].N        = N;
    args[i].B        = B;
    args[i].nb       = &nb;
    args[i].E        = E;
    args[i].H        = H;
    args[i].calc_neg = calc_neg;
    args[i].next     = &next;
    args[i].mtx      = &mtx;
    args[i].ok       = 1;
  }

  for (i = 0; i < nthreads; i++)
    if (pthread_create(&th[i], NULL, batch_thread, &args[i]) == 0) started++;
    else break;

  /* if threads could not be created, do the rest on this thread */
  if (started < nthreads) batch_thread(&args[0]);

  for (i = 0; i < started; i++) pthread_join(th[i], NULL);

  pthread_mutex_destroy(&mtx);
  if (aptr) { free(aptr); free(aidx); }

  for (i = 0; i < nthreads; i++)
    if (!args[i].ok) mexErrMsgTxt("Memory allocation error in worker thread.");
}
