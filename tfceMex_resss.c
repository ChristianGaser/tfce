/*
 * Fused residual sum of squares for the GLM.
 *
 *   ResSS = tfceMex_resss(Y, Beta, X [, n_threads])
 *
 * Y     - n_vox x n_data data matrix (single)
 * Beta  - n_vox x p parameter estimates (single), i.e. Y*pinv(W*X)'
 * X     - n_data x p design matrix (double)
 *
 * ResSS - n_vox x 1 residual sum of squares (double)
 *
 * Computes, per voxel v,
 *
 *     ResSS(v) = sum_j ( Y(v,j) - sum_k Beta(v,k)*X(j,k) )^2
 *
 * which is exactly what
 *
 *     res0 = Beta*single(X'); res0 = Y - res0; res0 = res0.^2; sum(res0,2)
 *
 * computes, but without materialising the n_vox x n_data residual matrix. The
 * original makes four trips through that array in DRAM; here the residual never
 * leaves the registers, so the only traffic is one read of Y and one of Beta.
 * The reduction is accumulated in double, so this is at least as accurate as
 * the expression it replaces (it does not use the y'y - b'X'X b identity, which
 * loses precision by cancellation when the design explains most of the variance).
 *
 * Voxels are blocked and distributed over threads.
 *
 * Christian Gaser
 */

#include "mex.h"
#include "math.h"
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#define MAX_THREADS 256
#define VBLOCK 4096          /* voxels per cache-resident block */

typedef struct {
  const float  *Y;           /* n_vox x n_data, column major */
  const float  *B;           /* n_vox x p,      column major */
  const double *X;           /* n_data x p,     column major */
  double *ResSS;
  int nv, nd, p;
  int *next;                 /* next block to claim */
  pthread_mutex_t *mtx;
  int nblocks;
} Args;

static void *worker(void *pa)
{
  Args *a = (Args *)pa;
  int nv = a->nv, nd = a->nd, p = a->p;
  double pred[VBLOCK];
  double ss[VBLOCK];

  for (;;) {
    int blk, v0, v1, n, i, j, k;

    pthread_mutex_lock(a->mtx);
    blk = (*a->next)++;
    pthread_mutex_unlock(a->mtx);
    if (blk >= a->nblocks) break;

    v0 = blk * VBLOCK;
    v1 = v0 + VBLOCK; if (v1 > nv) v1 = nv;
    n  = v1 - v0;

    for (i = 0; i < n; i++) ss[i] = 0.0;

    /* j outer, voxels inner: both Y and Beta are read down their columns */
    for (j = 0; j < nd; j++) {
      const float *Yj = a->Y + (size_t)j * nv + v0;

      for (i = 0; i < n; i++) pred[i] = 0.0;

      for (k = 0; k < p; k++) {
        double xjk = a->X[j + (size_t)k * nd];
        const float *Bk = a->B + (size_t)k * nv + v0;
        if (xjk == 0.0) continue;
        for (i = 0; i < n; i++) pred[i] += (double)Bk[i] * xjk;
      }

      for (i = 0; i < n; i++) {
        double d = (double)Yj[i] - pred[i];
        ss[i] += d * d;
      }
    }

    for (i = 0; i < n; i++) a->ResSS[v0 + i] = ss[i];
  }

  return NULL;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int nv, nd, p, i, nthreads, next = 0, nblocks, started = 0;
  Args args[MAX_THREADS];
  pthread_t th[MAX_THREADS];
  pthread_mutex_t mtx;

  if (nrhs < 3)
    mexErrMsgTxt("Usage: ResSS = tfceMex_resss(Y, Beta, X [, n_threads])");
  if (!mxIsSingle(prhs[0]) || !mxIsSingle(prhs[1]))
    mexErrMsgTxt("Y and Beta must be single.");
  if (!mxIsDouble(prhs[2]))
    mexErrMsgTxt("X must be double.");

  nv = (int) mxGetM(prhs[0]);
  nd = (int) mxGetN(prhs[0]);
  p  = (int) mxGetN(prhs[1]);

  if ((int) mxGetM(prhs[1]) != nv) mexErrMsgTxt("Beta must have as many rows as Y.");
  if ((int) mxGetM(prhs[2]) != nd) mexErrMsgTxt("X must have as many rows as Y has columns.");
  if ((int) mxGetN(prhs[2]) != p)  mexErrMsgTxt("X must have as many columns as Beta.");

  plhs[0] = mxCreateDoubleMatrix((mwSize)nv, 1, mxREAL);

  nblocks = (nv + VBLOCK - 1) / VBLOCK;

  nthreads = nblocks;
  if (nrhs > 3 && !mxIsEmpty(prhs[3])) nthreads = (int) mxGetScalar(prhs[3]);
  if (nthreads < 1) nthreads = 1;
  if (nthreads > nblocks) nthreads = nblocks;
  if (nthreads > MAX_THREADS) nthreads = MAX_THREADS;

  if (pthread_mutex_init(&mtx, NULL) != 0) mexErrMsgTxt("Mutex init failed.");

  for (i = 0; i < nthreads; i++) {
    args[i].Y       = (const float *)  mxGetData(prhs[0]);
    args[i].B       = (const float *)  mxGetData(prhs[1]);
    args[i].X       = (const double *) mxGetPr(prhs[2]);
    args[i].ResSS   = mxGetPr(plhs[0]);
    args[i].nv      = nv;
    args[i].nd      = nd;
    args[i].p       = p;
    args[i].next    = &next;
    args[i].mtx     = &mtx;
    args[i].nblocks = nblocks;
  }

  for (i = 0; i < nthreads; i++)
    if (pthread_create(&th[i], NULL, worker, &args[i]) == 0) started++;
    else break;

  if (started < nthreads) worker(&args[0]);

  for (i = 0; i < started; i++) pthread_join(th[i], NULL);

  pthread_mutex_destroy(&mtx);
}
