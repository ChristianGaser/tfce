/*
 * Exact TFCE via a max-tree (component tree) built with union-find.
 *
 *   tfce = tfceMex_maxtree(t, E, H, calc_neg)          % volume, t is 3D
 *   tfce = tfceMex_maxtree(t, E, H, calc_neg, faces)   % surface, t is Nx1
 *
 * In contrast to tfceMex_pthread, the TFCE integral is evaluated exactly, so
 * there is no step size dh and no discretisation error. See tfce_maxtree.h for
 * the algorithm. To compute many maps at once (permutations), use
 * tfceMex_maxtree_batch, which runs one map per thread.
 *
 * The positive and negative passes touch disjoint elements and run on two
 * threads without locking.
 *
 * Christian Gaser
 */

#include "mex.h"
#include <pthread.h>
#include "tfce_maxtree.h"

typedef struct {
  Workspace ws;
  const double *d;
  double *out;
  const Neigh *nb;
  double E, H, sign;
  int ok;
} PassArgs;

static void *pass_thread(void *pa)
{
  PassArgs *A = (PassArgs *)pa;
  maxtree_pass(&A->ws, A->d, A->out, A->nb, A->E, A->H, A->sign);
  return NULL;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *inData, *outData, E, H;
  int calc_neg = 1, i, N;
  Neigh nb;
  PassArgs ap, an;
  pthread_t th;
  int threaded = 0;
  int *aptr = NULL, *aidx = NULL;

  if (nrhs < 3)
    mexErrMsgTxt("Usage: tfce = tfceMex_maxtree(t, E, H [, calc_neg, faces])");
  if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
    mexErrMsgTxt("First argument must be real double.");

  inData = mxGetPr(prhs[0]);
  N      = (int) mxGetNumberOfElements(prhs[0]);
  E      = mxGetScalar(prhs[1]);
  H      = mxGetScalar(prhs[2]);
  if (nrhs > 3) calc_neg = (int) mxGetScalar(prhs[3]);

  memset(&nb, 0, sizeof(nb));

  if (nrhs > 4 && !mxIsEmpty(prhs[4])) {          /* surface */
    if (mxGetN(prhs[4]) != 3) mexErrMsgTxt("faces must be F x 3.");
    build_adjacency(mxGetPr(prhs[4]), (int) mxGetM(prhs[4]), N, &aptr, &aidx);
    nb.is_mesh = 1; nb.ptr = aptr; nb.idx = aidx;
  } else {                                        /* volume */
    const mwSize *dims;
    if (mxGetNumberOfDimensions(prhs[0]) != 3)
      mexErrMsgTxt("Volume input must have 3 dimensions (or pass faces).");
    dims = mxGetDimensions(prhs[0]);
    nb.is_mesh = 0;
    nb.nx = (int)dims[0]; nb.ny = (int)dims[1]; nb.nz = (int)dims[2];
  }

  plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
                                 mxGetDimensions(prhs[0]), mxDOUBLE_CLASS, mxREAL);
  outData = mxGetPr(plhs[0]);
  for (i = 0; i < N; i++) outData[i] = 0.0;

  if (!ws_alloc(&ap.ws, N)) mexErrMsgTxt("Memory allocation error.");
  ap.d = inData; ap.out = outData; ap.nb = &nb;
  ap.E = E; ap.H = H; ap.sign = 1.0;

  if (calc_neg) {
    if (!ws_alloc(&an.ws, N)) mexErrMsgTxt("Memory allocation error.");
    for (i = 0; i < N; i++) an.ws.neg[i] = -inData[i];
    an.d = an.ws.neg; an.out = outData; an.nb = &nb;
    an.E = E; an.H = H; an.sign = -1.0;

    /* the two passes write to disjoint elements: no lock needed */
    threaded = (pthread_create(&th, NULL, pass_thread, &an) == 0);
    maxtree_pass(&ap.ws, ap.d, ap.out, ap.nb, E, H, 1.0);
    if (threaded) pthread_join(th, NULL);
    else          maxtree_pass(&an.ws, an.d, an.out, an.nb, E, H, -1.0);

    ws_free(&an.ws);
  } else {
    maxtree_pass(&ap.ws, ap.d, ap.out, ap.nb, E, H, 1.0);
  }

  ws_free(&ap.ws);

  if (aptr) { free(aptr); free(aidx); }
}
