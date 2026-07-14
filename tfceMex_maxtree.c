/*
 * Exact TFCE via a max-tree (component tree) built with union-find.
 *
 *   tfce = tfceMex_maxtree(t, E, H, calc_neg)          % volume, t is 3D
 *   tfce = tfceMex_maxtree(t, E, H, calc_neg, faces)   % surface, t is Nx1
 *
 * The last argument selects the neighbourhood: if faces are omitted or empty, t
 * is treated as a 3D volume with 26-connectivity, otherwise as surface data on
 * the mesh defined by these faces.
 *
 * The TFCE integral is evaluated exactly, so there is no step size dh and no
 * discretisation error. See tfce_maxtree.h for the algorithm. To compute many
 * maps at once (permutations), use tfceMex_maxtree_batch, which runs one map
 * per thread.
 *
 * The positive and negative passes touch disjoint elements and run on two
 * threads without locking.
 *
 * Christian Gaser
 */

#include "mex.h"
#include "tfce_threads.h"
#include "tfce_maxtree.h"

typedef struct {
  Workspace ws;
  const tfce_val *d;
  tfce_val *out;
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
  const tfce_val *inData;
  tfce_val *outData, *inBuf = NULL, *outBuf = NULL;
  double E, H;
  int calc_neg = 1, i, N, is_single;
  Neigh nb;
  PassArgs ap, an;
  tfce_thread_t th;
  int threaded = 0;
  int *aptr = NULL, *aidx = NULL;

  if (nrhs < 3)
    mexErrMsgTxt("Usage: tfce = tfceMex_maxtree(t, E, H [, calc_neg, faces])");

  /* The core works in single precision. A single map is used as it stands, a
     double one is converted, and the result comes back in the class it went in
     as, so that a caller who has not moved to single sees no change. */
  is_single = mxIsSingle(prhs[0]);
  if ((!mxIsDouble(prhs[0]) && !is_single) || mxIsComplex(prhs[0]))
    mexErrMsgTxt("First argument must be a real single or double array.");

  N = (int) mxGetNumberOfElements(prhs[0]);
  E = mxGetScalar(prhs[1]);
  H = mxGetScalar(prhs[2]);
  if (nrhs > 3) calc_neg = (int) mxGetScalar(prhs[3]);

  if (is_single) {
    inData = (const tfce_val *) mxGetData(prhs[0]);
  } else {
    const double *dd = mxGetPr(prhs[0]);
    inBuf = (tfce_val *) malloc((size_t)N * sizeof(tfce_val));
    if (!inBuf) mexErrMsgTxt("Memory allocation error.");
    for (i = 0; i < N; i++) inBuf[i] = (tfce_val) dd[i];
    inData = inBuf;
  }

  memset(&nb, 0, sizeof(nb));

  if (nrhs > 4 && !mxIsEmpty(prhs[4])) {          /* surface */
    int nF = (int) mxGetM(prhs[4]), *fidx;
    const char *err;
    if (mxGetN(prhs[4]) != 3) mexErrMsgTxt("faces must be F x 3.");
    fidx = faces_to_int(prhs[4], nF, &err);
    if (!fidx) { if (inBuf) free(inBuf); mexErrMsgTxt(err); }
    if (!build_adjacency(fidx, nF, N, &aptr, &aidx)) {
      free(fidx);
      if (inBuf) free(inBuf);
      mexErrMsgTxt("faces name a vertex that is not in the data.");
    }
    free(fidx);
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
                                 mxGetDimensions(prhs[0]),
                                 is_single ? mxSINGLE_CLASS : mxDOUBLE_CLASS,
                                 mxREAL);

  if (is_single) {
    outData = (tfce_val *) mxGetData(plhs[0]);
  } else {
    outBuf = (tfce_val *) malloc((size_t)N * sizeof(tfce_val));
    if (!outBuf) mexErrMsgTxt("Memory allocation error.");
    outData = outBuf;
  }
  for (i = 0; i < N; i++) outData[i] = 0.0f;

  if (!ws_alloc(&ap.ws, N)) mexErrMsgTxt("Memory allocation error.");
  ap.d = inData; ap.out = outData; ap.nb = &nb;
  ap.E = E; ap.H = H; ap.sign = 1.0;

  if (calc_neg) {
    if (!ws_alloc(&an.ws, N)) mexErrMsgTxt("Memory allocation error.");
    for (i = 0; i < N; i++) an.ws.neg[i] = -inData[i];
    an.d = an.ws.neg; an.out = outData; an.nb = &nb;
    an.E = E; an.H = H; an.sign = -1.0;

    /* the two passes write to disjoint elements: no lock needed */
    threaded = (tfce_thread_create(&th, pass_thread, &an) == 0);
    maxtree_pass(&ap.ws, ap.d, ap.out, ap.nb, E, H, 1.0);
    if (threaded) tfce_thread_join(th);
    else          maxtree_pass(&an.ws, an.d, an.out, an.nb, E, H, -1.0);

    ws_free(&an.ws);
  } else {
    maxtree_pass(&ap.ws, ap.d, ap.out, ap.nb, E, H, 1.0);
  }

  ws_free(&ap.ws);

  if (!is_single) {                      /* hand the result back as double */
    double *od = mxGetPr(plhs[0]);
    for (i = 0; i < N; i++) od[i] = (double) outBuf[i];
  }

  if (inBuf)  free(inBuf);
  if (outBuf) free(outBuf);
  if (aptr) { free(aptr); free(aidx); }
}
