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
 * See tfce_maxtree.h for the algorithm, and tfce_batch.h for the threading. The
 * batch driver is shared with the C API that the Python binding sits on, so both
 * run the very same code rather than two copies of it that could drift apart.
 *
 * Christian Gaser
 */

#include "mex.h"
#include "tfce_maxtree.h"
#include "tfce_batch.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double E, H;
  int calc_neg = 1, N, B, nthreads, is_single, ok;
  size_t NB;
  const tfce_val *Tin;
  tfce_val *outData, *inBuf = NULL, *outBuf = NULL;
  Neigh nb;
  int *aptr = NULL, *aidx = NULL;

  if (nrhs < 5)
    mexErrMsgTxt("Usage: tfce = tfceMex_maxtree_batch(T, E, H, calc_neg, geom [, n_threads])");

  /* The core works in single precision. A single block of maps is used as it
     stands, which is the point of the exercise: it is half the memory traffic,
     and the max-tree is bound by bandwidth. A double block is converted, and the
     result comes back in the class it went in as. */
  is_single = mxIsSingle(prhs[0]);
  if ((!mxIsDouble(prhs[0]) && !is_single) || mxIsComplex(prhs[0]))
    mexErrMsgTxt("T must be a real single or double matrix.");
  if (mxGetNumberOfDimensions(prhs[0]) != 2)
    mexErrMsgTxt("T must be an N x B matrix, one map per column.");

  N = (int) mxGetM(prhs[0]);
  B = (int) mxGetN(prhs[0]);
  E = mxGetScalar(prhs[1]);
  H = mxGetScalar(prhs[2]);
  calc_neg = (int) mxGetScalar(prhs[3]);
  NB = (size_t)N * (size_t)B;

  memset(&nb, 0, sizeof(nb));

  /* geom is either [nx ny nz] or an F x 3 face list */
  if (mxGetNumberOfElements(prhs[4]) == 3 && mxGetM(prhs[4]) == 1) {
    const double *g = mxGetPr(prhs[4]);
    nb.is_mesh = 0;
    nb.nx = (int)g[0]; nb.ny = (int)g[1]; nb.nz = (int)g[2];
    if (nb.nx * nb.ny * nb.nz != N)
      mexErrMsgTxt("prod(geom) does not match the number of rows of T.");
  } else if (mxGetN(prhs[4]) == 3) {
    int nF = (int) mxGetM(prhs[4]), *fidx;
    const char *err;
    fidx = faces_to_int(prhs[4], nF, &err);
    if (!fidx) mexErrMsgTxt(err);
    if (!build_adjacency(fidx, nF, N, &aptr, &aidx)) {
      free(fidx);
      mexErrMsgTxt("faces name a vertex that is not in the data.");
    }
    free(fidx);
    nb.is_mesh = 1; nb.ptr = aptr; nb.idx = aidx;
  } else {
    mexErrMsgTxt("geom must be [nx ny nz] or an F x 3 face list.");
  }

  nthreads = B;   /* tfce_batch_run caps this; <= 0 would mean one per map */
  if (nrhs > 5 && !mxIsEmpty(prhs[5])) nthreads = (int) mxGetScalar(prhs[5]);

  if (is_single) {
    Tin = (const tfce_val *) mxGetData(prhs[0]);
  } else {
    const double *dd = mxGetPr(prhs[0]);
    size_t k;
    inBuf = (tfce_val *) malloc(NB * sizeof(tfce_val));
    if (!inBuf) mexErrMsgTxt("Memory allocation error.");
    for (k = 0; k < NB; k++) inBuf[k] = (tfce_val) dd[k];
    Tin = inBuf;
  }

  plhs[0] = mxCreateNumericMatrix((mwSize)N, (mwSize)B,
                                  is_single ? mxSINGLE_CLASS : mxDOUBLE_CLASS,
                                  mxREAL);

  if (is_single) {
    outData = (tfce_val *) mxGetData(plhs[0]);
  } else {
    outBuf = (tfce_val *) malloc(NB * sizeof(tfce_val));
    if (!outBuf) { free(inBuf); mexErrMsgTxt("Memory allocation error."); }
    outData = outBuf;
  }

  /* the driver itself is shared with the C API, see tfce_batch.h */
  ok = tfce_batch_run(Tin, outData, N, B, &nb, E, H, calc_neg, nthreads);

  if (aptr) { free(aptr); free(aidx); }

  if (!is_single) {                      /* hand the result back as double */
    double *od = mxGetPr(plhs[0]);
    size_t k;
    for (k = 0; k < NB; k++) od[k] = (double) outBuf[k];
  }

  if (inBuf)  free(inBuf);
  if (outBuf) free(outBuf);

  if (!ok) mexErrMsgTxt("Memory allocation error in worker thread.");
}
