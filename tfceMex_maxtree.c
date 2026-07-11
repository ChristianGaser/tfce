/*
 * Exact TFCE via a max-tree (component tree) built with union-find.
 *
 *   tfce = tfceMex_maxtree(t, E, H, calc_neg)          % volume, t is 3D
 *   tfce = tfceMex_maxtree(t, E, H, calc_neg, faces)   % surface, t is Nx1
 *
 * Computes  TFCE(v) = int_0^{t_v} e_v(h)^E h^H dh  exactly.  The extent
 * function e_v(h) is piecewise constant in h, so the integral is evaluated in
 * closed form over each constant piece rather than sampled on a dh grid.
 *
 * Voxels/vertices are swept in descending order and merged with their active
 * neighbours.  Each max-tree node is a component of fixed size s alive over
 * [death, birth], contributing  s^E * (birth^(H+1) - death^(H+1)) / (H+1)
 * to every element beneath it.  A single root->leaf accumulation distributes
 * these.  Cost is O(N log N) for the sort plus O(N a(N)) for the sweep, and is
 * independent of any precision parameter.
 *
 * Volume neighbourhood is 26-connected, matching tfceMex_pthread.c.
 * The positive and negative passes touch disjoint elements and run on two
 * threads without locking.
 *
 * Christian Gaser
 */

#include "mex.h"
#include "math.h"
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

/* ------------------------------------------------------------------ *
 * Neighbourhood: either an implicit 26-connected grid or an explicit
 * CSR adjacency (surface mesh).  The max-tree core sees only this.
 * ------------------------------------------------------------------ */
typedef struct {
  int is_mesh;
  int nx, ny, nz;          /* grid  */
  const int *ptr, *idx;    /* CSR   */
} Neigh;

#define MAX_NEIGH 26

static int neigh_get(const Neigh *nb, int u, int *buf)
{
  int n = 0;

  if (nb->is_mesh) {
    int s = nb->ptr[u], e = nb->ptr[u + 1], i;
    for (i = s; i < e; i++) buf[n++] = nb->idx[i];
    return n;
  } else {
    int nx = nb->nx, ny = nb->ny, nz = nb->nz, nxy = nx * ny;
    int uz  = u / nxy, rem = u - uz * nxy;
    int uy  = rem / nx, ux = rem - uy * nx;
    int x0 = ux > 0 ? ux - 1 : 0, x1 = ux < nx - 1 ? ux + 1 : nx - 1;
    int y0 = uy > 0 ? uy - 1 : 0, y1 = uy < ny - 1 ? uy + 1 : ny - 1;
    int z0 = uz > 0 ? uz - 1 : 0, z1 = uz < nz - 1 ? uz + 1 : nz - 1;
    int tx, ty, tz;
    for (tz = z0; tz <= z1; tz++)
      for (ty = y0; ty <= y1; ty++)
        for (tx = x0; tx <= x1; tx++) {
          int w = tz * nxy + ty * nx + tx;
          if (w != u) buf[n++] = w;
        }
    return n;
  }
}

/* ------------------------------------------------------------------ *
 * Sort
 * ------------------------------------------------------------------ */
typedef struct { double val; int idx; } VI;

static int cmp_desc(const void *a, const void *b)
{
  double x = ((const VI *)a)->val, y = ((const VI *)b)->val;
  if (x < y) return  1;
  if (x > y) return -1;
  return 0;
}

/* ------------------------------------------------------------------ *
 * Union-find, path halving + union by size.  Each component root also
 * carries a singly-linked list of the max-tree nodes that have merged
 * into it at the current level, spliced in O(1) on union.
 * ------------------------------------------------------------------ */
typedef struct {
  int *parent, *sz, *head, *tail;   /* per element */
  int *nd_next;                     /* per node    */
} UF;

static int uf_find(int *parent, int x)
{
  while (parent[x] != x) {
    parent[x] = parent[parent[x]];
    x = parent[x];
  }
  return x;
}

static void uf_union(UF *uf, int a, int b)
{
  a = uf_find(uf->parent, a);
  b = uf_find(uf->parent, b);
  if (a == b) return;
  if (uf->sz[a] < uf->sz[b]) { int t = a; a = b; b = t; }
  uf->parent[b] = a;
  uf->sz[a] += uf->sz[b];

  if (uf->head[b] >= 0) {                       /* splice b's node list into a */
    if (uf->head[a] < 0) {
      uf->head[a] = uf->head[b];
      uf->tail[a] = uf->tail[b];
    } else {
      uf->nd_next[uf->tail[a]] = uf->head[b];
      uf->tail[a] = uf->tail[b];
    }
    uf->head[b] = uf->tail[b] = -1;
  }
}

/* ------------------------------------------------------------------ *
 * Core: one signed pass over the positive part of `d`.
 * ------------------------------------------------------------------ */
typedef struct {
  const double *d;
  double *out;
  const Neigh *nb;
  int N;
  double E, H, sign;
} PassArgs;

static void *maxtree_pass(void *pa)
{
  PassArgs *A = (PassArgs *)pa;
  const double *d = A->d;
  const Neigh *nb = A->nb;
  int N = A->N;
  double E = A->E, Hp1 = A->H + 1.0, sign = A->sign;

  int i, g, gs, ge, nPos = 0, nNodes = 0, level = 0;
  int nbuf[MAX_NEIGH];
  VI *ord;
  UF uf;
  int *comp_node, *elem_node, *root_stamp;
  int *nd_parent, *nd_size;
  double *nd_birth, *nd_death, *cum;
  char *active;

  for (i = 0; i < N; i++) if (d[i] > 0.0) nPos++;
  if (nPos == 0) return NULL;

  ord = (VI *) malloc((size_t)nPos * sizeof(VI));
  for (i = 0, g = 0; i < N; i++)
    if (d[i] > 0.0) { ord[g].val = d[i]; ord[g].idx = i; g++; }
  qsort(ord, (size_t)nPos, sizeof(VI), cmp_desc);

  uf.parent  = (int *) malloc((size_t)N * sizeof(int));
  uf.sz      = (int *) malloc((size_t)N * sizeof(int));
  uf.head    = (int *) malloc((size_t)N * sizeof(int));
  uf.tail    = (int *) malloc((size_t)N * sizeof(int));
  uf.nd_next = (int *) malloc((size_t)nPos * sizeof(int));

  comp_node  = (int *) malloc((size_t)N * sizeof(int));
  elem_node  = (int *) malloc((size_t)N * sizeof(int));
  root_stamp = (int *) malloc((size_t)N * sizeof(int));
  active     = (char *)calloc((size_t)N, sizeof(char));

  nd_parent = (int *)    malloc((size_t)nPos * sizeof(int));
  nd_size   = (int *)    malloc((size_t)nPos * sizeof(int));
  nd_birth  = (double *) malloc((size_t)nPos * sizeof(double));
  nd_death  = (double *) malloc((size_t)nPos * sizeof(double));
  cum       = (double *) malloc((size_t)nPos * sizeof(double));

  for (i = 0; i < N; i++) root_stamp[i] = -1;

  gs = 0;
  while (gs < nPos) {
    double h = ord[gs].val;

    ge = gs;
    while (ge < nPos && ord[ge].val == h) ge++;   /* ties share one level */
    level++;

    /* activate this level as singletons */
    for (g = gs; g < ge; g++) {
      int u = ord[g].idx;
      active[u]   = 1;
      uf.parent[u] = u;
      uf.sz[u]     = 1;
      uf.head[u]   = uf.tail[u] = -1;
    }

    /* merge with active neighbours; old node lists splice along for free */
    for (g = gs; g < ge; g++) {
      int u = ord[g].idx;
      int n = neigh_get(nb, u, nbuf), j;
      for (j = 0; j < n; j++)
        if (active[nbuf[j]]) uf_union(&uf, u, nbuf[j]);
    }

    /* one fresh node per component that changed at this level */
    for (g = gs; g < ge; g++) {
      int r = uf_find(uf.parent, ord[g].idx);
      if (root_stamp[r] != level) {
        int n = nNodes++, o;
        root_stamp[r] = level;
        nd_parent[n]  = -1;
        nd_size[n]    = uf.sz[r];
        nd_birth[n]   = h;
        nd_death[n]   = 0.0;         /* overwritten if it later merges */

        /* every old component in this root dies here, under n */
        for (o = uf.head[r]; o >= 0; o = uf.nd_next[o]) {
          nd_parent[o] = n;
          nd_death[o]  = h;
        }
        uf.head[r] = uf.tail[r] = n;
        uf.nd_next[n] = -1;
        comp_node[r]  = n;
      }
    }

    for (g = gs; g < ge; g++) {
      int u = ord[g].idx;
      elem_node[u] = comp_node[uf_find(uf.parent, u)];
    }

    gs = ge;
  }

  /* root->leaf accumulation: a node's parent always has a larger index */
  for (i = nNodes - 1; i >= 0; i--) {
    double c = pow((double)nd_size[i], E) *
               (pow(nd_birth[i], Hp1) - pow(nd_death[i], Hp1)) / Hp1;
    cum[i] = c + (nd_parent[i] >= 0 ? cum[nd_parent[i]] : 0.0);
  }

  for (g = 0; g < nPos; g++) {
    int u = ord[g].idx;
    A->out[u] += sign * cum[elem_node[u]];
  }

  free(ord); free(uf.parent); free(uf.sz); free(uf.head); free(uf.tail);
  free(uf.nd_next); free(comp_node); free(elem_node); free(root_stamp);
  free(active); free(nd_parent); free(nd_size); free(nd_birth); free(nd_death);
  free(cum);
  return NULL;
}

/* ------------------------------------------------------------------ *
 * CSR adjacency from a face list (1-based, F x 3)
 * ------------------------------------------------------------------ */
static void build_adjacency(const double *faces, int nF, int N, int **ptr, int **idx)
{
  int *deg = (int *) calloc((size_t)N + 1, sizeof(int));
  int *p, *ix, *fill;
  int f, i, j, e, total;
  static const int ea[3] = {0, 1, 2}, eb[3] = {1, 2, 0};

  for (f = 0; f < nF; f++)
    for (e = 0; e < 3; e++) {
      int a = (int)faces[ea[e] * nF + f] - 1;
      int b = (int)faces[eb[e] * nF + f] - 1;
      deg[a]++; deg[b]++;
    }

  p = (int *) malloc(((size_t)N + 1) * sizeof(int));
  p[0] = 0;
  for (i = 0; i < N; i++) p[i + 1] = p[i] + deg[i];
  total = p[N];

  ix   = (int *) malloc((size_t)total * sizeof(int));
  fill = (int *) malloc((size_t)N * sizeof(int));
  memcpy(fill, p, (size_t)N * sizeof(int));

  for (f = 0; f < nF; f++)
    for (e = 0; e < 3; e++) {
      int a = (int)faces[ea[e] * nF + f] - 1;
      int b = (int)faces[eb[e] * nF + f] - 1;
      ix[fill[a]++] = b;
      ix[fill[b]++] = a;
    }

  /* dedupe each row (insertion sort: vertex degrees are small) */
  for (i = 0; i < N; i++) {
    int s = p[i], n = fill[i] - s, k = 0;
    for (j = 1; j < n; j++) {
      int v = ix[s + j], m = j - 1;
      while (m >= 0 && ix[s + m] > v) { ix[s + m + 1] = ix[s + m]; m--; }
      ix[s + m + 1] = v;
    }
    for (j = 0; j < n; j++)
      if (j == 0 || ix[s + j] != ix[s + j - 1]) ix[s + k++] = ix[s + j];
    deg[i] = k;
  }
  /* compact */
  {
    int w = 0;
    int *p2 = (int *) malloc(((size_t)N + 1) * sizeof(int));
    p2[0] = 0;
    for (i = 0; i < N; i++) {
      int s = p[i];
      for (j = 0; j < deg[i]; j++) ix[w++] = ix[s + j];
      p2[i + 1] = w;
    }
    free(p);
    p = p2;
  }

  free(deg); free(fill);
  *ptr = p; *idx = ix;
}

/* ------------------------------------------------------------------ */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *inData, *outData, *neg = NULL, E, H;
  int calc_neg = 1, i, N;
  Neigh nb;
  PassArgs ap, an;
  pthread_t th;
  int *aptr = NULL, *aidx = NULL;

  if (nrhs < 3)
    mexErrMsgTxt("Usage: tfce = tfceMex_maxtree(t, E, H [, calc_neg, faces])");
  if (!mxIsDouble(prhs[0]))
    mexErrMsgTxt("First argument must be double.");

  inData = mxGetPr(prhs[0]);
  N      = (int) mxGetNumberOfElements(prhs[0]);
  E      = mxGetScalar(prhs[1]);
  H      = mxGetScalar(prhs[2]);
  if (nrhs > 3) calc_neg = (int) mxGetScalar(prhs[3]);

  memset(&nb, 0, sizeof(nb));

  if (nrhs > 4 && !mxIsEmpty(prhs[4])) {          /* surface */
    const double *faces = mxGetPr(prhs[4]);
    int nF = (int) mxGetM(prhs[4]);
    if (mxGetN(prhs[4]) != 3) mexErrMsgTxt("faces must be F x 3.");
    build_adjacency(faces, nF, N, &aptr, &aidx);
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

  ap.d = inData; ap.out = outData; ap.nb = &nb; ap.N = N;
  ap.E = E; ap.H = H; ap.sign = 1.0;

  if (calc_neg) {
    neg = (double *) mxMalloc((size_t)N * sizeof(double));
    for (i = 0; i < N; i++) neg[i] = -inData[i];
    an = ap; an.d = neg; an.sign = -1.0;

    /* the two passes write to disjoint elements: no lock needed */
    if (pthread_create(&th, NULL, maxtree_pass, &an) == 0) {
      maxtree_pass(&ap);
      pthread_join(th, NULL);
    } else {
      maxtree_pass(&ap);
      maxtree_pass(&an);
    }
    mxFree(neg);
  } else {
    maxtree_pass(&ap);
  }

  if (aptr) { free(aptr); free(aidx); }
}
