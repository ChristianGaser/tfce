/*
 * Exact TFCE via a max-tree (component tree) built with union-find.
 *
 * Shared core for tfceMex_maxtree (one map) and tfceMex_maxtree_batch (many
 * maps in parallel). The TFCE integral
 *
 *     TFCE(v) = int_0^{t_v} e_v(h)^E h^H dh
 *
 * is evaluated exactly rather than sampled on a dh grid: the extent function
 * e_v(h) is piecewise constant, so the integral is computed in closed form over
 * each constant piece. Elements are swept in descending order and merged with
 * their active neighbours; every max-tree node is a component of fixed size s
 * alive over [death, birth] and contributes
 *
 *     s^E * (birth^(H+1) - death^(H+1)) / (H+1)
 *
 * to every element beneath it. A single root->leaf accumulation distributes
 * these. Cost is O(N log N) for the sort plus O(N a(N)) for the sweep, and is
 * independent of any precision parameter.
 *
 * The neighbourhood is abstracted, so the same core serves 26-connected volumes
 * and surface meshes.
 *
 * Christian Gaser
 */

#ifndef TFCE_MAXTREE_H
#define TFCE_MAXTREE_H

#include "math.h"
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ *
 * Neighbourhood: implicit 26-connected grid, or explicit CSR adjacency
 * ------------------------------------------------------------------ */
typedef struct {
  int is_mesh;
  int nx, ny, nz;          /* grid  */
  const int *ptr, *idx;    /* CSR   */
} Neigh;

/* a 3x3x3 neighbourhood has 26 elements around its centre, and that is the most
   the grid case can ever produce. It says nothing about a mesh, where the degree
   of a vertex is not bounded by anything. */
#define MAX_NEIGH 26

/*
 * Neighbours of element u. Returns how many there are and, in *out, where they
 * are.
 *
 * A mesh is NOT copied into buf: the adjacency is already a contiguous run in the
 * CSR index array, so *out simply points into it. Copying it would need a buffer
 * as large as the largest vertex degree in the mesh, which nothing bounds -- an
 * earlier version copied into a buffer of MAX_NEIGH and smashed the stack on any
 * vertex with more than 26 neighbours.
 */
static int neigh_get(const Neigh *nb, int u, int *buf, const int **out)
{
  int n = 0;

  if (nb->is_mesh) {
    *out = nb->idx + nb->ptr[u];
    return nb->ptr[u + 1] - nb->ptr[u];
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
    *out = buf;
    return n;
  }
}

/* ------------------------------------------------------------------ *
 * Precision
 *
 * The statistic itself, and everything derived from it that is N elements long,
 * is carried as float. The map is a permuted statistic, single precision is all
 * it is ever worth, and the max-tree is bound by memory bandwidth rather than by
 * arithmetic: the sort array halves, the working set shrinks by about a third,
 * and the whole transform speeds up accordingly.
 *
 * The integral itself is NOT computed in float. `cum` accumulates a node's
 * contribution along a root-to-leaf path, which can be long, and pow(birth, H+1)
 * with H = 2 spans a wide dynamic range. Those stay double, so only the values
 * that are stored per element lose precision, not the arithmetic that combines
 * them.
 * ------------------------------------------------------------------ */
typedef float tfce_val;

/* ------------------------------------------------------------------ *
 * Sort
 * ------------------------------------------------------------------ */
typedef struct { tfce_val val; int idx; } VI;

/*
 * Sort descending by value.
 *
 * This is the most expensive single step of the transform, and it is a sort of
 * hundreds of thousands of elements by a float key -- exactly the case a radix
 * sort is for. Every value reaching here is strictly positive, and for positive
 * IEEE-754 floats the bit pattern read as a uint32 increases monotonically with
 * the value, so the key needs no transformation at all. Four counting passes over
 * the array replace n log n comparisons through a function pointer.
 *
 * The sort is stable, so ties come out in a fixed order. Which order does not
 * matter: ties are consumed as one level, and the union-find partition they
 * produce does not depend on the order they are merged in.
 */
static void sort_desc(VI *a, VI *tmp, int n)
{
  int pass, i;
  unsigned int cnt[256], pos[256], k;
  VI *src = a, *dst = tmp, *sw;

  for (pass = 0; pass < 4; pass++) {
    int shift = pass * 8;

    memset(cnt, 0, sizeof(cnt));
    for (i = 0; i < n; i++) {
      unsigned int key;
      memcpy(&key, &src[i].val, sizeof(key));
      cnt[(key >> shift) & 0xFFu]++;
    }

    pos[0] = 0;
    for (k = 1; k < 256; k++) pos[k] = pos[k - 1] + cnt[k - 1];

    for (i = 0; i < n; i++) {
      unsigned int key;
      memcpy(&key, &src[i].val, sizeof(key));
      dst[pos[(key >> shift) & 0xFFu]++] = src[i];
    }

    sw = src; src = dst; dst = sw;
  }

  /* four passes is an even number of swaps, so the result is back in `a`,
     ascending. The sweep wants it descending. */
  for (i = 0; i < n / 2; i++) {
    VI t = a[i]; a[i] = a[n - 1 - i]; a[n - 1 - i] = t;
  }
}

/* ------------------------------------------------------------------ *
 * Workspace. Allocated once per thread and reused across maps, so a
 * permutation run does not malloc/free a dozen N-sized arrays per map.
 *
 * `level` is never reset between passes: it only ever increases, so stale
 * root_stamp entries from an earlier map always compare as older and the
 * O(N) reset can be skipped. `active` is cleared at the end of each pass,
 * touching only the elements the pass actually set.
 * ------------------------------------------------------------------ */
typedef struct {
  int N;
  int level;
  VI *ord, *ord_tmp;                 /* ord_tmp is the radix sort's scratch */
  int *parent, *sz, *head, *tail;   /* union-find, per element */
  int *nd_next;                     /* merge lists, per node   */
  int *elem_node, *root_stamp;
  int *nd_parent, *nd_size;
  tfce_val *nd_birth, *nd_death;
  double *cum;                      /* the integral is accumulated in double */
  char *active;
  tfce_val *neg;                    /* buffer for the negative pass */
} Workspace;

static int ws_alloc(Workspace *w, int N)
{
  int i;
  memset(w, 0, sizeof(*w));
  w->N     = N;
  w->level = 0;

  w->ord        = (VI *)     malloc((size_t)N * sizeof(VI));
  w->ord_tmp    = (VI *)     malloc((size_t)N * sizeof(VI));
  w->parent     = (int *)    malloc((size_t)N * sizeof(int));
  w->sz         = (int *)    malloc((size_t)N * sizeof(int));
  w->head       = (int *)    malloc((size_t)N * sizeof(int));
  w->tail       = (int *)    malloc((size_t)N * sizeof(int));
  w->nd_next    = (int *)    malloc((size_t)N * sizeof(int));
  w->elem_node  = (int *)    malloc((size_t)N * sizeof(int));
  w->root_stamp = (int *)    malloc((size_t)N * sizeof(int));
  w->nd_parent  = (int *)    malloc((size_t)N * sizeof(int));
  w->nd_size    = (int *)    malloc((size_t)N * sizeof(int));
  w->nd_birth   = (tfce_val *) malloc((size_t)N * sizeof(tfce_val));
  w->nd_death   = (tfce_val *) malloc((size_t)N * sizeof(tfce_val));
  w->cum        = (double *)   malloc((size_t)N * sizeof(double));
  w->active     = (char *)     calloc((size_t)N, sizeof(char));
  w->neg        = (tfce_val *) malloc((size_t)N * sizeof(tfce_val));

  if (!w->ord || !w->ord_tmp || !w->parent || !w->sz || !w->head || !w->tail || !w->nd_next ||
      !w->elem_node || !w->root_stamp || !w->nd_parent || !w->nd_size ||
      !w->nd_birth || !w->nd_death || !w->cum || !w->active || !w->neg)
    return 0;

  for (i = 0; i < N; i++) w->root_stamp[i] = -1;
  return 1;
}

static void ws_free(Workspace *w)
{
  free(w->ord); free(w->ord_tmp);
  free(w->parent); free(w->sz); free(w->head); free(w->tail);
  free(w->nd_next); free(w->elem_node); free(w->root_stamp);
  free(w->nd_parent); free(w->nd_size); free(w->nd_birth); free(w->nd_death);
  free(w->cum); free(w->active); free(w->neg);
  memset(w, 0, sizeof(*w));
}

/* ------------------------------------------------------------------ *
 * Union-find, path halving + union by size. Each root carries a singly
 * linked list of the max-tree nodes merged into it at the current level,
 * spliced in O(1) on union.
 * ------------------------------------------------------------------ */
static int uf_find(int *parent, int x)
{
  while (parent[x] != x) {
    parent[x] = parent[parent[x]];
    x = parent[x];
  }
  return x;
}

static void uf_union(Workspace *w, int a, int b)
{
  a = uf_find(w->parent, a);
  b = uf_find(w->parent, b);
  if (a == b) return;
  if (w->sz[a] < w->sz[b]) { int t = a; a = b; b = t; }
  w->parent[b] = a;
  w->sz[a] += w->sz[b];

  if (w->head[b] >= 0) {                       /* splice b's node list into a */
    if (w->head[a] < 0) {
      w->head[a] = w->head[b];
      w->tail[a] = w->tail[b];
    } else {
      w->nd_next[w->tail[a]] = w->head[b];
      w->tail[a] = w->tail[b];
    }
    w->head[b] = w->tail[b] = -1;
  }
}

/* ------------------------------------------------------------------ *
 * One signed pass over the positive part of `d`, accumulated into `out`
 * with weight `sign`.
 * ------------------------------------------------------------------ */
static void maxtree_pass(Workspace *w, const tfce_val *d, tfce_val *out,
                         const Neigh *nb, double E, double H, double sign)
{
  int N = w->N;
  int i, g, gs, ge, nPos = 0, nNodes = 0;
  int nbuf[MAX_NEIGH];
  double Hp1 = H + 1.0;

  for (i = 0, g = 0; i < N; i++)
    if (d[i] > 0.0f) { w->ord[g].val = d[i]; w->ord[g].idx = i; g++; }
  nPos = g;
  if (nPos == 0) return;

  sort_desc(w->ord, w->ord_tmp, nPos);

  gs = 0;
  while (gs < nPos) {
    tfce_val h = w->ord[gs].val;

    ge = gs;
    while (ge < nPos && w->ord[ge].val == h) ge++;   /* ties share one level */
    w->level++;

    /* activate this level as singletons */
    for (g = gs; g < ge; g++) {
      int u = w->ord[g].idx;
      w->active[u] = 1;
      w->parent[u] = u;
      w->sz[u]     = 1;
      w->head[u]   = w->tail[u] = -1;
    }

    /* merge with active neighbours; old node lists splice along for free */
    for (g = gs; g < ge; g++) {
      int u = w->ord[g].idx;
      const int *nl;
      int n = neigh_get(nb, u, nbuf, &nl), j;
      for (j = 0; j < n; j++)
        if (w->active[nl[j]]) uf_union(w, u, nl[j]);
    }

    /* one fresh node per component that changed at this level */
    for (g = gs; g < ge; g++) {
      int r = uf_find(w->parent, w->ord[g].idx);
      if (w->root_stamp[r] != w->level) {
        int n = nNodes++, o;
        w->root_stamp[r] = w->level;
        w->nd_parent[n]  = -1;
        w->nd_size[n]    = w->sz[r];
        w->nd_birth[n]   = h;
        w->nd_death[n]   = 0.0f;        /* overwritten if it later merges */

        /* every old component in this root dies here, under n */
        for (o = w->head[r]; o >= 0; o = w->nd_next[o]) {
          w->nd_parent[o] = n;
          w->nd_death[o]  = h;
        }
        w->head[r] = w->tail[r] = n;
        w->nd_next[n] = -1;
      }
    }

    for (g = gs; g < ge; g++) {
      int u = w->ord[g].idx;
      w->elem_node[u] = w->head[uf_find(w->parent, u)];
    }

    gs = ge;
  }

  /* root->leaf accumulation: a node's parent always has a larger index. The
     heights are stored as float, but everything done with them here is double */
  for (i = nNodes - 1; i >= 0; i--) {
    double c = pow((double)w->nd_size[i], E) *
               (pow((double)w->nd_birth[i], Hp1) -
                pow((double)w->nd_death[i], Hp1)) / Hp1;
    w->cum[i] = c + (w->nd_parent[i] >= 0 ? w->cum[w->nd_parent[i]] : 0.0);
  }

  for (g = 0; g < nPos; g++) {
    int u = w->ord[g].idx;
    out[u] += (tfce_val)(sign * w->cum[w->elem_node[u]]);
    w->active[u] = 0;                  /* clear only what this pass touched */
  }
}

/* TFCE of one map, both signs. `out` is overwritten. */
static void maxtree_map(Workspace *w, const tfce_val *d, tfce_val *out,
                        const Neigh *nb, double E, double H, int calc_neg)
{
  int i, N = w->N;

  for (i = 0; i < N; i++) out[i] = 0.0f;

  maxtree_pass(w, d, out, nb, E, H, 1.0);

  if (calc_neg) {
    for (i = 0; i < N; i++) w->neg[i] = -d[i];
    maxtree_pass(w, w->neg, out, nb, E, H, -1.0);
  }
}

/* ------------------------------------------------------------------ *
 * CSR adjacency from a face list (1-based, F x 3, column-major)
 *
 * Returns 0 and builds nothing if any index does not name a vertex that exists.
 * That check is not a formality: an out-of-range index is written straight past
 * the end of the degree array below, which corrupts the heap and takes the
 * process down somewhere else entirely, long after the fact.
 * ------------------------------------------------------------------ */
static int build_adjacency(const int *faces, int nF, int N, int **ptr, int **idx)
{
  int *deg;
  int *p, *ix, *fill;
  int f, i, j, e, total;
  static const int ea[3] = {0, 1, 2}, eb[3] = {1, 2, 0};

  for (i = 0; i < nF * 3; i++)
    if (faces[i] < 1 || faces[i] > N) return 0;

  deg = (int *) calloc((size_t)N + 1, sizeof(int));
  if (!deg) return 0;

  for (f = 0; f < nF; f++)
    for (e = 0; e < 3; e++) {
      int a = faces[ea[e] * nF + f] - 1;
      int b = faces[eb[e] * nF + f] - 1;
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
      int a = faces[ea[e] * nF + f] - 1;
      int b = faces[eb[e] * nF + f] - 1;
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
    int wpos = 0;
    int *p2 = (int *) malloc(((size_t)N + 1) * sizeof(int));
    p2[0] = 0;
    for (i = 0; i < N; i++) {
      int s = p[i];
      for (j = 0; j < deg[i]; j++) ix[wpos++] = ix[s + j];
      p2[i + 1] = wpos;
    }
    free(p);
    p = p2;
  }

  free(deg); free(fill);
  *ptr = p; *idx = ix;
  return 1;
}

#ifdef MATLAB_MEX_FILE
/* ------------------------------------------------------------------ *
 * Read a face list out of MATLAB, whatever numeric class it arrived in.
 *
 * This has to honour the class. GIFTI stores triangles as int32, so SPM hands
 * SPM.xVol.G.faces over as an int32 array -- and reading that with mxGetPr, which
 * assumes double, reinterprets the integer bit patterns as floating point and
 * produces vertex indices that are pure noise. Those indices then get written
 * past the ends of the adjacency arrays. That was a heap corruption on every
 * surface analysis, and it took MATLAB down wherever it happened to notice.
 * ------------------------------------------------------------------ */
static int *faces_to_int(const mxArray *A, int nF, const char **err)
{
  size_t n = (size_t)nF * 3, i;
  int *f;

  *err = NULL;

  if (mxIsComplex(A) || !mxIsNumeric(A)) {
    *err = "faces must be a real numeric F x 3 array.";
    return NULL;
  }

  f = (int *) malloc(n * sizeof(int));
  if (!f) { *err = "Memory allocation error."; return NULL; }

  switch (mxGetClassID(A)) {
    case mxDOUBLE_CLASS:
      { const double   *p = (const double   *)mxGetData(A);
        for (i = 0; i < n; i++) f[i] = (int)p[i]; break; }
    case mxSINGLE_CLASS:
      { const float    *p = (const float    *)mxGetData(A);
        for (i = 0; i < n; i++) f[i] = (int)p[i]; break; }
    case mxINT32_CLASS:
      { const int32_T  *p = (const int32_T  *)mxGetData(A);
        for (i = 0; i < n; i++) f[i] = (int)p[i]; break; }
    case mxUINT32_CLASS:
      { const uint32_T *p = (const uint32_T *)mxGetData(A);
        for (i = 0; i < n; i++) f[i] = (int)p[i]; break; }
    case mxINT16_CLASS:
      { const int16_T  *p = (const int16_T  *)mxGetData(A);
        for (i = 0; i < n; i++) f[i] = (int)p[i]; break; }
    case mxUINT16_CLASS:
      { const uint16_T *p = (const uint16_T *)mxGetData(A);
        for (i = 0; i < n; i++) f[i] = (int)p[i]; break; }
    default:
      free(f);
      *err = "faces must be double, single, int16/uint16 or int32/uint32.";
      return NULL;
  }

  return f;
}
#endif /* MATLAB_MEX_FILE */

#endif /* TFCE_MAXTREE_H */
