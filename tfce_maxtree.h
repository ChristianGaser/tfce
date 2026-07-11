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
  VI *ord;
  int *parent, *sz, *head, *tail;   /* union-find, per element */
  int *nd_next;                     /* merge lists, per node   */
  int *elem_node, *root_stamp;
  int *nd_parent, *nd_size;
  double *nd_birth, *nd_death, *cum;
  char *active;
  double *neg;                      /* buffer for the negative pass */
} Workspace;

static int ws_alloc(Workspace *w, int N)
{
  int i;
  memset(w, 0, sizeof(*w));
  w->N     = N;
  w->level = 0;

  w->ord        = (VI *)     malloc((size_t)N * sizeof(VI));
  w->parent     = (int *)    malloc((size_t)N * sizeof(int));
  w->sz         = (int *)    malloc((size_t)N * sizeof(int));
  w->head       = (int *)    malloc((size_t)N * sizeof(int));
  w->tail       = (int *)    malloc((size_t)N * sizeof(int));
  w->nd_next    = (int *)    malloc((size_t)N * sizeof(int));
  w->elem_node  = (int *)    malloc((size_t)N * sizeof(int));
  w->root_stamp = (int *)    malloc((size_t)N * sizeof(int));
  w->nd_parent  = (int *)    malloc((size_t)N * sizeof(int));
  w->nd_size    = (int *)    malloc((size_t)N * sizeof(int));
  w->nd_birth   = (double *) malloc((size_t)N * sizeof(double));
  w->nd_death   = (double *) malloc((size_t)N * sizeof(double));
  w->cum        = (double *) malloc((size_t)N * sizeof(double));
  w->active     = (char *)   calloc((size_t)N, sizeof(char));
  w->neg        = (double *) malloc((size_t)N * sizeof(double));

  if (!w->ord || !w->parent || !w->sz || !w->head || !w->tail || !w->nd_next ||
      !w->elem_node || !w->root_stamp || !w->nd_parent || !w->nd_size ||
      !w->nd_birth || !w->nd_death || !w->cum || !w->active || !w->neg)
    return 0;

  for (i = 0; i < N; i++) w->root_stamp[i] = -1;
  return 1;
}

static void ws_free(Workspace *w)
{
  free(w->ord); free(w->parent); free(w->sz); free(w->head); free(w->tail);
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
static void maxtree_pass(Workspace *w, const double *d, double *out,
                         const Neigh *nb, double E, double H, double sign)
{
  int N = w->N;
  int i, g, gs, ge, nPos = 0, nNodes = 0;
  int nbuf[MAX_NEIGH];
  double Hp1 = H + 1.0;

  for (i = 0, g = 0; i < N; i++)
    if (d[i] > 0.0) { w->ord[g].val = d[i]; w->ord[g].idx = i; g++; }
  nPos = g;
  if (nPos == 0) return;

  qsort(w->ord, (size_t)nPos, sizeof(VI), cmp_desc);

  gs = 0;
  while (gs < nPos) {
    double h = w->ord[gs].val;

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
      int n = neigh_get(nb, u, nbuf), j;
      for (j = 0; j < n; j++)
        if (w->active[nbuf[j]]) uf_union(w, u, nbuf[j]);
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
        w->nd_death[n]   = 0.0;         /* overwritten if it later merges */

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

  /* root->leaf accumulation: a node's parent always has a larger index */
  for (i = nNodes - 1; i >= 0; i--) {
    double c = pow((double)w->nd_size[i], E) *
               (pow(w->nd_birth[i], Hp1) - pow(w->nd_death[i], Hp1)) / Hp1;
    w->cum[i] = c + (w->nd_parent[i] >= 0 ? w->cum[w->nd_parent[i]] : 0.0);
  }

  for (g = 0; g < nPos; g++) {
    int u = w->ord[g].idx;
    out[u] += sign * w->cum[w->elem_node[u]];
    w->active[u] = 0;                  /* clear only what this pass touched */
  }
}

/* TFCE of one map, both signs. `out` is overwritten. */
static void maxtree_map(Workspace *w, const double *d, double *out,
                        const Neigh *nb, double E, double H, int calc_neg)
{
  int i, N = w->N;

  for (i = 0; i < N; i++) out[i] = 0.0;

  maxtree_pass(w, d, out, nb, E, H, 1.0);

  if (calc_neg) {
    for (i = 0; i < N; i++) w->neg[i] = -d[i];
    maxtree_pass(w, w->neg, out, nb, E, H, -1.0);
  }
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
}

#endif /* TFCE_MAXTREE_H */
