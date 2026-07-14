/*
 * Flat C API for the exact TFCE max-tree. See tfce_capi.h.
 *
 * Christian Gaser
 */

#include <stdlib.h>
#include <string.h>

#include "tfce_maxtree.h"
#include "tfce_batch.h"
#include "tfce_capi.h"

static void neigh_volume(Neigh *nb, int nx, int ny, int nz, int connectivity)
{
  memset(nb, 0, sizeof(*nb));
  nb->is_mesh = 0;
  nb->nx = nx; nb->ny = ny; nb->nz = nz;
  nb->conn = connectivity;      /* 0 is taken as 26 by neigh_get */
}

static void neigh_mesh(Neigh *nb, const int *ptr, const int *idx)
{
  memset(nb, 0, sizeof(*nb));
  nb->is_mesh = 1;
  nb->ptr = ptr; nb->idx = idx;
}

static int run_single(const float *data, int N, const Neigh *nb,
                      double E, double H, int calc_neg, float *out)
{
  Workspace ws;

  if (!ws_alloc(&ws, N)) return TFCE_ENOMEM;
  maxtree_map(&ws, data, out, nb, E, H, calc_neg);
  ws_free(&ws);

  return TFCE_OK;
}

int tfce_volume(const float *data, int nx, int ny, int nz, int connectivity,
                double E, double H, int calc_neg, float *out)
{
  Neigh nb;
  neigh_volume(&nb, nx, ny, nz, connectivity);
  return run_single(data, nx*ny*nz, &nb, E, H, calc_neg, out);
}

int tfce_mesh(const float *data, int N,
              const int *adj_ptr, const int *adj_idx,
              double E, double H, int calc_neg, float *out)
{
  Neigh nb;
  neigh_mesh(&nb, adj_ptr, adj_idx);
  return run_single(data, N, &nb, E, H, calc_neg, out);
}

int tfce_volume_batch(const float *data, int nx, int ny, int nz, int connectivity,
                      int B, double E, double H, int calc_neg, int n_threads,
                      float *out)
{
  Neigh nb;
  neigh_volume(&nb, nx, ny, nz, connectivity);
  return tfce_batch_run(data, out, nx*ny*nz, B, &nb, E, H, calc_neg, n_threads)
         ? TFCE_OK : TFCE_ENOMEM;
}

int tfce_mesh_batch(const float *data, int N, int B,
                    const int *adj_ptr, const int *adj_idx,
                    double E, double H, int calc_neg, int n_threads,
                    float *out)
{
  Neigh nb;
  neigh_mesh(&nb, adj_ptr, adj_idx);
  return tfce_batch_run(data, out, N, B, &nb, E, H, calc_neg, n_threads)
         ? TFCE_OK : TFCE_ENOMEM;
}

int tfce_adjacency_from_faces(const int *faces, int nF, int N,
                              int **ptr, int **idx)
{
  return build_adjacency(faces, nF, N, ptr, idx) ? TFCE_OK : TFCE_EFACES;
}

void tfce_free(void *p) { free(p); }
