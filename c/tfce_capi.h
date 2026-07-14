/* SPDX-License-Identifier: BSD-3-Clause */
/* Copyright (c) 2020-2026, Christian Gaser. See LICENSE. */
/*
 * Flat C API for the exact TFCE max-tree.
 *
 * The transform itself is plain C and always was; what tied it to MATLAB was the
 * mexFunction glue wrapped around it. This header exposes it as an ordinary C
 * library instead, so that a binding in any language -- Cython, ctypes, R, Julia
 * -- can call it without going anywhere near mex.h.
 *
 * Everything here works on float32 (the map is a permuted statistic; single
 * precision is all it is ever worth) and returns 0 on success, non-zero on
 * failure. Nothing here allocates anything the caller has to know about, except
 * tfce_adjacency_from_faces, whose two arrays are released with tfce_free.
 *
 * The TFCE integral
 *
 *     TFCE(v) = int_0^{t_v} e_v(h)^E h^H dh
 *
 * is evaluated exactly rather than sampled on a grid of thresholds: the extent
 * function e_v(h) is piecewise constant, so each piece is integrated in closed
 * form. There is no step size dh, and so no discretisation error and no accuracy
 * parameter to get wrong.
 *
 * Christian Gaser
 */

#ifndef TFCE_CAPI_H
#define TFCE_CAPI_H

#ifdef __cplusplus
extern "C" {
#endif

/* return codes */
#define TFCE_OK          0
#define TFCE_ENOMEM      1
#define TFCE_EFACES      2   /* a face names a vertex that is not in the data */

/*
 * TFCE of a single map on a 3D grid.
 *
 * data, out    : nx*ny*nz floats, column-major (Fortran order), as MATLAB and
 *                numpy with order='F' both lay a volume out
 * connectivity : 6 (shared face), 18 (+ shared edge) or 26 (+ shared corner).
 *                0 means 26, which is what the MATLAB toolbox and fslmaths use.
 *                nilearn's TFCE uses 6, so a binding that wants to reproduce it
 *                asks for 6.
 * calc_neg     : also enhance the negative part of the map, giving it a negative
 *                TFCE value. Asking for it on a map that has no negative values
 *                costs a little time and changes nothing.
 */
int tfce_volume(const float *data, int nx, int ny, int nz, int connectivity,
                double E, double H, int calc_neg, float *out);

/*
 * TFCE of a single map on a mesh.
 *
 * adj_ptr, adj_idx : CSR adjacency of N vertices, 0-based. adj_ptr has N+1
 *                    entries. Build it with tfce_adjacency_from_faces, or hand
 *                    over the indptr/indices of any scipy CSR matrix.
 */
int tfce_mesh(const float *data, int N,
              const int *adj_ptr, const int *adj_idx,
              double E, double H, int calc_neg, float *out);

/*
 * TFCE of B maps at once, one map per thread. This is what a permutation test
 * wants: the permutations are independent of one another, so parallelising
 * across them keeps every thread on its own data and never locks on the hot path.
 *
 * data, out  : N x B floats, column-major, one map per column
 * n_threads  : <= 0 means one per map
 */
int tfce_volume_batch(const float *data, int nx, int ny, int nz, int connectivity,
                      int B, double E, double H, int calc_neg, int n_threads,
                      float *out);

int tfce_mesh_batch(const float *data, int N, int B,
                    const int *adj_ptr, const int *adj_idx,
                    double E, double H, int calc_neg, int n_threads,
                    float *out);

/*
 * CSR adjacency from a triangle list.
 *
 * faces  : nF x 3, column-major, 1-BASED vertex indices (this is what GIFTI and
 *          MATLAB hand over). Every index must name a vertex in [1, N]; one that
 *          does not is refused with TFCE_EFACES rather than written past the end
 *          of the arrays.
 * ptr,idx: allocated here, released with tfce_free
 */
int tfce_adjacency_from_faces(const int *faces, int nF, int N,
                              int **ptr, int **idx);

void tfce_free(void *p);

#ifdef __cplusplus
}
#endif

#endif /* TFCE_CAPI_H */
