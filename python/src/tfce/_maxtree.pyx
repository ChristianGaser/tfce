# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
#
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2020-2026, Christian Gaser. See LICENSE.
"""
Cython binding to the exact TFCE max-tree.

This is a thin layer and is meant to stay thin: it hands numpy arrays to the C
API in tfce_capi.h and gets numpy arrays back. Everything that is actually TFCE
lives in the C, and is the very same C the MATLAB toolbox runs -- not a
reimplementation of it that could drift.

The GIL is released around every call, so the batched transform really does use
its threads, and a caller running several of these in a thread pool is not
serialised on the interpreter.
"""

import numpy as np
cimport numpy as cnp

cnp.import_array()


cdef extern from "tfce_capi.h":
    int TFCE_OK
    int TFCE_ENOMEM
    int TFCE_EFACES

    int tfce_volume(const float *data, int nx, int ny, int nz, int connectivity,
                    double E, double H, int calc_neg, float *out) nogil

    int tfce_mesh(const float *data, int N,
                  const int *adj_ptr, const int *adj_idx,
                  double E, double H, int calc_neg, float *out) nogil

    int tfce_volume_batch(const float *data, int nx, int ny, int nz,
                          int connectivity, int B,
                          double E, double H, int calc_neg, int n_threads,
                          float *out) nogil

    int tfce_mesh_batch(const float *data, int N, int B,
                        const int *adj_ptr, const int *adj_idx,
                        double E, double H, int calc_neg, int n_threads,
                        float *out) nogil

    int tfce_adjacency_from_faces(const int *faces, int nF, int N,
                                  int **ptr, int **idx) nogil

    void tfce_free(void *p) nogil


class TfceError(RuntimeError):
    """The C core refused the call."""


cdef _check(int rc):
    if rc == TFCE_OK:
        return
    if rc == TFCE_ENOMEM:
        raise MemoryError("TFCE: out of memory")
    if rc == TFCE_EFACES:
        raise TfceError("faces name a vertex that is not in the data")
    raise TfceError("TFCE: error %d" % rc)


def volume_batch(cnp.ndarray data, int nx, int ny, int nz, int connectivity,
                 double E, double H, bint calc_neg, int n_threads):
    """TFCE of B volume maps.

    data : (nx*ny*nz, B) float32, Fortran-ordered (one map per column)
    """
    cdef cnp.ndarray[float, ndim=2, mode="fortran"] d = \
        np.asfortranarray(data, dtype=np.float32)
    cdef int N = d.shape[0]
    cdef int B = d.shape[1]
    cdef cnp.ndarray[float, ndim=2, mode="fortran"] out = \
        np.zeros((N, B), dtype=np.float32, order="F")
    cdef int rc

    if N != nx * ny * nz:
        raise ValueError("data has %d rows, but nx*ny*nz = %d"
                         % (N, nx * ny * nz))

    with nogil:
        rc = tfce_volume_batch(&d[0, 0], nx, ny, nz, connectivity, B,
                               E, H, calc_neg, n_threads, &out[0, 0])
    _check(rc)
    return out


def mesh_batch(cnp.ndarray data,
               cnp.ndarray adj_ptr, cnp.ndarray adj_idx,
               double E, double H, bint calc_neg, int n_threads):
    """TFCE of B surface maps.

    data    : (n_vertices, B) float32, Fortran-ordered (one map per column)
    adj_ptr : (n_vertices + 1,) int32   CSR indptr, 0-based
    adj_idx : (nnz,) int32              CSR indices, 0-based
    """
    cdef cnp.ndarray[float, ndim=2, mode="fortran"] d = \
        np.asfortranarray(data, dtype=np.float32)
    cdef cnp.ndarray[int, ndim=1, mode="c"] p = \
        np.ascontiguousarray(adj_ptr, dtype=np.int32)
    cdef cnp.ndarray[int, ndim=1, mode="c"] ix = \
        np.ascontiguousarray(adj_idx, dtype=np.int32)
    cdef int N = d.shape[0]
    cdef int B = d.shape[1]
    cdef cnp.ndarray[float, ndim=2, mode="fortran"] out = \
        np.zeros((N, B), dtype=np.float32, order="F")
    cdef int rc

    if p.shape[0] != N + 1:
        raise ValueError("adj_ptr has %d entries, expected n_vertices + 1 = %d"
                         % (p.shape[0], N + 1))

    with nogil:
        rc = tfce_mesh_batch(&d[0, 0], N, B, &p[0], &ix[0],
                             E, H, calc_neg, n_threads, &out[0, 0])
    _check(rc)
    return out


def adjacency_from_faces(cnp.ndarray faces, int n_vertices):
    """CSR adjacency from a triangle list.

    faces : (n_faces, 3) integer array of 1-BASED vertex indices, as GIFTI and
            MATLAB store them. Any index outside [1, n_vertices] is refused.

    Returns (indptr, indices), both int32 and 0-based, ready to hand to
    scipy.sparse.csr_matrix.
    """
    # column-major, which is what the C expects of an nF x 3 list
    cdef cnp.ndarray[int, ndim=2, mode="fortran"] f = \
        np.asfortranarray(faces, dtype=np.int32)
    cdef int nF = f.shape[0]
    cdef int *ptr = NULL
    cdef int *idx = NULL
    cdef int rc
    cdef int nnz

    if f.shape[1] != 3:
        raise ValueError("faces must be n_faces x 3")

    with nogil:
        rc = tfce_adjacency_from_faces(&f[0, 0], nF, n_vertices, &ptr, &idx)
    _check(rc)

    try:
        indptr = np.empty(n_vertices + 1, dtype=np.int32)
        for i in range(n_vertices + 1):
            indptr[i] = ptr[i]
        nnz = ptr[n_vertices]
        indices = np.empty(nnz, dtype=np.int32)
        for i in range(nnz):
            indices[i] = idx[i]
    finally:
        tfce_free(ptr)
        tfce_free(idx)

    return indptr, indices
