"""
Exact threshold-free cluster enhancement.

The TFCE of an element is the integral

    TFCE(v) = int_0^{t_v} e_v(h)^E h^H dh

over the extent e_v(h) of the cluster containing v at height h. Implementations
normally approximate it by stepping h over a grid and summing, which costs a step
size dh, a discretisation error that depends on it, and an accuracy parameter that
the caller has to guess. This one does not: it builds the max-tree (the component
tree) of the map with union-find, and because e_v(h) is piecewise constant it
integrates each piece in closed form. The answer is the integral, not a sample of
it, and there is nothing to tune.

Everything here is numpy in, numpy out. There is no file I/O, no image object and
no framework: the functions take arrays and a neighbourhood, which is all TFCE
needs and all that a caller -- nilearn, fitlins, a bare script -- should have to
supply.
"""

import numpy as np

from . import _maxtree

__all__ = [
    "tfce",
    "adjacency_from_faces",
    "TfceError",
]

TfceError = _maxtree.TfceError

#: Neighbourhoods a volume can be enhanced under. 6 shares a face, 18 also shares
#: an edge, 26 also shares a corner. The MATLAB toolbox and fslmaths use 26;
#: nilearn uses 6.
VALID_CONNECTIVITY = (6, 18, 26)


def _as_csr(adjacency, n_elements):
    """Accept a scipy sparse matrix, or an (indptr, indices) pair, either way."""
    if isinstance(adjacency, tuple):
        indptr, indices = adjacency
    else:  # anything scipy-sparse-like, without importing scipy for it
        csr = adjacency.tocsr() if hasattr(adjacency, "tocsr") else adjacency
        indptr, indices = csr.indptr, csr.indices

    indptr = np.ascontiguousarray(indptr, dtype=np.int32)
    indices = np.ascontiguousarray(indices, dtype=np.int32)

    if indptr.size != n_elements + 1:
        raise ValueError(
            "adjacency describes %d elements, but the data has %d"
            % (indptr.size - 1, n_elements)
        )
    return indptr, indices


def tfce(
    data,
    *,
    adjacency=None,
    connectivity=26,
    E=0.5,
    H=2.0,
    two_sided=True,
    n_jobs=1,
):
    """Threshold-free cluster enhancement of one map or of many.

    Parameters
    ----------
    data : ndarray
        One of

        * ``(nx, ny, nz)``      -- a single volume
        * ``(nx, ny, nz, B)``   -- B volumes, e.g. B permutations
        * ``(n_vertices,)``     -- a single surface map, needs ``adjacency``
        * ``(n_vertices, B)``   -- B surface maps, needs ``adjacency``

        Any dtype; it is used as float32, which is all a permuted statistic is
        ever worth.

    adjacency : scipy.sparse matrix or (indptr, indices), optional
        Element adjacency, for surface data. Omit it for volumes, where the
        neighbourhood is the grid itself. Build one from a triangle list with
        :func:`adjacency_from_faces`.

    connectivity : {6, 18, 26}, default 26
        Volume neighbourhood. 26 is what the MATLAB TFCE toolbox and fslmaths
        use; nilearn uses 6. Ignored for surface data, where the mesh says what
        the neighbours are.

    E, H : float
        Weight of extent and of height. ``E=0.5, H=2`` are the values Smith and
        Nichols established for volumes; surfaces and TBSS use ``E=1, H=2``.

    two_sided : bool, default True
        Enhance the negative part of the map as well, giving it a negative TFCE
        value. A map with no negative values is unaffected by this.

    n_jobs : int, default 1
        Threads to spread B maps over. ``-1`` means one per map. Only the batched
        forms can use this -- the permutations are independent of one another, so
        that is where the parallelism belongs. The GIL is released for the whole
        call.

    Returns
    -------
    ndarray
        TFCE values, float32, in the shape the data came in.
    """
    if connectivity not in VALID_CONNECTIVITY:
        raise ValueError(
            "connectivity must be one of %r, got %r"
            % (VALID_CONNECTIVITY, connectivity)
        )

    arr = np.asarray(data)
    n_threads = 0 if n_jobs is None or n_jobs < 0 else int(n_jobs)

    is_surface = adjacency is not None

    if is_surface:
        if arr.ndim not in (1, 2):
            raise ValueError(
                "surface data must be (n_vertices,) or (n_vertices, B), got %r"
                % (arr.shape,)
            )
        squeeze = arr.ndim == 1
        flat = arr.reshape(arr.shape[0], -1)

        indptr, indices = _as_csr(adjacency, flat.shape[0])
        out = _maxtree.mesh_batch(
            flat, indptr, indices, E, H, bool(two_sided), n_threads
        )
    else:
        if arr.ndim not in (3, 4):
            raise ValueError(
                "volume data must be (nx, ny, nz) or (nx, ny, nz, B), got %r"
                % (arr.shape,)
            )
        squeeze = arr.ndim == 3
        nx, ny, nz = arr.shape[:3]

        # the core indexes a volume column-major, the way MATLAB and Fortran do
        flat = np.reshape(arr, (nx * ny * nz, -1), order="F")

        out = _maxtree.volume_batch(
            flat, nx, ny, nz, connectivity, E, H, bool(two_sided), n_threads
        )
        out = np.reshape(out, (nx, ny, nz, -1), order="F")

    if squeeze:
        out = out[..., 0]
    return out


def adjacency_from_faces(faces, n_vertices):
    """Element adjacency of a triangulated surface.

    Parameters
    ----------
    faces : (n_faces, 3) array of int
        Triangles, as **1-based** vertex indices -- which is how GIFTI stores
        them, and therefore what ``nibabel`` and ``SPM.xVol.G.faces`` hand over.
        An index that names no vertex is refused rather than read past the end of
        the arrays.
    n_vertices : int

    Returns
    -------
    (indptr, indices) : int32 arrays
        CSR, 0-based, ready for ``scipy.sparse.csr_matrix`` or to pass straight
        back into :func:`tfce` as ``adjacency``.
    """
    faces = np.asarray(faces)
    if faces.ndim != 2 or faces.shape[1] != 3:
        raise ValueError("faces must be (n_faces, 3), got %r" % (faces.shape,))
    return _maxtree.adjacency_from_faces(faces, int(n_vertices))
