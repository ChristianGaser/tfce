"""The array API, the surface adjacency, and the things that must be refused."""

import numpy as np
import pytest

from tfce import TfceError, adjacency_from_faces, tfce


def tetrahedron():
    """The smallest closed mesh there is: 4 vertices, 4 faces, 1-based."""
    faces = np.array([[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]], dtype=np.int32)
    return faces, 4


class TestShapes:
    def test_volume_3d(self):
        x = np.zeros((6, 7, 5))
        x[3, 3, 3] = 2.0
        out = tfce(x)
        assert out.shape == (6, 7, 5)
        assert out.dtype == np.float32
        assert out[3, 3, 3] > 0

    def test_volume_4d(self):
        x = np.zeros((6, 7, 5, 3))
        x[3, 3, 3, :] = 2.0
        out = tfce(x)
        assert out.shape == (6, 7, 5, 3)

    def test_surface_1d(self):
        faces, nv = tetrahedron()
        adj = adjacency_from_faces(faces, nv)
        out = tfce(np.array([3.0, 1.0, 0.0, 0.0]), adjacency=adj, E=1.0)
        assert out.shape == (4,)

    def test_surface_2d(self):
        faces, nv = tetrahedron()
        adj = adjacency_from_faces(faces, nv)
        data = np.tile(np.array([3.0, 1.0, 0.0, 0.0])[:, None], (1, 5))
        out = tfce(data, adjacency=adj, E=1.0)
        assert out.shape == (4, 5)
        # every column is the same map, so every column must give the same answer
        assert np.array_equal(out[:, 0], out[:, 4])


class TestSemantics:
    def test_axes_are_not_transposed(self):
        """The core indexes a volume column-major; a mistake here is silent."""
        x = np.zeros((4, 6, 8))
        x[0, 5, 7] = 5.0                     # a corner that is unambiguous
        out = tfce(x, two_sided=False)
        assert out[0, 5, 7] > 0
        assert out.sum() == pytest.approx(out[0, 5, 7]), (
            "the enhanced voxel is not where it was put"
        )

    def test_two_sided(self):
        x = np.zeros((5, 5, 5))
        x[2, 2, 2] = -3.0

        assert tfce(x, two_sided=True)[2, 2, 2] < 0
        assert tfce(x, two_sided=False)[2, 2, 2] == 0

    def test_positive_map_is_unaffected_by_two_sided(self):
        rng = np.random.default_rng(0)
        x = np.abs(rng.standard_normal((6, 6, 6)))
        assert np.array_equal(tfce(x, two_sided=True), tfce(x, two_sided=False))

    def test_all_zero(self):
        out = tfce(np.zeros((5, 5, 5)))
        assert np.all(out == 0)


class TestAdjacency:
    def test_round_trip(self):
        faces, nv = tetrahedron()
        indptr, indices = adjacency_from_faces(faces, nv)

        assert indptr.shape == (nv + 1,)
        assert indptr[0] == 0
        # every vertex of a tetrahedron touches the other three
        assert np.array_equal(np.diff(indptr), [3, 3, 3, 3])
        for v in range(nv):
            nb = set(indices[indptr[v] : indptr[v + 1]])
            assert nb == set(range(nv)) - {v}

    def test_accepts_scipy_sparse(self):
        from scipy.sparse import csr_matrix

        faces, nv = tetrahedron()
        indptr, indices = adjacency_from_faces(faces, nv)
        m = csr_matrix(
            (np.ones(indices.size), indices, indptr), shape=(nv, nv)
        )

        data = np.array([3.0, 1.0, 0.0, 0.0])
        assert np.array_equal(
            tfce(data, adjacency=m, E=1.0),
            tfce(data, adjacency=(indptr, indices), E=1.0),
        )

    def test_out_of_range_vertex_is_refused(self):
        """Not a formality: this used to be written past the end of the arrays."""
        faces, nv = tetrahedron()
        bad = faces.copy()
        bad[0, 0] = 99
        with pytest.raises(TfceError, match="not in the data"):
            adjacency_from_faces(bad, nv)

    def test_zero_based_faces_are_refused(self):
        """GIFTI is 1-based. Handing 0-based faces over is a mistake, not a mode."""
        faces, nv = tetrahedron()
        with pytest.raises(TfceError):
            adjacency_from_faces(faces - 1, nv)


class TestValidation:
    def test_bad_connectivity(self):
        with pytest.raises(ValueError, match="connectivity"):
            tfce(np.zeros((4, 4, 4)), connectivity=7)

    def test_bad_ndim(self):
        with pytest.raises(ValueError, match="volume data"):
            tfce(np.zeros((4, 4)))

    def test_adjacency_size_mismatch(self):
        faces, nv = tetrahedron()
        adj = adjacency_from_faces(faces, nv)
        with pytest.raises(ValueError, match="adjacency describes"):
            tfce(np.zeros(7), adjacency=adj)
