"""
Is it really a drop-in for nilearn?

nilearn computes TFCE by stepping the height (``calculate_tfce`` in
``nilearn.mass_univariate._utils``). The function here has the same signature and
returns the same thing, computed exactly. The checks are therefore about the
*contract*, not about agreeing numerically with nilearn -- the whole point is that
the numbers are better.
"""

import numpy as np
import pytest
from scipy.ndimage import gaussian_filter, generate_binary_structure

from tfce.nilearn_compat import calculate_tfce, connectivity_from_bin_struct

nilearn = pytest.importorskip("nilearn")


class TestConnectivity:
    @pytest.mark.parametrize("rank,expected", [(1, 6), (2, 18), (3, 26)])
    def test_read_out_of_bin_struct(self, rank, expected):
        bs = generate_binary_structure(3, rank)
        assert connectivity_from_bin_struct(bs) == expected

    def test_nilearn_uses_six(self):
        """Pinning the assumption: nilearn builds its structure with rank 1."""
        bs = generate_binary_structure(3, 1)
        assert connectivity_from_bin_struct(bs) == 6

    def test_rejects_nonsense(self):
        with pytest.raises(ValueError):
            connectivity_from_bin_struct(np.ones((3, 3)))


class TestDropIn:
    def test_signature_and_shape(self):
        rng = np.random.default_rng(0)
        arr4d = np.stack(
            [gaussian_filter(rng.standard_normal((10, 11, 9)), 1.3) for _ in range(3)],
            axis=-1,
        )
        bs = generate_binary_structure(3, 1)

        out = calculate_tfce(arr4d, bs, E=0.5, H=2, dh="auto", two_sided_test=True)

        assert out.shape == arr4d.shape

    def test_dh_is_accepted_and_irrelevant(self):
        """There is no step size. Passing one must not be an error, and must not
        change anything -- that is the entire selling point."""
        rng = np.random.default_rng(1)
        arr4d = gaussian_filter(rng.standard_normal((8, 8, 8)), 1.2)[..., None]
        bs = generate_binary_structure(3, 1)

        a = calculate_tfce(arr4d, bs, dh="auto")
        b = calculate_tfce(arr4d, bs, dh=0.1)
        c = calculate_tfce(arr4d, bs, dh=0.001)

        assert np.array_equal(a, b)
        assert np.array_equal(a, c)

    def test_one_sided(self):
        rng = np.random.default_rng(2)
        arr4d = gaussian_filter(rng.standard_normal((8, 8, 8)), 1.2)[..., None]
        bs = generate_binary_structure(3, 1)

        two = calculate_tfce(arr4d, bs, two_sided_test=True)
        one = calculate_tfce(arr4d, bs, two_sided_test=False)

        assert np.any(two < 0)
        assert np.all(one >= 0)

    def test_monkeypatching_nilearn_works(self, monkeypatch):
        """The actual use: swap our exact transform into nilearn and run its GLM."""
        from nilearn.mass_univariate import _utils

        monkeypatch.setattr(_utils, "calculate_tfce", calculate_tfce)

        rng = np.random.default_rng(3)
        arr4d = gaussian_filter(rng.standard_normal((6, 6, 6)), 1.0)[..., None]
        bs = generate_binary_structure(3, 1)

        out = _utils.calculate_tfce(arr4d, bs, E=0.5, H=2)
        assert out.shape == arr4d.shape
        assert np.isfinite(out).all()
