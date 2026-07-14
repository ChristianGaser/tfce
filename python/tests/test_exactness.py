"""
Is the max-tree really the exact TFCE integral?

This is the claim the whole transform rests on, so it is checked against an
*independent* implementation: the obvious one, which steps the height over a grid
and labels connected components with ``scipy.ndimage.label`` at each step. That
reference has a discretisation error of order dh, so it must converge onto the
max-tree at first order as the step size shrinks. If the max-tree were merely
another approximation, the two would converge onto each other at some other rate,
or onto nothing.

The criterion is the slope of log(error) against log(n_steps), which must be about
-1. Testing the ratio of successive errors instead is far too noisy, because the
error is a maximum over elements.
"""

import numpy as np
import pytest
from scipy.ndimage import gaussian_filter, generate_binary_structure, label

from tfce import tfce


def tfce_reference(arr3d, E, H, n_steps, bin_struct):
    """TFCE by stepping the threshold. Deliberately naive; deliberately not ours."""
    out = np.zeros_like(arr3d, dtype=float)

    for sign in (1, -1):
        a = arr3d * sign
        top = a.max()
        if top <= 0:
            continue

        dh = top / n_steps
        # midpoints, so the sum is a midpoint rule and converges at first order
        for i in range(n_steps):
            h = (i + 0.5) * dh
            lab, n = label(a >= h, bin_struct)
            if n == 0:
                continue
            sizes = np.bincount(lab.ravel())
            extent = sizes[lab]
            contrib = (extent.astype(float) ** E) * (h**H) * dh
            out += sign * np.where(lab > 0, contrib, 0.0)

    return out


@pytest.mark.parametrize("connectivity", [6, 26])
def test_volume_stepping_converges_onto_maxtree(connectivity):
    rng = np.random.default_rng(1)
    arr = gaussian_filter(rng.standard_normal((16, 18, 14)), 1.6)
    arr = arr / arr.std()

    E, H = 0.5, 2.0
    conn_rank = {6: 1, 18: 2, 26: 3}[connectivity]
    bin_struct = generate_binary_structure(3, conn_rank)

    exact = tfce(arr, connectivity=connectivity, E=E, H=H, two_sided=True)

    n_steps = np.array([50, 100, 200, 400])
    errs = np.array(
        [
            np.abs(
                tfce_reference(arr, E, H, int(n), bin_struct) - exact
            ).max()
            / np.abs(exact).max()
            for n in n_steps
        ]
    )

    slope = np.polyfit(np.log(n_steps), np.log(errs), 1)[0]

    # first order: halving dh halves the error
    assert -1.25 < slope < -0.75, (
        "dh-stepping does not converge onto the max-tree at first order "
        f"(slope {slope:.2f}, errors {errs})"
    )


def test_batch_equals_single_map():
    """A block of maps must give exactly what the maps give one at a time."""
    rng = np.random.default_rng(2)
    vols = np.stack(
        [gaussian_filter(rng.standard_normal((12, 13, 11)), 1.4) for _ in range(5)],
        axis=-1,
    )
    # half the maps have no negative values at all, which is the case that the
    # shared two_sided flag has to leave alone
    vols[..., 1] = np.abs(vols[..., 1])
    vols[..., 3] = np.abs(vols[..., 3])

    batched = tfce(vols, n_jobs=-1)
    one_by_one = np.stack(
        [tfce(vols[..., b]) for b in range(vols.shape[-1])], axis=-1
    )

    assert np.array_equal(batched, one_by_one)


def test_threads_do_not_change_the_answer():
    rng = np.random.default_rng(3)
    vols = np.stack(
        [gaussian_filter(rng.standard_normal((10, 11, 9)), 1.2) for _ in range(8)],
        axis=-1,
    )
    a = tfce(vols, n_jobs=1)
    b = tfce(vols, n_jobs=8)
    assert np.array_equal(a, b)


def test_connectivity_changes_the_answer():
    """A sanity check on the sanity check: 6 and 26 must not agree by accident."""
    rng = np.random.default_rng(4)
    arr = gaussian_filter(rng.standard_normal((12, 12, 12)), 1.2)

    t6 = tfce(arr, connectivity=6)
    t26 = tfce(arr, connectivity=26)

    assert not np.allclose(t6, t26)
    # a looser neighbourhood merges clusters, so extents -- and TFCE -- grow
    assert np.abs(t26).max() > np.abs(t6).max()
