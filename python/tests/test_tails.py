"""
Do the tail approximations do what they are for?

The point of both fits is to resolve a p-value that counting cannot reach, so the
question is never "is the fit exact" -- extrapolating a tail from a hundred points
cannot be -- but "does it beat what it replaces". Both checks below therefore
build a large reference by counting, show the fit only a fraction of it, and score
the fit *and* plain counting against the truth. The fit has to win, or at least
not lose.
"""

import numpy as np
import pytest

from tfce.tails import gamma_pvalue, gpd_fit_pwm, gpd_sf, pareto_pvalue


def gumbel_max(rng, n, n_elem=200, shift=0.0):
    """The maximum of many independent statistics is Gumbel-like."""
    return np.max(-np.log(-np.log(rng.random((n_elem, n)))), axis=0) + shift


class TestGamma:
    def test_matches_counting_in_the_tail(self):
        """Where counting is reliable, the Gamma fit must agree with it."""
        rng = np.random.default_rng(10)
        null = gumbel_max(rng, 20000)

        # thresholds whose true p-value we know from the big sample
        ref = np.sort(null)[::-1]
        for q in (0.05, 0.02, 0.01):
            x = ref[int(q * ref.size)]
            p = gamma_pvalue(np.array([x]), null[:1000])[0]
            # 3 standard errors of the counting estimate the fit is compared to
            se = np.sqrt(q * (1 - q) / 1000)
            assert abs(p - q) < 3 * se + 0.25 * q, (
                f"Gamma p={p:.4f} vs true {q:.4f}"
            )

    def test_resolves_below_the_counting_floor(self):
        """From 500 permutations, counting cannot go below 1/500 = 0.002."""
        rng = np.random.default_rng(11)
        null = gumbel_max(rng, 500)

        x = null.max() + 2.0                # far out in the tail
        p = gamma_pvalue(np.array([x]), null)[0]

        assert 0 < p <= 1 / 500, "the Gamma fit should reach past the floor"
        assert np.isfinite(p)

    def test_never_returns_zero(self):
        rng = np.random.default_rng(12)
        null = gumbel_max(rng, 1000)
        x = np.array([null.max() + 10.0])
        p = gamma_pvalue(x, null)[0]
        assert p > 0


class TestPareto:
    @pytest.fixture(scope="class")
    def null_and_ref(self):
        """A real permutation distribution: a one-sample t under sign-flipping."""
        rng = np.random.default_rng(13)
        n, n_elem, n_ref = 40, 400, 20000
        Y = rng.standard_normal((n_elem, n))

        null = np.empty((n_ref, n_elem))
        for p in range(n_ref):
            s = rng.choice([-1.0, 1.0], size=n)
            Ys = Y * s
            m = Ys.mean(axis=1)
            sd = Ys.std(axis=1, ddof=1)
            null[p] = m / (sd / np.sqrt(n))
        return null

    @pytest.mark.parametrize("q", [3e-3, 1e-3, 3e-4])
    def test_beats_counting_below_the_floor(self, null_and_ref, q):
        null = null_and_ref
        n_ref, n_elem = null.shape
        n_acc = 1000                       # what the fit is allowed to see

        ref = -np.sort(-null, axis=0)
        x = ref[int(q * n_ref)]            # threshold with true p-value q

        seen = null[:n_acc]
        cnt = (seen >= x[None, :]).sum(axis=0)
        p_cnt = cnt / n_acc

        # the fit only ever needs the tail, but the COUNT has to come from every
        # permutation -- a truncated tail cannot say how often an element was
        # exceeded, which is the whole of the p-value
        tail = -np.sort(-seen, axis=0)[:100]
        p_par = pareto_pvalue(x, tail, cnt, n_acc)

        # unbiased in the middle
        med = np.median(p_par / q)
        assert 0.5 < med < 2.0, f"median p_hat/p = {med:.2f}"

        # never zero, which counting is, constantly, down here
        assert np.all(p_par > 0)

        # and at least as often within a factor of two of the truth as counting
        def within2(p):
            r = p / q
            return np.mean((r >= 0.5) & (r <= 2.0))

        assert within2(p_par) >= within2(p_cnt) - 0.05, (
            f"Pareto {within2(p_par):.2f} vs counting {within2(p_cnt):.2f}"
        )

    def test_well_counted_elements_are_left_alone(self, null_and_ref):
        null = null_and_ref
        n_acc = 1000
        seen = null[:n_acc]

        # a threshold every element clears often: p ~ 0.2, so ~200 exceedances
        x = -np.sort(-null, axis=0)[int(0.2 * null.shape[0])]

        cnt = (seen >= x[None, :]).sum(axis=0)
        tail = -np.sort(-seen, axis=0)[:100]

        p = pareto_pvalue(x, tail, cnt, n_acc, n_exc_min=25)

        assert np.array_equal(p, cnt / n_acc), (
            "the fit touched elements counting resolves"
        )


class TestGPD:
    def test_pwm_recovers_known_parameters(self):
        """The moment estimator must find the shape and scale it was given."""
        from scipy import stats

        for k_true, sigma_true in [(0.2, 1.5), (-0.3, 2.0), (0.0, 1.0)]:
            rng = np.random.default_rng(20)
            # scipy's c is the negative of the Hosking-Wallis k
            y = stats.genpareto.rvs(
                c=-k_true, scale=sigma_true, size=(4000, 200), random_state=rng
            )
            y = np.sort(y, axis=0)

            k, sigma = gpd_fit_pwm(y)
            assert abs(np.median(k) - k_true) < 0.05, f"k {np.median(k)} vs {k_true}"
            assert abs(np.median(sigma) - sigma_true) < 0.1

    def test_sf_is_a_survival_function(self):
        z = np.linspace(0, 5, 50)
        for k in (-0.3, 0.0, 0.3):
            s = gpd_sf(z, np.full(z.shape, k), np.ones(z.shape))
            assert s[0] == pytest.approx(1.0)
            assert np.all(np.diff(s) <= 1e-12), "survival must be non-increasing"
            assert np.all(s >= 0) and np.all(s <= 1)

    def test_positive_shape_has_a_finite_end_point(self):
        """k > 0 bounds the distribution above -- the trap the fit has to survive."""
        k, sigma = 0.5, 1.0
        end = sigma / k
        s = gpd_sf(np.array([end + 1.0]), np.array([k]), np.array([sigma]))
        assert s[0] == 0.0, "beyond the end point the GPD really does say zero"
