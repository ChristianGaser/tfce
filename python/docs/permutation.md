# A permutation test, end to end

The package gives you three pieces - the GLM, the transform, the tail approximations - and leaves you
to assemble them, because the right way to permute depends on your design and nobody can guess that
for you.

This page assembles them for the two commonest designs. **The code here is run as part of the test
suite**, so it works.

---

## The shape of it

1. **Fit the unpermuted model** → observed statistic `t0`.
2. **Permute** many times → a null distribution of the statistic.
3. **TFCE** both.
4. Turn the nulls into p-values:
   - **FWE-corrected** ← the distribution of the *maximum* TFCE over the image (Gamma fit).
   - **Uncorrected** ← each element's *own* distribution (Pareto tail fit).

The two nulls are different objects and answer different questions. The maximum is what controls the
family-wise error; an element's own distribution is what gives it an uncorrected p.

---

## One-sample test (sign-flipping)

A one-sample design is not permuted by shuffling subjects - there is nothing to shuffle them with.
It is permuted by **flipping the sign** of each subject's map, which is exchangeable under the null
of "mean zero".

```python
import numpy as np
import tfce
from tfce.glm import PermutedGLM
from tfce.tails import gamma_pvalue, pareto_pvalue

rng = np.random.default_rng(0)

# Y: (n_elements, n_subjects). One row per voxel, one column per subject.
n_subj, shape = 20, (16, 16, 16)
Y = data.reshape(-1, n_subj)          # your data, elements x subjects
n_elem = Y.shape[0]

X = np.ones((n_subj, 1))              # intercept only
c = np.array([1.0])                   # test the mean

model = PermutedGLM(Y, X, c)

# --- observed ------------------------------------------------------------
t0 = model.fit(None).reshape(shape)
tfce0 = tfce.tfce(t0, E=0.5, H=2.0)

# --- permutations --------------------------------------------------------
n_perm = 1000
n_tail = 100                          # how much of each element's tail to keep

sgn = np.sign(tfce0.ravel())          # both sides become an upper tail
sgn[sgn == 0] = 1
obs = np.abs(tfce0).ravel()

null_max = np.empty(n_perm)                   # for FWE
tail = np.full((n_tail, n_elem), -np.inf)     # for the uncorrected p-values
cnt = np.zeros(n_elem)                        # exceedances, over ALL permutations

for p in range(n_perm):
    flips = rng.choice([-1.0, 1.0], size=n_subj)
    tp = model.fit_signs(flips).reshape(shape)
    tp = tfce.tfce(tp, E=0.5, H=2.0)

    null_max[p] = np.abs(tp).max()

    v = tp.ravel() * sgn                      # an upper tail on both sides

    cnt += v >= obs                           # exact, and cheap

    # keep only the largest n_tail values per element: that is all the fit needs
    tail[0] = v
    tail.sort(axis=0)                         # smallest first; row 0 is the weakest
```

> **Two accumulators, and you need both.** Storing the whole `n_perm × n_elem` null is out of the
> question - at 1000 permutations and 400k voxels that is 1.6 GB - so only the **tail** is kept. But a
> truncated tail cannot tell you *how often* an element was exceeded, and that count is what the
> p-value is made of. So `cnt` is accumulated separately, exactly, over every permutation. It costs
> one comparison per element and is the reason `pareto_pvalue` takes it as an argument rather than
> working it out.
>
> The loop keeps the tail the simple way (overwrite the weakest, re-sort). The MATLAB toolbox tracks a
> running minimum instead, which is cheaper; either is fine at these sizes.

### The p-values

```python
# FWE-corrected: how extreme is each element against the null of the MAXIMUM?
p_fwe = gamma_pvalue(obs, null_max).reshape(shape)

# Uncorrected: how extreme is each element against its OWN null?
p_unc = pareto_pvalue(obs, tail, cnt, n_perm).reshape(shape)

significant = p_fwe < 0.05
```

`gamma_pvalue` fits a Gamma to the maximum distribution from its first three moments, rather than
counting exceedances - so a corrected p-value is **not floored at `1/n_perm`** and converges with far
fewer permutations.

`pareto_pvalue` fits a Generalised Pareto to the tail of each element's own distribution, which does
the same for the uncorrected p-values, and with them for FDR (which is computed from them). Elements
that were exceeded often enough are left to plain counting - the count is already precise there and
the fit has nothing to add.

---

## Two-sample test (permuting group labels)

Here the exchangeable thing is the group label. Permute the *rows of the data*, not the design.

```python
n_per_group = 15
n_subj = 2 * n_per_group

g = np.repeat([1.0, 0.0], n_per_group)
age = rng.standard_normal(n_subj)          # a nuisance covariate

X = np.column_stack([g, 1 - g, age])       # group A, group B, age
c = np.array([1.0, -1.0, 0.0])             # A > B, controlling for age

model = PermutedGLM(Y, X, c)
t0 = model.fit(None).reshape(shape)

for p in range(n_perm):
    perm = rng.permutation(n_subj)
    tp = model.fit(perm).reshape(shape)
    ...                                     # exactly as above
```

`model.fit(perm)` never forms the permuted dataset. See [theory](theory.md#the-permuted-glm) for why
it does not have to.

---

## Exchangeability - the part that will bite you

Everything above rests on one assumption: **that the thing you are permuting is exchangeable under
the null.** Get that wrong and the test does not fail, it just gives you a confidently wrong answer.

The rules are not the library's to enforce, because they depend on your design:

- **Repeated measures / longitudinal**: sessions from the same subject are *not* exchangeable with
  sessions from another. You must permute **within blocks**.
- **Unequal variances between groups**: permutation assumes exchangeability, which unequal variance
  violates.
- **A nuisance variable correlated with the effect of interest**: how you permute matters a great
  deal. The MATLAB toolbox offers three schemes (Draper-Stoneman, Freedman-Lane, Smith); this package
  gives you the GLM and lets you implement the one you need. See
  [Winkler et al. 2014](https://doi.org/10.1016/j.neuroimage.2014.01.060).

The MATLAB toolbox checks itself for this at runtime (it verifies the permutation against the
exchangeability blocks, and warns if the permutation null comes out too narrow). **This package does
not** - it hands you the pieces and trusts you. That is the trade for being a library rather than a
pipeline, and it is the single biggest thing to be careful about.

A cheap sanity check you can run yourself, on the **upper** tail of the p-values:

```python
# The p-value above is ONE-SIDED and conditioned on the sign of the observed
# effect, so under a valid permutation it is uniform on (0, 0.5] -- not on (0, 1].
# Roughly 10% of the elements should therefore exceed 0.45.
#
# A permutation null that is too narrow makes the test anti-conservative, and shows
# up here as far FEWER. True effects only ever produce small p-values and cannot
# inflate the upper tail, which is what makes this a clean check.
frac = np.mean(p_unc > 0.45)
if frac < 0.05:
    print(f"WARNING: only {frac:.1%} of p exceed 0.45, expected ~10%. "
          "The permutation null looks too narrow - check exchangeability.")
```

Watch the bound: **0.5, not 1.0**. Because the p-value is conditioned on the sign of the effect that
was actually observed, it can never come out near 1 - a large positive statistic cannot also be
unusually *small*. Testing it against 0.95 as if it were a two-sided p-value is a mistake that looks
perfectly reasonable and fires on every analysis.

---

## How many permutations?

Fewer than you think, because of the tail fits.

- **FWE**: the Gamma fit converges quickly. **500–1000** is usually plenty.
- **Uncorrected / FDR**: the Pareto fit recovers `p ≈ 1e-4` from **1000** permutations with the median
  right, where plain counting returns zero for 91% of elements and tells you nothing at all.

Without the fits you would need 10,000+ just to get below the `1/n_perm` floor. That floor - not the
statistic - is what makes permutation testing expensive. See [theory](theory.md#fewer-permutations).
