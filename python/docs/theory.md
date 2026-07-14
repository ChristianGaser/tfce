# Theory

Three ideas, each of which removes something you would otherwise have to guess at.

---

## The exact integral

The TFCE of an element `v` is

```text
TFCE(v) = ∫₀^{t_v} e_v(h)^E · h^H dh
```

where `e_v(h)` is the **extent** of the connected component containing `v` when the map is
thresholded at height `h`. It combines focal effects (large `t_v`, so a long integral) with broad
ones (large `e_v`), and — the point of the whole method — needs no cluster-forming threshold.

### How it is normally computed, and why that costs you

Step `h` over a grid, label the connected components at each step, and sum:

```python
for h in np.linspace(0, t.max(), 100):
    labels, _ = label(t >= h)
    sizes = np.bincount(labels.ravel())
    out += sizes[labels]**E * h**H * dh
```

This is a Riemann sum. It has a discretisation error of order `dh`, and `dh` is a parameter you have
to choose without knowing what it costs you. At 100 steps the error is about **1.7%**; halving `dh`
halves it, forever, and it never reaches zero.

### How it is computed here

Notice that `e_v(h)` is **piecewise constant**. As `h` falls, the component containing `v` only ever
changes at the finitely many heights where two components merge. Between merges the integrand is
`s^E · h^H` with `s` fixed — and *that* integrates in closed form:

```text
s^E · (birth^(H+1) − death^(H+1)) / (H+1)
```

So the whole integral is a finite sum over the components that ever contain `v`, and there is nothing
to approximate.

Finding those components is the **max-tree** (or component tree): sweep the elements in descending
order of value, union-find them with their already-active neighbours, and every time a component
changes, record a node with its size, its birth height and its death height. Each element's TFCE is
then the sum of the contributions of its ancestors in that tree — one root-to-leaf accumulation.

Cost is `O(N log N)` for the sort plus `O(N α(N))` for the sweep, and it is **independent of any
precision parameter**, because there isn't one.

> The sort is the expensive part, and every value it sorts is a positive float. For positive IEEE-754
> floats the bit pattern read as an integer is monotonic, so the sort key needs no transformation and
> four counting passes replace `n log n` comparisons. That is where most of the speed comes from.

### Two-sided

The negative part of the map is enhanced by the same procedure applied to `-t`, and given a negative
TFCE value. The two passes touch disjoint elements — an element cannot be both positive and negative —
so a map with no negative values is completely unaffected by asking for them.

---

## The permuted GLM

A permutation test refits the same model thousands of times. The obvious way to do it is to build the
permuted dataset and regress:

```python
beta = pinv(X) @ Y[:, perm].T
```

That costs an `n_elements × n × n` multiply per permutation and is by far the dominant term. It is
also unnecessary.

**First**, for any design, `R = I − X·pinv(X)` is a projector, so

```text
ResSS = ‖y‖² − βᵀ (XᵀX) β
```

and the residual matrix never has to exist.

**Second**, `‖y‖²` does not depend on the permutation at all — a permutation matrix only reorders the
observations — so it is computed once.

**Third**, the permutation can be moved off the data and onto the design. Writing out
`β = pinv(X) @ Y[:, perm].T` and substituting `k = perm[j]` shows that permuting the columns of `Y` is
the same as permuting the columns of `pinv(X)` by the **inverse** permutation. Since `pinv(X)` is
`rank(X) × n` — tiny — that costs nothing.

So the whole permutation collapses into one small projection, and the permuted dataset is never
formed. That is what [`PermutedGLM`](api.md#tfceglm) does.

> The *inverse* is the whole trick, and getting it the wrong way round is silent: a wrong permutation
> is still a permutation, and the statistic still looks perfectly reasonable. (I got it wrong first
> time; the test suite caught it.)

Under **Freedman-Lane** there is a fourth simplification: because the nuisance columns are part of the
full design, the full-model residual-former already annihilates them, so the nuisance projection drops
out of the residuals entirely.

---

## Fewer permutations

Counting exceedances cannot report a p-value below `1/n_perm`. If you run 1000 permutations, the
smallest p-value in the universe is 0.001 — and for 91% of the elements that actually matter, the
count is **zero**, which tells you nothing at all.

**That floor is what makes permutation testing expensive.** Not the statistic, not the number of
elements. And it caps FDR too, since FDR is computed from the uncorrected p-values.

Two fits remove it.

### Gamma, for the FWE-corrected p-values

The family-wise error is controlled by the distribution of the **maximum** statistic over the image.
That distribution is right-skewed and well described by a Gamma, so rather than counting exceedances,
[`gamma_pvalue`](api.md#gamma_pvalue) fits one from its first three moments:

```text
mean = μ + kθ      var = kθ²      skew = 2/√k
```

Three moments, three parameters, closed form. The corrected p-value is then read off the fitted
survival function and is not floored at anything.

### Generalised Pareto, for the uncorrected ones

Each element has its own permutation distribution, and what you need is its far tail.
[`pareto_pvalue`](api.md#pareto_pvalue) fits a Generalised Pareto to it, by probability-weighted
moments (closed form again, so every element is fitted at once).

Two details make it work.

**The shape is pooled across elements.** Every element carries the *same* statistic under the *same*
design, so they differ in **scale**, not in **shape**. Fitting the shape separately for each element
spends a hundred tail values on a number they all share, and the error in it dominates the error of
the extrapolation. Pooling it lifts the fraction of elements landing within a factor of two of the
truth at `p = 1e-4` from 45% to **60%**.

**A positive shape has a finite upper endpoint.** Probability-weighted moments put the shape `k` above
zero for about half the elements *by chance*, which gives the fitted distribution a hard maximum at
`u + σ/k`. When the observed statistic lies beyond it, the fit returns exactly **zero** — and a
p-value of zero is unusable (it is `Inf` on a `-log10` map). But the observation is itself proof that
the endpoint is wrong, so such a fit is rejected rather than believed, and the exponential limit
(`k = 0`, infinite support) answers instead.

### What it buys

From **1000** permutations, against a reference counted over 50,000:

| true p | Pareto (within a factor of 2) | counting (within a factor of 2) |
| --- | --- | --- |
| 3e-3 | 88% | 78% |
| 1e-3 | 73% | 57% |
| 3e-4 | **56%** | **0%** |
| 1e-4 | **42%** | **0%** |

Median-unbiased throughout, and **never zero** — where counting returns zero for 91% of elements and
says nothing.

Be honest about the size of this, though. It is not "100× fewer permutations". Extrapolating a tail
from a hundred points has real variance. What it buys is roughly a **5× cut** for the same accuracy
around `p ≈ 1e-3`, and it makes `p ≈ 1e-4` reachable *at all*.

---

## References

- Smith & Nichols (2009), *Threshold-free cluster enhancement*, NeuroImage 44:83–98.
  [doi](https://doi.org/10.1016/j.neuroimage.2008.03.061)
- Winkler et al. (2014), *Permutation inference for the general linear model*, NeuroImage 92:381–397.
  [doi](https://doi.org/10.1016/j.neuroimage.2014.01.060)
- Winkler et al. (2016), *Faster permutation inference in brain imaging*, NeuroImage 141:502–516.
  [doi](https://doi.org/10.1016/j.neuroimage.2016.05.068)
- Hosking & Wallis (1987), *Parameter and quantile estimation for the generalized Pareto
  distribution*, Technometrics 29:339–349.
