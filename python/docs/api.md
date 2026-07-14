# API reference

Everything public, with what it takes and what it gives back.

```python
import tfce
from tfce.glm import PermutedGLM, partition_design
from tfce.tails import gamma_pvalue, pareto_pvalue, gpd_fit_pwm, gpd_sf
from tfce.nilearn_compat import calculate_tfce, connectivity_from_bin_struct
```

---

## `tfce`

### `tfce(data, *, adjacency=None, connectivity=26, E=0.5, H=2.0, two_sided=True, n_jobs=1)`

Threshold-free cluster enhancement of one map or of many.

| argument | | |
| --- | --- | --- |
| `data` | ndarray | `(nx, ny, nz)` or `(nx, ny, nz, B)` for volumes; `(n_vertices,)` or `(n_vertices, B)` for surfaces. Any dtype — used as `float32`. |
| `adjacency` | sparse or `(indptr, indices)` | The mesh, for surface data. Omit for volumes. |
| `connectivity` | `6`, `18`, `26` | Volume neighbourhood. Ignored for surfaces. |
| `E`, `H` | float | Extent and height weights. `0.5, 2` for volumes; `1, 2` for surfaces and TBSS. |
| `two_sided` | bool | Also enhance the negative part, giving it a negative TFCE value. |
| `n_jobs` | int | Threads across the `B` maps. `-1` = one per map. Releases the GIL. |

**Returns** `float32` array in the shape the data came in.

**Raises** `ValueError` on a bad shape or connectivity; `TfceError` if the C core refuses the call.

### `adjacency_from_faces(faces, n_vertices)`

CSR adjacency of a triangulated surface.

| argument | | |
| --- | --- | --- |
| `faces` | `(n_faces, 3)` int | **1-based** vertex indices, as GIFTI stores them. |
| `n_vertices` | int | |

**Returns** `(indptr, indices)`, both `int32` and 0-based — ready for `scipy.sparse.csr_matrix`, or to
pass straight back to `tfce()` as `adjacency`.

**Raises** `TfceError` if any index names a vertex that does not exist. This is a real check, not a
formality: an out-of-range index used to be written past the end of the adjacency arrays.

### `TfceError`

Raised when the C core refuses a call. Subclasses `RuntimeError`.

---

## `tfce.tails`

### `gamma_pvalue(stat, null_max, floor=None)`

FWE-corrected p-values, from a Gamma fitted to the distribution of the **maximum**.

| argument | | |
| --- | --- | --- |
| `stat` | ndarray | Observed statistic, one value per element. |
| `null_max` | `(n_perm,)` | The maximum of the statistic over the image, one per permutation. |
| `floor` | float | Smallest p-value to report. Defaults to `1/n_perm`. |

**Returns** p-values, same shape as `stat`.

The Gamma is fitted by matching the first three moments, which inverts in closed form. If the sample
skewness comes out non-positive there is nothing for a Gamma to hold on to, and the empirical count is
returned instead.

### `pareto_pvalue(stat, tail, cnt, n_perm, *, n_exc_min=25, n_shape=20000, ad_max=5.0)`

Uncorrected p-values, resolved **below** the `1/n_perm` floor of counting.

| argument | | |
| --- | --- | --- |
| `stat` | `(n_elements,)` | Observed statistic, **already an upper tail** — flip the sign of the negative elements first. |
| `tail` | `(K, n_elements)` | The largest `K` permuted values of each element. `K ≈ 100` is plenty. |
| `cnt` | `(n_elements,)` | How often the permutations reached or exceeded `stat`, counted over **all** `n_perm` — not just the ones in `tail`. |
| `n_perm` | int | How many permutations were drawn. |
| `n_exc_min` | int | Elements exceeded at least this often keep their count; the fit has nothing to add. |
| `n_shape` | int | How many elements the shape is pooled over. |
| `ad_max` | float | A tail no candidate fit describes this well is counted, not extrapolated. |

**Returns** p-values, `(n_elements,)`.

> **Why `cnt` and `n_perm` are arguments.** Storing the whole `n_perm × n_elements` null is impossible
> at any real size, so only the tail is kept — but a truncated tail cannot say *how often* an element
> was exceeded, and that count is what the p-value is made of. Accumulate it in your permutation loop:
> it costs one comparison per element. See [permutation.md](permutation.md).

### `gpd_fit_pwm(y)`

Generalised Pareto parameters by probability-weighted moments (Hosking & Wallis).

`y` is `(m, n)` — exceedances above a threshold, **ascending** down each column, one column per
element. Returns `(k, sigma)`, each `(n,)`, parameterised as `F(y) = 1 − (1 − k·y/σ)^(1/k)`.

### `gpd_sf(z, k, sigma)`

Survival function of that distribution. Handles the `k → 0` limit (exponential) and clamps to zero
beyond the upper endpoint of a `k > 0` fit.

---

## `tfce.glm`

### `PermutedGLM(Y, design, contrast)`

Fit `Y = X·B`, and re-fit it under permutations **without ever forming the permuted data**.

| argument | | |
| --- | --- | --- |
| `Y` | `(n_elements, n_samples)` | One **row** per element, one column per observation. Note the orientation: elements are the long axis and this keeps them contiguous. |
| `design` | `(n_samples, n_regressors)` | |
| `contrast` | `(n_regressors,)` or `(n_regressors, q)` | A vector gives a t-contrast, a matrix an F-contrast. |

**Attributes**: `df`, `rank`, `pinvX`, `XtX`, `ssy`, `is_t`.

**Raises** `ValueError` if the design leaves no residual degrees of freedom, or if `Y` and the design
disagree about the number of observations.

#### `.fit(perm=None)`

Statistic under one permutation. `perm` is read the way numpy reads a fancy index: the permuted data
is `Y[:, perm]`. `None` is the identity — the unpermuted fit.

**Returns** `(n_elements,)`.

#### `.fit_signs(signs)`

Statistic under one sign-flip, for a one-sample design. `signs` is `(n_samples,)` of ±1.

A one-sample design is not permuted by shuffling — there is nothing to shuffle it with. It is permuted
by flipping signs, which is exchangeable under the null of "mean zero".

### `partition_design(design, contrast)`

Guttman partition: which columns the contrast tests, and which are nuisance.

**Returns** `(idx_interest, idx_nuisance)`.

---

## `tfce.nilearn_compat`

### `calculate_tfce(arr4d, bin_struct, E=0.5, H=2, dh="auto", two_sided_test=True)`

Drop-in for `nilearn.mass_univariate._utils.calculate_tfce`, with the same signature and an exact
transform. See [nilearn.md](nilearn.md).

`dh` is **accepted and ignored** — there is no step size. Passing one is not an error; it simply has
nothing to change.

### `connectivity_from_bin_struct(bin_struct)`

Read `6`, `18` or `26` out of a `scipy.ndimage.generate_binary_structure` array, by counting the
neighbours it actually marks.

---

## Module-level

| | |
| --- | --- |
| `tfce.__version__` | the installed version |
| `tfce.core.VALID_CONNECTIVITY` | `(6, 18, 26)` |
