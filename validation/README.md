# Validation suite

Calibration and correctness checks for the TFCE toolbox. Run everything with

```matlab
addpath /path/to/spm12
spm('defaults','fmri')
cd /path/to/tfce
compile          % the mex-files must be built
validation/run_all
```

Each script prints a `[PASS]` / `[FAIL]` line per check and a summary.

The scripts do **not** re-implement the functions they test. Where a function is
a subfunction of `tfce_estimate_stat.m`, `val_util('extract', ...)` pulls it out
of the source file at run time, so the suite always exercises the code that
ships rather than a copy of it that could drift.

| script | what it establishes |
|---|---|
| `val_tfce_exactness` | The max-tree really is the exact TFCE integral, and the batched transform is identical to the sequential one — including under the single `calc_neg` flag that the permutation loop now shares across a whole block of permutations. |
| `val_gamma` | The Gamma fit to the maximum distribution gives calibrated FWE p-values. This is the most consequential check: it enters *every* corrected p-value the toolbox reports. |
| `val_pareto` | The Generalised Pareto fit to the tail of each element's permutation distribution recovers uncorrected p-values *below* the 1/n_perm floor that counting cannot reach, is unbiased, never returns zero, and is never less accurate than the counting it overrides. |
| `val_sequential` | Stopping early once the image is decisively null reaches the same answer as running every permutation, keeps the false-positive rate at alpha, and never cuts short an image that is significant or borderline. |
| `val_half_permutations` | The half-permutation shortcut is lossless, and correctly refuses unbalanced designs. |
| `val_nuisance` | Draper-Stoneman, Freedman-Lane and Smith all control the false-positive rate, including when the nuisance variable is correlated with the effect of interest. |
| `val_glm_fast` | The accelerated GLM returns exactly the statistic `calc_GLM` returns, for t- and F-contrasts under all three nuisance methods, and refuses itself when its algebraic precondition does not hold. |
| `val_voxel_covariate` | The voxel-/vertex-wise covariate path controls the false-positive rate for t- and F-contrasts, both when the covariate is the effect of interest and when it is a nuisance. Its batched solve returns what the per-voxel pseudoinverse loop returns, hands rank-deficient elements back to `pinv`, shrinks `ResMS` towards the maximum over all elements, and permutes a multi-column covariate whenever any of its columns is tested. Also reports its cost. |

## Notes on the checks

**Exactness.** The reference is an *independent* implementation in plain MATLAB
that steps through thresholds and labels components with `bwlabeln` (volume) or
the mesh adjacency (surface). It must converge onto the max-tree result at first
order as the step size goes to zero. The criterion is the slope of
log(error) vs log(n_steps), which should be −1; testing the ratio of successive
errors instead is too noisy, because the error is a maximum over elements.

**Gamma approximation.** Counting exceedances cannot resolve p-values below
1/n_perm, which is why the maximum distribution is approximated by a Gamma fitted
to its first three moments. The check compares the fit against counting *in the
tail only* — the bulk is irrelevant for inference — with a tolerance of three
standard errors of the counting estimator, `sqrt(p(1-p)/n_perm)`. A tighter
tolerance would be testing the noise of the reference rather than the fit.

**Batched permutations.** The TFCE transform dominates the cost of the
permutation loop, and permutations are independent of one another, so the loop
works through them a block at a time and hands the whole block to
`tfceMex_maxtree_batch`, which runs one permutation per thread. Only the order of
the work changes; the permutations are consumed in exactly the order they would
have been, and a block size of 1 reproduces the old permutation-by-permutation
path. One thing does genuinely change, and is checked: the loop used to decide
`calc_neg` per map, and now shares one flag across the block. That is the same
thing only if asking for negative TFCE values on a map that has none leaves it
unchanged — otherwise a block mixing signed and unsigned maps (an F-statistic, or
a t-statistic that happens to be all-positive) would come out differently.

Do not expect the speed-up to scale with the core count. The max-tree is bound by
memory rather than by arithmetic, and on an 8-core machine it saturates at about
2x from 6 threads onwards, getting *worse* beyond 8. The default block size is one
permutation per computational thread, which sits on that plateau.

**Single precision and the sort.** The max-tree carries the statistic, and
everything derived from it that is one value per element, as `float`. The map is a
permuted statistic and single precision is all it is ever worth. The integral
itself is *not* computed in float: `cum` accumulates a node's contribution along a
root-to-leaf path and `pow(birth, H+1)` spans a wide dynamic range, so those stay
double. Only the stored values lose precision, not the arithmetic combining them.

On its own this bought nothing — the transform is bound by random access (the sort
and the union-find), not by streaming bandwidth, and it measured the same to
within noise. What it *enables* is the sort: every value reaching the sort is
strictly positive, and for positive IEEE-754 floats the bit pattern read as a
uint32 increases monotonically with the value, so the key needs no transformation
and four counting passes replace `qsort`'s n log n comparisons through a function
pointer. That is where the time was. The two together take the transform from
0.128 to 0.081 s per permutation sequentially, and to 0.051 batched.

Both mex files accept `single` or `double` and return the class they were given,
so a caller that has not moved to single sees no change other than the precision
of the values themselves. `val_tfce_exactness` still recovers the first-order
convergence of dh-stepping onto the max-tree, which is what establishes that it is
still the exact TFCE integral.

**Sequential stopping.** The loop no longer always runs `n_perm` permutations. It
stops once the observed *global maximum* has been exceeded often enough by the
permutation maxima (Besag & Clifford, 1991; the negative binomial method of
Winkler et al., 2016). The global maximum is the right thing to watch because no
element rests on fewer exceedances than it does — an element with a smaller
statistic is exceeded at least as often — so once it is settled, every corrected
p-value in the image is settled with it.

Counting exceedances alone is *not* a safe rule, and the suite shows why. An image
whose true corrected p-value sits right at alpha reaches any fixed number of
exceedances quickly, and stopping there leaves the decision resting on an estimate
far too coarse to place it on one side of alpha or the other: the truncated run
and the full run then disagree about 40% of the time, in exactly the case where
the answer matters most. So the rule is not "enough exceedances" but "enough
exceedances **and** the estimate is 3 standard errors clear of the largest alpha
asked about" — which is what bounds the risk of the sequential decision differing
from the full run's (Gandy, 2009). A decisively null image clears that within the
floor of 500 permutations; anything significant, or anywhere near the boundary,
never clears it and runs the full `n_perm`.

The floor is not optional: the Gamma fit to the maximum distribution and the
Pareto fits to the element-wise tails need a few hundred permutations before they
have anything to work with, however quickly the FWE question itself settles. The
check also sits *after* the permutation-null calibration check in the loop, so
stopping early can never skip it.

**Pareto tail.** Counting exceedances cannot report a p-value below `1/n_perm`,
and that floor — not the FWE-corrected null, which the Gamma fit already handles
and which converges far sooner — is what forces a permutation test to run many
thousands of permutations. It also caps FDR, which is computed from the
uncorrected p-values. A Generalised Pareto distribution is therefore fitted to
the tail of each element's permutation distribution and the p-value read off it.

Extrapolating a tail from a hundred points cannot be exact, so the check is not
"is it right" but "does it beat what it replaces". The suite builds a real
permutation distribution from 50000 permutations, shows the fit only the first
1000 or 5000, and scores both the fit and plain counting against the truth. The
fit has to be unbiased, and has to land within a factor of two of the truth at
least as often as counting does. Two structural properties matter as much as the
accuracy. The fit must never return zero — counting does so constantly once the
true p drops below the floor, and a zero becomes `Inf` in the `-log10` map the
toolbox writes out. And it must leave alone the elements counting already
resolves. Both are checked.

The zero case is the subtle one. Probability-weighted moments put the shape `k`
above zero for roughly half the elements simply by chance, which gives the fitted
distribution a *finite upper end point*; where the observed statistic lies beyond
it, the fit returns exactly zero. But the observation is itself proof that the end
point is wrong, so such a fit is rejected rather than believed, and the
exponential limit (`k = 0`, infinite support) answers instead.

**Pooling the shape.** The elements all carry the permutation distribution of the
same statistic under the same design, so they differ in *scale*, not in *shape*.
Estimating `k` separately for each element therefore spends a hundred tail values
on a number they all share, and the error in it is what dominates the error of the
extrapolation. It is pooled instead — estimated once over elements sampled across
the whole image, with only the scale fitted per element, which then follows in
closed form from the mean of the exceedances. The suite checks the premise
directly: the per-element estimates of `k` scatter by ~0.15, but their median
agrees to ~0.02 between one half of the image and the other, which is exactly the
signature of one shared number seen through noise. Pooling lifts the share of
elements landing within a factor of two of the truth at p = 1e-4 from 45% to 60%
(5000 permutations) without moving the median, and the accuracy floors in the
script are set to catch a regression below that.

**Accelerated GLM.** `calc_GLM_fast` is not an approximation but an algebraic
rearrangement, so the tolerance is the single-precision rounding floor (1e-5
relative) rather than a statistical one. Two identities carry it. For any design,
`R = I - X*pinv(X)` is a projector, so the residual sum of squares is
`||y||^2 - Beta'*(X'*X)*Beta` and never needs the residual matrix. And under
Freedman-Lane the nuisance columns Z are part of the full design X, so `R*Rz = R`
and Rz drops out of the residuals altogether — which is what allows the
`n_vox x n_data x n_data` multiply `Y*(Pset'*Rz)` to disappear from the loop. The
second identity requires Z to lie in the space spanned by the design that is
actually fitted, and `calc_GLM` fits `W*X` while Rz is built from the unwhitened
Z, so it fails for a non-identity whitening matrix. `prepare_GLM_fast` checks
this rather than assuming it, and the suite confirms that it refuses the fast
path in exactly that case.

**Voxel-wise covariates.** The design of an element differs from `X` only in the
columns `xC.cols`, so `Xi'*Xi` and `Xi'*y` can be assembled for every element at
once and the small systems batch-solved, instead of taking a pseudoinverse of a
fresh design matrix per element. Like the accelerated GLM this is a rearrangement
rather than an approximation, so the tolerance is the single-precision rounding
floor. Two things make it subtle. Working from `Xi'*Xi` squares the condition
number, and a covariate that barely varies across subjects is nearly collinear
with the intercept, so the products are accumulated in double even though the
data are single. And an element whose covariate is *exactly* constant makes `Xi`
rank-deficient; an LU solve does not fail there but returns a huge finite answer
through a tiny pivot, so degeneracy is detected from the Schur complement of the
fixed block (the variance the covariate has left once the fixed columns are
projected out of it) and those elements are handed back to `pinv`.

Three things changed behaviour here rather than just speed, and each has its own
check. `ResMS` is now shrunk towards the largest residual variance in the image,
as `calc_GLM` does; the per-element loop applied that line to a *scalar*, where
`max(ResMS(isfinite(ResMS)))` is that same scalar, so it collapsed to
`ResMS*1.001` and the path was in effect not shrinking at all. A covariate spread
over several columns is now permuted whenever *any* of them is tested; the flag
that used to decide this was reassigned on every pass of its loop, so the last
column decided on its own and the covariate images could stay with their original
subjects while the design was permuted around them. And F-contrasts are now
supported, which the path used to refuse outright: replacing the contrast by an
orthonormal basis of its column space leaves ESS unchanged and makes the middle
matrix invertible, so the `pinv` of the textbook formula becomes a batched solve.
The F of a 1-df contrast is checked against the square of its t.

**Calibration.** All false-positive rates are for the *uncorrected* permutation
p-value at α = 0.05 unless stated otherwise, over independent elements. The
tolerance is ±0.015, which is generous relative to the Monte Carlo error of the
estimate itself.

## What is NOT covered

- Real data. Everything here is on synthetic phantoms with known ground truth.
- Comparison against PALM and FSL `randomise`.
- FDR correction (`snpm_P_FDR`) is used as inherited from SnPM and is not
  re-validated here.
- The design-recognition audit and the diagnostics sensitivity analysis were run
  during development but are not yet part of this suite.
