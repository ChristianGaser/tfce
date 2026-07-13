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
| `val_tfce_exactness` | The max-tree really is the exact TFCE integral, and the batched transform is identical to the sequential one. |
| `val_gamma` | The Gamma fit to the maximum distribution gives calibrated FWE p-values. This is the most consequential check: it enters *every* corrected p-value the toolbox reports. |
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
