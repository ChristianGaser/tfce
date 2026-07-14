# tfce - documentation

Exact threshold-free cluster enhancement, and the permutation inference around it.

The [package README](../README.md) is the two-minute version. These pages are the rest.

## Start here

| | |
| --- | --- |
| **[Installation](installation.md)** | `pip install tfce`, wheels, building from source, optional extras. |
| **[Usage](usage.md)** | Volumes, surfaces, batches of permutations, connectivity, `E` and `H`. |
| **[A permutation test, end to end](permutation.md)** | The whole thing: GLM → TFCE → FWE and uncorrected p-values. A worked example you can run. |

## Going deeper

| | |
| --- | --- |
| **[Theory](theory.md)** | Why the integral is exact, what a max-tree is, and why the tail approximations mean you can run far fewer permutations. |
| **[Using it with nilearn](nilearn.md)** | The drop-in, the connectivity mismatch, and a bug in nilearn's stepped TFCE worth knowing about. |
| **[API reference](api.md)** | Every public function, with its arguments and what they do. |
| **[Benchmarks](benchmarks.md)** | Accuracy and speed, with the scripts to reproduce them. |
| **[Validation](validation.md)** | What the test suite establishes, and how it establishes it. |

## The one-paragraph summary

The TFCE of an element is an integral over the extent of the cluster containing it, taken across all
heights. Implementations normally approximate that integral by stepping the height over a grid and
summing, which costs a step size, a discretisation error that depends on it, and an accuracy
parameter the caller has to guess. This one builds the **max-tree** (the component tree) with
union-find and integrates each piece of the extent function in closed form, because that function is
piecewise constant. The result is the integral, not a sample of it, and there is nothing to tune.

It is the same C the [MATLAB TFCE toolbox](https://github.com/ChristianGaser/tfce) runs - not a
reimplementation - so the two are bit-identical, and the test suite holds them to it.
