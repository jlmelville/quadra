[![Last Commit](https://img.shields.io/github/last-commit/jlmelville/quadra)](https://github.com/jlmelville/quadra)
[![R-CMD-check](https://github.com/jlmelville/quadra/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jlmelville/quadra/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/github/jlmelville/quadra/graph/badge.svg?token=HxTB6T6KmW)](https://codecov.io/github/jlmelville/quadra)

# Quadra: QUality Assessment of Dimensionality Reduction Algorithms

Quadra is an R package for measuring how well a dimensionality-reduction
embedding preserves structure from the input data.

## Installation

Install the development version from GitHub with
[pak](https://pak.r-lib.org/):

```R
pak::pak("jlmelville/quadra")
```

Quadra contains C++ code, so installing from source requires the appropriate
build tools. On Windows, install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/). On macOS, see the
[R for macOS FAQ](https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html#Installation-of-source-packages)
if compilation or linking fails.

## Quick Start

Nearest-neighbor preservation asks how many of each observation's input-space
neighbors remain neighbors in the embedding. Here is a small PCA example:

```R
library(quadra)

iris_x <- as.matrix(iris[, -5])
iris_pca <- stats::prcomp(iris_x, rank. = 2)$x

nn_preservation(
  iris_x, iris_pca, k = c(5, 15),
  nn_method_in = "brute"
)
```

Here `nnp5` is 0.64. On average, the PCA plot keeps just over three of each flower's five nearest
neighbors. Five neighbors chosen at random from the other 149 flowers would give an expected score
of about 0.034. Widen the neighborhood to 15 flowers and the fraction retained rises to about 0.80.

## Learn More

The [metric overview](https://jlmelville.github.io/quadra/articles/metrics-overview.html) helps
you choose a measure and decide whether to score the whole dataset or a sample. For individual
metrics, see the guides to
[local preservation](https://jlmelville.github.io/quadra/articles/local-preservation.html),
[global preservation](https://jlmelville.github.io/quadra/articles/global-preservation.html),
or [label retrieval](https://jlmelville.github.io/quadra/articles/label-retrieval.html).

For a worked example, pick the question closest to yours.

- Can you trust the digit next to yours in a UMAP plot?
  [Compare PCA and UMAP on MNIST](https://jlmelville.github.io/quadra/articles/mnist-pca-umap.html).
- Should a map preserve distances across a fold or along the sheet?
  [Unroll a Swiss roll](https://jlmelville.github.io/quadra/articles/reference-geometry.html).
- How can two clouds merge while their short distances still look perfect?
  [Compare the bands](https://jlmelville.github.io/quadra/articles/distance-band-fidelity.html).

## License

[GPLv3 or later](https://www.gnu.org/licenses/gpl-3.0.txt).
