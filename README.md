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
iris_pca <- stats::prcomp(iris_x, retx = TRUE, rank. = 2)$x

nn_preservation(iris_x, iris_pca, k = c(5, 15), nn_method_in = "brute")
```

The result contains one score for each `k`. A score of 1 means every requested
neighbor was preserved; for unrelated neighborhoods, the expected score is
approximately `k / (n - 1)`.

## Learn More

Start with the [metric overview](https://jlmelville.github.io/quadra/articles/metrics-overview.html),
then use the focused guides for
[local preservation](https://jlmelville.github.io/quadra/articles/local-preservation.html),
[global preservation](https://jlmelville.github.io/quadra/articles/global-preservation.html),
or [label retrieval](https://jlmelville.github.io/quadra/articles/label-retrieval.html).

## License

[GPLv3 or later](https://www.gnu.org/licenses/gpl-3.0.txt).
