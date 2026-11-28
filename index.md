# Quadra: QUality Assessment of Dimensionality Reduction Algorithms

Quadra is an R package for measuring how well a dimensionality-reduction
embedding preserves structure from the input data.

## Installation

Install the development version from GitHub with
[pak](https://pak.r-lib.org/):

``` r

pak::pak("jlmelville/quadra")
```

Quadra requires rnndescent 0.3.0 or later. If you are upgrading a source
installation, rebuild Quadra after updating rnndescent so its compiled
code uses the new distance functions.

Quadra contains C++ code, so installing from source requires the
appropriate build tools. On Windows, install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/). On macOS, see
the [R for macOS
FAQ](https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html#Installation-of-source-packages)
if compilation or linking fails.

## Quick Start

Nearest-neighbor preservation asks how many of each observation’s
input-space neighbors remain neighbors in the embedding. Here is a small
PCA example:

``` r

library(quadra)

iris_x <- as.matrix(iris[, -5])
iris_pca <- stats::prcomp(iris_x, rank. = 2)$x

nn_preservation(
  iris_x, iris_pca, k = c(5, 15),
  nn_method_in = "brute"
)
```

Here `nnp5` is 0.64. On average, the PCA plot keeps just over three of
each flower’s five nearest neighbors. Five neighbors chosen at random
from the other 149 flowers would give an expected score of about 0.034.
Widen the neighborhood to 15 flowers and the fraction retained rises to
about 0.80.

## Learn More

The [metric
overview](https://jlmelville.github.io/quadra/articles/metrics-overview.html)
helps you choose a measure and decide whether to score the whole dataset
or a sample. For individual metrics, see the guides to [local
preservation](https://jlmelville.github.io/quadra/articles/local-preservation.html),
[global
preservation](https://jlmelville.github.io/quadra/articles/global-preservation.html),
or [label
retrieval](https://jlmelville.github.io/quadra/articles/label-retrieval.html).

For a worked example, pick the question closest to yours.

- Can you trust the digit next to yours in a UMAP plot? [Compare PCA,
  UMAP and t-SNE on
  MNIST](https://jlmelville.github.io/quadra/articles/mnist-pca-umap.html).
- What happens to the scores when two digits set the scale of a map?
  [Compare two LTSA
  fits](https://jlmelville.github.io/quadra/articles/ltsa-map-spread.html).
- Should a map preserve distances across a fold or along the sheet?
  [Unroll a Swiss
  roll](https://jlmelville.github.io/quadra/articles/reference-geometry.html).
- How can two clouds merge while their short distances still look
  perfect? [Compare the
  bands](https://jlmelville.github.io/quadra/articles/distance-band-fidelity.html).

## Acknowledgments

[ZADU](https://github.com/hj-n/zadu), the Python library by Hyeon Jeon
and colleagues, pointed me to several measures, including the [Gap
Index](https://arxiv.org/abs/2607.28324) of Ros, Arleo and Paulovich.
For Steadiness and Cohesiveness, Quadra follows ZADU 0.5.4 in how it
weights shared neighbors and how it combines losses in the pointwise
maps.

- [drnb](https://github.com/jlmelville/drnb) inspired three comparisons
  in Quadra. I compare the distance ranks of sampled pairs, estimate
  stress from sampled pairs, and compare counts of reciprocal neighbors
  between spaces.
- The distance-band profile comes from [Yoshio Takaeda’s
  benchmark](https://github.com/toorpia/dr-fidelity-benchmark). It asks
  how well distances are preserved among pairs that fall within a given
  radius in the input.
- For the Class Angular Distortion Index, Quadra follows [Gunaratne,
  Kobourov and Miller’s code](https://github.com/KaviruGunaratne/cadi).
  The cosine of an angle can change by at most two, so dividing the mean
  squared change by four puts the score on a scale from zero to one.
- The idea of comparing distance distributions with Earth Mover’s
  Distance comes from [Heiser and Lau’s work on single-cell
  data](https://github.com/KenLauLab/DR-structure-preservation). Quadra
  estimates those distributions by sampling pairs of observations.

## License

[GPLv3 or later](https://www.gnu.org/licenses/gpl-3.0.txt).
