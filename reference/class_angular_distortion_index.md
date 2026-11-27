# Class Angular Distortion Index

Compare how a class appears when viewed from another class. For each
triplet, measure the cosine of the angle at the anchor between the rays
to two observations from a different class.

## Usage

``` r
class_angular_distortion_index(
  Xin,
  Xout,
  labels,
  n_triplets = NULL,
  triplets = NULL,
  ret_extra = FALSE
)
```

## Arguments

- Xin:

  Reference coordinates, as a dense numeric matrix or data frame with
  observations in rows.

- Xout:

  Embedding coordinates for the same observations in the same order.
  Nonnumeric data-frame columns are ignored in both inputs.

- labels:

  A vector without missing values, with one label per observation and at
  least two distinct labels. At least one class must have two members to
  supply endpoints.

- n_triplets:

  Total number of triplets to sample. The default is ten times the
  number of observations.

- triplets:

  Optional one-based plan such as
  [`sample_class_triplets()`](https://jlmelville.github.io/quadra/reference/sample_class_triplets.md).
  The three indices must be distinct, and the endpoints must share a
  class different from the anchor's. When supplied, `n_triplets` is
  ignored.

- ret_extra:

  Whether to return the plan, cosines, and directed class-pair
  summaries.

## Value

The scalar index. With `ret_extra = TRUE`, returns a list with `cadi`,
`triplets`, `cosine_in`, `cosine_out`, and `pairs`. `pairs` contains the
observed directed class pairs, with columns `anchor_class`,
`endpoint_class`, `n_triplets`, and `cadi`. Weighting the pair scores by
their counts recovers the scalar.

## Details

The index is the mean squared change in cosine divided by four, as in
the [authors' code](https://github.com/KaviruGunaratne/cadi). It ranges
from 0 (agreement) to 1. A zero-length ray gives a cosine of zero.
Quadra's sampler gives each eligible triplet equal probability, so class
frequencies affect their representation in the overall score. Reuse a
triplet plan for paired comparisons of embeddings.

The two directions of a class pair can have different distortions.
Uniform scaling, translation, and orthogonal transformations preserve
the angles. Angle agreement emphasizes relative angular spread. For
example, compact classes can subtend small angles under several
arrangements of their centers.

## Scaling

Use the full coordinates and control the cost with `n_triplets`, the
total number of sampled triplets (default `10 * N`). The calculation
processes coordinates in blocks and stores the triplet plan and cosines.
Check the per-pair counts before interpreting a rare class direction.

## Numerical failures

Nonfinite intermediate calculations cause an error. Extreme coordinate
scales may affect accuracy even when the calculations remain finite.

## References

Gunaratne, K., Kobourov, S., & Miller, J. (2026). Class Angular
Distortion Index for Dimensionality Reduction.
<https://arxiv.org/abs/2605.00637>

## See also

[`label_trustworthiness_continuity()`](https://jlmelville.github.io/quadra/reference/label_trustworthiness_continuity.md),
[`roc_auc()`](https://jlmelville.github.io/quadra/reference/roc_auc.md),
[`pr_auc()`](https://jlmelville.github.io/quadra/reference/pr_auc.md)

## Examples

``` r
x <- as.matrix(iris[, -5])
y <- prcomp(x, rank. = 2)$x
set.seed(42)
triplets <- sample_class_triplets(iris$Species, 1000)
class_angular_distortion_index(x, y, iris$Species, triplets = triplets)
#> [1] 0.002392586
```
