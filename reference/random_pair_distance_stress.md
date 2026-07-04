# Random Pair Distance Stress

Evaluates global distance preservation with a sampled stress summary.

## Usage

``` r
random_pair_distance_stress(
  Xin,
  Xout,
  n_pairs = 1000,
  metric_in = "sqeuclidean",
  metric_out = "sqeuclidean",
  range_scale = TRUE,
  is_transposed = FALSE,
  n_threads = 0
)
```

## Arguments

- Xin:

  the input data (usually high-dimensional), a matrix or data frame with
  one observation per row, or if `is_transposed = TRUE`, one observation
  per column.

- Xout:

  the output data (usually lower dimensional than `Xin`), a matrix or
  data frame with one observation per row, or if `is_transposed = TRUE`,
  one observation per column.

- n_pairs:

  the number of random pairs of observations to calculate distances for
  in the input and output space.

- metric_in:

  the distance calculation to apply to `Xin`. One of `"euclidean"`,
  `"sqeuclidean"` (squared Euclidean), `"cosine"`, `"manhattan"`,
  `"correlation"` (1 minus the Pearson correlation), or `"hamming"`.

- metric_out:

  the distance metric to apply to `Xout`. See `metric_in` for details.

- range_scale:

  if `TRUE` (the default) then scale each sampled distance vector to the
  range 0-1 before calculating stress.

- is_transposed:

  if `TRUE` then `Xin` and `Xout` are assumed to have been passed in
  transposed format, i.e. with one observation per column. Otherwise,
  `Xin` and `Xout` will be transposed. For large datasets, transposing
  can be slow, so if this function will be called multiple times with
  the same input data, it is more efficient to transpose the input data
  once outside of this function and set `is_transposed = TRUE`.

- n_threads:

  the maximum number of threads to use. `0` or `1` runs serially.

## Value

The sampled stress between the matched distances in the input and output
space.

## Details

This function repeatedly samples random pairs of observations and
calculates the distance between the points in both the original data and
the embedding space. The returned value is the root mean squared
difference between the matched sampled distances.

By default, each sampled distance vector is scaled to the range 0-1
before stress is calculated. This makes the result comparable across
embeddings with different distance scales, but it also means the value
is mainly useful for comparing methods under identical sampling and
scaling settings.

## See also

[`random_pair_distance_correlation()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_correlation.md),
[`random_pair_distance_emd()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_emd.md),
and
[`random_triplet_accuracy()`](https://jlmelville.github.io/quadra/reference/random_triplet_accuracy.md)
for other measures of global structure preservation.

## Examples

``` r
iris_pca2 <- stats::prcomp(iris[, -5], rank. = 2, scale = FALSE, retx = TRUE)$x
random_pair_distance_stress(iris, iris_pca2)
#> [1] 0.005466072

# If you plan on comparing the results of multiple output methods, then
# pre-transposing the input data can save time
tiris <- t(iris[, -5])
iris_pca1 <- stats::prcomp(iris[, -5], rank. = 1, scale = FALSE, retx = TRUE)$x
iris_pca3 <- stats::prcomp(iris[, -5], rank. = 3, scale = FALSE, retx = TRUE)$x
random_pair_distance_stress(tiris, t(iris_pca1), is_transposed = TRUE)
#> [1] 0.01893072
random_pair_distance_stress(tiris, t(iris_pca2), is_transposed = TRUE)
#> [1] 0.006364007
random_pair_distance_stress(tiris, t(iris_pca3), is_transposed = TRUE)
#> [1] 0.001815861
```
