# Random Triplet Accuracy

Evaluates the preservation of global structure of dimensionality
reduction results using the random triplet accuracy method of Wang and
co-workers (2020).

## Usage

``` r
random_triplet_accuracy(
  Xin,
  Xout,
  n_triplets = 5,
  metric_in = "sqeuclidean",
  metric_out = "sqeuclidean",
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

- n_triplets:

  the number of triplets per observation to generate.

- metric_in:

  the distance calculation to apply to `Xin`. One of `"euclidean"`,
  `"sqeuclidean"` (squared Euclidean), `"cosine"`, `"manhattan"`,
  `"correlation"` (1 minus the Pearson correlation), or `"hamming"`.

- metric_out:

  the distance metric to apply to `Xout`. See `metric_in` for details.

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

The triplet accuracy, ranging from 0 (no relative distances agree) to 1
(all of them agree). For randomly distributed `Xout`, the accuracy will
be 0.5.

## Details

The random triplet accuracy is calculated by randomly selecting three
points in the input data and calculating the distances for two sides of
the resulting triangle. This is repeated for the output data, and the
relative ordering of the distances are compared. The returned accuracy
is the proportion of triangles where the relative distances agree
between the input and output data. Triplets with tied input-space
distances are excluded from the denominator because they do not define a
relative ordering. If no sampled input triplets define an ordering, the
result is `NA_real_`.

## References

Wang, Y., Huang, H., Rudin, C., & Shaposhnik, Y. (2021). Understanding
how dimension reduction tools work: an empirical approach to deciphering
t-SNE, UMAP, TriMAP, and PaCMAP for data visualization. *J Mach. Learn.
Res*, *22*, 1-73. <https://jmlr.org/papers/v22/20-1061.html>.

## See also

[`random_pair_distance_emd()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_emd.md)
and
[`random_pair_distance_correlation()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_correlation.md)
for another measure of global structure preservation.

## Examples

``` r
iris_pca2 <- stats::prcomp(iris[, -5], rank. = 2, scale = FALSE, retx = TRUE)$x
random_triplet_accuracy(iris, iris_pca2)
#> [1] 0.9826667

# If you plan on comparing the results of multiple output methods, then
# pre-transposing the input data can save time
tiris <- t(iris[, -5])
iris_pca1 <- stats::prcomp(iris[, -5], rank. = 1, scale = FALSE, retx = TRUE)$x
iris_pca3 <- stats::prcomp(iris[, -5], rank. = 3, scale = FALSE, retx = TRUE)$x
random_triplet_accuracy(tiris, t(iris_pca1), is_transposed = TRUE)
#> [1] 0.928
random_triplet_accuracy(tiris, t(iris_pca2), is_transposed = TRUE)
#> [1] 0.9813333
random_triplet_accuracy(tiris, t(iris_pca3), is_transposed = TRUE)
#> [1] 0.9933333
```
