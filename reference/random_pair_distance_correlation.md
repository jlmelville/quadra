# Random Pair Distance Correlation

Evaluates the preservation of global structure of dimensionality
reduction results using the Pearson correlation coefficient between
randomly selected distances, similar to the method of Becht and
co-workers (2019).

## Usage

``` r
random_pair_distance_correlation(
  Xin,
  Xout,
  n_pairs = 1000,
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

- n_pairs:

  the number of random pairs of observations to calculate distances for
  in the input and output space.

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

The Pearson correlation between the distances in the input and output
space. For randomly distributed data, the expected value is 0.

## Details

This function repeatedly samples random pairs of observation and
calculates the distance between the points in both the original data and
the embedding space. The Pearson correlation coefficient between the two
sets of distances is reported. This differs slightly from the procedure
in the Becht paper which randomly samples a subset of observations and
then exhaustively calculates all pair-wise distances within that subset.

## References

Becht, E., McInnes, L., Healy, J., Dutertre, C. A., Kwok, I. W., Ng, L.
G., ... & Newell, E. W. (2019). Dimensionality reduction for visualizing
single-cell data using UMAP. *Nature biotechnology*, *37*(1), 38-44.

## See also

[`random_pair_distance_emd()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_emd.md)
and
[`random_triplet_accuracy()`](https://jlmelville.github.io/quadra/reference/random_triplet_accuracy.md)
for another measure of global structure preservation.

## Examples

``` r
iris_pca2 <- stats::prcomp(iris[, -5], rank. = 2, scale = FALSE, retx = TRUE)$x
random_pair_distance_correlation(iris, iris_pca2)
#> [1] 0.9997105

# If you plan on comparing the results of multiple output methods, then
# pre-transposing the input data can save time
tiris <- t(iris[, -5])
iris_pca1 <- stats::prcomp(iris[, -5], rank. = 1, scale = FALSE, retx = TRUE)$x
iris_pca3 <- stats::prcomp(iris[, -5], rank. = 3, scale = FALSE, retx = TRUE)$x
random_pair_distance_correlation(tiris, t(iris_pca1), is_transposed = TRUE)
#> [1] 0.9977412
random_pair_distance_correlation(tiris, t(iris_pca2), is_transposed = TRUE)
#> [1] 0.9997525
random_pair_distance_correlation(tiris, t(iris_pca3), is_transposed = TRUE)
#> [1] 0.9999674
```
