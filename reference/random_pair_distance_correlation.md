# Random Pair Distance Correlation

Correlates distances for sampled pairs of observations in `Xin` and
`Xout`. Pearson correlation measures linear association; Spearman
correlation measures rank association.

## Usage

``` r
random_pair_distance_correlation(
  Xin,
  Xout,
  n_pairs = 1000,
  metric_in = "sqeuclidean",
  metric_out = "sqeuclidean",
  method = c("pearson", "spearman"),
  is_transposed = FALSE,
  n_threads = 0,
  pairs = NULL,
  ret_extra = FALSE
)
```

## Arguments

- Xin:

  Input data, with observations in rows by default. Data must be dense
  and finite; nonnumeric data-frame columns are ignored.

- Xout:

  Output data, with the same input conventions and number of
  observations as `Xin`.

- n_pairs:

  Number of pairs to sample.

- metric_in:

  Distance metric for `Xin`. One of `"euclidean"`, `"sqeuclidean"`
  (squared Euclidean), `"cosine"`, `"manhattan"`, `"correlation"` (1
  minus the Pearson correlation), or `"hamming"`.

- metric_out:

  Distance metric for `Xout`. See `metric_in` for details.

- method:

  Correlation method, either `"pearson"` or `"spearman"`.

- is_transposed:

  Whether observations are stored in columns rather than rows.

- n_threads:

  Maximum number of threads to use. `0` or `1` runs serially.

- pairs:

  Optional one-based matrix with one pair per row. The endpoints in each
  row must differ. If supplied, `n_pairs` is ignored.

- ret_extra:

  Whether to return the pairs and their raw distances.

## Value

The correlation, or a list with `correlation`, `pairs`, `distance_in`,
and `distance_out` when `ret_extra = TRUE`.

## Details

`Xin` and `metric_in` define the reference geometry. Supply `pairs` to
reuse exact comparisons, or reset the R seed and keep `n_threads` fixed.

## References

Becht, E., McInnes, L., Healy, J., Dutertre, C. A., Kwok, I. W., Ng, L.
G., ... & Newell, E. W. (2019). Dimensionality reduction for visualizing
single-cell data using UMAP. *Nature biotechnology*, *37*(1), 38-44.

## See also

[`random_pair_distance_emd()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_emd.md),
[`random_pair_distance_stress()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_stress.md),
and
[`random_triplet_accuracy()`](https://jlmelville.github.io/quadra/reference/random_triplet_accuracy.md)
for other global preservation measures.

## Examples

``` r
iris_x <- iris[, -5]
iris_pca2 <- stats::prcomp(iris_x, rank. = 2, scale = FALSE, retx = TRUE)$x
random_pair_distance_correlation(iris_x, iris_pca2)
#> [1] 0.9997105
```
