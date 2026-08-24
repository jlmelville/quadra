# Random Triplet Accuracy

Returns the proportion of sampled anchored triplets whose relative
distance ordering in `Xin` is preserved in `Xout`. Input-space ties are
excluded; the result is `NA_real_` if no triplet defines an ordering.

## Usage

``` r
random_triplet_accuracy(
  Xin,
  Xout,
  n_triplets = 5,
  metric_in = "sqeuclidean",
  metric_out = "sqeuclidean",
  is_transposed = FALSE,
  n_threads = 0,
  triplets = NULL,
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

- n_triplets:

  Number of triplets to sample per observation, or a zero-based endpoint
  matrix. A matrix has one column per anchor and a pair of rows per
  triplet. Its endpoints must lie in `[0, n - 1]` and differ from each
  other and the anchor.

- metric_in:

  Distance metric for `Xin`. One of `"euclidean"`, `"sqeuclidean"`
  (squared Euclidean), `"cosine"`, `"manhattan"`, `"correlation"` (1
  minus the Pearson correlation), or `"hamming"`.

- metric_out:

  Distance metric for `Xout`. See `metric_in` for details.

- is_transposed:

  Whether observations are stored in columns rather than rows.

- n_threads:

  Maximum number of threads to use. `0` or `1` runs serially.

- triplets:

  Optional one-based matrix with columns for the anchor and two
  endpoints. All three indices in a row must differ. If supplied,
  numeric `n_triplets` is ignored.

- ret_extra:

  Whether to return the triplets and their row-aligned agreement
  outcomes.

## Value

Triplet accuracy in `[0, 1]`, or `NA_real_` if every input comparison is
tied. With `ret_extra = TRUE`, returns a list with `accuracy`,
`triplets`, and `agreement`. Agreement is `NA` for input-distance ties,
`TRUE` for matching strict orderings, and `FALSE` otherwise.

## Details

`Xin` and `metric_in` define the reference geometry. Supply `triplets`
to reuse exact comparisons, or reset the R seed and keep `n_threads`
fixed. A matrix-valued `n_triplets` uses its documented zero-based
column layout.

## References

Wang, Y., Huang, H., Rudin, C., & Shaposhnik, Y. (2021). Understanding
how dimension reduction tools work: an empirical approach to deciphering
t-SNE, UMAP, TriMAP, and PaCMAP for data visualization. *J Mach. Learn.
Res*, *22*, 1-73. <https://jmlr.org/papers/v22/20-1061.html>.

## See also

[`random_pair_distance_emd()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_emd.md)
and
[`random_pair_distance_correlation()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_correlation.md)
for other global preservation measures.

## Examples

``` r
iris_x <- iris[, -5]
iris_pca2 <- stats::prcomp(iris_x, rank. = 2, scale = FALSE, retx = TRUE)$x
random_triplet_accuracy(iris_x, iris_pca2)
#> [1] 0.976
```
