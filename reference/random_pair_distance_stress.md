# Random Pair Distance Stress

Returns the root mean squared difference between matched sampled input
and output distances. Each distance vector is scaled to `[0, 1]` by
default.

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

- range_scale:

  Whether to scale each sampled distance vector to `[0, 1]`.

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

The stress, or a list with `stress`, `pairs`, `distance_in`, and
`distance_out` when `ret_extra = TRUE`.

## Details

`Xin` and `metric_in` define the reference geometry. Supply `pairs` to
reuse exact comparisons, or reset the R seed and keep `n_threads` fixed.

## See also

[`random_pair_distance_correlation()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_correlation.md),
[`random_pair_distance_emd()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_emd.md),
and
[`random_triplet_accuracy()`](https://jlmelville.github.io/quadra/reference/random_triplet_accuracy.md)
for other measures of global structure preservation.

## Examples

``` r
iris_x <- iris[, -5]
iris_pca2 <- stats::prcomp(iris_x, rank. = 2, scale = FALSE, retx = TRUE)$x
random_pair_distance_stress(iris_x, iris_pca2)
#> [1] 0.005397666
```
