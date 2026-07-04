# Local Radius Correlation

Compares the local scale around each observation in the input data and
output embedding. The local scale is summarized from nearest-neighbor
distances, either as the distance to the `k`th nearest non-self neighbor
or as the mean distance to the first `k` nearest non-self neighbors.

## Usage

``` r
local_radius_correlation(
  Xin,
  Xout,
  k = 15,
  statistic = c("radius", "mean"),
  method = c("spearman", "pearson"),
  log = FALSE,
  nn_method_in = "nnd",
  metric_in = "sqeuclidean",
  nn_method_out = "brute",
  metric_out = "sqeuclidean",
  is_transposed = FALSE,
  n_threads = 0,
  verbose = FALSE,
  ret_extra = FALSE,
  nn_args_in = list(),
  nn_args_out = list()
)
```

## Arguments

- Xin:

  the input data (usually high-dimensional), a matrix or data frame with
  one observation per row, or if `is_transposed = TRUE`, one observation
  per column. Alternatively, it can be a pre-computed nearest neighbor
  graph with `idx` and `dist` matrix elements.

- Xout:

  the output data (usually lower dimensional than `Xin`), a matrix or
  data frame with one observation per row, or if `is_transposed = TRUE`,
  one observation per column. Alternatively, it can be a pre-computed
  nearest neighbor graph with `idx` and `dist` matrix elements.

- k:

  the number of nearest neighbors to use. Can be a numeric vector, in
  which case the local radius correlation is calculated for each value
  separately.

- statistic:

  the local scale statistic. `"radius"` uses the distance to the `k`th
  nearest neighbor. `"mean"` uses the mean distance to the first `k`
  nearest neighbors.

- method:

  correlation method, either `"spearman"` or `"pearson"`.

- log:

  if `TRUE`, calculate the correlation on log local radii.

- nn_method_in:

  the nearest neighbor method to calculate the neighbors of `Xin`. Can
  be one of `"brute"` (brute force calculation) or `"nnd"`, the nearest
  neighbor descent method of Dong and co-workers (2011).

- metric_in:

  the distance calculation to apply to `Xin`. One of `"euclidean"`,
  `"sqeuclidean"` (squared Euclidean), `"cosine"`, `"manhattan"`,
  `"correlation"` (1 minus the Pearson correlation), or `"hamming"`.

- nn_method_out:

  the nearest neighbor method to calculate the neighbors of `Xout`. See
  `nn_method_in` for details.

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

  the maximum number of threads to use.

- verbose:

  if `TRUE`, log information about the calculation to the console.

- ret_extra:

  if `TRUE`, additionally return the nearest neighbor graphs for `Xin`
  and `Xout`.

- nn_args_in:

  list of extra arguments to pass to the nearest neighbor methods,
  [`rnndescent::brute_force_knn()`](https://jlmelville.github.io/rnndescent/reference/brute_force_knn.html)
  or
  [`rnndescent::nnd_knn()`](https://jlmelville.github.io/rnndescent/reference/nnd_knn.html),
  depending on the value of `nn_method_in`.

- nn_args_out:

  list of extra arguments to pass to the nearest neighbor methods,
  [`rnndescent::brute_force_knn()`](https://jlmelville.github.io/rnndescent/reference/brute_force_knn.html)
  or
  [`rnndescent::nnd_knn()`](https://jlmelville.github.io/rnndescent/reference/nnd_knn.html),
  depending on the value of `nn_method_out`.

## Value

A named numeric vector of local radius correlations, one value for each
`k`. Items are named `lrc<k>`, where `<k>` refers to the values provided
in the `k` parameter. If `ret_extra = TRUE`, then a list is returned
containing:

- `lrc`: the vector of local radius correlations.

- `nn_in`: the nearest neighbor graph for `Xin`.

- `nn_out`: the nearest neighbor graph for `Xout`.

- `scale_in`: a matrix of input local scale values, one column per `k`.

- `scale_out`: a matrix of output local scale values, one column per
  `k`.

## Details

This is a local scale diagnostic, not a density estimator. A high value
means observations with small or large local radii in the input data
tend to have small or large local radii in the output embedding. It
complements
[`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md),
which measures whether neighbor identities are preserved.

`Xin` and `Xout` can be raw observation matrices or pre-computed nearest
neighbor graphs. Unlike
[`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md),
supplied nearest-neighbor graphs must contain a `dist` matrix as well as
an `idx` matrix because this metric uses neighbor distances. Graphs
calculated by this function are self-excluded. Older cached
self-inclusive graphs are detected and stripped with a warning.

If either local-radius vector is constant, the correlation is undefined
and the corresponding result is `NA_real_`. If `log = TRUE`, all
selected local radius values must be positive, so duplicate points that
produce zero local radii should be handled before requesting a log-scale
comparison.

## See also

[`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md)
for neighbor-identity preservation and
[`trustworthiness()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)
for exact rank-penalty neighborhood preservation on distance matrices.

## Examples

``` r
iris_pca <- stats::prcomp(iris[, -5], rank. = 2, scale = FALSE, retx = TRUE)$x
local_radius_correlation(
  iris[, -5],
  iris_pca,
  k = 15,
  nn_method_in = "brute"
)
#>     lrc15 
#> 0.8930213 

cached <- local_radius_correlation(
  iris[, -5],
  iris_pca,
  k = c(15, 30),
  nn_method_in = "brute",
  ret_extra = TRUE
)
local_radius_correlation(cached$nn_in, cached$nn_out, k = c(15, 30))
#>     lrc15     lrc30 
#> 0.8930213 0.9397574 
```
