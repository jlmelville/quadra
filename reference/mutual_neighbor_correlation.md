# Mutual Neighbor Correlation

Compares the mutual-neighbor count pattern in the input data and output
embedding. For each observation, its mutual-neighbor count is the number
of its first `k` nearest non-self neighbors that also include the
observation among their first `k` nearest non-self neighbors.

## Usage

``` r
mutual_neighbor_correlation(
  Xin,
  Xout,
  k = 15,
  method = c("pearson", "spearman"),
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
  graph. In the latter case, `nn_method_in`, `metric_in` and
  `nn_args_in` are ignored. If `Xin` is a data-frame, non-numeric
  columns are ignored.

- Xout:

  the output data (usually lower dimensional than `Xin`), a matrix or
  data frame with one observation per row, or if `is_transposed = TRUE`,
  one observation per column. Alternatively, it can be a pre-computed
  nearest neighbor graph. In the latter case, `nn_method_out`,
  `metric_out` and `nn_args_out` are ignored. If `Xout` is a data-frame,
  non-numeric columns are ignored.

- k:

  the number of nearest neighbors to use. Can be a numeric vector, in
  which case the mutual neighbor correlation is calculated for each
  value separately.

- method:

  correlation method, either `"pearson"` or `"spearman"`.

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

A named numeric vector of mutual neighbor correlations, one value for
each `k`. Items are named `mnc<k>`, where `<k>` refers to the values
provided in the `k` parameter. If `ret_extra = TRUE`, then a list is
returned containing:

- `mnc`: the vector of mutual neighbor correlations.

- `nn_in`: the nearest neighbor graph for `Xin`.

- `nn_out`: the nearest neighbor graph for `Xout`.

- `mutual_neighbor_in`: a matrix of input mutual-neighbor counts, one
  column per `k`.

- `mutual_neighbor_out`: a matrix of output mutual-neighbor counts, one
  column per `k`.

## Details

This is a local graph diagnostic, not a replacement for nearest-neighbor
preservation. A high value means observations with many or few
mutual-neighbor relationships in the input graph tend to have many or
few mutual-neighbor relationships in the output graph. It complements
[`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md),
which measures whether neighbor identities are preserved directly.

`Xin` and `Xout` can be raw observation matrices or pre-computed nearest
neighbor graphs. Unlike
[`local_radius_correlation()`](https://jlmelville.github.io/quadra/reference/local_radius_correlation.md),
supplied nearest-neighbor graphs only need an `idx` matrix because this
metric uses neighbor identities, not distances. Graphs calculated by
this function are self-excluded. Older cached self-inclusive graphs are
detected and stripped with a warning.

If either mutual-neighbor count vector is constant, the correlation is
undefined and the corresponding result is `NA_real_`.

## See also

[`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md)
for neighbor-identity preservation and
[`local_radius_correlation()`](https://jlmelville.github.io/quadra/reference/local_radius_correlation.md)
for local scale preservation.

## Examples

``` r
iris_pca <- stats::prcomp(iris[, -5], rank. = 2, scale = FALSE, retx = TRUE)$x
mutual_neighbor_correlation(
  iris[, -5],
  iris_pca,
  k = 15,
  nn_method_in = "brute"
)
#>     mnc15 
#> 0.7653062 

cached <- mutual_neighbor_correlation(
  iris[, -5],
  iris_pca,
  k = c(15, 30),
  nn_method_in = "brute",
  ret_extra = TRUE
)
mutual_neighbor_correlation(cached$nn_in, cached$nn_out, k = c(15, 30))
#>     mnc15     mnc30 
#> 0.7680180 0.9102906 
```
