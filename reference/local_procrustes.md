# Local Procrustes Error

Measure how closely each reference neighborhood retains its geometry in
an embedding after allowing translation, rotation, and reflection. Each
patch contains its anchor and the same `k` reference neighbors in both
spaces.

## Usage

``` r
local_procrustes(
  Xin,
  Xout,
  k = 15,
  nn_method_in = c("brute", "nnd"),
  is_transposed = FALSE,
  n_threads = 0,
  nn_in = NULL,
  ret_extra = FALSE
)
```

## Arguments

- Xin:

  Reference coordinates, as a dense numeric matrix or data frame.

- Xout:

  Embedding coordinates for the same observations in the same order.
  Nonnumeric data-frame columns are ignored in both inputs.

- k:

  Number of reference neighbors per anchor. The anchor is included in
  addition to these neighbors.

- nn_method_in:

  Reference neighbor search method. Use `"brute"` for exact search or
  `"nnd"` for approximate neighbor descent.

- is_transposed:

  Whether observations are stored in columns rather than rows.

- n_threads:

  Maximum number of threads to use. `0` or `1` runs serially.

- nn_in:

  Optional reference neighbor graph, such as `nn_in` from a detailed
  [`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md)
  result. Its first `k` non-self indices are used in each row. Supplying
  it skips neighbor search.

- ret_extra:

  Whether to also return the row-aligned residuals and graph.

## Value

The mean residual over defined reference patches, or `NA_real_` when all
are undefined. With `ret_extra = TRUE`, a list with `error`,
`pointwise_error`, `n_defined`, and `nn_in`. The pointwise vector
follows the rows and names of the reference observations.

## Details

Coordinates are centered within each patch. The alignment uses
orthogonal least squares, with zero padding to compare different
dimensions. Scale is fixed. To account for arbitrary embedding units,
multiply the entire embedding by a single calibration factor before
calling this function.

The residual is the sum of squared aligned coordinate differences
divided by the reference patch's centered sum of squares. Zero indicates
exact agreement up to numerical precision. Values can exceed one. A
constant reference patch has an undefined residual (`NA_real_`). Moving
an intact group can leave its internal patch scores unchanged.

Automatic neighbor search uses Euclidean ordering. Supplied graphs
select the patch members in their given order. The Procrustes fit always
uses the coordinate geometry of `Xin` and `Xout`.

## Scaling

Each alignment uses one patch of `k + 1` observations. With `k` fixed,
the number of fits grows linearly with `N`. The default brute-force
reference search compares all pairs. Supply `nn_in` to reuse a reference
graph across embeddings, or use `nn_method_in = "nnd"` to find
approximate neighbors.

## References

Goldberg, Y., & Ritov, Y. (2009). Local Procrustes for Manifold
Embedding: A Measure of Embedding Quality and Embedding Algorithms.
<https://arxiv.org/abs/0806.2669> Quadra follows their construction by
matching the same observations in both patches. The fit keeps scale
fixed within each patch. [ZADU
0.5.4](https://github.com/hj-n/zadu/tree/v0.5.4) selects neighbors
separately in each space.

## See also

[`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md),
[`local_radius_correlation()`](https://jlmelville.github.io/quadra/reference/local_radius_correlation.md)

## Examples

``` r
x <- as.matrix(iris[, -5])
embedding <- prcomp(x, rank. = 2)$x
result <- local_procrustes(x, embedding, k = 10, ret_extra = TRUE)
result$error
#> [1] 0.304925
head(result$pointwise_error)
#> [1] 0.4039810 0.2655268 0.4390545 0.2287525 0.3931438 0.2965453
```
