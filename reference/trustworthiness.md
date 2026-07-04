# Trustworthiness and Continuity Between Distance Matrices

`trustworthiness()` penalizes observations that appear among the `k`
nearest neighbors in `dout` but have input-space rank greater than `k`
in `din`. `continuity()` applies the dual penalty to input-space
neighbors that are no longer among the `k` nearest neighbors in `dout`.

## Usage

``` r
trustworthiness(din, dout, k)

continuity(din, dout, k)
```

## Arguments

- din:

  Input distance matrix. The "ground truth" or reference distances.

- dout:

  Output distance matrix. A set of distances to compare to the reference
  distances.

- k:

  The size of the neighborhood. Must be a positive integer less than
  half the number of observations so the standard 0-1 normalization
  remains bounded.

## Value

A scalar score. A value of 1 indicates no rank-penalty errors at
neighborhood size `k`; lower values indicate worse preservation.

## Details

Both functions use exact ranks from the supplied distance matrices and
exclude the diagonal self-neighbor from each row. Tied distances are
ranked in their original column order after self-neighbor exclusion,
matching `rank(ties.method = "first")`.

Because these functions require full `n` by `n` distance matrices, they
are practical only for small datasets. For larger datasets, use
nearest-neighbor preservation metrics such as
[`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md)
or
[`nbr_pres_knn()`](https://jlmelville.github.io/quadra/reference/nbr_pres_knn.md).

Unlike
[`nbr_pres()`](https://jlmelville.github.io/quadra/reference/nbr_pres.md),
which only counts shared neighbors, these metrics weight each unexpected
or missing neighbor by how far its rank lies outside the
`k`-neighborhood.
[`rnx_auc()`](https://jlmelville.github.io/quadra/reference/rnx_auc.md)
also uses rank-based neighborhood agreement, but aggregates across
neighborhood sizes; these functions report the standard trustworthiness
or continuity score at one `k`.

## References

Venna, J., & Kaski, S. (2001). Neighborhood preservation in nonlinear
projection methods: An experimental study. In *Artificial Neural
Networks - ICANN 2001* (pp. 485-491).

## Examples

``` r
iris_pca <- stats::prcomp(iris[, -5], rank. = 2, scale = FALSE, retx = TRUE)$x
din <- as.matrix(stats::dist(iris[, -5]))
dout <- as.matrix(stats::dist(iris_pca))
trustworthiness(din, dout, k = 15)
#> [1] 0.9859563
continuity(din, dout, k = 15)
#> [1] 0.9927384
```
