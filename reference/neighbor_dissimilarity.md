# Neighbor Dissimilarity Between Two Graphs

Compare how many neighbors each pair of observations shares in two
directed nearest-neighbor graphs on the same observations.

## Usage

``` r
neighbor_dissimilarity(nn_in, nn_out, k = 15, ret_extra = FALSE)
```

## Arguments

- nn_in, nn_out:

  Reference and output graphs with the same rows in the same order.
  Supply lists with a nearest-first `idx` matrix as returned by
  [`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md).
  Distances are optional. Self indices are stripped.

- k:

  Number of non-self neighbors to use from each row.

- ret_extra:

  Return the separate scores for gains and losses and the row pairs
  whose shared-neighbor counts changed?

## Value

A scalar, with smaller values indicating closer shared-neighbor
structure. With `ret_extra = TRUE`, a list containing `nd`, `gain` and
`loss` (the nonnegative directional norms), and `pairs`. The data frame
contains one row per changed unordered pair (`i < j`), with one-based
row positions and neighbor counts `shared_in` and `shared_out`.

## Details

For each graph, divide the number of shared neighbors of each pair by
`k`. Subtract the reference fractions from the output fractions. Compute
the Frobenius norms of the positive and negative parts of this
difference matrix. The larger norm is Neighbor Dissimilarity (ND). Both
orders of each pair enter the norms. Zero means all shared counts agree,
which can happen even when neighbor identities change.

ND depends on graph size and `k`. Compare embeddings on the same rows
and at the same `k`.

## Scaling

Each input graph needs only `k` neighbors per row, but the sparse matrix
product represents every pair of rows sharing a neighbor. One neighbor
listed by all `N` rows creates `N * (N - 1) / 2` unordered pairs. Even
the scalar score can therefore need quadratic memory. `ret_extra = TRUE`
also returns a table of changed pairs. The calculation uses the `Matrix`
package.

## References

Fujiwara, T., Kuo, Y.-H., Ynnerman, A., and Ma, K.-L. (2023). Feature
Learning for Nonlinear Dimensionality Reduction toward Maximal
Extraction of Hidden Patterns, equation 3.
<https://arxiv.org/abs/2206.13891> The authors implement the graph
comparison in [FEALM](https://github.com/takanori-fujiwara/fealm).
Quadra also reports which observation pairs gained or lost shared
neighbors.

## See also

[`mutual_neighbor_correlation()`](https://jlmelville.github.io/quadra/reference/mutual_neighbor_correlation.md),
[`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md)

## Examples

``` r
x <- as.matrix(iris[, -5])
g <- nn_preservation(x, prcomp(x)$x[, 1:2], k = 5,
                     nn_method_in = "brute", ret_extra = TRUE)
neighbor_dissimilarity(g$nn_in, g$nn_out, k = 5)
#> [1] 9.887366
```
