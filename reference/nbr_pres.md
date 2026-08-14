# Neighborhood Preservation Between Distance Matrices

Calculates the neighborhood preservation for each observation in a
dataset, represented by two distance matrices. The first matrix is the
"ground truth", the second being the estimation or approximation. The
neighborhood preservation is calculated for each row where each element
`d[i, j]` is taken to be the distance between observation `i` and `j`.

## Usage

``` r
nbr_pres(din, dout, k)
```

## Arguments

- din:

  Distance matrix. The "ground truth" or reference distances.

- dout:

  Distance matrix. A set of distances to compare to the reference
  distances.

- k:

  The size of the neighborhood, where k is the number of neighbors to
  include in the neighborhood.

## Value

Vector of preservation values, one for each row of the distance matrix.

## Details

The neighborhood preservation can vary between 0 (no neighbors in
common) and 1 (perfect preservation). However, random performance gives
an approximate value of k / (n - 1), where k is the size of the
neighborhood and n is the number of observations or items in the
dataset.

Self-neighbors on the diagonal are excluded from each row before the
neighborhood overlap is calculated.

If the `k`th-smallest distance is tied, every observation at that
distance is treated as a neighbor. The shared count is capped at `k`
before it is divided by `k`. This differs from exact-rank metrics such
as
[`trustworthiness()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)
and
[`rnx_auc()`](https://jlmelville.github.io/quadra/reference/rnx_auc.md),
which break ties by original column order, and from supplied-graph
metrics, which use the first `k` indices.

Distance matrices are used as supplied without a finite-value check.

## Note

This is not a very efficient way to calculate the preservation if you
want to calculate the value for multiple values of `k`. For more global
measures of preservation, see
[`rnx_auc()`](https://jlmelville.github.io/quadra/reference/rnx_auc.md).
