# Neighborhood Preservation Between Nearest Neighbor Matrices

Returns the overlap between reference and comparison nearest-neighbor
matrices, separately for each observation.

## Usage

``` r
nbr_pres_knn(kin, kout, k = ncol(kin), n_threads = 0)
```

## Arguments

- kin:

  Reference nearest-neighbor index matrix, with observations in rows and
  neighbors in nearest-first order.

- kout:

  Nearest-neighbor index matrix to compare with `kin`.

- k:

  Neighborhood size.

- n_threads:

  Maximum number of threads to use. `0` or `1` runs serially.

## Value

Per-observation neighborhood overlap in `[0, 1]`.

## Details

Rows contain distinct one-based indices in nearest-first order.
Self-indices are removed, so self-inclusive matrices need at least
`k + 1` columns. Supplied order resolves ties.

## Scaling

With `k` stored neighbors per observation, each input uses `N * k`
indices. Reuse these lists across comparisons to avoid rebuilding
neighbor graphs.

## References

Chen and Buja (2009) count the neighbors shared by two neighborhoods and
divide by the neighborhood size.
[doi:10.1198/jasa.2009.0111](https://doi.org/10.1198/jasa.2009.0111)

## See also

[`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md)
for computing the graphs from coordinates.
