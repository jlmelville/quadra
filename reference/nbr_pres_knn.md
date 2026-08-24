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
