# Neighborhood Preservation Between Nearest Neighbor Matrices

Calculates the neighborhood preservation for each observation in a
dataset, represented by two matrices of the indices of the nearest
neighbors. The first matrix is the "ground truth", the second being the
estimation or approximation. The neighborhood preservation is calculated
for each row where each element `d[i, k]` is taken to be the index of
the kth nearest neighbor of `i`.

## Usage

``` r
nbr_pres_knn(kin, kout, k = ncol(kin), n_threads = 0)
```

## Arguments

- kin:

  Nearest-neighbor index matrix. The "ground truth" or reference
  indices, with observations in rows and neighbors in nearest-first
  order.

- kout:

  Nearest-neighbor index matrix to compare with `kin`, using the same
  row and ordering conventions.

- k:

  The size of the neighborhood, where k is the number of neighbors to
  include in the neighborhood.

- n_threads:

  the maximum number of threads to use. `0` or `1` runs serially.

## Value

Vector of preservation values, one for each row of `kin`.

## Details

Approximate nearest neighbor methods, e.g.
[RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy), can find
k-nearest neighbors quite efficiently and so makes calculating
preservation values for larger datasets feasible.

The neighborhood preservation can vary between 0 (no neighbors in
common) and 1 (perfect preservation). For nearest-neighbor matrices that
exclude self-neighbors, random performance gives an approximate value of
k / (n - 1), where k is the size of the neighborhood and n is the number
of observations or items in the dataset.

Rows contain distinct one-based indices in nearest-first order. All
supplied columns are checked before self-indices are removed and the
first `k` non-self neighbors are retained. Self-inclusive inputs
therefore need at least `k + 1` columns. Supplied order resolves ties.
