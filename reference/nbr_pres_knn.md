# Neighborhood Preservation Between Nearest Neighbor Matrices

Calculates the neighborhood preservation for each observation in a
dataset, represented by two matrices of the indices of the nearest
neighbors. The first matrix is the "ground truth", the second being the
estimation or approximation. The neighborhood preservation is calculated
for each row where each element `d[i, k]` is taken to be the index of
the kth nearest neighbor of `i`.

## Usage

``` r
nbr_pres_knn(kin, kout, k = ncol(kin))
```

## Arguments

- kin:

  Nearest neighbor matrix. The "ground truth" or reference indices.

- kout:

  Nearest neighbor matrix. A set of distances to compare to the
  reference indices.

- k:

  The size of the neighborhood, where k is the number of neighbors to
  include in the neighborhood.

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
