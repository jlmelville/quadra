# Neighborhood Preservation Between Distance Matrices

Returns the overlap between the `k`-neighborhoods defined by two
distance matrices, separately for each observation.

## Usage

``` r
nbr_pres(din, dout, k)
```

## Arguments

- din:

  Reference distance matrix.

- dout:

  Distance matrix to compare with `din`.

- k:

  Neighborhood size.

## Value

Per-observation neighborhood overlap in `[0, 1]`.

## Details

Diagonal entries are excluded; off-diagonal entries must be finite. Ties
at the `k`th distance are included, with each score capped at 1. A
collapsed output can therefore score 1. Inspect output-distance
magnitudes or coordinate spread alongside neighborhood scores when
distances are tied.

## Scaling

Supply two dense `N` by `N` distance matrices. For large datasets, use
[`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md)
or
[`nbr_pres_knn()`](https://jlmelville.github.io/quadra/reference/nbr_pres_knn.md)
to score overlap from neighbor lists. The graph functions use the stored
order to resolve boundary ties.

## References

Chen and Buja (2009) count the neighbors shared by two neighborhoods and
divide by the neighborhood size.
[doi:10.1198/jasa.2009.0111](https://doi.org/10.1198/jasa.2009.0111)

## See also

[`pointwise_rank_preservation()`](https://jlmelville.github.io/quadra/reference/pointwise_rank_preservation.md)
for rank penalties per observation.
