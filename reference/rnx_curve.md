# RNX Curve

Returns exact rank-based neighborhood agreement at selected neighborhood
sizes. A value of 1 is perfect, 0 is the random-neighborhood baseline,
and values may be negative. Diagonal entries are excluded; off-diagonal
entries must be finite. Ties follow column order. The curve diagnoses
scale-specific preservation; it does not establish a universally
preferred embedding.

## Usage

``` r
rnx_curve(din, dout, k = NULL)
```

## Arguments

- din:

  Reference distance matrix.

- dout:

  Distance matrix to compare with `din`.

- k:

  Optional unique neighborhood sizes from 1 through `n - 2`, where `n`
  is the number of observations. By default, return all of them.

## Value

A named numeric vector. Items are named `rnx<k>`.

## Scaling

Supply two dense `N` by `N` distance matrices. At 5,000 observations
their double entries use about 381.5 MiB together, before working
storage. For Euclidean coordinates,
[`exact_rank_profiles()`](https://jlmelville.github.io/quadra/reference/exact_rank_profiles.md)
ranks one query at a time. Set `query_rows` to score a sample against
all other observations.

## References

Lee, Renard, Bernard, Dupont and Verleysen (2013) define the RNX curve
by adjusting neighborhood overlap for chance.
[doi:10.1016/j.neucom.2012.12.036](https://doi.org/10.1016/j.neucom.2012.12.036)

## See also

[`rnx_auc()`](https://jlmelville.github.io/quadra/reference/rnx_auc.md)
for a weighted summary of the complete curve.

## Examples

``` r
iris_pca <- stats::prcomp(iris[, -5], rank. = 2, scale = FALSE, retx = TRUE)$x
din <- as.matrix(stats::dist(iris[, -5]))
dout <- as.matrix(stats::dist(iris_pca))
rnx_curve(din, dout, k = c(5, 15, 30))
#>      rnx5     rnx15     rnx30 
#> 0.6219815 0.7771177 0.8942670 
```
