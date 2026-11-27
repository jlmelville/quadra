# Pointwise RNX Agreement Across Neighborhood Sizes

For each observation, compare which neighbors are shared by the
reference and the embedding across neighborhood sizes.

## Usage

``` r
pointwise_rnx(din, dout, k = NULL)
```

## Arguments

- din:

  Reference distance matrix.

- dout:

  Distance matrix to compare with `din`.

- k:

  Optional unique neighborhood sizes from 1 through `n - 2`, where `n`
  is the number of observations. Supply these sizes to return curve
  columns alongside the per-observation AUC. The default returns AUC
  only.

## Value

A list with `auc`, a numeric vector in observation order, and `curve`, a
matrix with one row per observation and columns `rnx<k>` in requested
order, or `NULL` when `k` is omitted. Both retain the reference row
names if present. The mean of `auc` equals
[`rnx_auc()`](https://jlmelville.github.io/quadra/reference/rnx_auc.md).
The column means of `curve` equal
[`rnx_curve()`](https://jlmelville.github.io/quadra/reference/rnx_curve.md)
at the requested sizes, up to numerical precision.

## Details

RNX adjusts nearest-neighbor overlap for chance. Its value at a single
`k` is an affine transformation of overlap for that observation. The AUC
averages RNX over every `k` from 1 through `n - 2`, weighting each value
by `1 / k`, regardless of the curve sizes requested.

Distances use the same exact ranking and column-order tie convention as
[`rnx_curve()`](https://jlmelville.github.io/quadra/reference/rnx_curve.md).
Diagonal entries are excluded. A value of 1 means perfect agreement, 0
is the random-neighborhood baseline, and values can be negative.

## Scaling

Supply two dense `N` by `N` distance matrices. At 5,000 observations
their double entries use about 381.5 MiB together, before working
storage. For Euclidean coordinates,
[`exact_rank_profiles()`](https://jlmelville.github.io/quadra/reference/exact_rank_profiles.md)
ranks one query at a time. Set `query_rows` to score a sample against
all other observations.

## References

Lee et al. (2013) define the RNX curve.
[doi:10.1016/j.neucom.2012.12.036](https://doi.org/10.1016/j.neucom.2012.12.036)

Lee, Peluffo-Ordóñez and Verleysen (2015) summarize the curve across
neighborhood sizes with weights proportional to `1 / k`.
[doi:10.1016/j.neucom.2014.12.095](https://doi.org/10.1016/j.neucom.2014.12.095)

## See also

[`pointwise_rank_preservation()`](https://jlmelville.github.io/quadra/reference/pointwise_rank_preservation.md),
[`rnx_curve()`](https://jlmelville.github.io/quadra/reference/rnx_curve.md),
[`rnx_auc()`](https://jlmelville.github.io/quadra/reference/rnx_auc.md)

## Examples

``` r
din <- as.matrix(dist(iris[, -5]))
dout <- as.matrix(dist(prcomp(iris[, -5], rank. = 2)$x))
scores <- pointwise_rnx(din, dout, k = c(5, 15, 40))
head(scores$auc)
#>         1         2         3         4         5         6 
#> 0.9004926 0.6447428 0.6861277 0.5879835 0.6954234 0.9432241 
colMeans(scores$curve)
#>      rnx5     rnx15     rnx40 
#> 0.6219815 0.7771177 0.9239052 
```
