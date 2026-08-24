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

## References

Lee, J. A., Peluffo-Ordo'nez, D. H., & Verleysen, M. (2015). Multi-scale
similarities in stochastic neighbour embedding: Reducing dimensionality
while preserving both local and global structure. *Neurocomputing*,
*169*, 246-261.

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
