# Area Under the RNX Curve

Summarizes rank-based neighborhood agreement across neighborhood sizes,
weighting smaller neighborhoods more heavily. A value of 1 is perfect, 0
is the random-neighborhood baseline, and values may be negative.
Diagonal entries are excluded; off-diagonal entries must be finite. Ties
follow column order.

## Usage

``` r
rnx_auc(din, dout)
```

## Arguments

- din:

  Reference distance matrix.

- dout:

  Distance matrix to compare with `din`.

## Value

Area under the RNX curve.

## References

Lee, J. A., Peluffo-Ordo'nez, D. H., & Verleysen, M. (2015). Multi-scale
similarities in stochastic neighbour embedding: Reducing dimensionality
while preserving both local and global structure. *Neurocomputing*,
*169*, 246-261.

## See also

[`rnx_curve()`](https://jlmelville.github.io/quadra/reference/rnx_curve.md)
for agreement at individual neighborhood sizes.
