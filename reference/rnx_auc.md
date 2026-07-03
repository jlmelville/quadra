# Area Under the RNX Curve

The RNX curve is formed by calculating the `rnx_crm` metric for
different sizes of neighborhood. Each value of RNX is scaled according
to the natural log of the neighborhood size, to give a higher weight to
smaller neighborhoods. An AUC of 1 indicates perfect neighborhood
preservation, an AUC of 0 is due to random results. Self-neighbors on
the distance-matrix diagonal are excluded before the co-ranking matrix
is calculated.

## Usage

``` r
rnx_auc(din, dout)
```

## Arguments

- din:

  Input distance matrix.

- dout:

  Output distance matrix.

## Value

Area under the RNX curve.

## References

Lee, J. A., Peluffo-Ordo'nez, D. H., & Verleysen, M. (2015). Multi-scale
similarities in stochastic neighbour embedding: Reducing dimensionality
while preserving both local and global structure. *Neurocomputing*,
*169*, 246-261.
