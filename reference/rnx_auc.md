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

[`rnx_curve()`](https://jlmelville.github.io/quadra/reference/rnx_curve.md)
for agreement at individual neighborhood sizes.
