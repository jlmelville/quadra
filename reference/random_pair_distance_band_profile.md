# Distance-Band Fidelity for Random Pairs

Compare the order and magnitude of matched distances within cumulative
bands of input-space distances. Each band includes pairs at or below its
input radius, including ties.

## Usage

``` r
random_pair_distance_band_profile(
  Xin,
  Xout,
  bands = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1),
  band_type = c("proportion", "radius"),
  n_pairs = 10000,
  metric_in = "euclidean",
  metric_out = "euclidean",
  is_transposed = FALSE,
  n_threads = 0,
  pairs = NULL,
  ret_extra = FALSE
)
```

## Arguments

- Xin:

  Input data, with observations in rows by default. Data must be dense
  and finite; nonnumeric data-frame columns are ignored.

- Xout:

  Output data, with the same input conventions and number of
  observations as `Xin`.

- bands:

  Strictly increasing, finite band selectors. Proportions must be in
  `(0, 1]`; radii must be nonnegative.

- band_type:

  Interpret `bands` as empirical proportions or as radii in the units of
  `metric_in`. Proportions use the inverse empirical distribution of the
  evaluated input distances.

- n_pairs:

  Number of pairs to sample.

- metric_in:

  Distance metric for `Xin`. One of `"euclidean"`, `"sqeuclidean"`
  (squared Euclidean), `"cosine"`, `"manhattan"`, `"correlation"` (1
  minus the Pearson correlation), or `"hamming"`.

- metric_out:

  Distance metric for `Xout`. See `metric_in` for details.

- is_transposed:

  Whether observations are stored in columns rather than rows.

- n_threads:

  Maximum number of threads to use. `0` or `1` runs serially.

- pairs:

  Optional one-based matrix with one pair per row. The endpoints in each
  row must differ. If supplied, `n_pairs` is ignored.

- ret_extra:

  Whether to return the profile and the evaluated pairs with their raw
  distances.

## Value

A data frame with columns `band_type`, `band`, `radius`, `n_pairs`,
`pair_fraction`, `spearman`, and `normalized_stress`. If
`ret_extra = TRUE`, a list containing that data frame as `profile` along
with `pairs`, `distance_in`, `distance_out`, and `fitted_scale`.
Multiply output distances by the corresponding entry of `fitted_scale`
to put them in input-distance units. The multiplier is undefined (NA)
for an empty band or all-zero output distances. A multiplier larger than
the floating-point range is returned as `Inf`.

## Details

Both distance metrics default to ordinary Euclidean distance. When you
supply all unordered pairs, the profile is exact for those observations.
Repeated pair rows each count toward the result. With proportion bands,
ties can make the realized fraction larger than the requested
proportion. Set `bands = 1` in proportion mode to use every evaluated
pair.

Each band reports Spearman correlation and normalized stress. Stress
divides the residual Euclidean norm by the input distance vector's norm,
after fitting one least-squares multiplier to the output distances
within that band. A uniformly rescaled band has zero stress. Separate
bands can fit different multipliers.

Spearman correlation needs at least three pairs and variation in both
distance vectors. Check `n_pairs` in the narrowest band and increase the
pair sample if few distances qualify. Stress needs nonzero norms in both
vectors. Undefined values are `NA_real_` with one warning per call.

## Numerical failures

A nonfinite computed distance causes an error, even with finite
coordinates. Finite distances alone do not guarantee accuracy at extreme
coordinate scales.

## Scaling

For large datasets, use the full coordinates and set `n_pairs` to a
total sample size. Distance calculations and stored distances scale with
the number of sampled pairs. Reuse `pairs` across embeddings. Supplying
every unordered pair makes both the pair count and distance work
quadratic.

## References

Takaeda (2026) compares distance ranks and fitted-scale stress within a
fixed input radius. <https://arxiv.org/abs/2608.21779> His
[benchmark](https://github.com/toorpia/dr-fidelity-benchmark) defines
the band scores. Quadra estimates them from sampled pairs by default.

Smelser, Miller and Kobourov (2024) fit the scale that minimizes
normalized stress. Takaeda fits it separately for each band.
<https://arxiv.org/abs/2408.07724>

## See also

[`sample_pairs()`](https://jlmelville.github.io/quadra/reference/sample_pairs.md),
[`random_pair_distance_correlation()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_correlation.md),
[`random_pair_distance_stress()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_stress.md)

## Examples

``` r
x <- as.matrix(iris[, -5])
y <- prcomp(x, rank. = 2)$x
set.seed(42)
pairs <- sample_pairs(nrow(x), n_pairs = 2000)
random_pair_distance_band_profile(x, y, bands = c(0.2, 1), pairs = pairs)
#>    band_type band   radius n_pairs pair_fraction  spearman normalized_stress
#> 1 proportion  0.2 0.900000     400           0.2 0.7479456        0.25987538
#> 2 proportion  1.0 6.857113    2000           1.0 0.9954018        0.04192187
```
