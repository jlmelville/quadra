# Sample Observation Pairs

Generates pairs for reuse with the random-pair distance metrics.

## Usage

``` r
sample_pairs(n_obs, n_pairs = 1000)
```

## Arguments

- n_obs:

  Number of observations.

- n_pairs:

  Number of pairs to sample.

## Value

An integer matrix with one row per pair and columns `endpoint1` and
`endpoint2`.

## Details

Pairs are sampled independently with replacement across rows, so
repeated or reversed pairs can occur. Endpoints within each row are
distinct.

## See also

[`random_pair_distance_correlation()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_correlation.md),
[`random_pair_distance_emd()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_emd.md),
and
[`random_pair_distance_stress()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_stress.md).

## Examples

``` r
set.seed(42)
sample_pairs(5, 3)
#>      endpoint1 endpoint2
#> [1,]         1         2
#> [2,]         5         2
#> [3,]         1         5
```
