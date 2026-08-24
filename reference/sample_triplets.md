# Sample Anchored Triplets

Generates anchored triplets for reuse with
[`random_triplet_accuracy()`](https://jlmelville.github.io/quadra/reference/random_triplet_accuracy.md).

## Usage

``` r
sample_triplets(n_obs, n_triplets = 5)
```

## Arguments

- n_obs:

  Number of observations.

- n_triplets:

  Number of triplets to sample per observation.

## Value

An integer matrix with columns `anchor`, `endpoint1`, and `endpoint2`,
ordered by anchor.

## Examples

``` r
set.seed(42)
sample_triplets(5, 2)
#>       anchor endpoint1 endpoint2
#>  [1,]      1         2         5
#>  [2,]      1         2         5
#>  [3,]      2         1         3
#>  [4,]      2         1         3
#>  [5,]      3         2         4
#>  [6,]      3         5         2
#>  [7,]      4         2         3
#>  [8,]      4         2         5
#>  [9,]      5         1         4
#> [10,]      5         4         1
```
