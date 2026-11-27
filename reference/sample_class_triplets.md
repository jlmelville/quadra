# Sample Triplets for Class Angular Distortion

The anchor comes from one class and two distinct endpoints come from
another. Draws are independent and may repeat. Each eligible triplet has
equal probability, so a directed class pair with `a` anchor members and
`b` endpoint members has weight proportional to `a * choose(b, 2)`.

## Usage

``` r
sample_class_triplets(labels, n_triplets = NULL)
```

## Arguments

- labels:

  A vector without missing values, with one label per observation and at
  least two distinct labels. At least one class must have two members to
  supply endpoints.

- n_triplets:

  Total number of triplets to sample. The default is ten times the
  number of observations.

## Value

An integer matrix with one-based columns `anchor`, `endpoint1`, and
`endpoint2`.

## References

Gunaratne, Kobourov and Miller (2026) define the class triplets used by
CADI. <https://arxiv.org/abs/2605.00637> Their
[implementation](https://github.com/KaviruGunaratne/cadi) samples
classes first. Quadra gives each eligible triplet equal probability.

## See also

[`class_angular_distortion_index()`](https://jlmelville.github.io/quadra/reference/class_angular_distortion_index.md)

## Examples

``` r
set.seed(42)
sample_class_triplets(c("A", "A", "A", "B", "B", "C"), 5)
#>      anchor endpoint1 endpoint2
#> [1,]      6         5         4
#> [2,]      1         5         4
#> [3,]      5         1         3
#> [4,]      1         5         4
#> [5,]      5         3         1
```
