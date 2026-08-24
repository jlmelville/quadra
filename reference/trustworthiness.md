# Trustworthiness and Continuity Between Distance Matrices

`trustworthiness()` penalizes observations that appear among the `k`
nearest neighbors in `dout` but have input-space rank greater than `k`
in `din`. `continuity()` applies the dual penalty to input-space
neighbors that are no longer among the `k` nearest neighbors in `dout`.

## Usage

``` r
trustworthiness(din, dout, k)

continuity(din, dout, k)
```

## Arguments

- din:

  Reference distance matrix.

- dout:

  Distance matrix to compare with `din`.

- k:

  One or more unique neighborhood sizes, each less than half the number
  of observations.

## Value

For scalar `k`, an unnamed rank-preservation score. For multiple values,
a named numeric vector with items `trustworthiness<k>` or
`continuity<k>`. A value of 1 indicates no penalty.

## Details

Diagonal entries are excluded; off-diagonal entries must be finite.
Unlike
[`nbr_pres()`](https://jlmelville.github.io/quadra/reference/nbr_pres.md),
the penalty reflects how far a misplaced neighbor falls beyond `k`.

## References

Venna, J., & Kaski, S. (2001). Neighborhood preservation in nonlinear
projection methods: An experimental study. In *Artificial Neural
Networks - ICANN 2001* (pp. 485-491).

## Examples

``` r
iris_pca <- stats::prcomp(iris[, -5], rank. = 2, scale = FALSE, retx = TRUE)$x
din <- as.matrix(stats::dist(iris[, -5]))
dout <- as.matrix(stats::dist(iris_pca))
trustworthiness(din, dout, k = 15)
#> [1] 0.9859563
continuity(din, dout, k = 15)
#> [1] 0.9927384
trustworthiness(din, dout, k = c(5, 15, 30))
#>  trustworthiness5 trustworthiness15 trustworthiness30 
#>         0.9786667         0.9859563         0.9935587 
```
