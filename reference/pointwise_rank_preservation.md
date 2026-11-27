# Pointwise Trustworthiness and Continuity

Locate false neighbors and missing neighbors for each observation.
Trustworthiness penalizes neighbors that are close in the embedding but
ranked beyond `k` in the reference. Continuity penalizes reference
neighbors that fall beyond `k` in the embedding.

## Usage

``` r
pointwise_rank_preservation(din, dout, k)
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

A list with matrices `trustworthiness` and `continuity`. Rows follow the
input observations and retain the row names of `din`, if present.
Columns are named `k<value>` in requested order. Column means reproduce
the corresponding aggregate scores up to numerical precision.

## Details

Like
[`trustworthiness()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)
and
[`continuity()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md),
this function uses exact distance matrices, excludes the diagonal, and
breaks ties by column order. Ranking is shared between the two measures.

## Scaling

Supply two dense `N` by `N` distance matrices. At 5,000 observations
their double entries use about 381.5 MiB together, before working
storage. For Euclidean coordinates,
[`exact_rank_profiles()`](https://jlmelville.github.io/quadra/reference/exact_rank_profiles.md)
ranks one query at a time. Set `query_rows` to score a sample against
all other observations.

## References

Venna, J., & Kaski, S. (2001). Neighborhood preservation in nonlinear
projection methods: An experimental study. In *Artificial Neural
Networks - ICANN 2001* (pp. 485-491).
[doi:10.1007/3-540-44668-0_68](https://doi.org/10.1007/3-540-44668-0_68)

## See also

[`trustworthiness()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md),
[`continuity()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)

## Examples

``` r
x <- as.matrix(iris[, -5])
y <- prcomp(x, rank. = 2)$x
pointwise <- pointwise_rank_preservation(
  as.matrix(dist(x)), as.matrix(dist(y)), k = c(5, 15)
)
head(pointwise$trustworthiness)
#>          k5       k15
#> 1 0.9929577 0.9963255
#> 2 0.9915493 0.9926509
#> 3 0.9943662 1.0000000
#> 4 0.9873239 1.0000000
#> 5 0.9929577 0.9947507
#> 6 0.9971831 0.9968504
```
