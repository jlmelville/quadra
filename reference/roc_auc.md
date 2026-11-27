# Average ROC AUC for Label Retrieval

Ranks other observations by distance from each query and treats matching
labels as relevant.

## Usage

``` r
roc_auc(dm = NULL, labels, ret_extra = FALSE, coords = NULL, query_rows = NULL)
```

## Arguments

- dm:

  Finite square numeric distance matrix. Supply exactly one of `dm` and
  `coords`.

- labels:

  A label vector without missing values, with one label per row of the
  supplied data.

- ret_extra:

  Whether to return query AUCs in the requested row order and the number
  of defined queries.

- coords:

  Dense Euclidean coordinates, one observation per row. Numeric columns
  are used for data frames.

- query_rows:

  Unique one-based row positions to score, in the order requested. The
  default scores every row. Only the query itself is excluded from its
  candidate pool.

## Value

A list whose first components are:

- `av_auc`: the mean AUC over defined queries.

- `label_av`: a named list of per-label means in first-appearance order.

With `ret_extra = TRUE`, the list also contains `query_auc`, with
`NA_real_` for undefined queries, and the integer `n_defined`. When
supplied, `query_rows` identifies the returned queries. Means use the
selected queries. Labels without a defined query have mean `NA_real_`.

## Details

`av_auc` is the mean ROC AUC over defined queries; `label_av` gives
per-label means in first-appearance order. Queries without both relevant
and irrelevant observations are excluded. An empty average is
`NA_real_`.

## Note

[`PRROC::roc.curve()`](https://rdrr.io/pkg/PRROC/man/roc.curve.html)
computes the area for each query.

## Scaling

For large datasets, pass `coords` and sample `query_rows`. Each query
computes one row of Euclidean distances and ranks all other
observations. The distance workspace grows linearly with the candidate
pool. Scoring all queries still takes quadratic distance work. If you
already have `dm`, you can sample the same queries from it. Subsetting
the coordinates or the matrix before the call also changes the
candidates.

## References

Bamber (1975) relates ROC area to ranking probability.
[doi:10.1016/0022-2496(75)90001-2](https://doi.org/10.1016/0022-2496%2875%2990001-2)

Venna et al. (2010) describe the retrieval view of dimensionality
reduction. <https://www.jmlr.org/papers/v11/venna10a.html>

Grau, Grosse and Keilwagen (2015) describe the PRROC implementation.
[doi:10.1093/bioinformatics/btv153](https://doi.org/10.1093/bioinformatics/btv153)

## See also

[`pr_auc()`](https://jlmelville.github.io/quadra/reference/pr_auc.md),
[`class_pair_auc()`](https://jlmelville.github.io/quadra/reference/class_pair_auc.md)

## Examples

``` r
x <- as.matrix(iris[, 1:4])
roc_auc(coords = x, labels = iris$Species, query_rows = c(1, 51, 101))
#> $av_auc
#> [1] 0.9172109
#> 
#> $label_av
#> $label_av$setosa
#> [1] 1
#> 
#> $label_av$versicolor
#> [1] 0.7717347
#> 
#> $label_av$virginica
#> [1] 0.979898
#> 
#> 
```
