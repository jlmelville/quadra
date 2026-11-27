# Average Precision-Recall AUC for Label Retrieval

Precision-recall counterpart to
[`roc_auc()`](https://jlmelville.github.io/quadra/reference/roc_auc.md),
often more informative for imbalanced labels.

## Usage

``` r
pr_auc(dm = NULL, labels, ret_extra = FALSE, coords = NULL, query_rows = NULL)
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

Aggregation and undefined-query handling are the same as for
[`roc_auc()`](https://jlmelville.github.io/quadra/reference/roc_auc.md).

## Note

Quadra uses the Davis-Goadrich area returned by
[`PRROC::pr.curve()`](https://rdrr.io/pkg/PRROC/man/pr.curve.html).

## Scaling

For large datasets, pass `coords` and sample `query_rows`. Each query
computes one row of Euclidean distances and ranks all other
observations. The distance workspace grows linearly with the candidate
pool. Scoring all queries still takes quadratic distance work. If you
already have `dm`, you can sample the same queries from it. Subsetting
the coordinates or the matrix before the call also changes the
candidates.

## References

Keilwagen, J., Grosse, I., & Grau, J. (2014). Area under
precision-recall curves for weighted and unweighted data. *PLOS ONE*,
*9*(3), e92209.
[doi:10.1371/journal.pone.0092209](https://doi.org/10.1371/journal.pone.0092209)

Davis, J., & Goadrich, M. (2006, June). The relationship between
Precision-Recall and ROC curves. In *Proceedings of the 23rd
international conference on Machine learning* (pp. 233-240). ACM.
[doi:10.1145/1143844.1143874](https://doi.org/10.1145/1143844.1143874)

Grau, Grosse and Keilwagen (2015) describe the PRROC implementation.
[doi:10.1093/bioinformatics/btv153](https://doi.org/10.1093/bioinformatics/btv153)
