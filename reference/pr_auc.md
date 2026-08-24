# Average Precision-Recall AUC for Label Retrieval

Precision-recall counterpart to
[`roc_auc()`](https://jlmelville.github.io/quadra/reference/roc_auc.md),
often more informative for imbalanced labels.

## Usage

``` r
pr_auc(dm, labels, ret_extra = FALSE)
```

## Arguments

- dm:

  Finite square numeric distance matrix. The query itself is excluded.

- labels:

  A label vector without missing values, with one label per row of `dm`.

- ret_extra:

  Whether to append row-aligned query AUCs and the number of defined
  queries.

## Value

A list whose first components are:

- `av_auc`: the mean PR AUC over defined queries.

- `label_av`: a named list of per-label means in first-appearance order.

With `ret_extra = TRUE`, the list also contains `query_auc`, with
`NA_real_` for undefined queries, and the integer `n_defined`.

## Details

Aggregation and undefined-query handling are the same as for
[`roc_auc()`](https://jlmelville.github.io/quadra/reference/roc_auc.md).

## Note

Requires the `PRROC` package.

## References

Keilwagen, J., Grosse, I., & Grau, J. (2014). Area under
precision-recall curves for weighted and unweighted data. *PloS One*,
*9*(3), e92209.

Davis, J., & Goadrich, M. (2006, June). The relationship between
Precision-Recall and ROC curves. In *Proceedings of the 23rd
international conference on Machine learning* (pp. 233-240). ACM.
