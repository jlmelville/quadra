# Average ROC AUC for Label Retrieval

Ranks other observations by distance from each query and treats matching
labels as relevant.

## Usage

``` r
roc_auc(dm, labels, ret_extra = FALSE)
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

- `av_auc`: the mean AUC over defined queries.

- `label_av`: a named list of per-label means in first-appearance order.

With `ret_extra = TRUE`, the list also contains `query_auc`, with
`NA_real_` for undefined queries, and the integer `n_defined`.

## Details

`av_auc` is the mean ROC AUC over defined queries; `label_av` gives
per-label means in first-appearance order. Queries without both relevant
and irrelevant observations are excluded. An empty average is
`NA_real_`.

## Note

Requires the `PRROC` package.
