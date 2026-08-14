# Average Area Under the ROC Curve

Embedding quality measure.

## Usage

``` r
roc_auc(dm, labels)
```

## Arguments

- dm:

  A finite square numeric distance matrix. The diagonal is excluded from
  each retrieval problem.

- labels:

  A label vector without missing values, with one label per row of `dm`.

## Value

A list with components in this order:

- `av_auc`: the AUC averaged over defined observation rows.

- `label_av`: a named list of per-label means in first-appearance order.

## Details

The ROC curve plots the true positive rate vs false positive rate. This
function calculates the curve N times, where N is the number of the
observations. The label of the Nth observation is set as the positive
class and then the other observations are ranked according to their
distance from the Nth observation in the output coordinates (lower
distances being better). Observations with the same label as the Nth
observation count as positive observations. The final reported result is
the average over all observations. Rows with undefined AUC values, such
as rows that cannot define both positive and negative examples, are
excluded from overall and per-label averages. If no rows remain for an
average, that average is `NA_real_`.

`av_auc` weights each defined query equally; `label_av` reports
per-label means in first-appearance order.

Perfect retrieval results in an AUC of 1. For random retrieval gives a
value of 0.5.

## Note

Use of this function requires that the `PRROC` package be installed.
