# Average Area Under the ROC Curve

Embedding quality measure.

## Usage

``` r
roc_auc(dm, labels)
```

## Arguments

- dm:

  Distance matrix of an embedding.

- labels:

  Vector of labels for each observation in the dataset in the same order
  as the observations in the distance matrix.

## Value

Area Under the ROC curve, averaged over each observation.

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

Perfect retrieval results in an AUC of 1. For random retrieval gives a
value of 0.5.

## Note

Use of this function requires that the `PRROC` package be installed.
