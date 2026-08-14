# Area Under the Precision-Recall Curve for an Embedding

Embedding quality measure.

## Usage

``` r
pr_auc(dm, labels)
```

## Arguments

- dm:

  A finite square numeric distance matrix. The diagonal is excluded from
  each retrieval problem.

- labels:

  A label vector without missing values, with one label per row of `dm`.

## Value

A list with components in this order:

- `av_auc`: the PR AUC averaged over defined observation rows.

- `label_av`: a named list of per-label means in first-appearance order.

## Details

The PR curve plots precision (also known as positive predictive value,
PPV) against recall (also known as the true positive rate). The area
under the curve provides similar information compared to the area under
the ROC curve, but may be more appropriate when classes are highly
imbalanced.

This function calculates the PR curve N times, where N is the number of
the observations. The label of the Nth observation is set as the
positive class and then the other observations are ranked according to
their distance from the Nth observation in the output coordinates (lower
distances being better). Observations with the same label as the Nth
observation count as positive observations. The final reported result is
the average over all observations. Rows with undefined AUC values, such
as rows that cannot define both positive and negative examples, are
excluded from overall and per-label averages. If no rows remain for an
average, that average is `NA_real_`.

`av_auc` weights each defined query equally; `label_av` reports
per-label means in first-appearance order.

Perfect retrieval results in an AUC of 1. For random retrieval, the
value is the proportion of the positive class labels for that curve.

## Note

Use of this function requires that the `PRROC` package be installed.

## References

Keilwagen, J., Grosse, I., & Grau, J. (2014). Area under
precision-recall curves for weighted and unweighted data. *PloS One*,
*9*(3), e92209.

Davis, J., & Goadrich, M. (2006, June). The relationship between
Precision-Recall and ROC curves. In *Proceedings of the 23rd
international conference on Machine learning* (pp. 233-240). ACM.
