# Label Trustworthiness and Continuity

For each pair of classes, what fraction of observations lie closer to
their own class centroid? Compare that fraction in the reference and the
embedding to see which classes have become harder or easier to tell
apart.

## Usage

``` r
label_trustworthiness_continuity(Xin, Xout, labels, ret_extra = FALSE)
```

## Arguments

- Xin:

  Reference coordinates, as a dense numeric matrix or data frame with
  observations in rows.

- Xout:

  Embedding coordinates for the same observations in the same order.
  Nonnumeric data-frame columns are ignored in both inputs.

- labels:

  A vector without missing values, with one label per row and at least
  two distinct labels.

- ret_extra:

  Whether to return the class-pair consistency scores and their signed
  changes.

## Value

A named numeric vector with `trustworthiness` and `continuity`. With
`ret_extra = TRUE`, a list containing `scores` and `pairs`. The data
frame `pairs` has columns `label_a`, `label_b`, `consistency_in`,
`consistency_out`, and `change` (output minus input). Pairs follow the
first appearance of the labels.

## Details

For each pair of labels and in each space, compute the two class
centroids. Distance consistency is the fraction of observations nearer
their own centroid in Euclidean distance. A tie counts as half a correct
assignment. This implementation uses the unrescaled fraction in `[0, 1]`
and weights observations equally within each pair. Centroids are
estimated on the observations being evaluated.

Label trustworthiness is one minus the mean decrease in consistency,
counting gains as zero. Label continuity is one minus the mean increase,
counting losses as zero. Each unordered pair has equal weight in these
means. A pair can gain consistency when the embedding exaggerates the
separation visible in the reference.

Distance consistency uses one centroid per class. Compare it with
[`roc_auc()`](https://jlmelville.github.io/quadra/reference/roc_auc.md)
or [`pr_auc()`](https://jlmelville.github.io/quadra/reference/pr_auc.md)
when classes have multiple clusters or irregular shapes. Nonfinite
computed centroid distances cause an error.

## Scaling

This is a full-data option when the number of classes is modest. Each
row is compared with its own centroid and each other class centroid.
With `C` classes the summaries cover `C * (C - 1) / 2` pairs. No
observation-level distance matrix is required.

## References

Jeon, H., Kuo, Y.-H., Aupetit, M., Ma, K.-L., & Seo, J. (2024). Classes
are not Clusters: Improving Label-based Evaluation of Dimensionality
Reduction.
[doi:10.1109/TVCG.2023.3327187](https://doi.org/10.1109/TVCG.2023.3327187)
Their [ltnc package](https://github.com/hj-n/ltnc) compares classes in
pairs. Quadra uses the fraction closer to their own centroid for each
pair, following the distance consistency criterion of Sips, Neubert,
Lewis and Hanrahan (2009).
[doi:10.1111/j.1467-8659.2009.01467.x](https://doi.org/10.1111/j.1467-8659.2009.01467.x)

## See also

[`roc_auc()`](https://jlmelville.github.io/quadra/reference/roc_auc.md),
[`pr_auc()`](https://jlmelville.github.io/quadra/reference/pr_auc.md)

## Examples

``` r
x <- as.matrix(iris[, -5])
y <- prcomp(x, rank. = 2)$x
label_trustworthiness_continuity(x, y, iris$Species, ret_extra = TRUE)
#> $scores
#> trustworthiness      continuity 
#>               1               1 
#> 
#> $pairs
#>      label_a    label_b consistency_in consistency_out change
#> 1     setosa versicolor           1.00            1.00      0
#> 2     setosa  virginica           1.00            1.00      0
#> 3 versicolor  virginica           0.89            0.89      0
#> 
```
