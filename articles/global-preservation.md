# Global preservation

Use this guide to choose a sampled global-distance metric and compare
embeddings on the same sampled relationships. For every metric here,
`Xin` and `metric_in` define the reference geometry.

| Question | Metric | Better direction | What is compared? |
|----|----|----|----|
| Are anchored distance orderings retained? | [`random_triplet_accuracy()`](https://jlmelville.github.io/quadra/reference/random_triplet_accuracy.md) | Higher; 1 is perfect | Which of two endpoints is closer to the same anchor |
| Are sampled distances associated? | [`random_pair_distance_correlation()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_correlation.md) | Higher; 1 is perfect positive association | Matched input/output pairs, linearly (Pearson) or by rank (Spearman) |
| Do sampled distance distributions have similar shapes? | [`random_pair_distance_emd()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_emd.md) | Lower; 0 is perfect | Marginal distributions, without preserving pair correspondence |
| Are matched sampled distance magnitudes close? | [`random_pair_distance_stress()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_stress.md) | Lower; 0 is perfect | The same input/output pairs after optional range scaling |

Create reusable samples with
[`sample_pairs()`](https://jlmelville.github.io/quadra/reference/sample_pairs.md)
or
[`sample_triplets()`](https://jlmelville.github.io/quadra/reference/sample_triplets.md)
and supply them with `pairs` or `triplets`. Set `ret_extra = TRUE` to
inspect the sampled comparisons behind a score.

For example, reuse one pair sample to compare two embeddings on exactly
the same input-space relationships:

``` r

library(quadra)

iris_x <- as.matrix(iris[, -5])
pca_unscaled <- stats::prcomp(iris_x, rank. = 2, scale. = FALSE)$x
pca_scaled <- stats::prcomp(iris_x, rank. = 2, scale. = TRUE)$x
```

``` r

set.seed(42)
pair_plan <- sample_pairs(nrow(iris_x), n_pairs = 1000)
distance_correlation <- c(
  pca_unscaled = random_pair_distance_correlation(
    iris_x, pca_unscaled, pairs = pair_plan, n_threads = 1
  ),
  pca_scaled = random_pair_distance_correlation(
    iris_x, pca_scaled, pairs = pair_plan, n_threads = 1
  )
)
signif(distance_correlation, 5)
```

    ## pca_unscaled   pca_scaled 
    ##      0.99967      0.94536

On this fixed Iris pair sample, the unscaled PCA has the higher Pearson
correlation, so its sampled output distances are more linearly
associated with the corresponding input distances. This is evidence
about these embeddings and this sample, not a universal ranking of
scaled and unscaled PCA.

## Random Triplet Accuracy

Random triplet accuracy is the proportion of sampled anchored triplets
whose input-space distance ordering is retained in the embedding.

Triplets with tied input-space distances are excluded from the
denominator because they do not define a relative ordering. If no
sampled input triplets define an ordering, the result is `NA_real_`.

Detailed results contain the evaluated triplets and a row-aligned
agreement vector; input-distance ties are `NA`. Reset the R seed and
keep `n_threads` fixed to repeat implicit sampling.

## Random Pair Distance Correlation

[`random_pair_distance_correlation()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_correlation.md)
measures Pearson correlation between matched input and output distances
by default. Set `method = "spearman"` to compare their ranks.

Euclidean and squared Euclidean have the same ordering, so their
Spearman interpretation is equivalent. Squaring changes magnitudes and
can change the default Pearson result.

## Earth-Mover’s Distance

[`random_pair_distance_emd()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_emd.md)
compares the empirical input and output distance distributions using
Earth Mover’s Distance (1D Wasserstein).

EMD compares marginal distributions rather than pairwise correspondence.

EMD and stress below are magnitude statistics. Squared Euclidean can
change both results relative to Euclidean distance, even when each
sampled distance vector is range-scaled.

## Random Pair Distance Stress

[`random_pair_distance_stress()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_stress.md)
compares the matched sampled distances with a root mean squared
difference, after scaling each distance vector to `[0, 1]` by default.

RNX is sometimes described as spanning local and global scales, but it
measures rank-based neighborhood agreement rather than sampled global
distances. See the [local
preservation](https://jlmelville.github.io/quadra/articles/local-preservation.html#multi-scale-rnx-agreement)
guide for its AUC, scale-specific curve, and interpretation.

## Further Reading

### Random Triplet Accuracy

Wang, Y., Huang, H., Rudin, C., & Shaposhnik, Y. (2021). Understanding
how dimension reduction tools work: an empirical approach to deciphering
t-SNE, UMAP, TriMAP, and PaCMAP for data visualization. *J Mach. Learn.
Res*, *22*, 1-73. <https://jmlr.org/papers/v22/20-1061.html>

### Random Pair Distance Correlation

Becht, E., McInnes, L., Healy, J., Dutertre, C. A., Kwok, I. W., Ng, L.
G., … & Newell, E. W. (2019). Dimensionality reduction for visualizing
single-cell data using UMAP. *Nature biotechnology*, *37*(1), 38-44.
<https://doi.org/10.1038/nbt.4314>

### Earth-Mover’s Distance

Heiser, C. N., & Lau, K. S. (2020). A quantitative framework for
evaluating single-cell data structure preservation by dimensionality
reduction techniques. *Cell reports*, *31*(5), 107576.
<https://doi.org/10.1016/j.celrep.2020.107576>
