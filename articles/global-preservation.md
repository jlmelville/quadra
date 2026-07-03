# Global preservation

## Random Triplet Accuracy

Random triplet accuracy is the proportion of triangle distances where
the relative ordering is retained. It was used by Wang et al. (2021) to
measure global structure preservation.

Triplets with tied input-space distances are excluded from the
denominator because they do not define a relative ordering. If no
sampled input triplets define an ordering, the result is `NA_real_`.

## Pearson Correlation for Distances

[`random_pair_distance_correlation()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_correlation.md)
measures Pearson correlation between equivalent input and output
distances, similar to the method used by Becht and co-workers.

## Earth-Mover’s Distance

[`random_pair_distance_emd()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_emd.md)
converts distances to an empirical distribution and compares them via
Earth Mover’s Distance (1D Wasserstein), based on the method used by
Heiser and Lau.

## Further Reading

### Random Triplet Accuracy

Wang, Y., Huang, H., Rudin, C., & Shaposhnik, Y. (2021). Understanding
how dimension reduction tools work: an empirical approach to deciphering
t-SNE, UMAP, TriMAP, and PaCMAP for data visualization. *J Mach. Learn.
Res*, *22*, 1-73. <https://jmlr.org/papers/v22/20-1061.html>

### Pearson Correlation for Distances

Becht, E., McInnes, L., Healy, J., Dutertre, C. A., Kwok, I. W., Ng, L.
G., … & Newell, E. W. (2019). Dimensionality reduction for visualizing
single-cell data using UMAP. *Nature biotechnology*, *37*(1), 38-44.
<https://doi.org/10.1038/nbt.4314>

### Earth-Mover’s Distance

Heiser, C. N., & Lau, K. S. (2020). A quantitative framework for
evaluating single-cell data structure preservation by dimensionality
reduction techniques. *Cell reports*, *31*(5), 107576.
<https://doi.org/10.1016/j.celrep.2020.107576>
