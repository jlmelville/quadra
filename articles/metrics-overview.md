# Metric overview

This package provides several ways to evaluate the performance of an
embedding. They fall into four broad families.

## Choosing a Starting Metric

| User question | Reasonable starting metric | Decision supported by the result |
|----|----|----|
| Are the same nearby observations retained? | [`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md) | Whether local-neighborhood identity is preserved well enough for the intended use. |
| Do dense and sparse regions keep their relative local scale? | [`local_radius_correlation()`](https://jlmelville.github.io/quadra/reference/local_radius_correlation.md) | Whether scale distortion needs attention even when neighbor identities look acceptable. |
| Are broad distance orderings retained? | [`random_triplet_accuracy()`](https://jlmelville.github.io/quadra/reference/random_triplet_accuracy.md) | Whether the embedding preserves enough global ordering for the intended use. |
| Do known labels remain retrievable? | [`roc_auc()`](https://jlmelville.github.io/quadra/reference/roc_auc.md) or [`pr_auc()`](https://jlmelville.github.io/quadra/reference/pr_auc.md) | Whether same-label observations rank ahead of other labels in the embedding. |

Choose metrics to match the relationships an embedding should preserve.
The table gives starting points.

## Local Neighborhood Preservation

[`nbr_pres()`](https://jlmelville.github.io/quadra/reference/nbr_pres.md)
compares neighborhoods from exact distance matrices.
[`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md)
and
[`nbr_pres_knn()`](https://jlmelville.github.io/quadra/reference/nbr_pres_knn.md)
work with nearest-neighbor graphs and scale to larger datasets.

[`trustworthiness()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)
and
[`continuity()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)
are also local-neighborhood metrics, but they use rank penalties rather
than only counting shared neighbors. They are exact distance-matrix
functions for small datasets.
[`rnx_auc()`](https://jlmelville.github.io/quadra/reference/rnx_auc.md)
summarizes rank-based neighborhood agreement over many neighborhood
sizes, while
[`rnx_curve()`](https://jlmelville.github.io/quadra/reference/rnx_curve.md)
retains the individual scales.

See the [local
preservation](https://jlmelville.github.io/quadra/articles/local-preservation.md)
article for examples and related local diagnostics.

## Local Scale Preservation

Neighbor identity does not say whether dense and sparse regions keep
their relative scale.
[`local_radius_correlation()`](https://jlmelville.github.io/quadra/reference/local_radius_correlation.md)
compares the local radius around each observation in the input data and
output embedding, using either the distance to the `k`th neighbor or the
mean distance to the first `k` neighbors.

## Local Graph Diagnostics

[`mutual_neighbor_correlation()`](https://jlmelville.github.io/quadra/reference/mutual_neighbor_correlation.md)
compares the mutual-neighbor count pattern in the input and output
nearest-neighbor graphs. It asks whether points with many or few
reciprocal nearest-neighbor relationships in the input graph have
similar counts in the output graph.

## Global Distance Preservation

`quadra` also contains random triplet and random pair distance methods
for global structure preservation:

- [`random_triplet_accuracy()`](https://jlmelville.github.io/quadra/reference/random_triplet_accuracy.md)
  checks whether sampled distance orderings are preserved.
- [`random_pair_distance_correlation()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_correlation.md)
  measures Pearson or Spearman correlation between sampled input and
  output distances.
- [`random_pair_distance_emd()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_emd.md)
  compares sampled distance distributions using Earth Mover’s Distance.
- [`random_pair_distance_stress()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_stress.md)
  compares matched sampled distances with a root mean squared
  difference.

For these metrics, `Xin` and `metric_in` define the reference geometry.

See the [global
preservation](https://jlmelville.github.io/quadra/articles/global-preservation.md)
article for those metrics.

## Input and Ordering Conventions

Sampled pair and triplet metrics require dense, finite inputs; graph
metrics accept the sparse inputs supported by `rnndescent`.

[`nbr_pres()`](https://jlmelville.github.io/quadra/reference/nbr_pres.md)
includes boundary ties, exact-rank metrics break ties by column order,
and supplied-graph metrics use their stored order.

Euclidean and squared Euclidean induce the same ordering but different
magnitudes. This distinction affects Pearson correlation, EMD, stress,
and local-scale summaries, but not neighbor, triplet, or Spearman
ordering.

## Label Retrieval

[`roc_auc()`](https://jlmelville.github.io/quadra/reference/roc_auc.md)
and
[`pr_auc()`](https://jlmelville.github.io/quadra/reference/pr_auc.md)
rank other observations by distance from each query and treat matching
labels as relevant. They require the [PRROC
package](https://cran.r-project.org/package=PRROC).

See the [label
retrieval](https://jlmelville.github.io/quadra/articles/label-retrieval.md)
article for examples and return value details.
