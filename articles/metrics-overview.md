# Metric overview

This package provides several ways to evaluate the performance of an
embedding. They fall into four broad families.

## Local Neighborhood Preservation

The most generic approach is to consider “neighborhood preservation”:

1.  Find the *k*-nearest neighbors of a point in the input space.
2.  Do the same in the output space.
3.  Count the overlap.

[`nbr_pres()`](https://jlmelville.github.io/quadra/reference/nbr_pres.md)
does this from exact distance matrices.
[`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md)
and
[`nbr_pres_knn()`](https://jlmelville.github.io/quadra/reference/nbr_pres_knn.md)
work with nearest-neighbor graphs, which is the more practical route for
larger datasets.

[`trustworthiness()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)
and
[`continuity()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)
are also local-neighborhood metrics, but they use rank penalties rather
than only counting shared neighbors. They are exact distance-matrix
functions for small datasets.
[`rnx_auc()`](https://jlmelville.github.io/quadra/reference/rnx_auc.md)
summarizes rank-based neighborhood agreement over many neighborhood
sizes.

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

See the global preservation article for those metrics.

## Label Retrieval

Alternatively, if some sort of labelling is applied to the points, each
point can be treated as the target in a retrieval procedure:

1.  Rank all the other points by distance to the target point.
2.  See how highly in the ranked list the points with the same label are
    found.
3.  Construct a Receiver Operating Characteristic (ROC) curve, or
    something like it (e.g. Precision-Recall curve)
4.  Calculate the Area Under the Curve (AUC).
5.  Average over all points.

This only needs the output distance matrix, but requires the sort of
labeling usually reserved for data intended for supervised
classification. Quadra can also provide some help with this, but
requires the [PRROC package](https://cran.r-project.org/package=PRROC)
to be installed.
