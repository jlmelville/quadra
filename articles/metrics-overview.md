# Metric overview

An embedding can put the right digit next to yours and still put two
whole groups much too far apart. Which of those would matter for the way
you intend to use the map?

## Choosing a Starting Metric

| What do you want to know? | Where to start |
|----|----|
| Who stayed nearby? | [`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md) counts how many of the original neighbors each point keeps. |
| Do dense regions stay dense? | [`local_radius_correlation()`](https://jlmelville.github.io/quadra/reference/local_radius_correlation.md) checks whether points with small neighborhood radii in the input still have small radii on the map. |
| Do distances keep their order? | [`random_triplet_accuracy()`](https://jlmelville.github.io/quadra/reference/random_triplet_accuracy.md) takes a point and two others. How often does the closer one stay closer? |
| Can I retrieve a label? | [`roc_auc()`](https://jlmelville.github.io/quadra/reference/roc_auc.md) and [`pr_auc()`](https://jlmelville.github.io/quadra/reference/pr_auc.md) tell you how well images of the same digit, for example, rank ahead of the other digits. |
| Did a patch change shape? | [`local_procrustes()`](https://jlmelville.github.io/quadra/reference/local_procrustes.md) aligns the same points in both spaces. See what happened to [one patch of the Swiss roll](https://jlmelville.github.io/quadra/articles/reference-geometry.html#compare-local-shape-after-alignment). |
| Can a few points dominate the map? | Compare coordinate spread with fitted-scale stress. See [the two digits that set the scale of an LTSA map](https://jlmelville.github.io/quadra/articles/ltsa-map-spread.md). |
| Did a gap open up? | [`gap_index()`](https://jlmelville.github.io/quadra/reference/gap_index.md) compares triangle areas. See the [gap between three digits](https://jlmelville.github.io/quadra/articles/mnist-pca-umap.html#inspect-the-gaps-between-digits). |
| Can I trust an island? | [`steadiness_cohesiveness()`](https://jlmelville.github.io/quadra/reference/steadiness_cohesiveness.md) checks groups [in both directions](https://jlmelville.github.io/quadra/articles/mnist-pca-umap.html#check-groups-without-digit-labels), without using the digit labels. |
| Did two classes become easier to tell apart? | Compare [`class_pair_auc()`](https://jlmelville.github.io/quadra/reference/class_pair_auc.md) with the fraction of points closer to their own class centroid. To check angular spread, look at [the 5s as seen from the 3s](https://jlmelville.github.io/quadra/articles/mnist-pca-umap.html#compare-angular-spread-from-each-digit). |

## Which metrics work on large datasets?

Around **5,000 observations**, check the memory budget before
constructing a distance matrix. One dense `N` by `N` double matrix takes
`8 * N^2` bytes for its entries. At 5,000 rows, budget 191 MiB for one
space, 381.5 MiB for a reference and one embedding, or 572 MiB for a
reference and two embeddings. A single matrix grows from 7.6 MiB at
1,000 rows to 763 MiB at 10,000 and 26.8 GiB at 60,000. Allow additional
memory for coordinates, copies and working storage. The practical limit
depends on which metric you run and how much memory is available.

- **Distances and their ordering.** Use the full coordinates and reuse a
  sample of pairs across embeddings for correlation, EMD, stress or a
  distance-band profile. Check how many pairs reach the narrowest band.
  For
  [`random_triplet_accuracy()`](https://jlmelville.github.io/quadra/reference/random_triplet_accuracy.md),
  the sample size is *per observation*. [Global distance
  metrics](#global-distance-preservation)
- **Neighbors.** A graph with `k` neighbors per row stores `N * k`
  indices. Keep the graphs from `nn_preservation(..., ret_extra = TRUE)`
  for overlap, local radii (which also need distances) and reciprocal
  neighbors. Brute-force search still compares every pair. For
  approximate search, [spot-check a few
  queries](https://jlmelville.github.io/quadra/articles/mnist-pca-umap.html#spot-check-the-pixel-neighbor-search).
- **Exact ranks.** With Euclidean coordinates,
  [`exact_rank_profiles()`](https://jlmelville.github.io/quadra/reference/exact_rank_profiles.md)
  can score a sample of queries against every candidate. It needs only
  one query’s ranks at a time, but scoring *all* queries still requires
  quadratic distance work. [Try a query
  sample](#sample-queries-when-the-full-neighbor-pool-matters).
- **Label retrieval.** Give
  [`roc_auc()`](https://jlmelville.github.io/quadra/reference/roc_auc.md)
  or
  [`pr_auc()`](https://jlmelville.github.io/quadra/reference/pr_auc.md)
  the map coordinates in `coords` and a sample of `query_rows`. Every
  other observation stays in the ranking. Use the same queries in each
  map. [Compare UMAP and t-SNE on the same
  digits](https://jlmelville.github.io/quadra/articles/mnist-pca-umap.html#umap-versus-tsne).
- **Class geometry.**
  [`label_trustworthiness_continuity()`](https://jlmelville.github.io/quadra/reference/label_trustworthiness_continuity.md)
  checks nearest-centroid consistency for every pair of classes.
  [`class_angular_distortion_index()`](https://jlmelville.github.io/quadra/reference/class_angular_distortion_index.md)
  estimates angular distortion from sampled triplets. Here `n_triplets`
  is the *total* number to draw. Both work on coordinates.
- **Patches and gaps.** Reuse a reference graph for
  [`local_procrustes()`](https://jlmelville.github.io/quadra/reference/local_procrustes.md)
  and fit one small alignment per point. For
  [`gap_index()`](https://jlmelville.github.io/quadra/reference/gap_index.md),
  only the edges of the 2D triangulation need distances.
- **Shared neighbors.**
  [`neighbor_dissimilarity()`](https://jlmelville.github.io/quadra/reference/neighbor_dissimilarity.md)
  starts with sparse graphs. One popular neighbor can occur in many
  lists, and every pair of those lists shares that neighbor. Even the
  scalar score can then require space for nearly every pair of
  observations.

**When to use a smaller dataset.**
[`nbr_pres()`](https://jlmelville.github.io/quadra/reference/nbr_pres.md),
[`trustworthiness()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md),
[`continuity()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md),
[`rnx_auc()`](https://jlmelville.github.io/quadra/reference/rnx_auc.md),
[`rnx_curve()`](https://jlmelville.github.io/quadra/reference/rnx_curve.md),
[`pointwise_rank_preservation()`](https://jlmelville.github.io/quadra/reference/pointwise_rank_preservation.md)
and
[`pointwise_rnx()`](https://jlmelville.github.io/quadra/reference/pointwise_rnx.md)
each take **two dense distance matrices**.
[`class_pair_auc()`](https://jlmelville.github.io/quadra/reference/class_pair_auc.md)
takes **one per space**. If you need class-pair AUC for just two
classes, select those classes before computing distances. If the full
matrices exceed your budget, use the same representative subset in every
space and report the candidate pool size.

For
[`steadiness_cohesiveness()`](https://jlmelville.github.io/quadra/reference/steadiness_cohesiveness.md),
start with a representative subset. It takes coordinates, then builds
dense distance and similarity matrices and repeatedly clusters sampled
groups. [Jeon and colleagues introduced the
measures](https://doi.org/10.1109/TVCG.2021.3114833). Quadra follows
[ZADU 0.5.4](https://github.com/hj-n/zadu) for shared-neighbor and map
weights. Asking for pointwise scores adds more square matrices. Fewer
iterations save CPU while the dense allocations remain. The MNIST
workflow ran exact coordinate ranks for all 5,000 images and used 1,000
for the group comparison. [MNIST group
audit](https://jlmelville.github.io/quadra/articles/mnist-pca-umap.html#check-groups-without-digit-labels)

### Which scores would I use to compare UMAP or t-SNE fits?

For a sweep over initializations or parameter settings, I’d start with
two scores on every fit.

1.  **Neighbor overlap.** I use
    [`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md)
    at two fixed sizes, say `k = c(15, 50)`. A map can keep the closest
    few images and lose its lead farther down the list. Build the input
    graph once for the larger size and reuse it for every fit, even when
    varying UMAP’s `n_neighbors` or t-SNE’s perplexity. Build a fresh
    graph for each map. At full MNIST size I’d supply an exact 2D graph
    from a k-d tree or use `nn_method_out = "nnd"` and spot-check the
    approximate search. The default output search is brute force.
2.  **Random triplet accuracy.** Two maps can keep the same neighbors
    and put whole groups in different places. Given an anchor and two
    other images, which image was closer to the anchor in the input, and
    is it still closer on the map? I’d start with five triplets per
    observation. For 60,000 rows that’s 300,000 comparisons. Draw them
    once with
    [`sample_triplets()`](https://jlmelville.github.io/quadra/reference/sample_triplets.md)
    and pass the same `triplets` to
    [`random_triplet_accuracy()`](https://jlmelville.github.io/quadra/reference/random_triplet_accuracy.md)
    for every fit. This calculation needs dense coordinates. If the
    input is too large to store densely, I’d stick with graph overlap.

For an initialization comparison, keep the input representation,
distance metric and optimizer settings fixed. For a parameter
comparison, vary the parameter of interest. Run stochastic settings over
the same set of seeds and compare the paired scores. I’d want to see a
difference hold up across seeds before choosing a setting on its
strength.

If finding other images of the same digit is the goal, I’d also report
**same-digit precision at 1, 5 and 15** from the map’s neighbor lists.
That answers a different question from recovering the same individual
images. In the MNIST comparison, t-SNE’s same-digit lead shrank at each
cutoff in all 12 pairings. [MNIST fit
comparison](https://jlmelville.github.io/quadra/articles/mnist-pca-umap.html#does-the-early-lead-depend-on-the-settings)

For a closer look at the finalists, I’d run
[`exact_rank_profiles()`](https://jlmelville.github.io/quadra/reference/exact_rank_profiles.md)
on the same few hundred queries against the full candidate pool.
Trustworthiness penalizes distant points that turned up nearby on the
map. Continuity penalizes close neighbors that moved far away. With
labels, pass the map coordinates and the same `query_rows` to
[`pr_auc()`](https://jlmelville.github.io/quadra/reference/pr_auc.md) to
check the whole ranking.

### Sample queries when the full neighbor pool matters

To estimate rank preservation for the full dataset, sample the
**queries** and keep every row of the reference and embeddings. With `x`
and `maps` from the MNIST workflow, these 500 queries each rank the
other 4,999 images. The mean gives each sampled query equal weight.

``` r

set.seed(42)
queries <- sample.int(nrow(x), 500)
profiles <- exact_rank_profiles(x, maps, k = 15, query_rows = queries)
vapply(profiles, function(p) mean(p$trustworthiness[, 1]), numeric(1))
```

Subsetting the coordinate matrices also removes candidates and changes
what the ranks mean. On the same MNIST queries at `k = 100`, 7s retained
neighbors better than 0s in the 1,000-image pool, while 0s did better in
the 5,000-image pool. At `k = 500` in the larger pool, 7s led again.
Report the number of queries, the size of the candidate pool and `k`.
[Full-pool MNIST
comparison](https://jlmelville.github.io/quadra/articles/mnist-pca-umap.html#expand-the-candidate-pool-to-all-5000-images)

## Local Neighborhood Preservation

[`nbr_pres()`](https://jlmelville.github.io/quadra/reference/nbr_pres.md)
compares neighborhoods from exact distance matrices.
[`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md)
and
[`nbr_pres_knn()`](https://jlmelville.github.io/quadra/reference/nbr_pres_knn.md)
compare nearest-neighbor graphs. All three count the fraction of
neighbors retained, as in [Chen and
Buja](https://doi.org/10.1198/jasa.2009.0111).

[`trustworthiness()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)
and
[`continuity()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)
score misplaced neighbors by their ranks.
[`rnx_auc()`](https://jlmelville.github.io/quadra/reference/rnx_auc.md)
summarizes rank-based neighborhood agreement over many neighborhood
sizes, while
[`rnx_curve()`](https://jlmelville.github.io/quadra/reference/rnx_curve.md)
retains the individual scales.
[`pointwise_rnx()`](https://jlmelville.github.io/quadra/reference/pointwise_rnx.md)
gives the comparison per observation. The rank penalties are from [Venna
and Kaski](https://doi.org/10.1007/3-540-44668-0_68). [Lee and
colleagues](https://doi.org/10.1016/j.neucom.2012.12.036) describe the
RNX curve and its [weighted
summary](https://doi.org/10.1016/j.neucom.2014.12.095).

To see how far the closest neighbors move in the ranking, use
[`exact_rank_profiles()`](https://jlmelville.github.io/quadra/reference/exact_rank_profiles.md)
with `mrre = TRUE`. A first neighbor that slips to fifth counts more
than a tenth neighbor that slips to fourteenth. You can follow the input
neighbors into the embedding or trace the embedding’s neighbors back to
the input. [Lee and Verleysen’s rank
errors](https://proceedings.mlr.press/v4/lee08a.html) [Worked
example](https://jlmelville.github.io/quadra/articles/local-preservation.html#can-the-same-neighbors-come-back-in-a-different-order)

See the [local
preservation](https://jlmelville.github.io/quadra/articles/local-preservation.md)
article for examples and related local diagnostics.

## Local Scale Preservation

Neighbor identity does not say whether dense and sparse regions keep
their relative scale.
[`local_radius_correlation()`](https://jlmelville.github.io/quadra/reference/local_radius_correlation.md)
compares the local radius around each observation in the input data and
output embedding, using either the distance to the `k`th neighbor or the
mean distance to the first `k` neighbors. [Narayan, Berger and
Cho](https://doi.org/10.1038/s41587-020-00801-7) use local radii to
assess density preservation in den-SNE and densMAP.

## Local Shape Preservation

[`local_procrustes()`](https://jlmelville.github.io/quadra/reference/local_procrustes.md)
checks how much the same reference neighborhood changes shape after
rigid alignment. It can find distortions even when all members stay
together. The matched patches follow [Goldberg and
Ritov](https://arxiv.org/abs/0806.2669). Quadra uses a fixed-scale
residual.

## Regional Area Distortion

[`gap_index()`](https://jlmelville.github.io/quadra/reference/gap_index.md)
locates expansion and compression in the empty triangles of a 2D
embedding. [Ros, Arleo and Paulovich introduced the
index](https://arxiv.org/abs/2607.28324). Their
[implementation](https://codeberg.org/jros/gap-index) is available.
[MNIST
example](https://jlmelville.github.io/quadra/articles/mnist-pca-umap.html#inspect-the-gaps-between-digits)

## Local Graph Diagnostics

[`mutual_neighbor_correlation()`](https://jlmelville.github.io/quadra/reference/mutual_neighbor_correlation.md)
compares the mutual-neighbor count pattern in the input and output
nearest-neighbor graphs. It asks whether points with many or few
reciprocal nearest-neighbor relationships in the input graph have
similar counts in the output graph. The comparison draws on
[drnb](https://github.com/jlmelville/drnb).

Use
[`neighbor_dissimilarity()`](https://jlmelville.github.io/quadra/reference/neighbor_dissimilarity.md)
when the question is which pairs of observations gain or lose shared
neighbors. [Fujiwara and colleagues](https://arxiv.org/abs/2206.13891)
define the measure. [Graph
comparison](https://jlmelville.github.io/quadra/articles/local-preservation.html#which-pairs-share-the-same-neighbors)

## Global Distance Preservation

`quadra` also contains random triplet and random pair distance methods
for global structure preservation:

- [`random_triplet_accuracy()`](https://jlmelville.github.io/quadra/reference/random_triplet_accuracy.md)
  checks whether sampled distance orderings are preserved, as in [Wang
  and colleagues](https://jmlr.org/papers/v22/20-1061.html).
- [`random_pair_distance_correlation()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_correlation.md)
  measures Pearson or Spearman correlation between sampled input and
  output distances. See [Becht and
  colleagues](https://doi.org/10.1038/nbt.4314) and
  [drnb](https://github.com/jlmelville/drnb).
- [`random_pair_distance_emd()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_emd.md)
  compares sampled distance distributions using Earth Mover’s Distance,
  following [Heiser and
  Lau](https://doi.org/10.1016/j.celrep.2020.107576) and their
  [code](https://github.com/KenLauLab/DR-structure-preservation).
- [`random_pair_distance_stress()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_stress.md)
  compares matched sampled distances with a root mean squared
  difference. This comparison draws on
  [drnb](https://github.com/jlmelville/drnb).

For these metrics, `Xin` and `metric_in` define the reference geometry.

The [identity-permutation
example](https://jlmelville.github.io/quadra/articles/global-preservation.html#identical-scatterplots-different-point-identities)
shows which metrics notice when the same plotted positions belong to
different observations. Use [distance
bands](https://jlmelville.github.io/quadra/articles/distance-band-fidelity.md)
to compare ordering and fitted-scale stress from nearby input pairs
through the entire pair sample. The rank and fitted-scale stress
comparisons follow [Takaeda’s
benchmark](https://arxiv.org/abs/2608.21779). Quadra samples pairs by
default.

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

For each class pair,
[`label_trustworthiness_continuity()`](https://jlmelville.github.io/quadra/reference/label_trustworthiness_continuity.md)
asks how many observations fall closer to their own class centroid and
how that fraction changes in the embedding.
[`class_angular_distortion_index()`](https://jlmelville.github.io/quadra/reference/class_angular_distortion_index.md)
measures changes in the angles that one class subtends from points in
another. These measures come from [Jeon and
colleagues](https://arxiv.org/abs/2308.00278) and [Gunaratne, Kobourov
and Miller](https://arxiv.org/abs/2605.00637), respectively. For the
former, Quadra uses the unrescaled distance consistency of [Sips and
colleagues](https://doi.org/10.1111/j.1467-8659.2009.01467.x).

See the [label
retrieval](https://jlmelville.github.io/quadra/articles/label-retrieval.md)
article for examples and return value details.
