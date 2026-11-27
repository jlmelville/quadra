# Changelog

## quadra 0.3.0

### New features and requirements

- Require rnndescent 0.3.0 or later. Its compiled cosine and correlation
  distances preserve the expected values when otherwise ordinary
  coordinates are multiplied by `1e100`.
- Use `coords` and `query_rows` in
  [`pr_auc()`](https://jlmelville.github.io/quadra/reference/pr_auc.md)
  or
  [`roc_auc()`](https://jlmelville.github.io/quadra/reference/roc_auc.md)
  to score selected observations against the full candidate pool without
  building a square distance matrix.
- Read the number of shared neighbors at each requested size from the
  `retained` matrix in each embedding’s
  [`exact_rank_profiles()`](https://jlmelville.github.io/quadra/reference/exact_rank_profiles.md)
  result.
- Add
  [`neighbor_dissimilarity()`](https://jlmelville.github.io/quadra/reference/neighbor_dissimilarity.md)
  to compare shared neighbors in two graphs and identify the observation
  pairs that gain or lose them.
- Add
  [`exact_rank_profiles()`](https://jlmelville.github.io/quadra/reference/exact_rank_profiles.md)
  to compare Euclidean neighborhood ranks across embeddings without
  storing a full distance matrix. Use `query_rows` to check selected
  observations against the full candidate pool. Use `mrre = TRUE` to see
  how far the closest neighbors move in either direction.
- Add
  [`steadiness_cohesiveness()`](https://jlmelville.github.io/quadra/reference/steadiness_cohesiveness.md)
  to check whether groups on the map hold together in the input and
  whether input groups hold together on the map. Neither check needs
  class labels.
- Return per-band fitted scale multipliers with detailed distance-band
  profiles.
- Add
  [`class_pair_auc()`](https://jlmelville.github.io/quadra/reference/class_pair_auc.md)
  to see which pairs of classes became easier or harder to retrieve. It
  reports ROC and PR AUC and the number of queries that contribute to
  each score.
- Add
  [`class_angular_distortion_index()`](https://jlmelville.github.io/quadra/reference/class_angular_distortion_index.md)
  for directed class-pair angular distortion.
- Add
  [`gap_index()`](https://jlmelville.github.io/quadra/reference/gap_index.md)
  for regional area distortion in two-dimensional embeddings.
- Add per-anchor counts and accuracy with
  `random_triplet_accuracy(..., by_anchor = TRUE)`.
- Add
  [`pointwise_rnx()`](https://jlmelville.github.io/quadra/reference/pointwise_rnx.md)
  for observation-level RNX AUC and selected neighborhood sizes.
- Add
  [`label_trustworthiness_continuity()`](https://jlmelville.github.io/quadra/reference/label_trustworthiness_continuity.md)
  to compare the fraction of points closer to their own centroid for
  each pair of classes in the input and the embedding.
- Add
  [`random_pair_distance_band_profile()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_band_profile.md)
  for Spearman correlation and fitted-scale stress across input-distance
  bands, with reusable pair samples.
- Add
  [`pointwise_rank_preservation()`](https://jlmelville.github.io/quadra/reference/pointwise_rank_preservation.md)
  for observation-level trustworthiness and continuity.
- Add
  [`local_procrustes()`](https://jlmelville.github.io/quadra/reference/local_procrustes.md)
  for fixed-scale shape error on matched reference neighborhoods.

### Worked examples

- Compare PCA, UMAP and t-SNE on the same MNIST images. See which pixel
  neighbors each map keeps and how many of its first few neighbors show
  the same digit as the query.
- Compare two LTSA maps whose neighbor overlap is nearly identical even
  though two images dominate the spread of one map. The Swiss-roll
  example shows how the choice of reference geometry changes which
  embedding best preserves distances.

### Bug fixes and other improvements

- Make random-triplet result construction safe during garbage
  collection, without changing sampled comparisons or random-number
  sequences.
- Neighbor-overlap calculations no longer create `.Random.seed` when it
  is absent.
- Sampled pair and triplet metrics now report nonfinite computed
  distances as numerical failures.
- [`random_pair_distance_stress()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_stress.md)
  avoids overflow and underflow when squaring large or small residuals.
- Exact-rank metrics can now be interrupted during calculation.
- [`mutual_neighbor_correlation()`](https://jlmelville.github.io/quadra/reference/mutual_neighbor_correlation.md)
  is faster.
- The global and local preservation guides now show how to detect
  collapsed embeddings whose scores appear perfect.

## quadra 0.2.1

### New features

- Random-pair distance metrics and
  [`random_triplet_accuracy()`](https://jlmelville.github.io/quadra/reference/random_triplet_accuracy.md)
  now accept reusable samples from
  [`sample_pairs()`](https://jlmelville.github.io/quadra/reference/sample_pairs.md)
  and
  [`sample_triplets()`](https://jlmelville.github.io/quadra/reference/sample_triplets.md).
  With `ret_extra = TRUE`, they return the sampled comparisons with pair
  distances or triplet agreement outcomes.
- New
  [`rnx_curve()`](https://jlmelville.github.io/quadra/reference/rnx_curve.md)
  returns exact rank agreement by neighborhood size.
  [`trustworthiness()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)
  and
  [`continuity()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)
  now accept multiple `k` values.
- [`roc_auc()`](https://jlmelville.github.io/quadra/reference/roc_auc.md)
  and
  [`pr_auc()`](https://jlmelville.github.io/quadra/reference/pr_auc.md)
  now accept `ret_extra = TRUE` to return row-aligned query AUCs and the
  number of defined queries.

### Bug fixes and minor improvements

- [`nbr_pres()`](https://jlmelville.github.io/quadra/reference/nbr_pres.md),
  [`trustworthiness()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md),
  [`continuity()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md),
  and
  [`rnx_auc()`](https://jlmelville.github.io/quadra/reference/rnx_auc.md)
  now reject nonfinite off-diagonal distances; diagonal values remain
  ignored.

## quadra 0.2.0

### New features

- New metric:
  [`trustworthiness()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)
  and
  [`continuity()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)
  for exact distance-matrix rank-penalty neighborhood preservation at a
  fixed `k`. Only useful for small datasets, though because it needs the
  full distance matrix.
- New metric:
  [`local_radius_correlation()`](https://jlmelville.github.io/quadra/reference/local_radius_correlation.md)
  for comparing local radius or scale preservation from raw observations
  or nearest-neighbor graphs with distances. Possibly useful for
  dimensionality reduction methods that attempt to retain the density of
  the input data (which standard t-SNE and UMAP do not).
- Added Spearman rank-correlation support to
  [`random_pair_distance_correlation()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_correlation.md)
  via `method = "spearman"`.
- New metric:
  [`random_pair_distance_stress()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_stress.md)
  for sampled root mean squared differences between matched random-pair
  distances. This is like a traditional STRESS metric (basically the
  root mean square error between the input and output distances), but
  only sampling a subset of distances, so this can scale to larger
  datasets.
- New metric:
  [`mutual_neighbor_correlation()`](https://jlmelville.github.io/quadra/reference/mutual_neighbor_correlation.md)
  for comparing mutual-neighbor count patterns between input and output
  nearest-neighbor graphs.
- [`random_triplet_accuracy()`](https://jlmelville.github.io/quadra/reference/random_triplet_accuracy.md),
  [`rnx_auc()`](https://jlmelville.github.io/quadra/reference/rnx_auc.md),
  [`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md)
  and
  [`nbr_pres_knn()`](https://jlmelville.github.io/quadra/reference/nbr_pres_knn.md)
  are now implemented with C++ so should be faster to calculate.

### Bug fixes and minor improvements

- Substantially improved input validation across the package.
- Corrected self-neighbor handling in
  [`nbr_pres_knn()`](https://jlmelville.github.io/quadra/reference/nbr_pres_knn.md)
  and supplied nearest-neighbor graphs.
- [`roc_auc()`](https://jlmelville.github.io/quadra/reference/roc_auc.md)
  and
  [`pr_auc()`](https://jlmelville.github.io/quadra/reference/pr_auc.md)
  now handle retrieval problems with undefined rows consistently.
- [`nbr_pres_knn()`](https://jlmelville.github.io/quadra/reference/nbr_pres_knn.md)
  now supports multithreaded calculations through `n_threads`.
- rnndescent is now on
  [CRAN](https://cran.r-project.org/package=rnndescent), so we no longer
  need to install that from github.
- `grain_size` for thread-handling has been removed from any API that
  exposed it.
- Numeric
  [`random_triplet_accuracy()`](https://jlmelville.github.io/quadra/reference/random_triplet_accuracy.md)
  sampling uses a new sampler, so the exact sample sequence can differ
  from earlier versions. Set the seed and thread count for reproducible
  results within this version. Matrix triplet inputs are unchanged.
- Fix for issue where index-only nearest neighbor graphs were not
  allowed for metrics that didn’t need the distances.
- [`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md)
  now excludes the “self”-neighbor (i.e. the item itself) from
  nearest-neighbor graphs. When calculating the graph internally,
  `k + 1` extra neighbors will be requested from `rnndescent`.

## quadra 0.1.0 (November 20 2023)

### New features

- New function: `random_triplet_accuracy`, as used in the [PaCMAP
  paper](https://jmlr.org/papers/v22/20-1061.html) for evaluating global
  preservation.
- New function: `nn_preservation` a nearest neighbor preservation
  function which can use approximate nearest neighbors and multiple
  threads for faster calculations. It also supports metrics other than
  Euclidean.
- New function: `random_pair_distance_correlation` for evaluating global
  preservation, similar to the method used by [Becht and
  co-workers](https://doi.org/10.1038/nbt.4314).
- New function: `random_pair_distance_emd` for evaluating global
  preservation, similar to the method of [Heiser and
  Lau](https://doi.org/10.1016/j.celrep.2020.107576).
- Changed license to GPL3+.
- Building this package now requires compiling C++.
- This package also now has a non-CRAN dependency:
  [rnndescent](https://github.com/jlmelville/rnndescent), which is on
  github.

### Bug fixes and minor improvements

- Added a `NEWS.md` file to track changes to the package.

## quadra 0.0.0.9000 (December 31 2021)

- Initial version. MIT Licensed.
