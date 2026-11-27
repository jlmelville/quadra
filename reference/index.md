# Package index

## Local neighborhood preservation

Compare neighbors or rank agreement between the input and embedding.

- [`nbr_pres()`](https://jlmelville.github.io/quadra/reference/nbr_pres.md)
  : Neighborhood Preservation Between Distance Matrices
- [`nbr_pres_knn()`](https://jlmelville.github.io/quadra/reference/nbr_pres_knn.md)
  : Neighborhood Preservation Between Nearest Neighbor Matrices
- [`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md)
  : Nearest Neighbor Preservation
- [`trustworthiness()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)
  [`continuity()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)
  : Trustworthiness and Continuity Between Distance Matrices
- [`pointwise_rank_preservation()`](https://jlmelville.github.io/quadra/reference/pointwise_rank_preservation.md)
  : Pointwise Trustworthiness and Continuity
- [`exact_rank_profiles()`](https://jlmelville.github.io/quadra/reference/exact_rank_profiles.md)
  : Exact Rank Profiles from Euclidean Coordinates
- [`rnx_curve()`](https://jlmelville.github.io/quadra/reference/rnx_curve.md)
  : RNX Curve
- [`rnx_auc()`](https://jlmelville.github.io/quadra/reference/rnx_auc.md)
  : Area Under the RNX Curve
- [`pointwise_rnx()`](https://jlmelville.github.io/quadra/reference/pointwise_rnx.md)
  : Pointwise RNX Agreement Across Neighborhood Sizes

## Local and regional geometry

Compare neighborhood shape, regional area, radii, and shared or
reciprocal neighbors.

- [`local_procrustes()`](https://jlmelville.github.io/quadra/reference/local_procrustes.md)
  : Local Procrustes Error
- [`gap_index()`](https://jlmelville.github.io/quadra/reference/gap_index.md)
  : Gap Index for Two-Dimensional Embeddings
- [`local_radius_correlation()`](https://jlmelville.github.io/quadra/reference/local_radius_correlation.md)
  : Local Radius Correlation
- [`mutual_neighbor_correlation()`](https://jlmelville.github.io/quadra/reference/mutual_neighbor_correlation.md)
  : Mutual Neighbor Correlation
- [`neighbor_dissimilarity()`](https://jlmelville.github.io/quadra/reference/neighbor_dissimilarity.md)
  : Neighbor Dissimilarity Between Two Graphs

## Global distance preservation

Compare sampled distances and distance orderings.

- [`random_triplet_accuracy()`](https://jlmelville.github.io/quadra/reference/random_triplet_accuracy.md)
  : Random Triplet Accuracy
- [`random_pair_distance_correlation()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_correlation.md)
  : Random Pair Distance Correlation
- [`random_pair_distance_emd()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_emd.md)
  : Earth Mover's Distance Between Random-Pair Distances
- [`random_pair_distance_stress()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_stress.md)
  : Random Pair Distance Stress
- [`random_pair_distance_band_profile()`](https://jlmelville.github.io/quadra/reference/random_pair_distance_band_profile.md)
  : Distance-Band Fidelity for Random Pairs

## Label retrieval and class geometry

Measure label retrieval, centroid consistency, and angular distortion.

- [`roc_auc()`](https://jlmelville.github.io/quadra/reference/roc_auc.md)
  : Average ROC AUC for Label Retrieval
- [`pr_auc()`](https://jlmelville.github.io/quadra/reference/pr_auc.md)
  : Average Precision-Recall AUC for Label Retrieval
- [`class_pair_auc()`](https://jlmelville.github.io/quadra/reference/class_pair_auc.md)
  : Retrieval AUC Within Each Pair of Classes
- [`label_trustworthiness_continuity()`](https://jlmelville.github.io/quadra/reference/label_trustworthiness_continuity.md)
  : Label Trustworthiness and Continuity
- [`class_angular_distortion_index()`](https://jlmelville.github.io/quadra/reference/class_angular_distortion_index.md)
  : Class Angular Distortion Index

## Group preservation

- [`steadiness_cohesiveness()`](https://jlmelville.github.io/quadra/reference/steadiness_cohesiveness.md)
  : Steadiness and Cohesiveness of Groups

## Sampling helpers

Sample observation pairs or anchored triplets.

- [`sample_pairs()`](https://jlmelville.github.io/quadra/reference/sample_pairs.md)
  : Sample Observation Pairs
- [`sample_triplets()`](https://jlmelville.github.io/quadra/reference/sample_triplets.md)
  : Sample Anchored Triplets
- [`sample_class_triplets()`](https://jlmelville.github.io/quadra/reference/sample_class_triplets.md)
  : Sample Triplets for Class Angular Distortion
