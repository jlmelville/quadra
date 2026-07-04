
# quadra 0.2.0

## New features

* New metric: `trustworthiness()` and `continuity()` for exact distance-matrix rank-penalty
neighborhood preservation at a fixed `k`. Only useful for small datasets, though because it needs
the full distance matrix.
* New metric: `local_radius_correlation()` for comparing local radius or scale preservation from
raw observations or nearest-neighbor graphs with distances. Possibly useful for dimensionality
reduction methods that attempt to retain the density of the input data (which standard t-SNE and
UMAP do not).
* Added Spearman rank-correlation support to `random_pair_distance_correlation()` via
`method = "spearman"`.
* New metric: `random_pair_distance_stress()` for sampled root mean squared differences between
matched random-pair distances. This is like a traditional STRESS metric (basically the root mean
square error between the input and output distances), but only sampling a subset of distances, so
this can scale to larger datasets.
* New metric: `mutual_neighbor_correlation()` for comparing mutual-neighbor count patterns between
input and output nearest-neighbor graphs.
* `random_triplet_accuracy()`, `rnx_auc()`, `nn_preservation()` and `nbr_pres_knn()` are now
implemented with C++ so should be faster to calculate.

## Bug fixes and minor improvements

* Many small fixes for corner cases and improved validation.
* rnndescent is now on [CRAN](https://cran.r-project.org/package=rnndescent), so we no longer need
to install that from github.
* `grain_size` for thread-handling has been removed from any API that exposed it.
* Numeric `random_triplet_accuracy()` sampling uses a new sampler, so the exact sample sequence can
differ from earlier versions. Set the seed and thread count for reproducible results within this
version. Matrix triplet inputs are unchanged.
* Fix for issue where index-only nearest neighbor graphs were not allowed for metrics that didn't
need the distances.
* `nn_preservation()` now excludes the "self"-neighbor (i.e. the item itself) from nearest-neighbor
graphs. When calculating the graph internally, `k + 1` extra neighbors will be requested from
`rnndescent`.

# quadra 0.1.0 (November 20 2023)

## New features

* New function: `random_triplet_accuracy`, as used in the
[PaCMAP paper](https://jmlr.org/papers/v22/20-1061.html) for evaluating global preservation.
* New function: `nn_preservation` a nearest neighbor preservation function which can use
approximate nearest neighbors and multiple threads for faster calculations. It also supports
metrics other than Euclidean.
* New function: `random_pair_distance_correlation` for evaluating global preservation, similar to
the method used by [Becht and co-workers](https://doi.org/10.1038/nbt.4314).
* New function: `random_pair_distance_emd` for evaluating global preservation, similar to the
method of [Heiser and Lau](https://doi.org/10.1016/j.celrep.2020.107576).
* Changed license to GPL3+.
* Building this package now requires compiling C++.
* This package also now has a non-CRAN dependency:
[rnndescent](https://github.com/jlmelville/rnndescent), which is on github.

## Bug fixes and minor improvements

* Added a `NEWS.md` file to track changes to the package.

# quadra 0.0.0.9000 (December 31 2021)

* Initial version. MIT Licensed.
