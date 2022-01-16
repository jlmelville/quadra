# quadra 0.1.0.9000

## New features

* New function: `random_triplet_accuracy`, as used in the [PaCMAP
paper](https://jmlr.org/papers/v22/20-1061.html) for evaluating global
preservation.
* New function: `nn_preservation` a nearest neighbor preservation function which
can use approximate nearest neighbors and multiple threads for faster
calculations. It also supports metrics other than Euclidean.
* New function: `random_distance_pair_correlation` for evaluating global
preservation, similar to the method used by [Becht and
co-workers](https://doi.org/10.1038/nbt.4314).
* New function: `random_pair_distance_emd` for evaluating global preservation,
similar to the method of [Heiser and
Lau](https://doi.org/10.1016/j.celrep.2020.107576).
* Changed license to GPL3+.
* Building this package now requires compiling C++.
* This package also now has a non-CRAN dependency:
[rnndescent](https://github.com/jlmelville/rnndescent), which is on github.

## Bug fixes and minor improvements

* Added a `NEWS.md` file to track changes to the package.

# quadra 0.0.0.9000

* Initial version. MIT Licensed.
