# Local preservation

## Neighborhood Preservation

For local preservation use nearest neighbor preservation. Pass a vector
to `k` to get back the preservation for different numbers of neighbors.

``` r

iris_x <- as.matrix(iris[, -5])
pca_iris <- stats::prcomp(iris_x, retx = TRUE, rank. = 2)$x
nn_preservation(iris_x, pca_iris, k = c(15, 30))
```

## Trustworthiness and Continuity

[`trustworthiness()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)
and
[`continuity()`](https://jlmelville.github.io/quadra/reference/trustworthiness.md)
are exact rank-penalty neighborhood metrics for distance matrices.
Trustworthiness penalizes observations that appear as low-dimensional
neighbors even though they are not close in the input space. Continuity
applies the dual penalty to input-space neighbors that are no longer
nearby in the embedding.

``` r

din <- as.matrix(stats::dist(iris_x))
dout <- as.matrix(stats::dist(pca_iris))
trustworthiness(din, dout, k = 15)
continuity(din, dout, k = 15)
```

These functions use full `n` by `n` distance matrices and require `k` to
be less than half the number of observations, so they are intended for
exact small-dataset checks. For larger datasets, use nearest-neighbor
preservation metrics such as
[`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md)
or
[`nbr_pres_knn()`](https://jlmelville.github.io/quadra/reference/nbr_pres_knn.md).

## Local Radius Correlation

Nearest neighbor preservation asks whether the same observations remain
nearby.
[`local_radius_correlation()`](https://jlmelville.github.io/quadra/reference/local_radius_correlation.md)
asks a different local-scale question: do points in dense or sparse
regions of the input data still have relatively small or large
nearest-neighbor radii in the embedding?

``` r

local_radius_correlation(
  iris_x,
  pca_iris,
  k = c(15, 30),
  nn_method_in = "brute"
)
```

By default, the metric uses Spearman correlation between the distance to
the `k`th nearest non-self neighbor in the input data and output
embedding. Set `statistic = "mean"` to use the mean distance to the
first `k` neighbors, or `log = TRUE` to compare log radii when all
selected radii are positive.

If you already have nearest-neighbor graphs with distance matrices, you
can reuse them:

``` r

cached <- local_radius_correlation(
  iris_x,
  pca_iris,
  k = c(15, 30),
  nn_method_in = "brute",
  ret_extra = TRUE
)

local_radius_correlation(cached$nn_in, cached$nn_out, k = c(15, 30))
```

This metric is a radius or scale diagnostic. It does not estimate local
density directly and should be interpreted alongside neighbor-identity
metrics such as
[`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md).

## Mutual Neighbor Correlation

[`mutual_neighbor_correlation()`](https://jlmelville.github.io/quadra/reference/mutual_neighbor_correlation.md)
compares mutual-neighbor count patterns. For each observation, the
mutual-neighbor count is the number of its first `k` nearest neighbors
that also include the observation among their first `k` nearest
neighbors.

``` r

mutual_neighbor_correlation(
  iris_x,
  pca_iris,
  k = c(15, 30),
  nn_method_in = "brute"
)
```

This is a graph diagnostic. It does not require neighbor distances, so
cached idx-only nearest-neighbor graphs are enough. It should be
interpreted alongside neighbor-identity metrics such as
[`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md)
and scale metrics such as
[`local_radius_correlation()`](https://jlmelville.github.io/quadra/reference/local_radius_correlation.md).

For larger datasets,
[`nn_preservation()`](https://jlmelville.github.io/quadra/reference/nn_preservation.md)
can use approximate nearest neighbors via
[rnndescent](https://github.com/jlmelville/rnndescent). If you have
nearest-neighbor matrices from another approximate nearest neighbor
package, e.g. [RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy),
use the `nbr_pres_knn` function:

``` r

# A function to wrap the RcppAnnoy API to return a matrix of the nearest
# neighbor indices
find_nn <- function(X, k = 10, n_trees = 50, search_k = k * n_trees) {
  nr <- nrow(X)
  nc <- ncol(X)

  ann <- methods::new(RcppAnnoy::AnnoyEuclidean, nc)
  for (i in 1:nr) {
    ann$addItem(i - 1, X[i, ])
  }
  ann$build(n_trees)

  idx <- matrix(nrow = nr, ncol = k)
  for (i in 1:nr) {
    res <- ann$getNNsByItemList(i - 1, k, search_k, FALSE)
    if (length(res$item) != k) {
      stop("search_k/n_trees settings were unable to find ", k,
           " neighbors for item ", i)
    }
    idx[i, ] <- res$item
  }
  idx + 1
}

kin <- find_nn(as.matrix(iris[, -5]), k = 5)
kout <- find_nn(pca_iris, k = 5)

# This should give very similar results to using nbr_pres on the distance
# matrices, subject to the approximations employed by Annoy
nbr_pres_knn(kin, kout, k = 5)
```

## Further Reading

France, S., & Carroll, D. (2007, July). Development of an agreement
metric based upon the RAND index for the evaluation of dimensionality
reduction techniques, with applications to mapping customer data. In
*International Workshop on Machine Learning and Data Mining in Pattern
Recognition* (pp. 499-517). Springer, Berlin, Heidelberg.
<https://doi.org/10.1007/978-3-540-73499-4_38>

Chen, L., & Buja, A. (2009). Local multidimensional scaling for
nonlinear dimension reduction, graph drawing, and proximity analysis.
*Journal of the American Statistical Association*, *104*(485), 209-219.
<http://dx.doi.org/10.1198/jasa.2009.0111>

Lee, J. A., & Verleysen, M. (2009). Quality assessment of dimensionality
reduction: Rank-based criteria. *Neurocomputing*, *72*(7), 1431-1443.
<https://dx.doi.org/10.1016/j.neucom.2008.12.017>

Lee, J. A., Peluffo-Ordonez, D. H., & Verleysen, M. (2015). Multi-scale
similarities in stochastic neighbour embedding: Reducing dimensionality
while preserving both local and global structure. *Neurocomputing*,
*169*, 246-261. <https://dx.doi.org/10.1016/j.neucom.2014.12.095>

Cooley, S. M., Hamilton, T., Deeds, E. J., & Ray, J. C. J. (2019). A
novel metric reveals previously unrecognized distortion in
dimensionality reduction of scRNA-Seq data. *BioRxiv*, 689851.
<https://www.biorxiv.org/content/10.1101/689851v6>

Narayan, A., Berger, B., & Cho, H. (2021). Assessing single-cell
transcriptomic variability through density-preserving data
visualization. *Nature Biotechnology*, *39*, 765-774.
<https://doi.org/10.1038/s41587-020-00801-7>
