# Quadra: QUantitative Assessment of Dimensionality Reduction Algorithms

An R Package for evaluating the success of embeddings from
dimensionality reduction methods (e.g. Principal Component Analysis,
Sammon Maps, t-Distributed Stochastic Neighbor Embedding).

## Description

This package provides two ways to evaluate the performance of an
embedding. The most generic is to consider “neighborhood preservation”:

1.  Find the *k*-nearest neighbors of a point in the input space.
2.  Do the same in the output space.
3.  Count the overlap.

This only requires calculating the Euclidean distance matrix in the
input and output spaces, without any other assumptions about the type of
dimensionality reduction carried out. This package provides some quality
measures based around this concept.

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

## Install

``` r

install.packages("pak")
pak::pak("jlmelville/quadra")
```

`quadra` makes use of C++ code which must be compiled. You may have to
carry out a few extra steps before being able to build this package:

**Windows**: install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) and ensure
`C:\Rtools\bin` is on your path.

**macOS**: using a custom `~/.R/Makevars` [may cause linking
errors](https://github.com/jlmelville/uwot/issues/1). This sort of thing
is a potential problem on all platforms but seems to bite Mac owners
more. [The R for Mac OS X
FAQ](https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html#Installation-of-source-packages)
may be helpful here to work out what you can get away with. To be on the
safe side, I would advise building `quadra` without a custom `Makevars`.

## Examples

``` r

# Embed with first two scores from PCA
pca_iris <- stats::prcomp(iris[, -5], retx = TRUE, rank. = 2)$x

# Random triplet accuracy: proportion of triangle distances where the relative
# ordering is retained. Used by Wang et al (2021) to measure global structure
# preservation
random_triplet_accuracy(iris, pca_iris)

# For large datasets, you can set the number of threads to use
# (pretend this is a big dataset)
random_triplet_accuracy(iris, pca_iris, n_threads = 2)

# If you plan to carry out multiple comparisons with the same (high dimensional)
# data, transpose once outside the function and set is_transposed = TRUE for a
# slight speed-up
tiris <- t(iris[, -5])

pca_iris2 <- t(pca_iris)
pca_iris1 <- t(stats::prcomp(iris[, -5], retx = TRUE, rank. = 1)$x)
pca_iris3 <- t(stats::prcomp(iris[, -5], retx = TRUE, rank. = 3)$x)
random_triplet_accuracy(tiris, pca_iris2, n_threads = 2, is_transposed = TRUE)
random_triplet_accuracy(tiris, pca_iris1, is_transposed = TRUE, n_threads = 2)
random_triplet_accuracy(tiris, pca_iris3, is_transposed = TRUE, n_threads = 2)

# Other ways to measure global preservation:

# measure Pearson correlation between equivalent input and output distances
# similar to the method used by Becht and co-workers
random_pair_distance_correlation(tiris, pca_iris3, is_transposed = TRUE, n_threads = 2)

# convert distances to an empirical distribution and compare them via Earth
# Mover's Distance (1D Wasserstein) based on the method used by Heiser and Lau
random_pair_distance_emd(tiris, pca_iris2, is_transposed = TRUE, n_threads = 2)

# For local preservation use nearest neighbor preservation
# Pass a vector to k to get back the preservation for different numbers of
# neighbors
nn_preservation(tiris, pca_iris2, k = c(15, 30), is_transposed = TRUE)

# slower methods not recommended for large datasets

din <- as.matrix(dist(iris[, -5]))
dout <- as.matrix(dist(pca_iris))

# Preservation of the 5-nearest neighbors for each point: returns a vector
nbr_pres(din, dout, k = 5)

# Area under the RNX curve. This is like a weighted average over neighborhood preservation for a range of k
# with a bias towards smaller k: returns a single value
rnx_auc(din, dout)

# Install the PRROC package
install.packages("PRROC")

# ROC AUC using Species as the label: returns a list with one AUC per factor level
# Also the average over all levels
roc_auc(dout, iris$Species)

# Precision-Recall AUC which some people prefer over ROC
pr_auc(dout, iris$Species)
```

## More Detail

Longer metric notes, larger-dataset examples, and references have moved
to the pkgdown articles:

- [Metric
  overview](https://jlmelville.github.io/quadra/articles/metrics-overview.html)
- [Global
  preservation](https://jlmelville.github.io/quadra/articles/global-preservation.html)
- [Local
  preservation](https://jlmelville.github.io/quadra/articles/local-preservation.html)
- [Label
  retrieval](https://jlmelville.github.io/quadra/articles/label-retrieval.html)

## License

[GPLv3 or later](https://www.gnu.org/licenses/gpl-3.0.txt).

## See Also

Some other R packages I maintain that might come in handy:

- [snedata](https://github.com/jlmelville/snedata) and
  [coil20](https://github.com/jlmelville/coil20) provide ways to access
  some datasets for embedding.
- [vizier](https://github.com/jlmelville/vizier) for a… let’s call it
  cheap and cheerful way to visualize the results of the embedding.
