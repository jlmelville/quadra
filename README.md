[![Last Commit](https://img.shields.io/github/last-commit/jlmelville/quadra)](https://github.com/jlmelville/quadra)

# Quadra: QUantitative Assessment of Dimensionality Reduction Algorithms

An R Package for evaluating the success of embeddings from dimensionality
reduction methods (e.g. Principal Component Analysis, Sammon Maps, 
t-Distributed Stochastic Neighbor Embedding).

## News

*January 2 2022*: To build this package you will now need to compile some C++.
Also, it depends on the [rnndescent](https://github.com/jlmelville/rnndescent)
package which is not on CRAN. In return, a new function is available: a
multi-threaded implementation of the `random_triplet_accuracy` technique for
evaluating global structure preservation (as used in the 
[PaCMAP paper](https://jmlr.org/papers/v22/20-1061.html)).

*December 31 2021*: Package is getting relicensed to GPL3+ so I can use some
other packages and code. Last MIT license version is release 0.0.0.9000.

## Description

This package provides two ways to evaluate the performance of an embedding.
The most generic is to consider "neighborhood preservation": 

1. Find the *k*-nearest neighbors of a point in the input space.
1. Do the same in the output space.
1. Count the overlap.

This only requires calculating the Euclidean distance matrix in the input and
output spaces, without any other assumptions about the type of dimensionality
reduction carried out. This package provides some quality measures based
around this concept.

Alternatively, if some sort of labelling is applied to the points, each point
can be treated as the target in a retrieval procedure: 

1. Rank all the other points by distance to the target point.
1. See how highly in the ranked list the points with the same label are found.
1. Construct a Receiver Operating Characteristic (ROC) curve, or something like
it (e.g. Precision-Recall curve)
1. Calculate the Area Under the Curve (AUC).
1. Average over all points.

This only needs the output distance matrix, but requires the sort of labeling
usually reserved for data intended for supervised classification. Quadra can
also provide some help with this, but requires the [PRROC package](https://cran.r-project.org/package=PRROC) to be installed.

## Installing

```R
install.packages("devtools")
devtools::install_github("jlmelville/quadra")
```

`quadra` makes use of C++ code which must be compiled. You may have to carry out
a few extra steps before being able to build this package:

**Windows**: install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) and ensure
`C:\Rtools\bin` is on your path.

**Mac OS X**: using a custom `~/.R/Makevars`
[may cause linking errors](https://github.com/jlmelville/uwot/issues/1).
This sort of thing is a potential problem on all platforms but seems to bite
Mac owners more.
[The R for Mac OS X FAQ](https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html#Installation-of-source-packages)
may be helpful here to work out what you can get away with. To be on the safe
side, I would advise building `quadra` without a custom `Makevars`.

## Examples

```R
# Embed with first two scores from PCA
pca_iris <- stats::prcomp(iris[, -5], retx = TRUE, rank. = 2)

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
random_triplet_accuracy(tiris, pca_iris2), n_threads = 2, is_transposed = TRUE)
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
dout <- as.matrix(dist(pca_iris$x))

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

## Limitations

Early days for this package. The interface could be more friendly, and the
implementation a lot more efficient. Currently, I don't recommend using this on
datasets larger than ~1000 observations.

For larger datasets, use an Approximate Nearest Neighbors package, e.g
[RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy) and the 
`nbr_pres_knn` function:

```R
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
kout <- find_nn(pca_iris$x, k = 5)

# This should give very similar results to using nbr_pres on the distance
# matrices, subject to the approximations employed by Annoy
nbr_pres_knn(kin, kout, k = 5)
```

## Further Reading

### Random Triplet Accuracy

Wang, Y., Huang, H., Rudin, C., & Shaposhnik, Y. (2021). 
Understanding how dimension reduction tools work: an empirical approach to 
deciphering t-SNE, UMAP, TriMAP, and PaCMAP for data visualization.
*J Mach. Learn. Res*, *22*, 1-73.
<https://jmlr.org/papers/v22/20-1061.html>

### Pearson Correlation for Distances

Becht, E., McInnes, L., Healy, J., Dutertre, C. A., Kwok, I. W., Ng, L. G., ... & Newell, E. W. (2019).
Dimensionality reduction for visualizing single-cell data using UMAP.
*Nature biotechnology*, *37*(1), 38-44.
<https://doi.org/10.1038/nbt.4314>

### Earth-Mover's Distance

Heiser, C. N., & Lau, K. S. (2020). 
A quantitative framework for evaluating single-cell data structure preservation by dimensionality reduction techniques. 
*Cell reports*, *31*(5), 107576.
<https://doi.org/10.1016/j.celrep.2020.107576>

### Neighborhood Preservation

France, S., & Carroll, D. (2007, July). 
Development of an agreement metric based upon the RAND index for the evaluation 
of dimensionality reduction techniques, with applications to mapping customer 
data.
In *International Workshop on Machine Learning and Data Mining in Pattern Recognition*
(pp. 499-517). Springer, Berlin, Heidelberg.
<https://doi.org/10.1007/978-3-540-73499-4_38>

Chen, L., & Buja, A. (2009). Local multidimensional scaling for nonlinear 
dimension reduction, graph drawing, and proximity analysis. 
*Journal of the American Statistical Association*, *104*(485), 209-219.
<http://dx.doi.org/10.1198/jasa.2009.0111>

Lee, J. A., & Verleysen, M. (2009).
Quality assessment of dimensionality reduction: Rank-based criteria.
*Neurocomputing*, *72*(7), 1431-1443.
<https://dx.doi.org/10.1016/j.neucom.2008.12.017>

Venna, J., Peltonen, J., Nybo, K., Aidos, H., & Kaski, S. (2010). 
Information retrieval perspective to nonlinear dimensionality reduction for data visualization. 
*Journal of Machine Learning Research*, *11*(Feb), 451-490.
<http://www.jmlr.org/papers/v11/venna10a.html>

Lee, J. A., Peluffo-Ordónez, D. H., & Verleysen, M. (2015).
Multi-scale similarities in stochastic neighbour embedding: Reducing
dimensionality while preserving both local and global structure.
*Neurocomputing*, *169*, 246-261.
<https://dx.doi.org/10.1016/j.neucom.2014.12.095>

Cooley, S. M., Hamilton, T., Deeds, E. J., & Ray, J. C. J. (2019).
A novel metric reveals previously unrecognized distortion in dimensionality reduction of scRNA-Seq data. 
*BioRxiv*, 689851.
<https://www.biorxiv.org/content/10.1101/689851v6>

### Precision-Recall AUC and Receiver Operating Characteristic Area Under the Curve

Davis, J., & Goadrich, M. (2006, June).
The relationship between Precision-Recall and ROC curves.
In *Proceedings of the 23rd international conference on Machine learning*
(pp. 233-240). ACM.
<http://pages.cs.wisc.edu/~jdavis/davisgoadrichcamera2.pdf>

Keilwagen, J., Grosse, I., & Grau, J. (2014).
Area under precision-recall curves for weighted and unweighted data.
*PloS One*, *9*(3), e92209.
<https://dx.doi.org/10.1371/journal.pone.0092209>

## License

[GPLv3 or later](https://www.gnu.org/licenses/gpl-3.0.txt).

## See Also

Some other R packages I maintain that might come in handy:

* [snedata](https://github.com/jlmelville/snedata) and 
(https://github.com/jlmelville/coil20) provide ways to access some datasets
for embedding.
* [vizier](https://github.com/jlmelville/vizier) for a... let's call it cheap
and cheerful way to visualize the results of the embedding.
