# Quadra: QUantitative Assessment of Dimensionality Reduction Algorithms

An R Package for evaluating the success of embeddings from dimensionality
reduction methods (e.g. Principal Component Analysis, Sammon Maps, 
t-Distributed Stochastic Neighbor Embedding).

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

## Examples

```R
# Embed with first two scores from PCA
pca_iris <- stats::prcomp(iris[, -5], retx = TRUE, rank. = 2)

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

### Neighborhood Preservation

France, S., & Carroll, D. (2007, July). 
Development of an agreement metric based upon the RAND index for the evaluation 
of dimensionality reduction techniques, with applications to mapping customer 
data.
In *International Workshop on Machine Learning and Data Mining in Pattern Recognition*
(pp. 499-517). Springer, Berlin, Heidelberg.
https://doi.org/10.1007/978-3-540-73499-4_38

Chen, L., & Buja, A. (2009). Local multidimensional scaling for nonlinear 
dimension reduction, graph drawing, and proximity analysis. 
*Journal of the American Statistical Association*, *104*(485), 209-219.
http://dx.doi.org/10.1198/jasa.2009.0111

Lee, J. A., & Verleysen, M. (2009).
Quality assessment of dimensionality reduction: Rank-based criteria.
*Neurocomputing*, *72*(7), 1431-1443.
https://dx.doi.org/10.1016/j.neucom.2008.12.017

Venna, J., Peltonen, J., Nybo, K., Aidos, H., & Kaski, S. (2010). 
Information retrieval perspective to nonlinear dimensionality reduction for data visualization. 
*Journal of Machine Learning Research*, *11*(Feb), 451-490.
http://www.jmlr.org/papers/v11/venna10a.html

Lee, J. A., Peluffo-Ordónez, D. H., & Verleysen, M. (2015).
Multi-scale similarities in stochastic neighbour embedding: Reducing
dimensionality while preserving both local and global structure.
*Neurocomputing*, *169*, 246-261.
https://dx.doi.org/10.1016/j.neucom.2014.12.095

### Precision-Recall AUC and Receiver Operating Characteristic Area Under the Curve

Davis, J., & Goadrich, M. (2006, June).
The relationship between Precision-Recall and ROC curves.
In *Proceedings of the 23rd international conference on Machine learning*
(pp. 233-240). ACM.
http://pages.cs.wisc.edu/~jdavis/davisgoadrichcamera2.pdf

Keilwagen, J., Grosse, I., & Grau, J. (2014).
Area under precision-recall curves for weighted and unweighted data.
*PloS One*, *9*(3), e92209.
https://dx.doi.org/10.1371/journal.pone.0092209.

## License

[MIT](https://opensource.org/licenses/MIT).
