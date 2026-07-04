# Neighborhood Retrieval Metrics

#' Neighborhood Preservation Between Distance Matrices
#'
#' Calculates the neighborhood preservation for each observation in a dataset,
#' represented by two distance matrices. The first matrix is the "ground truth",
#' the second being the estimation or approximation.
#' The neighborhood preservation is calculated for each row where each element
#' `d[i, j]` is taken to be the distance between observation `i` and `j`.
#'
#' The neighborhood preservation can vary between 0 (no neighbors in common)
#' and 1 (perfect preservation). However, random performance gives an
#' approximate value of k / (n - 1), where k is the size of the neighborhood and
#' n is the number of observations or items in the dataset.
#'
#' Self-neighbors on the diagonal are excluded from each row before the
#' neighborhood overlap is calculated.
#'
#' @note This is not a very efficient way to calculate the preservation if you
#'  want to calculate the value for multiple values of `k`. For more global
#'  measures of preservation, see [rnx_auc()].
#'
#' @param din Distance matrix. The "ground truth" or reference distances.
#' @param dout Distance matrix. A set of distances to compare to the reference
#'  distances.
#' @param k The size of the neighborhood, where k is the number of neighbors to
#'  include in the neighborhood.
#' @return Vector of preservation values, one for each row of the distance
#'  matrix.
#' @export
nbr_pres <- function(din, dout, k) {
  validate_distance_matrix_pair(din, dout)
  k <- validate_positive_integer(k, "k")
  max_k <- ncol(din) - 1L
  if (k > max_k) {
    stop(
      "k cannot be larger than the number of non-self observations",
      call. = FALSE
    )
  }
  preservations <- vector(mode = "numeric", length = nrow(din))
  for (i in seq_len(nrow(din))) {
    di <- din[i, ]
    dj <- dout[i, ]
    di[i] <- Inf
    dj[i] <- Inf
    preservations[i] <- nbr_pres_i(di, dj, k)
  }
  preservations
}


#' Neighborhood Preservation Between Nearest Neighbor Matrices
#'
#' Calculates the neighborhood preservation for each observation in a dataset,
#' represented by two matrices of the indices of the nearest neighbors. The
#' first matrix is the "ground truth", the second being the estimation or
#' approximation. The neighborhood preservation is calculated for each row where
#' each element `d[i, k]` is taken to be the index of the kth nearest neighbor
#' of `i`.
#'
#' Approximate nearest neighbor methods, e.g.
#' [RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy), can find
#' k-nearest neighbors quite efficiently and so makes calculating preservation
#' values for larger datasets feasible.
#'
#' The neighborhood preservation can vary between 0 (no neighbors in common)
#' and 1 (perfect preservation). For nearest-neighbor matrices that exclude
#' self-neighbors, random performance gives an approximate value of k / (n - 1),
#' where k is the size of the neighborhood and n is the number of observations
#' or items in the dataset.
#'
#' @param kin Nearest neighbor matrix. The "ground truth" or reference indices.
#' @param kout Nearest neighbor matrix. A set of distances to compare to the reference
#'  indices.
#' @param k The size of the neighborhood, where k is the number of neighbors to
#'  include in the neighborhood.
#' @return Vector of preservation values, one for each row of `kin`.
#' @export
nbr_pres_knn <- function(kin, kout, k = ncol(kin)) {
  if (!methods::is(kin, "matrix")) {
    stop("kin must be a matrix", call. = FALSE)
  }
  if (!methods::is(kout, "matrix")) {
    stop("kout must be a matrix", call. = FALSE)
  }
  k <- validate_positive_integer(k, "k")
  if (k > ncol(kin) || k > ncol(kout)) {
    stop("k cannot be larger than the number of columns in kin or kout")
  }
  if (nrow(kin) != nrow(kout)) {
    stop("kin and kout must have the same number of rows", call. = FALSE)
  }
  counts <- neighbor_overlap_counts(kin, kout, k)
  counts[, 1] * (1 / k)
}

#' Trustworthiness and Continuity Between Distance Matrices
#'
#' `trustworthiness()` penalizes observations that appear among the `k` nearest
#' neighbors in `dout` but have input-space rank greater than `k` in `din`.
#' `continuity()` applies the dual penalty to input-space neighbors that are no
#' longer among the `k` nearest neighbors in `dout`.
#'
#' Both functions use exact ranks from the supplied distance matrices and
#' exclude the diagonal self-neighbor from each row. Tied distances are ranked in
#' their original column order after self-neighbor exclusion, matching
#' `rank(ties.method = "first")`.
#'
#' Because these functions require full `n` by `n` distance matrices, they are
#' practical only for small datasets. For larger datasets, use nearest-neighbor
#' preservation metrics such as [nn_preservation()] or [nbr_pres_knn()].
#'
#' Unlike [nbr_pres()], which only counts shared neighbors, these metrics weight
#' each unexpected or missing neighbor by how far its rank lies outside the
#' `k`-neighborhood. [rnx_auc()] also uses rank-based neighborhood agreement,
#' but aggregates across neighborhood sizes; these functions report the standard
#' trustworthiness or continuity score at one `k`.
#'
#' @param din Input distance matrix. The "ground truth" or reference distances.
#' @param dout Output distance matrix. A set of distances to compare to the
#'   reference distances.
#' @param k The size of the neighborhood. Must be a positive integer less than
#'   half the number of observations so the standard 0-1 normalization remains
#'   bounded.
#' @return A scalar score. A value of 1 indicates no rank-penalty errors at
#'   neighborhood size `k`; lower values indicate worse preservation.
#' @references
#' Venna, J., & Kaski, S. (2001). Neighborhood preservation in nonlinear
#' projection methods: An experimental study. In *Artificial Neural Networks -
#' ICANN 2001* (pp. 485-491).
#' @export
trustworthiness <- function(din, dout, k) {
  validate_distance_matrix_pair(din, dout)
  k <- validate_rank_penalty_k(k, nrow(din))
  trustworthiness_exact(din, dout, k)
}

#' @rdname trustworthiness
#' @export
continuity <- function(din, dout, k) {
  validate_distance_matrix_pair(din, dout)
  k <- validate_rank_penalty_k(k, nrow(din))
  continuity_exact(din, dout, k)
}

#' Area Under the RNX Curve
#'
#' The RNX curve is formed by calculating the `rnx_crm` metric for
#' different sizes of neighborhood. Each value of RNX is scaled according to
#' the natural log of the neighborhood size, to give a higher weight to smaller
#' neighborhoods. An AUC of 1 indicates perfect neighborhood preservation, an
#' AUC of 0 is due to random results. Self-neighbors on the distance-matrix
#' diagonal are excluded before the co-ranking matrix is calculated.
#'
#' @param din Input distance matrix.
#' @param dout Output distance matrix.
#' @return Area under the RNX curve.
#' @references
#' Lee, J. A., Peluffo-Ordo'nez, D. H., & Verleysen, M. (2015).
#' Multi-scale similarities in stochastic neighbour embedding: Reducing
#' dimensionality while preserving both local and global structure.
#' *Neurocomputing*, *169*, 246-261.
#' @export
rnx_auc <- function(din, dout) {
  validate_distance_matrix_pair(din, dout)
  if (nrow(din) < 3L) {
    stop("RNX AUC requires at least three observations", call. = FALSE)
  }
  rnx_auc_direct(din, dout)
}

# Co-ranking Matrix
#
# Calculates the co-ranking matrix for an embedding.
#
# The co-ranking matrix is the basic data structure used for calculating
# various quality metrics, such as \code{qnx_crm},
# \code{rnx_crm} and \code{bnx_crm}.
#
# The co-ranking matrix is an (N - 1) x (N - 1) matrix where N is the number of
# observations. The diagonal self-neighbor is excluded. The element (i, j) is the
# number of times an ith-nearest neighbor of an observation in the input
# distance matrix was the jth-nearest neighbor in the output space.
#
# The lower diagonal represents "intrusions". This is when observations
# have a larger rank in the input space than in the output space,
# i.e. non-neighbors are falsely marked as neighbors in the output space.
#
# The upper diagonal represents "extrusions". This occurs when observations
# have a smaller rank in the input space than in the output space,
# i.e. true neighbors are falsely marked as non-neighbors in the output space.
#
# @param din Input distance matrix.
# @param dout Output distance matrix.
# @return Co-ranking matrix.
# @references
# Lee, J. A., & Verleysen, M. (2009).
# Quality assessment of dimensionality reduction: Rank-based criteria.
# \emph{Neurocomputing}, \emph{72(7)}, 1431-1443.
coranking_matrix <- function(din, dout) {
  validate_distance_matrix_pair(din, dout)
  n <- nrow(din)
  crm <- matrix(0, nrow = n - 1L, ncol = n - 1L)
  for (i in seq_len(nrow(din))) {
    rin <- rank(din[i, -i], ties.method = "first")
    rout <- rank(dout[i, -i], ties.method = "first")
    for (j in seq_along(rin)) {
      crm[rin[j], rout[j]] <- crm[rin[j], rout[j]] + 1
    }
  }
  crm
}

# Area Under the RNX Curve
#
# The RNX curve is formed by calculating the \code{rnx_crm} metric for
# different sizes of neighborhood. Each value of RNX is scaled according to
# the natural log of the neighborhood size, to give a higher weight to smaller
# neighborhoods. An AUC of 1 indicates perfect neighborhood preservation, an
# AUC of 0 is due to random results.
#
# @param crm Co-ranking matrix.
# @return Area under the curve.
# @references
# Lee, J. A., Peluffo-Ordo'nez, D. H., & Verleysen, M. (2015).
# Multi-scale similarities in stochastic neighbour embedding: Reducing
# dimensionality while preserving both local and global structure.
# \emph{Neurocomputing}, \emph{169}, 246-261.
rnx_auc_crm <- function(crm) {
  n_ranks <- nrow(crm)
  if (n_ranks < 2L) {
    return(NA_real_)
  }
  k <- seq_len(n_ranks - 1L)
  top_left <- diag(cumulative_top_left_sums(crm))[k]
  qnx <- top_left / (k * (n_ranks + 1L))
  rnx <- ((qnx * n_ranks) - k) / (n_ranks - k)
  sum(rnx / k) / sum(1 / k)
}

cumulative_top_left_sums <- function(crm) {
  sums <- crm + 0
  nr <- nrow(sums)
  nc <- ncol(sums)

  if (nc > 1L) {
    for (j in 2:nc) {
      sums[, j] <- sums[, j] + sums[, j - 1L]
    }
  }
  if (nr > 1L) {
    for (i in 2:nr) {
      sums[i, ] <- sums[i, ] + sums[i - 1L, ]
    }
  }
  sums
}

# Rescaled Agreement Between K-ary Neighborhoods (RNX)
#
# RNX is a scaled version of QNX which measures the agreement between two
# embeddings in terms of the shared number of k-nearest neighbors for each
# observation. RNX gives a value of 1 if the neighbors are all preserved
# perfectly and a value of 0 for a random embedding.
#
# @param crm Co-ranking matrix. Create from a pair of distance matrices with
# \code{coranking_matrix}.
# @param k Neighborhood size.
# @return RNX for \code{k}.
# @references
# Lee, J. A., Renard, E., Bernard, G., Dupont, P., & Verleysen, M. (2013).
# Type 1 and 2 mixtures of Kullback-Leibler divergences as cost functions in
# dimensionality reduction based on similarity preservation.
# \emph{Neurocomputing}, \emph{112}, 92-108.
rnx_crm <- function(crm, k) {
  n_ranks <- nrow(crm)
  ((qnx_crm(crm, k) * n_ranks) - k) / (n_ranks - k)
}

# Average Normalized Agreement Between K-ary Neighborhoods (QNX)
#
# QNX measures the degree to which an embedding preserves the local
# neighborhood around each observation. For a value of K, the K closest
# neighbors of each observation are retrieved in the input and output space.
# For each observation, the number of shared neighbors can vary between 0
# and K. QNX is simply the average value of the number of shared neighbors,
# normalized by K, so that if the neighborhoods are perfectly preserved, QNX
# is 1, and if there is no neighborhood preservation, QNX is 0.
#
# For a random embedding, the expected value of QNX is approximately
# K / (N - 1) where N is the number of observations. Using RNX
# (\code{rnx_crm}) removes this dependency on K and the number of
# observations.
#
# @param crm Co-ranking matrix. Create from a pair of distance matrices with
# \code{coranking_matrix}.
# @param k Neighborhood size.
# @return QNX for \code{k}.
# @references
# Lee, J. A., & Verleysen, M. (2009).
# Quality assessment of dimensionality reduction: Rank-based criteria.
# \emph{Neurocomputing}, \emph{72(7)}, 1431-1443.
qnx_crm <- function(crm, k) {
  n_obs <- nrow(crm) + 1L
  sum(crm[1:k, 1:k]) / (k * n_obs)
}


# Intrusions and Extrusions for K-ary Neighborhoods (BNX)
#
# BNX measures the degree of intrusions versus extrusions that contributes
# to the QNX measure of embedding error. If BNX > 0 this means that intrusions
# dominate over extrusions: i.e. non-neighbors in the input space are neighbors
# in the output space. BNX < 0 means that extrusions dominate over intrusions:
# neighbors in the input space tend to be non-neighbors in the output space.
#
# @param crm Co-ranking matrix. Create from a pair of distance matrices with
# \code{coranking_matrix}.
# @param k Neighborhood size.
# @return BNX for \code{k}.
bnx_crm <- function(crm, k) {
  kcrm <- crm[1:k, 1:k]
  intrusions <- sum(kcrm[lower.tri(kcrm)])
  extrusions <- sum(kcrm[upper.tri(kcrm)])
  (intrusions - extrusions) / (k * (nrow(crm) + 1L))
}

# Indexes of the k-largest numbers.
#
# Given a vector of numbers, return the indexes of the k-largest
# values.
#
# @param x Vector of numbers.
# @param k Top k results to return
# @return Vector of the indexes of the \code{k} largest values in \code{x}.
k_largest_ind <- function(x, k) {
  which(x >= sort(x, decreasing = TRUE)[k])
}

# Indexes of the k-smallest numbers in a vector.
#
# Given a vector of numbers, return the indexes of the k-smallest
# values.
#
# @param x Vector of numbers.
# @param k Top k results to return
# @return Vector of the indexes of the \code{k} smallest values in \code{x}.
k_smallest_ind <- function(x, k) {
  k_largest_ind(-x, k)
}

# Indexes of the shared neighbors between two distance vectors
#
# Return the indexes of shared k-closest neighbors in two lists of distances.
#
# @param di list of distances
# @param dj list of distances
# @param k The size of the shared neighborhood
# @return Vector of the indexes of the elements which are among both the
# \code{k}-smallest values of \code{di} and the \code{k}-smallest
# values of \code{dj}. If there aren't exactly k values (i.e. because of ties),
# more than k results will be returned.
k_shared_nbrs_ind <- function(di, dj, k) {
  nindi <- k_smallest_ind(di, k)
  nindj <- k_smallest_ind(dj, k)

  Reduce(intersect, list(nindi, nindj))
}

# Neighborhood Preservation
#
# For the K nearest neighbors in one set of distances, returns the number of
# those neighbors which are also K nearest neighbors in another list,
# normalized with respect to K.
#
# The neighborhood preservation can vary between 0 (no neighbors in common)
# and 1 (perfect preservation). With self-neighbors excluded, random performance
# gives an approximate value of K / (N - 1), where N is the number of distances.
#
# @param di Vector of distances.
# @param dj Vector of distances.
# @param k Size of the neighborhood to consider.
# @return The number of shared neighbors in the equivalent neighbor lists of
# \code{di} and \code{dj}.
nbr_pres_i <- function(di, dj, k) {
  base::min(k, length(k_shared_nbrs_ind(di, dj, k))) / k
}

validate_distance_matrix_pair <- function(din, dout) {
  if (!is.matrix(din) || !is.numeric(din)) {
    stop("din must be a numeric matrix", call. = FALSE)
  }
  if (!is.matrix(dout) || !is.numeric(dout)) {
    stop("dout must be a numeric matrix", call. = FALSE)
  }
  if (!all(dim(din) == dim(dout))) {
    stop("din and dout must have the same dimensions", call. = FALSE)
  }
  if (nrow(din) != ncol(din)) {
    stop("din and dout must be square distance matrices", call. = FALSE)
  }
  invisible(TRUE)
}

validate_rank_penalty_k <- function(k, n_obs) {
  k <- validate_positive_integer(k, "k")
  if (n_obs < 3L) {
    stop(
      "trustworthiness and continuity require at least three observations",
      call. = FALSE
    )
  }
  if ((2L * k) >= n_obs) {
    stop(
      "k must be less than half the number of observations",
      call. = FALSE
    )
  }
  k
}


distance2_matrix <- function(X) {
  X <- as.matrix(X)
  D2 <- rowSums(X * X)
  D2 <- D2 + sweep(-2 * X %*% t(X), 2, -t(D2))
  D2[D2 < 0] <- 0
  D2
}

distance_matrix <- function(X) {
  sqrt(distance2_matrix(X))
}
