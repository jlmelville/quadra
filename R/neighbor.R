# Neighborhood Retrieval Metrics

#' Neighborhood Preservation Between Distance Matrices
#'
#' Returns the overlap between the `k`-neighborhoods defined by two distance
#' matrices, separately for each observation.
#'
#' Diagonal entries are excluded. Ties at the `k`th distance are included, with
#' each score capped at 1.
#' Nonfinite entries are not rejected.
#'
#' @param din Reference distance matrix.
#' @param dout Distance matrix to compare with `din`.
#' @param k Neighborhood size.
#' @return Per-observation neighborhood overlap in `[0, 1]`.
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
    preservations[i] <- neighbor_preservation_for_row(di, dj, k)
  }
  preservations
}


#' Neighborhood Preservation Between Nearest Neighbor Matrices
#'
#' Returns the overlap between reference and comparison nearest-neighbor
#' matrices, separately for each observation.
#'
#' Rows contain distinct one-based indices in nearest-first order. Self-indices
#' are removed, so self-inclusive matrices need at least `k + 1` columns.
#' Supplied order resolves ties.
#'
#' @param kin Reference nearest-neighbor index matrix, with observations in rows
#'   and neighbors in nearest-first order.
#' @param kout Nearest-neighbor index matrix to compare with `kin`.
#' @param k Neighborhood size.
#' @param n_threads Maximum number of threads to use. `0` or `1` runs
#'   serially.
#' @return Per-observation neighborhood overlap in `[0, 1]`.
#' @export
nbr_pres_knn <- function(kin, kout, k = ncol(kin), n_threads = 0) {
  if (!methods::is(kin, "matrix")) {
    stop("kin must be a matrix", call. = FALSE)
  }
  if (!methods::is(kout, "matrix")) {
    stop("kout must be a matrix", call. = FALSE)
  }
  if (nrow(kin) != nrow(kout)) {
    stop("kin and kout must have the same number of rows", call. = FALSE)
  }
  k <- validate_positive_integer(k, "k")
  n_threads <- validate_n_threads(n_threads)
  kin <- prepare_supplied_nn_graph(
    list(idx = kin),
    k = k,
    name = "kin",
    warn_self = FALSE
  )$idx
  kout <- prepare_supplied_nn_graph(
    list(idx = kout),
    k = k,
    name = "kout",
    warn_self = FALSE
  )$idx
  counts <- neighbor_overlap_counts(kin, kout, k, n_threads)
  counts[, 1] * (1 / k)
}

#' Trustworthiness and Continuity Between Distance Matrices
#'
#' `trustworthiness()` penalizes observations that appear among the `k` nearest
#' neighbors in `dout` but have input-space rank greater than `k` in `din`.
#' `continuity()` applies the dual penalty to input-space neighbors that are no
#' longer among the `k` nearest neighbors in `dout`.
#'
#' Diagonal entries are excluded. Unlike [nbr_pres()], the penalty reflects how
#' far a misplaced neighbor falls beyond `k`.
#' Nonfinite entries are ranked rather than rejected; ties follow column order.
#'
#' @param din Reference distance matrix.
#' @param dout Distance matrix to compare with `din`.
#' @param k Neighborhood size. Must be less than half the number of
#'   observations.
#' @return Scalar rank-preservation score; 1 indicates no penalty.
#' @references
#' Venna, J., & Kaski, S. (2001). Neighborhood preservation in nonlinear
#' projection methods: An experimental study. In *Artificial Neural Networks -
#' ICANN 2001* (pp. 485-491).
#' @examples
#' iris_pca <- stats::prcomp(iris[, -5], rank. = 2, scale = FALSE, retx = TRUE)$x
#' din <- as.matrix(stats::dist(iris[, -5]))
#' dout <- as.matrix(stats::dist(iris_pca))
#' trustworthiness(din, dout, k = 15)
#' continuity(din, dout, k = 15)
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
#' Summarizes rank-based neighborhood agreement across neighborhood sizes,
#' weighting smaller neighborhoods more heavily. A value of 1 is perfect, 0 is
#' the random-neighborhood baseline, and values may be negative. Diagonal
#' entries are excluded. Nonfinite entries are ranked rather than rejected;
#' ties follow column order.
#'
#' @param din Reference distance matrix.
#' @param dout Distance matrix to compare with `din`.
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

# Co-ranking utilities.
#
# Co-ranking Matrix
#
# Calculates the co-ranking matrix for an embedding.
#
# The co-ranking matrix is the basic data structure used for calculating
# various quality metrics, such as \code{qnx_crm} and \code{rnx_crm}.
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
shared_neighbor_indices <- function(di, dj, k) {
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
neighbor_preservation_for_row <- function(di, dj, k) {
  base::min(k, length(shared_neighbor_indices(di, dj, k))) / k
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
