#' Local Procrustes Error
#'
#' Measure how closely each reference neighborhood retains its geometry
#' in an embedding after allowing translation, rotation, and reflection.
#' Each patch contains its anchor and the same `k` reference neighbors in
#' both spaces.
#'
#' @param Xin Reference coordinates, as a dense numeric matrix or data frame.
#' @param Xout Embedding coordinates for the same observations in the same
#'   order. Nonnumeric data-frame columns are ignored in both inputs.
#' @param k Number of reference neighbors per anchor. The anchor is included
#'   in addition to these neighbors.
#' @param nn_method_in Reference neighbor search method. Use `"brute"` for
#'   exact search or `"nnd"` for approximate neighbor descent.
#' @param nn_in Optional reference neighbor graph, such as `nn_in` from
#'   a detailed [nn_preservation()] result. Its first `k` non-self indices
#'   are used in each row. Supplying it skips neighbor search.
#' @inheritParams random_pair_distance_correlation
#' @param ret_extra Whether to also return the row-aligned residuals and graph.
#' @details Coordinates are centered within each patch. The alignment uses
#'   orthogonal least squares, with zero padding to compare different
#'   dimensions. Scale is fixed. To account for arbitrary embedding units,
#'   multiply the entire embedding by a single calibration factor before
#'   calling this function.
#'
#'   The residual is the sum of squared aligned coordinate differences
#'   divided by the reference patch's centered sum of squares. Zero indicates
#'   exact agreement up to numerical precision. Values can exceed one.
#'   A constant reference patch has an undefined residual (`NA_real_`).
#'   Moving an intact group can leave its internal patch scores unchanged.
#'
#'   Automatic neighbor search uses Euclidean ordering. Supplied graphs
#'   select the patch members in their given order. The Procrustes fit
#'   always uses the coordinate geometry of `Xin` and `Xout`.
#' @section Scaling:
#' Each alignment uses one patch of `k + 1` observations. With `k` fixed, the
#' number of fits grows linearly with `N`. The default brute-force reference
#' search compares all pairs. Supply `nn_in` to reuse a reference graph across
#' embeddings, or use `nn_method_in = "nnd"` to find approximate neighbors.
#' @return The mean residual over defined reference patches, or `NA_real_`
#'   when all are undefined. With `ret_extra = TRUE`, a list with `error`,
#'   `pointwise_error`, `n_defined`, and `nn_in`. The pointwise vector
#'   follows the rows and names of the reference observations.
#' @references Goldberg, Y., & Ritov, Y. (2009). Local Procrustes for Manifold
#'   Embedding: A Measure of Embedding Quality and Embedding Algorithms.
#'   <https://arxiv.org/abs/0806.2669>
#'   Quadra follows their construction by matching the same observations
#'   in both patches. The fit keeps scale fixed within each patch.
#'   [ZADU 0.5.4](https://github.com/hj-n/zadu/tree/v0.5.4) selects
#'   neighbors separately in each space.
#' @seealso [nn_preservation()], [local_radius_correlation()]
#' @examples
#' x <- as.matrix(iris[, -5])
#' embedding <- prcomp(x, rank. = 2)$x
#' result <- local_procrustes(x, embedding, k = 10, ret_extra = TRUE)
#' result$error
#' head(result$pointwise_error)
#' @export
local_procrustes <- function(
  Xin,
  Xout,
  k = 15,
  nn_method_in = c("brute", "nnd"),
  is_transposed = FALSE,
  n_threads = 0,
  nn_in = NULL,
  ret_extra = FALSE
) {
  k <- validate_positive_integer(k, "k")
  is_transposed <- validate_scalar_logical(is_transposed, "is_transposed")
  n_threads <- validate_n_threads(n_threads)
  ret_extra <- validate_scalar_logical(ret_extra, "ret_extra")
  Xin <- prepare_input_matrix(Xin, name = "Xin", require_finite = TRUE)
  Xout <- prepare_input_matrix(Xout, name = "Xout", require_finite = TRUE)
  if (is_transposed) {
    Xin <- t(Xin)
    Xout <- t(Xout)
  }
  n_obs <- nrow(Xin)
  if (nrow(Xout) != n_obs) {
    stop(
      "Xin and Xout must have the same number of observations",
      call. = FALSE
    )
  }
  check_k_for_n_obs(k, n_obs)

  if (is.null(nn_in)) {
    nn_in <- resolve_nn_graph(
      Xin,
      k = k,
      nn_method = match.arg(nn_method_in),
      metric = "euclidean",
      is_transposed = FALSE,
      n_threads = n_threads,
      verbose = FALSE,
      nn_args = list()
    )
  } else {
    nn_in <- prepare_supplied_nn_graph(nn_in, k, name = "nn_in")
    if (nrow(nn_in$idx) != n_obs) {
      stop("nn_in must have one row per reference observation", call. = FALSE)
    }
  }

  errors <- vapply(
    seq_len(n_obs),
    function(i) {
      ids <- c(i, nn_in$idx[i, ])
      procrustes_patch_error(
        Xin[ids, , drop = FALSE],
        Xout[ids, , drop = FALSE]
      )
    },
    numeric(1)
  )
  n_defined <- sum(!is.na(errors))
  error <- if (n_defined > 0L) mean(errors, na.rm = TRUE) else NA_real_
  if (!ret_extra) {
    return(error)
  }
  names(errors) <- rownames(Xin)
  list(
    error = error,
    pointwise_error = errors,
    n_defined = n_defined,
    nn_in = nn_in
  )
}

procrustes_patch_error <- function(a, b) {
  a <- scale(a, scale = FALSE)
  if (all(a == 0)) {
    return(NA_real_)
  }
  b <- scale(b, scale = FALSE)
  magnitude <- max(abs(a), abs(b))
  if (!is.finite(magnitude)) {
    stop(
      "Nonfinite centered coordinates in local Procrustes alignment",
      call. = FALSE
    )
  }

  a <- procrustes_patch_basis(a / magnitude)
  b <- procrustes_patch_basis(b / magnitude)
  width <- max(ncol(a), ncol(b))
  a <- cbind(a, matrix(0, nrow(a), width - ncol(a)))
  b <- cbind(b, matrix(0, nrow(b), width - ncol(b)))
  rotation <- svd(crossprod(b, a))
  aligned <- b %*% rotation$u %*% t(rotation$v)
  error <- sum((a - aligned)^2) / sum(a^2)
  if (!is.finite(error)) {
    stop("Nonfinite local Procrustes residual", call. = FALSE)
  }
  error
}

procrustes_patch_basis <- function(x) {
  if (ncol(x) <= nrow(x)) {
    return(x)
  }
  # Use an orthogonal basis for wide patches so the SVD stays patch-sized.
  decomposition <- qr(t(x), LAPACK = TRUE)
  t(qr.R(decomposition)[, order(decomposition$pivot), drop = FALSE])
}
