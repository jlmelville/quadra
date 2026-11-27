#' Exact Rank Profiles from Euclidean Coordinates
#'
#' Compare several embeddings against the same candidate neighbors. For each
#' query, compute its reference ranks once, score every embedding, then
#' release the ranks before moving to the next query.
#'
#' @param Xin Dense reference coordinates with observations in rows.
#' @param embeddings A nonempty list of dense coordinate matrices or data
#'   frames with the same observations in the same order as `Xin`. Pass
#'   `list(y)` for a single embedding. Nonnumeric data-frame columns are
#'   ignored. All coordinate values must be finite.
#' @param k Unique neighborhood sizes, each positive and less than half the
#'   number of observations. This range supports both rank penalties.
#' @param mrre Logical. Include mean relative rank errors (MRRE) for input
#'   neighbors and output neighbors? Defaults to `FALSE`.
#' @param query_rows Optional vector of unique, one-based row positions in
#'   `Xin` to score, in the requested order. The default (`NULL`) scores all
#'   rows. Each query still ranks all other observations, including rows
#'   outside the selection.
#' @details Each row's candidate pool contains all other observations.
#'   Tied Euclidean distances follow observation order, as with
#'   [pointwise_rank_preservation()] and [pointwise_rnx()] on matrices from
#'   [stats::dist()]. RNX AUC uses every size from one through `n - 2`.
#'   Only the selected queries incur distance and ranking work. Temporary
#'   workspace grows linearly with the number of observations. Returned
#'   profiles use space proportional to scored queries, requested sizes, and
#'   embeddings.
#'
#' For MRRE, select the first `k` neighbors in one space. For each of these
#' neighbors, divide the absolute change in rank by its rank in the selecting
#' space. Sum and divide by `sum((n - 2 * r) / r)` over ranks 1 through `k`.
#' Identical rankings give zero and reversed rankings give one. For large
#' neighborhoods, other permutations can exceed one under this normalization.
#' @references Lee, J. A. and Verleysen, M. (2008).
#' Quality assessment of nonlinear dimensionality reduction based on K-ary
#' neighborhoods. <https://proceedings.mlr.press/v4/lee08a.html>
#' @return A list in embedding order, retaining the names of `embeddings`.
#'   Each entry contains pointwise `trustworthiness` and `continuity`
#'   matrices, the RNX `auc` vector, and the RNX `curve` matrix at the
#'   requested sizes. Rows retain the selected reference row names. Column
#'   means summarize the selected queries. With `mrre = TRUE`, each
#'   entry also contains `mrre_input` and `mrre_output` matrices at the
#'   requested sizes. Smaller errors indicate better rank preservation.
#' @examples
#' x <- as.matrix(iris[, -5])
#' y <- prcomp(x)$x[, 1:2]
#' profile <- exact_rank_profiles(x, list(PCA = y), k = c(5, 15), mrre = TRUE)
#' colMeans(profile$PCA$mrre_input)
#' spot <- exact_rank_profiles(x, list(PCA = y), k = 5, query_rows = c(80, 3))
#' spot$PCA$curve
#' @export
exact_rank_profiles <- function(
  Xin,
  embeddings,
  k,
  mrre = FALSE,
  query_rows = NULL
) {
  Xin <- prepare_input_matrix(Xin, name = "Xin", require_finite = TRUE)
  if (
    !is.list(embeddings) || is.data.frame(embeddings) || !length(embeddings)
  ) {
    stop(
      "embeddings must be a nonempty list of coordinate matrices",
      call. = FALSE
    )
  }
  mrre <- validate_scalar_logical(mrre, "mrre")
  k <- validate_rank_penalty_k(k, nrow(Xin))
  if (is.null(query_rows)) {
    query_rows <- seq_len(nrow(Xin))
  } else {
    if (!is.null(dim(query_rows))) {
      stop("query_rows must be a vector of row positions", call. = FALSE)
    }
    query_rows <- validate_positive_integer_vector(query_rows, "query_rows")
    if (any(query_rows > nrow(Xin))) {
      stop("query_rows contains a position outside Xin", call. = FALSE)
    }
  }
  embeddings <- lapply(
    embeddings,
    prepare_input_matrix,
    name = "embedding",
    require_finite = TRUE
  )
  result <- coordinate_rank_profiles(Xin, embeddings, k, mrre, query_rows)
  names(result) <- names(embeddings)
  rows <- rownames(Xin)[query_rows]
  lapply(result, function(profile) {
    dimnames(profile$trustworthiness) <- dimnames(profile$continuity) <- list(
      rows,
      paste0("k", k)
    )
    if (mrre) {
      dimnames(profile$mrre_input) <- dimnames(profile$trustworthiness)
      dimnames(profile$mrre_output) <- dimnames(profile$trustworthiness)
    }
    dimnames(profile$curve) <- list(rows, paste0("rnx", k))
    names(profile$auc) <- rows
    profile
  })
}
