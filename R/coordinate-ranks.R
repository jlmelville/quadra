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
#' @details Each row's candidate pool contains all other observations.
#'   Tied Euclidean distances follow observation order, as with
#'   [pointwise_rank_preservation()] and [pointwise_rnx()] on matrices from
#'   [stats::dist()]. RNX AUC uses every size from one through `n - 2`.
#'   The method evaluates all query-candidate distances. Distance and rank
#'   workspace grows linearly with the number of observations. Returned
#'   profiles use space proportional to observations, requested sizes, and
#'   embeddings.
#' @return A list in embedding order, retaining the names of `embeddings`.
#'   Each entry contains pointwise `trustworthiness` and `continuity`
#'   matrices, the RNX `auc` vector, and the RNX `curve` matrix at the
#'   requested sizes. Rows retain the reference row names. Average over
#'   rows for the corresponding aggregate scores.
#' @examples
#' x <- as.matrix(iris[, -5])
#' y <- prcomp(x)$x[, 1:2]
#' profile <- exact_rank_profiles(x, list(PCA = y), k = c(5, 15))
#' colMeans(profile$PCA$trustworthiness)
#' @export
exact_rank_profiles <- function(Xin, embeddings, k) {
  Xin <- prepare_input_matrix(Xin, name = "Xin", require_finite = TRUE)
  if (
    !is.list(embeddings) || is.data.frame(embeddings) || !length(embeddings)
  ) {
    stop(
      "embeddings must be a nonempty list of coordinate matrices",
      call. = FALSE
    )
  }
  k <- validate_rank_penalty_k(k, nrow(Xin))
  embeddings <- lapply(
    embeddings,
    prepare_input_matrix,
    name = "embedding",
    require_finite = TRUE
  )
  result <- coordinate_rank_profiles(Xin, embeddings, k)
  names(result) <- names(embeddings)
  rows <- rownames(Xin)
  lapply(result, function(profile) {
    dimnames(profile$trustworthiness) <- dimnames(profile$continuity) <- list(
      rows,
      paste0("k", k)
    )
    dimnames(profile$curve) <- list(rows, paste0("rnx", k))
    names(profile$auc) <- rows
    profile
  })
}
