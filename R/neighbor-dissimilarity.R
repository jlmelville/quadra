#' Neighbor Dissimilarity Between Two Graphs
#'
#' Compare how many neighbors each pair of observations shares in two
#' directed nearest-neighbor graphs on the same observations.
#'
#' @param nn_in,nn_out Reference and output graphs with the same rows in the
#'   same order. Supply lists with a nearest-first `idx` matrix as returned by
#'   [nn_preservation()]. Distances are optional. Self indices are stripped.
#' @param k Number of non-self neighbors to use from each row.
#' @param ret_extra Return the separate scores for gains and losses and the
#'   row pairs whose shared-neighbor counts changed?
#' @details For each graph, divide the number of shared neighbors of each
#'   pair by `k`. Subtract the reference fractions from the output fractions.
#'   Compute the Frobenius norms of the positive and negative parts of this
#'   difference matrix. The larger norm is Neighbor Dissimilarity (ND).
#'   Both orders of each pair enter the norms. Zero means all shared counts
#'   agree, which can happen even when neighbor identities change.
#'
#'   ND depends on graph size and `k`. Compare embeddings on the same rows
#'   and at the same `k`.
#' @section Scaling:
#' Each input graph needs only `k` neighbors per row, but the sparse matrix
#' product represents every pair of rows sharing a neighbor. One neighbor
#' listed by all `N` rows creates `N * (N - 1) / 2` unordered pairs. Even the
#' scalar score can therefore need quadratic memory. `ret_extra = TRUE` also
#' returns a table of changed pairs. The calculation uses the `Matrix` package.
#' @return A scalar, with smaller values indicating closer shared-neighbor
#'   structure. With `ret_extra = TRUE`, a list containing `nd`, `gain` and
#'   `loss` (the nonnegative directional norms), and `pairs`. The data frame
#'   contains one row per changed unordered pair (`i < j`), with one-based
#'   row positions and neighbor counts `shared_in` and `shared_out`.
#' @references Fujiwara, T., Kuo, Y.-H., Ynnerman, A., and Ma, K.-L.
#'   (2023). Feature Learning for Nonlinear Dimensionality Reduction toward
#'   Maximal Extraction of Hidden Patterns, equation 3.
#'   <https://arxiv.org/abs/2206.13891>
#'   The authors implement the graph comparison in
#'   [FEALM](https://github.com/takanori-fujiwara/fealm). Quadra also
#'   reports which observation pairs gained or lost shared neighbors.
#' @seealso [mutual_neighbor_correlation()], [nn_preservation()]
#' @examplesIf requireNamespace("Matrix", quietly = TRUE)
#' x <- as.matrix(iris[, -5])
#' g <- nn_preservation(x, prcomp(x)$x[, 1:2], k = 5,
#'                      nn_method_in = "brute", ret_extra = TRUE)
#' neighbor_dissimilarity(g$nn_in, g$nn_out, k = 5)
#' @export
neighbor_dissimilarity <- function(nn_in, nn_out, k = 15, ret_extra = FALSE) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("neighbor_dissimilarity requires the Matrix package", call. = FALSE)
  }
  k <- validate_positive_integer(k, "k")
  ret_extra <- validate_scalar_logical(ret_extra, "ret_extra")
  input <- prepare_supplied_nn_graph(nn_in, k, name = "nn_in")$idx
  output <- prepare_supplied_nn_graph(nn_out, k, name = "nn_out")$idx
  if (nrow(input) != nrow(output)) {
    stop("nn_in and nn_out must have the same number of rows", call. = FALSE)
  }

  shared_counts <- function(idx) {
    n <- nrow(idx)
    adjacency <- Matrix::sparseMatrix(
      i = rep(seq_len(n), k),
      j = as.vector(idx),
      x = 1,
      dims = c(n, n)
    )
    Matrix::tcrossprod(adjacency)
  }
  reference <- shared_counts(input)
  delta <- Matrix::triu(shared_counts(output) - reference, k = 1L)
  changed <- Matrix::summary(Matrix::drop0(delta))
  gain <- sqrt(2 * sum(pmax(changed$x, 0)^2)) / k
  loss <- sqrt(2 * sum(pmin(changed$x, 0)^2)) / k
  nd <- max(gain, loss)
  if (!ret_extra) return(nd)

  before <- as.integer(reference[cbind(changed$i, changed$j)])
  pairs <- data.frame(
    i = changed$i,
    j = changed$j,
    shared_in = before,
    shared_out = before + as.integer(changed$x)
  )
  list(nd = nd, gain = gain, loss = loss, pairs = pairs)
}
