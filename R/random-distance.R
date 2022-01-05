#' Random Pair Distance Correlation
#'
#' Evaluates the preservation of global structure of dimensionality reduction
#' results using the Pearson correlation coefficient between randomly selected
#' distances, similar to the method of Becht and co-workers (2019).
#'
#' This function repeatedly samples random pairs of observation and calculates
#' the distance between the points in both the original data and the embedding
#' space. The Pearson correlation coefficient between the two sets of distances
#' is reported. This differs slightly from the procedure in the Becht paper
#' which randomly samples a subset of observations and then exhaustively
#' calculates all pair-wise distances within that subset.
#'
#' @param Xin the input data (usually high-dimensional), a matrix or data frame
#'   with one observation per row, or if `is_transposed = TRUE`, one observation
#'   per column.
#' @param Xout the output data (usually lower dimensional than `Xin`), a matrix
#'   or data frame with one observation per row, or if `is_transposed = TRUE`,
#'   one observation per column.
#' @param n_pairs the number of random pairs of observations to calculate
#'   distances for in the input and output space.
#' @param metric_in the distance calculation to apply to `Xin`. One of
#'   `"euclidean"`, `"l2sqr"` (squared Euclidean), `"cosine"`, `"manhattan"`,
#'   `"correlation"` (1 minus the Pearson correlation), or `"hamming"`.
#' @param metric_out the distance metric to apply to `Xout`. See `metric_in` for
#'   details.
#' @param is_transposed if `TRUE` then `Xin` and `Xout` are assumed to have been
#'   passed in transposed format, i.e. with one observation per column.
#'   Otherwise, `Xin` and `Xout` will be transposed. For large datasets,
#'   transposing can be slow, so if this function will be called multiple times
#'   with the same input data, it is more efficient to transpose the input data
#'   once outside of this function and set `is_transposed = TRUE`.
#' @param n_threads the maximum number of threads to use.
#' @return The Pearson correlation between the distances in the input and output
#' space.
#' @references Becht, E., McInnes, L., Healy, J., Dutertre, C. A., Kwok, I. W.,
#' Ng, L. G., ... & Newell, E. W. (2019).
#' Dimensionality reduction for visualizing single-cell data using UMAP.
#' *Nature biotechnology*, *37*(1), 38-44.
#' @seealso [random_triplet_accuracy()] for another measure of global structure
#' preservation.
#' @md
#' @examples
#' iris_pca2 <- stats::prcomp(iris[, -5], rank. = 2, scale = FALSE, retx = TRUE)$x
#' random_pair_distance_correlation(iris, iris_pca2)
#'
#' # If you plan on comparing the results of multiple output methods, then
#' # pre-transposing the input data can save time
#' tiris <- t(iris[, -5])
#' iris_pca1 <- stats::prcomp(iris[, -5], rank. = 1, scale = FALSE, retx = TRUE)$x
#' iris_pca3 <- stats::prcomp(iris[, -5], rank. = 3, scale = FALSE, retx = TRUE)$x
#' random_pair_distance_correlation(tiris, t(iris_pca1), is_transposed = TRUE)
#' random_pair_distance_correlation(tiris, t(iris_pca2), is_transposed = TRUE)
#' random_pair_distance_correlation(tiris, t(iris_pca3), is_transposed = TRUE)
#' @export
random_pair_distance_correlation <- function(Xin,
                                             Xout,
                                             n_pairs = 1000,
                                             metric_in = "l2sqr",
                                             metric_out = "l2sqr",
                                             is_transposed = FALSE,
                                             n_threads = 0) {
  metric_in <- validate_distance(metric_in)
  metric_out <- validate_distance(metric_out)

  Xin <- x2m(Xin)
  Xout <- x2m(Xout)
  if (!is_transposed) {
    Xin <- t(Xin)
    Xout <- t(Xout)
  }
  n_obs <- ncol(Xin)
  if (n_obs != ncol(Xout)) {
    stop("Xin and Xout must have the same number of observations")
  }
  randlist <-
    random_distances(
      Xin,
      Xout,
      n_pairs = n_pairs,
      metric_in = metric_in,
      metric_out = metric_out,
      n_threads = n_threads
    )

  stats::cor(
    x = randlist$din,
    y = randlist$dout,
    method = "pearson"
  )
}
