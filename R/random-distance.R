#' Random Pair Distance Correlation
#'
#' Evaluates the preservation of global structure of dimensionality reduction
#' results using the correlation coefficient between randomly selected
#' distances. The default Pearson correlation is similar to the method of Becht
#' and co-workers (2019).
#'
#' This function repeatedly samples random pairs of observations and calculates
#' the distance between the points in both the original data and the embedding
#' space. The correlation coefficient between the two sets of distances is
#' reported. Pearson correlation measures linear agreement in the sampled
#' distances, while Spearman correlation measures rank agreement. This differs
#' slightly from the procedure in the Becht paper which randomly samples a
#' subset of observations and then exhaustively calculates all pair-wise
#' distances within that subset.
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
#'   `"euclidean"`, `"sqeuclidean"` (squared Euclidean), `"cosine"`, `"manhattan"`,
#'   `"correlation"` (1 minus the Pearson correlation), or `"hamming"`.
#' @param metric_out the distance metric to apply to `Xout`. See `metric_in` for
#'   details.
#' @param method correlation method, either `"pearson"` or `"spearman"`.
#' @param is_transposed if `TRUE` then `Xin` and `Xout` are assumed to have been
#'   passed in transposed format, i.e. with one observation per column.
#'   Otherwise, `Xin` and `Xout` will be transposed. For large datasets,
#'   transposing can be slow, so if this function will be called multiple times
#'   with the same input data, it is more efficient to transpose the input data
#'   once outside of this function and set `is_transposed = TRUE`.
#' @param n_threads the maximum number of threads to use. `0` or `1` runs
#'   serially.
#' @return The correlation between the distances in the input and output space.
#'   For randomly distributed data, the expected value is 0.
#' @references Becht, E., McInnes, L., Healy, J., Dutertre, C. A., Kwok, I. W.,
#' Ng, L. G., ... & Newell, E. W. (2019).
#' Dimensionality reduction for visualizing single-cell data using UMAP.
#' *Nature biotechnology*, *37*(1), 38-44.
#' @seealso [random_pair_distance_emd()], [random_pair_distance_stress()], and
#'   [random_triplet_accuracy()] for
#'   another measure of global structure preservation.
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
random_pair_distance_correlation <- function(
  Xin,
  Xout,
  n_pairs = 1000,
  metric_in = "sqeuclidean",
  metric_out = "sqeuclidean",
  method = c("pearson", "spearman"),
  is_transposed = FALSE,
  n_threads = 0
) {
  method <- match.arg(method)
  randlist <-
    random_pair_distances(
      Xin,
      Xout,
      n_pairs = n_pairs,
      metric_in = metric_in,
      metric_out = metric_out,
      n_threads = n_threads,
      is_transposed = is_transposed
    )

  stats::cor(
    x = randlist$din,
    y = randlist$dout,
    method = method
  )
}

#' Random Pair Distance Earth Mover's Distance
#'
#' Evaluates the preservation of global structure of dimensionality reduction
#' results using the Earth Mover's Distances between the distribution of
#' randomly selected pairs of distances, similar to the method of Heiser and Lau
#' (2020).
#'
#' This function repeatedly samples random pairs of observations and calculates
#' the distance between the points in both the original data and the embedding
#' space. Each set of distances are scaled between `(0, 1)`, converted into an
#' empirical distribution and the Earth Mover's (or Wasserstein) distance is
#' calculated.
#'
#' This function differs slightly from the procedure in the Heiser and Lau paper
#' which considers all pair-wise distances.
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
#'   `"euclidean"`, `"sqeuclidean"` (squared Euclidean), `"cosine"`, `"manhattan"`,
#'   `"correlation"` (1 minus the Pearson correlation), or `"hamming"`.
#' @param metric_out the distance metric to apply to `Xout`. See `metric_in` for
#'   details.
#' @param range_scale if `TRUE` (the default) then scale each distance vector
#' to the range 0-1 before converting to a distribution.
#' @param is_transposed if `TRUE` then `Xin` and `Xout` are assumed to have been
#'   passed in transposed format, i.e. with one observation per column.
#'   Otherwise, `Xin` and `Xout` will be transposed. For large datasets,
#'   transposing can be slow, so if this function will be called multiple times
#'   with the same input data, it is more efficient to transpose the input data
#'   once outside of this function and set `is_transposed = TRUE`.
#' @param n_threads the maximum number of threads to use. `0` or `1` runs
#'   serially.
#' @return The Earth Mover's distance between the empirical distributions formed
#'   from the distances in the input and output space.
#' @references Heiser, C. N., & Lau, K. S. (2020).
#' A quantitative framework for evaluating single-cell data structure preservation by dimensionality reduction techniques.
#' *Cell reports*, *31*(5), 107576.
#' <https://github.com/KenLauLab/DR-structure-preservation>
#' @seealso [random_pair_distance_correlation()],
#'   [random_pair_distance_stress()], and [random_triplet_accuracy()] for
#'   another measure of global structure preservation.
#' @examples
#' iris_pca2 <- stats::prcomp(iris[, -5], rank. = 2, scale = FALSE, retx = TRUE)$x
#' random_pair_distance_emd(iris, iris_pca2)
#'
#' # If you plan on comparing the results of multiple output methods, then
#' # pre-transposing the input data can save time
#' tiris <- t(iris[, -5])
#' iris_pca1 <- stats::prcomp(iris[, -5], rank. = 1, scale = FALSE, retx = TRUE)$x
#' iris_pca3 <- stats::prcomp(iris[, -5], rank. = 3, scale = FALSE, retx = TRUE)$x
#' random_pair_distance_emd(tiris, t(iris_pca1), is_transposed = TRUE)
#' random_pair_distance_emd(tiris, t(iris_pca2), is_transposed = TRUE)
#' random_pair_distance_emd(tiris, t(iris_pca3), is_transposed = TRUE)
#' @export
random_pair_distance_emd <- function(
  Xin,
  Xout,
  n_pairs = 1000,
  metric_in = "sqeuclidean",
  metric_out = "sqeuclidean",
  range_scale = TRUE,
  is_transposed = FALSE,
  n_threads = 0
) {
  randlist <-
    random_pair_distances(
      Xin,
      Xout,
      n_pairs = n_pairs,
      metric_in = metric_in,
      metric_out = metric_out,
      n_threads = n_threads,
      is_transposed = is_transposed
    )

  x <- randlist$din
  y <- randlist$dout
  if (range_scale) {
    x <- scale01(x)
    y <- scale01(y)
  }

  emd(
    x = x,
    y = y
  )
}

# Earth-Mover's distance (equivalent to 1D Wasserstein with p = 1)
emd <- function(x, y) {
  mean(abs(sort(x) - sort(y)))
}

#' Random Pair Distance Stress
#'
#' Evaluates global distance preservation with a sampled stress summary.
#'
#' This function repeatedly samples random pairs of observations and calculates
#' the distance between the points in both the original data and the embedding
#' space. The returned value is the root mean squared difference between the
#' matched sampled distances.
#'
#' By default, each sampled distance vector is scaled to the range 0-1 before
#' stress is calculated. This makes the result comparable across embeddings with
#' different distance scales, but it also means the value is mainly useful for
#' comparing methods under identical sampling and scaling settings.
#'
#' @inheritParams random_pair_distance_emd
#' @param range_scale if `TRUE` (the default) then scale each sampled distance
#'   vector to the range 0-1 before calculating stress.
#' @return The sampled stress between the matched distances in the input and
#'   output space.
#' @seealso [random_pair_distance_correlation()],
#'   [random_pair_distance_emd()], and [random_triplet_accuracy()] for other
#'   measures of global structure preservation.
#' @examples
#' iris_pca2 <- stats::prcomp(iris[, -5], rank. = 2, scale = FALSE, retx = TRUE)$x
#' random_pair_distance_stress(iris, iris_pca2)
#'
#' # If you plan on comparing the results of multiple output methods, then
#' # pre-transposing the input data can save time
#' tiris <- t(iris[, -5])
#' iris_pca1 <- stats::prcomp(iris[, -5], rank. = 1, scale = FALSE, retx = TRUE)$x
#' iris_pca3 <- stats::prcomp(iris[, -5], rank. = 3, scale = FALSE, retx = TRUE)$x
#' random_pair_distance_stress(tiris, t(iris_pca1), is_transposed = TRUE)
#' random_pair_distance_stress(tiris, t(iris_pca2), is_transposed = TRUE)
#' random_pair_distance_stress(tiris, t(iris_pca3), is_transposed = TRUE)
#' @export
random_pair_distance_stress <- function(
  Xin,
  Xout,
  n_pairs = 1000,
  metric_in = "sqeuclidean",
  metric_out = "sqeuclidean",
  range_scale = TRUE,
  is_transposed = FALSE,
  n_threads = 0
) {
  randlist <-
    random_pair_distances(
      Xin,
      Xout,
      n_pairs = n_pairs,
      metric_in = metric_in,
      metric_out = metric_out,
      n_threads = n_threads,
      is_transposed = is_transposed
    )

  x <- randlist$din
  y <- randlist$dout
  if (range_scale) {
    x <- scale01(x)
    y <- scale01(y)
  }

  sqrt(mean((x - y)^2))
}

random_pair_distances <- function(
  Xin,
  Xout,
  n_pairs = 1000,
  metric_in = "sqeuclidean",
  metric_out = "sqeuclidean",
  is_transposed = FALSE,
  n_threads = 0
) {
  n_pairs <- validate_positive_integer(n_pairs, "n_pairs")
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
    stop(
      "Xin and Xout must have the same number of observations",
      call. = FALSE
    )
  }
  if (n_obs < 2) {
    stop(
      "Xin and Xout must contain at least 2 observations for random-pair distances",
      call. = FALSE
    )
  }
  random_distances(
    Xin,
    Xout,
    n_pairs = n_pairs,
    metric_in = metric_in,
    metric_out = metric_out,
    n_threads = n_threads
  )
}

scale01 <- function(x) {
  x <- x - min(x)
  xmax <- max(x)
  if (xmax == 0) {
    return(rep(0, length(x)))
  }
  x / xmax
}
