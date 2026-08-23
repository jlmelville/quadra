#' Random Pair Distance Correlation
#'
#' Correlates distances for sampled pairs of observations in `Xin` and `Xout`.
#' Pearson correlation measures linear association; Spearman correlation
#' measures rank association.
#'
#' `Xin` and `metric_in` define the reference geometry. Reset the R seed and
#' keep `n_threads` fixed to reuse the same pairs across calls.
#'
#' @param Xin Input data, with observations in rows by default. Data must be
#'   dense and finite; nonnumeric data-frame columns are ignored.
#' @param Xout Output data, with the same input conventions and number of
#'   observations as `Xin`.
#' @param n_pairs Number of pairs to sample.
#' @param metric_in Distance metric for `Xin`. One of
#'   `"euclidean"`, `"sqeuclidean"` (squared Euclidean), `"cosine"`, `"manhattan"`,
#'   `"correlation"` (1 minus the Pearson correlation), or `"hamming"`.
#' @param metric_out Distance metric for `Xout`. See `metric_in` for
#'   details.
#' @param method Correlation method, either `"pearson"` or `"spearman"`.
#' @param is_transposed Whether observations are stored in columns rather than
#'   rows.
#' @param n_threads Maximum number of threads to use. `0` or `1` runs
#'   serially.
#' @return Correlation between the sampled input and output distances.
#' @references Becht, E., McInnes, L., Healy, J., Dutertre, C. A., Kwok, I. W.,
#' Ng, L. G., ... & Newell, E. W. (2019).
#' Dimensionality reduction for visualizing single-cell data using UMAP.
#' *Nature biotechnology*, *37*(1), 38-44.
#' @seealso [random_pair_distance_emd()], [random_pair_distance_stress()], and
#'   [random_triplet_accuracy()] for other global preservation measures.
#' @examples
#' iris_x <- iris[, -5]
#' iris_pca2 <- stats::prcomp(iris_x, rank. = 2, scale = FALSE, retx = TRUE)$x
#' random_pair_distance_correlation(iris_x, iris_pca2)
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
  is_transposed <- validate_scalar_logical(is_transposed, "is_transposed")
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

#' Earth Mover's Distance Between Random-Pair Distances
#'
#' Compares the empirical distributions of sampled input and output distances
#' using Earth Mover's Distance. Each distance vector is scaled to `[0, 1]` by
#' default.
#'
#' EMD compares marginal distributions, not corresponding pairs. `Xin` and
#' `metric_in` define the reference geometry. Reset the R seed and keep
#' `n_threads` fixed to reuse the same pairs across calls.
#'
#' @param Xin Input data, with observations in rows by default. Data must be
#'   dense and finite; nonnumeric data-frame columns are ignored.
#' @param Xout Output data, with the same input conventions and number of
#'   observations as `Xin`.
#' @param n_pairs Number of pairs to sample.
#' @param metric_in Distance metric for `Xin`. One of
#'   `"euclidean"`, `"sqeuclidean"` (squared Euclidean), `"cosine"`, `"manhattan"`,
#'   `"correlation"` (1 minus the Pearson correlation), or `"hamming"`.
#' @param metric_out Distance metric for `Xout`. See `metric_in` for
#'   details.
#' @param range_scale Whether to scale each distance vector to `[0, 1]`.
#' @param is_transposed Whether observations are stored in columns rather than
#'   rows.
#' @param n_threads Maximum number of threads to use. `0` or `1` runs
#'   serially.
#' @return Earth Mover's Distance between the sampled distance distributions.
#' @references Heiser, C. N., & Lau, K. S. (2020).
#' A quantitative framework for evaluating single-cell data structure preservation by dimensionality reduction techniques.
#' *Cell reports*, *31*(5), 107576.
#' <https://github.com/KenLauLab/DR-structure-preservation>
#' @seealso [random_pair_distance_correlation()],
#'   [random_pair_distance_stress()], and [random_triplet_accuracy()] for other
#'   global preservation measures.
#' @examples
#' iris_x <- iris[, -5]
#' iris_pca2 <- stats::prcomp(iris_x, rank. = 2, scale = FALSE, retx = TRUE)$x
#' random_pair_distance_emd(iris_x, iris_pca2)
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
  range_scale <- validate_scalar_logical(range_scale, "range_scale")
  is_transposed <- validate_scalar_logical(is_transposed, "is_transposed")
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
    x <- rescale_unit_interval(x)
    y <- rescale_unit_interval(y)
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
#' Returns the root mean squared difference between matched sampled input and
#' output distances. Each distance vector is scaled to `[0, 1]` by default.
#'
#' `Xin` and `metric_in` define the reference geometry. Reset the R seed and
#' keep `n_threads` fixed to reuse the same pairs across calls.
#'
#' @inheritParams random_pair_distance_emd
#' @param range_scale Whether to scale each sampled distance vector to `[0, 1]`.
#' @return Root mean squared difference between the sampled distances.
#' @seealso [random_pair_distance_correlation()],
#'   [random_pair_distance_emd()], and [random_triplet_accuracy()] for other
#'   measures of global structure preservation.
#' @examples
#' iris_x <- iris[, -5]
#' iris_pca2 <- stats::prcomp(iris_x, rank. = 2, scale = FALSE, retx = TRUE)$x
#' random_pair_distance_stress(iris_x, iris_pca2)
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
  range_scale <- validate_scalar_logical(range_scale, "range_scale")
  is_transposed <- validate_scalar_logical(is_transposed, "is_transposed")
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
    x <- rescale_unit_interval(x)
    y <- rescale_unit_interval(y)
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
  n_threads <- validate_n_threads(n_threads)
  metric_in <- validate_distance(metric_in)
  metric_out <- validate_distance(metric_out)

  Xin <- prepare_input_matrix(Xin, name = "Xin", require_finite = TRUE)
  Xout <- prepare_input_matrix(Xout, name = "Xout", require_finite = TRUE)
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

rescale_unit_interval <- function(x) {
  x <- x - min(x)
  xmax <- max(x)
  if (xmax == 0) {
    return(rep(0, length(x)))
  }
  x / xmax
}
