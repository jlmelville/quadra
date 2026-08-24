#' Sample Observation Pairs
#'
#' Generates pairs for reuse with the random-pair distance metrics.
#'
#' @param n_obs Number of observations.
#' @param n_pairs Number of pairs to sample.
#' @return An integer matrix with one row per pair and columns `endpoint1` and
#'   `endpoint2`.
#' @seealso [random_pair_distance_correlation()], [random_pair_distance_emd()],
#'   and [random_pair_distance_stress()].
#' @examples
#' set.seed(42)
#' sample_pairs(5, 3)
#' @export
sample_pairs <- function(n_obs, n_pairs = 1000) {
  n_obs <- validate_positive_integer(n_obs, "n_obs")
  if (n_obs < 2L) {
    stop("n_obs must describe at least 2 observations", call. = FALSE)
  }
  n_pairs <- validate_positive_integer(n_pairs, "n_pairs")

  endpoint1 <- sample.int(n_obs, n_pairs, replace = TRUE)
  endpoint2 <- sample.int(n_obs - 1L, n_pairs, replace = TRUE)
  endpoint2 <- endpoint2 + (endpoint2 >= endpoint1)
  cbind(endpoint1 = endpoint1, endpoint2 = endpoint2)
}

#' Random Pair Distance Correlation
#'
#' Correlates distances for sampled pairs of observations in `Xin` and `Xout`.
#' Pearson correlation measures linear association; Spearman correlation
#' measures rank association.
#'
#' `Xin` and `metric_in` define the reference geometry. Supply `pairs` to reuse
#' exact comparisons, or reset the R seed and keep `n_threads` fixed.
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
#' @param pairs Optional one-based matrix with one pair per row. The endpoints
#'   in each row must differ. If supplied, `n_pairs` is ignored.
#' @param ret_extra Whether to return the pairs and their raw distances.
#' @return The correlation, or a list with `correlation`, `pairs`,
#'   `distance_in`, and `distance_out` when `ret_extra = TRUE`.
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
  n_threads = 0,
  pairs = NULL,
  ret_extra = FALSE
) {
  method <- match.arg(method)
  is_transposed <- validate_scalar_logical(is_transposed, "is_transposed")
  ret_extra <- validate_scalar_logical(ret_extra, "ret_extra")
  randlist <-
    random_pair_distances(
      Xin,
      Xout,
      n_pairs = n_pairs,
      metric_in = metric_in,
      metric_out = metric_out,
      n_threads = n_threads,
      is_transposed = is_transposed,
      pairs = pairs,
      ret_pairs = ret_extra
    )

  correlation <- stats::cor(
    x = randlist$din,
    y = randlist$dout,
    method = method
  )
  random_pair_result(correlation, "correlation", randlist, ret_extra)
}

#' Earth Mover's Distance Between Random-Pair Distances
#'
#' Compares the empirical distributions of sampled input and output distances
#' using Earth Mover's Distance. Each distance vector is scaled to `[0, 1]` by
#' default.
#'
#' EMD compares marginal distributions, not corresponding pairs. `Xin` and
#' `metric_in` define the reference geometry. Supply `pairs` to reuse exact
#' comparisons, or reset the R seed and keep `n_threads` fixed.
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
#' @param pairs Optional one-based matrix with one pair per row. The endpoints
#'   in each row must differ. If supplied, `n_pairs` is ignored.
#' @param ret_extra Whether to return the pairs and their raw distances.
#' @return The EMD, or a list with `emd`, `pairs`, `distance_in`, and
#'   `distance_out` when `ret_extra = TRUE`.
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
  n_threads = 0,
  pairs = NULL,
  ret_extra = FALSE
) {
  range_scale <- validate_scalar_logical(range_scale, "range_scale")
  is_transposed <- validate_scalar_logical(is_transposed, "is_transposed")
  ret_extra <- validate_scalar_logical(ret_extra, "ret_extra")
  randlist <-
    random_pair_distances(
      Xin,
      Xout,
      n_pairs = n_pairs,
      metric_in = metric_in,
      metric_out = metric_out,
      n_threads = n_threads,
      is_transposed = is_transposed,
      pairs = pairs,
      ret_pairs = ret_extra
    )

  x <- randlist$din
  y <- randlist$dout
  if (range_scale) {
    x <- rescale_unit_interval(x)
    y <- rescale_unit_interval(y)
  }

  value <- emd(
    x = x,
    y = y
  )
  random_pair_result(value, "emd", randlist, ret_extra)
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
#' `Xin` and `metric_in` define the reference geometry. Supply `pairs` to reuse
#' exact comparisons, or reset the R seed and keep `n_threads` fixed.
#'
#' @inheritParams random_pair_distance_emd
#' @param range_scale Whether to scale each sampled distance vector to `[0, 1]`.
#' @return The stress, or a list with `stress`, `pairs`, `distance_in`, and
#'   `distance_out` when `ret_extra = TRUE`.
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
  n_threads = 0,
  pairs = NULL,
  ret_extra = FALSE
) {
  range_scale <- validate_scalar_logical(range_scale, "range_scale")
  is_transposed <- validate_scalar_logical(is_transposed, "is_transposed")
  ret_extra <- validate_scalar_logical(ret_extra, "ret_extra")
  randlist <-
    random_pair_distances(
      Xin,
      Xout,
      n_pairs = n_pairs,
      metric_in = metric_in,
      metric_out = metric_out,
      n_threads = n_threads,
      is_transposed = is_transposed,
      pairs = pairs,
      ret_pairs = ret_extra
    )

  x <- randlist$din
  y <- randlist$dout
  if (range_scale) {
    x <- rescale_unit_interval(x)
    y <- rescale_unit_interval(y)
  }

  value <- sqrt(mean((x - y)^2))
  random_pair_result(value, "stress", randlist, ret_extra)
}

random_pair_distances <- function(
  Xin,
  Xout,
  n_pairs = 1000,
  metric_in = "sqeuclidean",
  metric_out = "sqeuclidean",
  is_transposed = FALSE,
  n_threads = 0,
  pairs = NULL,
  ret_pairs = FALSE
) {
  if (is.null(pairs)) {
    n_pairs <- validate_positive_integer(n_pairs, "n_pairs")
  } else {
    n_pairs <- 0L
  }
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
  if (!is.null(pairs)) {
    pairs <- validate_pair_sample(pairs, n_obs)
  }
  random_distances(
    Xin,
    Xout,
    n_pairs = n_pairs,
    metric_in = metric_in,
    metric_out = metric_out,
    n_threads = n_threads,
    pairs = pairs,
    ret_pairs = ret_pairs
  )
}

validate_pair_sample <- function(pairs, n_obs) {
  if (
    !is.matrix(pairs) ||
      !is.numeric(pairs) ||
      is.complex(pairs) ||
      nrow(pairs) < 1L ||
      ncol(pairs) != 2L
  ) {
    stop(
      "pairs must be a nonempty numeric matrix with 2 columns",
      call. = FALSE
    )
  }
  if (
    anyNA(pairs) ||
      any(!is.finite(pairs)) ||
      any(pairs != floor(pairs))
  ) {
    stop("pairs must contain finite whole-number indices", call. = FALSE)
  }
  if (any(pairs < 1 | pairs > n_obs)) {
    stop(
      "pairs indices must be between 1 and the number of observations",
      call. = FALSE
    )
  }
  if (any(pairs[, 1L] == pairs[, 2L])) {
    stop("pairs must contain distinct endpoints within each row", call. = FALSE)
  }

  matrix(
    as.integer(pairs),
    nrow = nrow(pairs),
    ncol = 2L,
    dimnames = list(NULL, c("endpoint1", "endpoint2"))
  )
}

random_pair_result <- function(value, name, randlist, ret_extra) {
  if (!ret_extra) {
    return(value)
  }

  result <- list(value, randlist$pairs, randlist$din, randlist$dout)
  names(result) <- c(name, "pairs", "distance_in", "distance_out")
  result
}

rescale_unit_interval <- function(x) {
  x <- x - min(x)
  xmax <- max(x)
  if (xmax == 0) {
    return(rep(0, length(x)))
  }
  x / xmax
}
