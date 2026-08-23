#' Random Triplet Accuracy
#'
#' Returns the proportion of sampled anchored triplets whose relative distance
#' ordering in `Xin` is preserved in `Xout`. Input-space ties are excluded; the
#' result is `NA_real_` if no triplet defines an ordering.
#'
#' `Xin` and `metric_in` define the reference geometry. Reset the R seed and
#' keep `n_threads` fixed to reuse the same triplets across calls. A
#' matrix-valued `n_triplets` bypasses sampling.
#'
#' @param Xin Input data, with observations in rows by default. Data must be
#'   dense and finite; nonnumeric data-frame columns are ignored.
#' @param Xout Output data, with the same input conventions and number of
#'   observations as `Xin`.
#' @param n_triplets Number of triplets to sample per observation, or a
#'   zero-based endpoint matrix. A matrix has one column per anchor and a pair of
#'   rows per triplet. Its endpoints must lie in `[0, n - 1]` and differ from
#'   each other and the anchor.
#' @param metric_in Distance metric for `Xin`. One of
#'   `"euclidean"`, `"sqeuclidean"` (squared Euclidean), `"cosine"`, `"manhattan"`,
#'   `"correlation"` (1 minus the Pearson correlation), or `"hamming"`.
#' @param metric_out Distance metric for `Xout`. See `metric_in` for
#'   details.
#' @param is_transposed Whether observations are stored in columns rather than
#'   rows.
#' @param n_threads Maximum number of threads to use. `0` or `1` runs
#'   serially.
#' @return Triplet accuracy in `[0, 1]`, or `NA_real_` if every input comparison
#'   is tied.
#' @references Wang, Y., Huang, H., Rudin, C., & Shaposhnik, Y. (2021).
#' Understanding how dimension reduction tools work: an empirical approach to
#' deciphering t-SNE, UMAP, TriMAP, and PaCMAP for data visualization.
#' *J Mach. Learn. Res*, *22*, 1-73. <https://jmlr.org/papers/v22/20-1061.html>.
#' @seealso [random_pair_distance_emd()] and
#'   [random_pair_distance_correlation()] for other global preservation
#'   measures.
#' @examples
#' iris_x <- iris[, -5]
#' iris_pca2 <- stats::prcomp(iris_x, rank. = 2, scale = FALSE, retx = TRUE)$x
#' random_triplet_accuracy(iris_x, iris_pca2)
#' @export
random_triplet_accuracy <-
  function(
    Xin,
    Xout,
    n_triplets = 5,
    metric_in = "sqeuclidean",
    metric_out = "sqeuclidean",
    is_transposed = FALSE,
    n_threads = 0
  ) {
    is_transposed <- validate_scalar_logical(is_transposed, "is_transposed")
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
    if (n_obs < 3) {
      stop(
        "Xin and Xout must contain at least 3 observations for triplet accuracy",
        call. = FALSE
      )
    }

    if (is.matrix(n_triplets)) {
      triplets <-
        validate_triplet_matrix(n_obs, n_triplets)

      return(
        triplet_sample(
          triplets,
          Xin,
          Xout,
          metric_in = metric_in,
          metric_out = metric_out,
          n_threads = n_threads
        )
      )
    }

    n_triplets <- validate_positive_integer(n_triplets, "n_triplets")
    random_triplet_sample(
      Xin,
      Xout,
      n_triplets = n_triplets,
      metric_in = metric_in,
      metric_out = metric_out,
      n_threads = n_threads
    )
  }

# Validate a pre-generated zero-indexed triplet matrix.
validate_triplet_matrix <- function(n_obs, triplets) {
  if (!is.numeric(triplets) || anyNA(triplets) || any(!is.finite(triplets))) {
    stop("Triplets matrix must contain finite numeric indices", call. = FALSE)
  }
  if (any(triplets != floor(triplets))) {
    stop("Triplets matrix must contain integer indices", call. = FALSE)
  }
  if (ncol(triplets) != n_obs) {
    stop("Triplets matrix must have ", n_obs, " columns", call. = FALSE)
  }
  if (nrow(triplets) %% 2 != 0) {
    stop("Triplets matrix must have even number of rows", call. = FALSE)
  }
  if (nrow(triplets) == 0L) {
    stop("Triplets matrix must be nonempty", call. = FALSE)
  }
  if (min(triplets) < 0) {
    stop("Triplet matrix must have non-negative values", call. = FALSE)
  }
  max_trip_idx <- max(triplets)
  if (max_trip_idx > n_obs - 1) {
    stop(
      "Triplet matrix elements must be between 0 and ",
      n_obs - 1,
      call. = FALSE
    )
  }
  anchors <- rep(seq_len(n_obs) - 1L, each = nrow(triplets))
  if (any(triplets == anchors)) {
    stop(
      "Each triplet endpoint must be distinct from its anchor",
      call. = FALSE
    )
  }
  first_endpoints <- seq.int(1L, nrow(triplets), by = 2L)
  second_endpoints <- first_endpoints + 1L
  if (any(triplets[first_endpoints, ] == triplets[second_endpoints, ])) {
    stop("Each triplet must contain distinct endpoints", call. = FALSE)
  }
  triplets
}
