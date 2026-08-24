#' Sample Anchored Triplets
#'
#' Generates anchored triplets for reuse with [random_triplet_accuracy()].
#'
#' @param n_obs Number of observations.
#' @param n_triplets Number of triplets to sample per observation.
#' @return An integer matrix with columns `anchor`, `endpoint1`, and
#'   `endpoint2`, ordered by anchor.
#' @examples
#' set.seed(42)
#' sample_triplets(5, 2)
#' @export
sample_triplets <- function(n_obs, n_triplets = 5) {
  n_obs <- validate_positive_integer(n_obs, "n_obs")
  if (n_obs < 3L) {
    stop("n_obs must describe at least 3 observations", call. = FALSE)
  }
  n_triplets <- validate_positive_integer(n_triplets, "n_triplets")
  if (n_obs > .Machine$integer.max %/% n_triplets) {
    stop(
      "n_obs * n_triplets exceeds the supported integer range",
      call. = FALSE
    )
  }

  n_comparisons <- n_obs * n_triplets
  anchor <- rep(seq_len(n_obs), each = n_triplets)

  endpoint1 <- sample.int(n_obs - 1L, n_comparisons, replace = TRUE)
  endpoint1 <- endpoint1 + (endpoint1 >= anchor)

  endpoint2 <- sample.int(n_obs - 2L, n_comparisons, replace = TRUE)
  excluded1 <- pmin(anchor, endpoint1)
  excluded2 <- pmax(anchor, endpoint1)
  endpoint2 <- endpoint2 + (endpoint2 >= excluded1)
  endpoint2 <- endpoint2 + (endpoint2 >= excluded2)

  cbind(anchor = anchor, endpoint1 = endpoint1, endpoint2 = endpoint2)
}

#' Random Triplet Accuracy
#'
#' Returns the proportion of sampled anchored triplets whose relative distance
#' ordering in `Xin` is preserved in `Xout`. Input-space ties are excluded; the
#' result is `NA_real_` if no triplet defines an ordering.
#'
#' `Xin` and `metric_in` define the reference geometry. Supply `triplets` to
#' reuse exact comparisons, or reset the R seed and keep `n_threads` fixed. A
#' matrix-valued `n_triplets` uses its documented zero-based column layout.
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
#' @param triplets Optional one-based matrix with columns for the anchor and two
#'   endpoints. All three indices in a row must differ. If supplied,
#'   numeric `n_triplets` is ignored.
#' @param ret_extra Whether to return the triplets and their row-aligned
#'   agreement outcomes.
#' @return Triplet accuracy in `[0, 1]`, or `NA_real_` if every input comparison
#'   is tied. With `ret_extra = TRUE`, returns a list with `accuracy`,
#'   `triplets`, and `agreement`. Agreement is `NA` for input-distance ties,
#'   `TRUE` for matching strict orderings, and `FALSE` otherwise.
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
    n_threads = 0,
    triplets = NULL,
    ret_extra = FALSE
  ) {
    is_transposed <- validate_scalar_logical(is_transposed, "is_transposed")
    ret_extra <- validate_scalar_logical(ret_extra, "ret_extra")
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

    if (!is.null(triplets)) {
      if (is.matrix(n_triplets)) {
        stop(
          "only one of n_triplets and triplets may supply a triplet sample",
          call. = FALSE
        )
      }
      triplets <- validate_triplet_sample(triplets, n_obs)
      return(
        triplet_plan_sample(
          triplets,
          Xin,
          Xout,
          metric_in = metric_in,
          metric_out = metric_out,
          n_threads = n_threads,
          ret_triplets = ret_extra
        )
      )
    }

    if (is.matrix(n_triplets)) {
      triplet_matrix <- validate_triplet_matrix(n_obs, n_triplets)
      if (!ret_extra) {
        return(
          triplet_sample(
            triplet_matrix,
            Xin,
            Xout,
            metric_in = metric_in,
            metric_out = metric_out,
            n_threads = n_threads
          )
        )
      }
      triplets <- canonicalize_triplet_matrix(n_obs, triplet_matrix)
      return(
        triplet_plan_sample(
          triplets,
          Xin,
          Xout,
          metric_in = metric_in,
          metric_out = metric_out,
          n_threads = n_threads,
          ret_triplets = ret_extra
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
      n_threads = n_threads,
      ret_triplets = ret_extra
    )
  }

validate_triplet_sample <- function(triplets, n_obs) {
  if (
    !is.matrix(triplets) ||
      !is.numeric(triplets) ||
      is.complex(triplets) ||
      nrow(triplets) < 1L ||
      ncol(triplets) != 3L
  ) {
    stop(
      "triplets must be a nonempty numeric matrix with 3 columns",
      call. = FALSE
    )
  }
  if (
    anyNA(triplets) ||
      any(!is.finite(triplets)) ||
      any(triplets != floor(triplets))
  ) {
    stop("triplets must contain finite whole-number indices", call. = FALSE)
  }
  if (any(triplets < 1 | triplets > n_obs)) {
    stop(
      "triplets indices must be between 1 and the number of observations",
      call. = FALSE
    )
  }
  if (
    any(triplets[, 1L] == triplets[, 2L]) ||
      any(triplets[, 1L] == triplets[, 3L]) ||
      any(triplets[, 2L] == triplets[, 3L])
  ) {
    stop("triplets must contain three distinct indices per row", call. = FALSE)
  }

  matrix(
    as.integer(triplets),
    nrow = nrow(triplets),
    ncol = 3L,
    dimnames = list(NULL, c("anchor", "endpoint1", "endpoint2"))
  )
}

canonicalize_triplet_matrix <- function(n_obs, triplets) {
  endpoint1_rows <- seq.int(1L, nrow(triplets), by = 2L)
  endpoint2_rows <- endpoint1_rows + 1L
  n_triplets <- length(endpoint1_rows)

  cbind(
    anchor = rep(seq_len(n_obs), each = n_triplets),
    endpoint1 = as.integer(triplets[endpoint1_rows, , drop = FALSE]) + 1L,
    endpoint2 = as.integer(triplets[endpoint2_rows, , drop = FALSE]) + 1L
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
