#' Distance-Band Fidelity for Random Pairs
#'
#' Compare the order and magnitude of matched distances within cumulative
#' bands of input-space distances. Each band includes pairs at or below its
#' input radius, including ties.
#'
#' @inheritParams random_pair_distance_correlation
#' @inheritSection random_pair_distance_correlation Numerical failures
#' @param bands Strictly increasing, finite band selectors. Proportions must
#'   be in `(0, 1]`; radii must be nonnegative.
#' @param band_type Interpret `bands` as empirical proportions or as radii in
#'   the units of `metric_in`. Proportions use the inverse empirical
#'   distribution of the evaluated input distances.
#' @param ret_extra Whether to return the profile and the evaluated pairs
#'   with their raw distances.
#' @details
#'   Both distance metrics default to ordinary Euclidean distance.
#'   When you supply all unordered pairs, the profile is exact for those
#'   observations. Repeated pair rows each count toward the result. With
#'   proportion bands, ties can make the realized fraction larger than the
#'   requested proportion. Set `bands = 1` in proportion mode to use every
#'   evaluated pair.
#'
#'   Each band reports Spearman correlation and normalized stress. Stress
#'   divides the residual Euclidean norm by the input distance vector's norm,
#'   after fitting one least-squares multiplier to the output distances
#'   within that band. A uniformly rescaled band has zero stress. Separate
#'   bands can fit different multipliers.
#'
#'   Spearman correlation needs at least three pairs and variation in both
#'   distance vectors. Stress needs nonzero norms in both vectors.
#'   Undefined values are `NA_real_` with one warning per call.
#' @return A data frame with columns `band_type`, `band`, `radius`,
#'   `n_pairs`, `pair_fraction`, `spearman`, and `normalized_stress`.
#'   If `ret_extra = TRUE`, a list containing that data frame as `profile`
#'   along with `pairs`, `distance_in`, `distance_out`, and `fitted_scale`.
#'   Multiply output distances by the corresponding entry of `fitted_scale`
#'   to put them in input-distance units. The multiplier is undefined (NA)
#'   for an empty band or all-zero output distances. A multiplier larger
#'   than the floating-point range is returned as `Inf`.
#' @seealso [sample_pairs()], [random_pair_distance_correlation()],
#'   [random_pair_distance_stress()]
#' @examples
#' x <- as.matrix(iris[, -5])
#' y <- prcomp(x, rank. = 2)$x
#' set.seed(42)
#' pairs <- sample_pairs(nrow(x), n_pairs = 2000)
#' random_pair_distance_band_profile(x, y, bands = c(0.2, 1), pairs = pairs)
#' @export
random_pair_distance_band_profile <- function(
  Xin,
  Xout,
  bands = c(0.05, 0.10, 0.20, 0.30, 0.50, 0.75, 1),
  band_type = c("proportion", "radius"),
  n_pairs = 10000,
  metric_in = "euclidean",
  metric_out = "euclidean",
  is_transposed = FALSE,
  n_threads = 0,
  pairs = NULL,
  ret_extra = FALSE
) {
  band_type <- match.arg(band_type)
  if (
    !is.numeric(bands) ||
      is.complex(bands) ||
      !is.null(dim(bands)) ||
      !length(bands) ||
      any(!is.finite(bands)) ||
      is.unsorted(bands, strictly = TRUE)
  ) {
    stop(
      "bands must be a strictly increasing finite numeric vector",
      call. = FALSE
    )
  }
  if (
    band_type == "proportion" && (bands[1L] <= 0 || bands[length(bands)] > 1)
  ) {
    stop("proportion bands must be in (0, 1]", call. = FALSE)
  }
  if (band_type == "radius" && bands[1L] < 0) {
    stop("radii must be nonnegative", call. = FALSE)
  }
  is_transposed <- validate_scalar_logical(is_transposed, "is_transposed")
  ret_extra <- validate_scalar_logical(ret_extra, "ret_extra")
  distances <- random_pair_distances(
    Xin,
    Xout,
    n_pairs = n_pairs,
    metric_in = metric_in,
    metric_out = metric_out,
    is_transposed = is_transposed,
    n_threads = n_threads,
    pairs = pairs,
    ret_pairs = ret_extra
  )

  radii <- bands
  if (band_type == "proportion") {
    radii <- stats::quantile(distances$din, bands, type = 1, names = FALSE)
  }
  summaries <- vapply(
    radii,
    function(radius) {
      inside <- distances$din <= radius
      distance_band_summary(distances$din[inside], distances$dout[inside])
    },
    numeric(4)
  )
  profile <- data.frame(
    band_type = band_type,
    band = as.numeric(bands),
    radius = as.numeric(radii),
    n_pairs = as.integer(summaries[1L, ]),
    pair_fraction = summaries[1L, ] / length(distances$din),
    spearman = summaries[2L, ],
    normalized_stress = summaries[3L, ]
  )
  if (anyNA(profile$spearman) || anyNA(profile$normalized_stress)) {
    warning(
      "Some band scores are undefined. Check pair counts and distance variation, or use wider bands.",
      call. = FALSE
    )
  }
  result <- random_pair_result(profile, "profile", distances, ret_extra)
  if (ret_extra) {
    result$fitted_scale <- summaries[4L, ]
  }
  result
}

distance_band_summary <- function(din, dout) {
  count <- length(din)
  spearman <- stress <- fitted_scale <- NA_real_
  if (count == 0L) {
    return(c(count, spearman, stress, fitted_scale))
  }
  if (count >= 3L && min(din) < max(din) && min(dout) < max(dout)) {
    spearman <- stats::cor(din, dout, method = "spearman")
  }

  scale_in <- max(abs(din))
  scale_out <- max(abs(dout))
  if (scale_in == 0 && scale_out > 0) {
    fitted_scale <- 0
  }
  if (scale_in > 0 && scale_out > 0) {
    # Normalizing before the fit keeps the sums of squares representable.
    din <- din / scale_in
    dout <- dout / scale_out
    multiplier <- sum(din * dout) / sum(dout^2)
    fitted_scale <- exp(log(multiplier) + (log(scale_in) - log(scale_out)))
    stress <- sqrt(sum((din - multiplier * dout)^2) / sum(din^2))
  }
  c(count, spearman, stress, fitted_scale)
}
