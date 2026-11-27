#' Sample Triplets for Class Angular Distortion
#'
#' The anchor comes from one class and two distinct endpoints come from
#' another. Draws are independent and may repeat. Each eligible triplet
#' has equal probability, so a directed class pair with `a` anchor members
#' and `b` endpoint members has weight proportional to `a * choose(b, 2)`.
#'
#' @param labels A vector without missing values, with one label per
#'   observation and at least two distinct labels. At least one class
#'   must have two members to supply endpoints.
#' @param n_triplets Total number of triplets to sample. The default is
#'   ten times the number of observations.
#' @return An integer matrix with one-based columns `anchor`,
#'   `endpoint1`, and `endpoint2`.
#' @seealso [class_angular_distortion_index()]
#' @examples
#' set.seed(42)
#' sample_class_triplets(c("A", "A", "A", "B", "B", "C"), 5)
#' @export
sample_class_triplets <- function(labels, n_triplets = NULL) {
  n_obs <- length(labels)
  validate_label_vector(labels, n_obs, "labels")
  groups <- split(seq_along(labels), match(labels, unique(labels)))
  sizes <- lengths(groups)
  weights <- (n_obs - sizes) * choose(sizes, 2)
  if (!any(weights > 0)) {
    stop(
      "Sampling requires two classes, with at least two members in one class",
      call. = FALSE
    )
  }
  if (is.null(n_triplets)) {
    n_triplets <- 10 * n_obs
  }
  n_triplets <- validate_positive_integer(n_triplets, "n_triplets")

  endpoint_classes <- sample.int(
    length(groups),
    n_triplets,
    replace = TRUE,
    prob = weights
  )
  by_class <- split(seq_len(n_triplets), endpoint_classes)
  # Place each class in a contiguous block so anchors can skip it directly.
  starts <- cumsum(sizes) - sizes
  ordered <- unlist(groups, use.names = FALSE)
  plan <- matrix(
    0L,
    n_triplets,
    3L,
    dimnames = list(NULL, c("anchor", "endpoint1", "endpoint2"))
  )
  for (class in names(by_class)) {
    rows <- by_class[[class]]
    b <- as.integer(class)
    count <- length(rows)
    size <- sizes[b]
    anchor <- sample.int(n_obs - size, count, replace = TRUE)
    anchor <- anchor + size * (anchor > starts[b])
    first <- sample.int(size, count, replace = TRUE)
    second <- sample.int(size - 1L, count, replace = TRUE)
    second <- second + (second >= first)
    plan[rows, ] <- cbind(
      ordered[anchor],
      groups[[b]][first],
      groups[[b]][second]
    )
  }
  plan
}

#' Class Angular Distortion Index
#'
#' Compare how a class appears when viewed from another class. For each
#' triplet, measure the cosine of the angle at the anchor between the
#' rays to two observations from a different class.
#'
#' @inheritParams sample_class_triplets
#' @param Xin Reference coordinates, as a dense numeric matrix or data frame
#'   with observations in rows.
#' @param Xout Embedding coordinates for the same observations in the same
#'   order. Nonnumeric data-frame columns are ignored in both inputs.
#' @param triplets Optional one-based plan such as [sample_class_triplets()].
#'   The three indices must be distinct, and the endpoints must share a
#'   class different from the anchor's. When supplied, `n_triplets` is ignored.
#' @param ret_extra Whether to return the plan, cosines, and directed
#'   class-pair summaries.
#' @details The index is the mean squared change in cosine divided by four,
#'   following the normalization in the authors' implementation. It ranges
#'   from 0 (agreement) to 1. A zero-length ray gives a cosine of zero.
#'   Reuse a triplet plan for paired comparisons of embeddings.
#'
#'   The two directions of a class pair can have different distortions.
#'   Uniform scaling, translation, and orthogonal transformations preserve
#'   the angles. Angle agreement emphasizes relative angular spread. For
#'   example, compact classes can subtend small angles under several
#'   arrangements of their centers.
#' @section Scaling:
#' Use the full coordinates and control the cost with `n_triplets`, the total
#' number of sampled triplets (default `10 * N`). The calculation processes
#' coordinates in blocks and stores the triplet plan and cosines. Check the
#' per-pair counts before interpreting a rare class direction.
#' @section Numerical failures:
#' Nonfinite intermediate calculations cause an error. Extreme coordinate
#' scales may affect accuracy even when the calculations remain finite.
#' @return The scalar index. With `ret_extra = TRUE`, returns a list with
#'   `cadi`, `triplets`, `cosine_in`, `cosine_out`, and `pairs`.
#'   `pairs` contains the observed directed class pairs, with columns
#'   `anchor_class`, `endpoint_class`, `n_triplets`, and `cadi`.
#'   Weighting the pair scores by their counts recovers the scalar.
#' @references Gunaratne, K., Kobourov, S., & Miller, J. (2026).
#'   Class Angular Distortion Index for Dimensionality Reduction.
#'   \url{https://arxiv.org/abs/2605.00637}
#' @seealso [label_trustworthiness_continuity()], [roc_auc()], [pr_auc()]
#' @examples
#' x <- as.matrix(iris[, -5])
#' y <- prcomp(x, rank. = 2)$x
#' set.seed(42)
#' triplets <- sample_class_triplets(iris$Species, 1000)
#' class_angular_distortion_index(x, y, iris$Species, triplets = triplets)
#' @export
class_angular_distortion_index <- function(
  Xin,
  Xout,
  labels,
  n_triplets = NULL,
  triplets = NULL,
  ret_extra = FALSE
) {
  ret_extra <- validate_scalar_logical(ret_extra, "ret_extra")
  Xin <- prepare_input_matrix(Xin, name = "Xin", require_finite = TRUE)
  Xout <- prepare_input_matrix(Xout, name = "Xout", require_finite = TRUE)
  n_obs <- nrow(Xin)
  if (nrow(Xout) != n_obs) {
    stop(
      "Xin and Xout must have the same number of observations",
      call. = FALSE
    )
  }
  validate_label_vector(labels, n_obs, "Xin")
  if (is.null(triplets)) {
    triplets <- sample_class_triplets(labels, n_triplets)
  } else {
    triplets <- validate_triplet_sample(triplets, n_obs)
  }
  classes <- unique(labels)
  codes <- match(labels, classes)
  anchor <- codes[triplets[, 1L]]
  endpoint <- codes[triplets[, 2L]]
  if (any(anchor == endpoint | endpoint != codes[triplets[, 3L]])) {
    stop(
      "Triplet endpoints must share a class different from the anchor's",
      call. = FALSE
    )
  }

  cosine_in <- class_triplet_cosines(Xin, triplets)
  cosine_out <- class_triplet_cosines(Xout, triplets)
  loss <- ((cosine_in - cosine_out) / 2)^2
  cadi <- mean(loss)
  if (!ret_extra) {
    return(cadi)
  }

  group <- paste(anchor, endpoint, sep = "/")
  totals <- rowsum(
    cbind(count = rep(1L, length(loss)), loss),
    group,
    reorder = FALSE
  )
  first <- match(rownames(totals), group)
  pairs <- data.frame(
    anchor_class = classes[anchor[first]],
    endpoint_class = classes[endpoint[first]],
    n_triplets = as.integer(totals[, "count"]),
    cadi = unname(totals[, "loss"] / totals[, "count"])
  )
  list(
    cadi = cadi,
    triplets = triplets,
    cosine_in = cosine_in,
    cosine_out = cosine_out,
    pairs = pairs
  )
}

class_triplet_cosines <- function(x, triplets) {
  if (is.integer(x)) {
    storage.mode(x) <- "double"
  }
  n <- nrow(triplets)
  # Limit the size of temporary displacement matrices for wide inputs.
  block_size <- max(1L, floor(1e6 / ncol(x)))
  cosines <- numeric(n)
  for (start in seq.int(1L, n, by = block_size)) {
    rows <- seq.int(start, min(n, start + block_size - 1L))
    anchor <- x[triplets[rows, 1L], , drop = FALSE]
    u <- x[triplets[rows, 2L], , drop = FALSE] - anchor
    v <- x[triplets[rows, 3L], , drop = FALSE] - anchor
    norm_u <- rowSums(u^2)
    norm_v <- rowSums(v^2)
    dot <- rowSums(u * v)
    if (any(!is.finite(norm_u) | !is.finite(norm_v) | !is.finite(dot))) {
      stop("Nonfinite computed triplet angle", call. = FALSE)
    }
    defined <- norm_u > 0 & norm_v > 0
    values <- (dot[defined] / sqrt(norm_u[defined])) / sqrt(norm_v[defined])
    cosines[rows[defined]] <- pmin(1, pmax(-1, values))
  }
  cosines
}
