#' Label Trustworthiness and Continuity
#'
#' Compare how well each pair of classes separates around its centroids in
#' the reference coordinates and the embedding.
#'
#' @param Xin Reference coordinates, as a dense numeric matrix or data
#'   frame with observations in rows.
#' @param Xout Embedding coordinates for the same observations in the same
#'   order. Nonnumeric data-frame columns are ignored in both inputs.
#' @param labels A vector without missing values, with one label per row
#'   and at least two distinct labels.
#' @param ret_extra Whether to return the class-pair consistency scores and
#'   their signed changes.
#' @details For each pair of labels and in each space, compute the two
#'   class centroids. Distance consistency is the fraction of observations
#'   nearer their own centroid in Euclidean distance. A tie counts as half
#'   a correct assignment. This implementation uses the unrescaled
#'   fraction in `[0, 1]` and weights observations equally within each
#'   pair. Centroids are estimated on the observations being evaluated.
#'
#'   Label trustworthiness is one minus the mean decrease in consistency,
#'   counting gains as zero. Label continuity is one minus the mean
#'   increase, counting losses as zero. Each unordered pair has equal
#'   weight in these means. A pair can gain consistency when the embedding
#'   exaggerates the separation visible in the reference.
#'
#'   Distance consistency uses one centroid per class. Compare it with
#'   [roc_auc()] or [pr_auc()] when classes have multiple clusters or
#'   irregular shapes. Nonfinite computed centroid distances cause an error.
#' @return A named numeric vector with `trustworthiness` and `continuity`.
#'   With `ret_extra = TRUE`, a list containing `scores` and `pairs`.
#'   The data frame `pairs` has columns `label_a`, `label_b`,
#'   `consistency_in`, `consistency_out`, and `change` (output minus
#'   input). Pairs follow the first appearance of the labels.
#' @references Jeon, H., Kuo, Y.-H., Aupetit, M., Ma, K.-L., & Seo, J.
#'   (2024). Classes are not Clusters: Improving Label-based Evaluation
#'   of Dimensionality Reduction. \doi{10.1109/TVCG.2023.3327187}
#' @seealso [roc_auc()], [pr_auc()]
#' @examples
#' x <- as.matrix(iris[, -5])
#' y <- prcomp(x, rank. = 2)$x
#' label_trustworthiness_continuity(x, y, iris$Species, ret_extra = TRUE)
#' @export
label_trustworthiness_continuity <- function(
  Xin,
  Xout,
  labels,
  ret_extra = FALSE
) {
  ret_extra <- validate_scalar_logical(ret_extra, "ret_extra")
  Xin <- prepare_input_matrix(Xin, name = "Xin", require_finite = TRUE)
  Xout <- prepare_input_matrix(Xout, name = "Xout", require_finite = TRUE)
  if (nrow(Xout) != nrow(Xin)) {
    stop(
      "Xin and Xout must have the same number of observations",
      call. = FALSE
    )
  }
  validate_label_vector(labels, nrow(Xin), "Xin")
  class_labels <- unique(labels)
  if (length(class_labels) < 2L) {
    stop("labels must contain at least two distinct labels", call. = FALSE)
  }

  groups <- split(seq_along(labels), match(labels, class_labels))
  pair_indices <- utils::combn(seq_along(class_labels), 2L)
  consistency_in <- label_pair_consistency(Xin, groups, pair_indices)
  consistency_out <- label_pair_consistency(Xout, groups, pair_indices)
  change <- consistency_out - consistency_in
  scores <- c(
    trustworthiness = 1 - mean(pmax(-change, 0)),
    continuity = 1 - mean(pmax(change, 0))
  )
  if (!ret_extra) {
    return(scores)
  }
  pairs <- data.frame(
    label_a = class_labels[pair_indices[1L, ]],
    label_b = class_labels[pair_indices[2L, ]],
    consistency_in = consistency_in,
    consistency_out = consistency_out,
    change = change
  )
  list(scores = scores, pairs = pairs)
}

label_pair_consistency <- function(x, groups, pair_indices) {
  members <- lapply(groups, function(rows) x[rows, , drop = FALSE])
  centroids <- lapply(members, colMeans)
  within <- Map(squared_centroid_distances, members, centroids)
  vapply(
    seq_len(ncol(pair_indices)),
    function(i) {
      a <- pair_indices[1L, i]
      b <- pair_indices[2L, i]
      own <- c(within[[a]], within[[b]])
      other <- c(
        squared_centroid_distances(members[[a]], centroids[[b]]),
        squared_centroid_distances(members[[b]], centroids[[a]])
      )
      mean((own < other) + 0.5 * (own == other))
    },
    numeric(1)
  )
}

squared_centroid_distances <- function(x, centroid) {
  distances <- rowSums(sweep(x, 2L, centroid, "-")^2)
  if (any(!is.finite(distances))) {
    stop("Nonfinite computed centroid distances", call. = FALSE)
  }
  distances
}
