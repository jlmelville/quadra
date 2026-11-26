#' Gap Index for Two-Dimensional Embeddings
#'
#' Measure how the relative areas of empty regions change in an embedding.
#' The Delaunay triangles of the embedding define the regions. Each triangle
#' is compared with the same three observations in the reference coordinates.
#'
#' @param Xin Reference coordinates, as a dense numeric matrix or data frame
#'   with observations in rows.
#' @param Xout Two-dimensional embedding coordinates for the same observations
#'   in the same row order. Use at least three non-collinear points with
#'   distinct coordinates. Nonnumeric data-frame columns are ignored in both
#'   inputs.
#' @param ret_extra Whether to return the triangles and their deformations.
#' @inheritParams random_pair_distance_correlation
#' @inheritSection random_pair_distance_correlation Numerical failures
#' @details Each area's relative share is computed by dividing by the sum of
#'   the areas of the same triangles in its own space. Both spaces use
#'   Euclidean geometry. Deformation is the output share minus the input
#'   share, divided by the larger share. Positive values indicate expansion,
#'   and negative values indicate compression. A triangle with zero area in
#'   both spaces has zero deformation. Positive total area is needed in
#'   each space.
#'
#'   The scalar is the mean absolute deformation weighted by each triangle's
#'   larger share across the two spaces. It ranges from 0 (no change in
#'   relative area) to 1. Relative areas are invariant to uniform
#'   rescaling of either space. A shear can preserve them while distorting shape. The triangles
#'   fill the convex hull, so regions at the boundary participate in the
#'   result. The triangulation can change with the embedding or the rows
#'   included. Requires the `geometry` package.
#' @return The scalar Gap Index. With `ret_extra = TRUE`, a list with
#'   `gap_index`, a matrix `triangles` of one-based row indices, and the
#'   corresponding vectors `relative_area_in`, `relative_area_out`, and
#'   signed `deformation`.
#' @references Ros, J., Arleo, A., & Paulovich, F. (2026).
#'   Measuring Distortion in the Empty Regions of Dimensionality Reduction
#'   Scatterplots with the Gap Index. \url{https://arxiv.org/abs/2607.28324}
#' @seealso [local_procrustes()], [pointwise_rnx()]
#' @examplesIf requireNamespace("geometry", quietly = TRUE)
#' x <- as.matrix(iris[, -5])
#' y <- prcomp(x, rank. = 2)$x
#' # Use one observation per distinct location in the embedding.
#' keep <- !duplicated(y)
#' gap_index(x[keep, ], y[keep, ])
#' @export
gap_index <- function(Xin, Xout, n_threads = 0, ret_extra = FALSE) {
  if (!requireNamespace("geometry", quietly = TRUE)) {
    stop("gap_index() requires the geometry package", call. = FALSE)
  }
  n_threads <- validate_n_threads(n_threads)
  ret_extra <- validate_scalar_logical(ret_extra, "ret_extra")
  Xin <- prepare_input_matrix(Xin, name = "Xin", require_finite = TRUE)
  Xout <- prepare_input_matrix(Xout, name = "Xout", require_finite = TRUE)
  if (nrow(Xout) != nrow(Xin)) {
    stop(
      "Xin and Xout must have the same number of observations",
      call. = FALSE
    )
  }
  if (ncol(Xout) != 2L || nrow(Xout) < 3L) {
    stop("Xout must have two columns and at least three rows", call. = FALSE)
  }
  if (anyDuplicated(Xout)) {
    stop("Xout must have distinct embedded coordinates", call. = FALSE)
  }

  triangles <- tryCatch(geometry::delaunayn(Xout), error = function(e) NULL)
  if (is.null(triangles) || nrow(triangles) == 0L) {
    stop(
      "Could not triangulate Xout. Supply non-collinear embedded points.",
      call. = FALSE
    )
  }
  pairs <- rbind(
    triangles[, c(1, 2), drop = FALSE],
    triangles[, c(1, 3), drop = FALSE],
    triangles[, c(2, 3), drop = FALSE]
  )
  distances <- random_pair_distances(
    Xin,
    Xout,
    pairs = pairs,
    metric_in = "euclidean",
    metric_out = "euclidean",
    n_threads = n_threads
  )
  area_in <- relative_triangle_areas(distances$din, nrow(triangles), "Xin")
  area_out <- relative_triangle_areas(distances$dout, nrow(triangles), "Xout")
  delta <- area_out - area_in
  weight <- pmax(area_in, area_out)
  score <- sum(abs(delta)) / sum(weight)

  if (!ret_extra) {
    return(score)
  }
  deformation <- numeric(length(weight))
  positive <- weight > 0
  deformation[positive] <- delta[positive] / weight[positive]
  list(
    gap_index = score,
    triangles = triangles,
    relative_area_in = area_in,
    relative_area_out = area_out,
    deformation = deformation
  )
}

relative_triangle_areas <- function(edges, n, name) {
  scale <- max(edges)
  if (scale > 0) {
    edges <- edges / scale
  }
  first <- edges[seq_len(n)]
  second <- edges[n + seq_len(n)]
  third <- edges[2L * n + seq_len(n)]
  longest <- pmax(first, second, third)
  shortest <- pmin(first, second, third)
  middle <- pmax(pmin(first, second), pmin(pmax(first, second), third))

  # Rearranged Heron formula. The common factor 1/4 cancels when normalized.
  gap <- shortest - (longest - middle)
  gap[gap <= 8 * .Machine$double.eps * longest] <- 0
  areas <- sqrt(
    (longest + (middle + shortest)) *
      gap *
      (shortest + (longest - middle)) *
      (longest + (middle - shortest))
  )
  total <- sum(areas)
  if (total == 0) {
    stop(name, " must have positive total triangle area", call. = FALSE)
  }
  areas / total
}
