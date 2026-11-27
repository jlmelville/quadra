#' Steadiness and Cohesiveness of Groups
#'
#' Steadiness checks whether groups sampled in the embedding are coherent
#' in the reference. Cohesiveness checks the reverse direction. Larger
#' scores indicate better preservation. Class labels are unnecessary.
#'
#' @inheritParams class_angular_distortion_index
#' @param k Number of other Euclidean neighbors. The default is the integer
#'   part of the square root of the number of rows. Distance ties follow
#'   row order.
#' @param n_iter Number of independently grown groups per direction.
#'   Use `set.seed()` to reproduce a draw.
#' @param ret_extra If true, also return scores for the reliability map.
#' @details The walk follows links with probability given by normalized shared
#'   neighbor similarity, using rank weights from k + 1 down to 2. A walk
#'   accepts about 0.3 times the number of observations in links, or stops
#'   when its queue is empty. The opposite space partitions the sampled
#'   group with [dbscan::hdbscan()] (minimum points 5). Each noise point
#'   forms a singleton. The dissimilarity between two subgroups is the
#'   reciprocal of the sum of their mean similarity and 0.1.
#'
#'   Each score subtracts the weighted mean positive directional distortion
#'   from one. Weights are the products of subgroup sizes. A direction
#'   with no observed positive distortion receives one. Both `dbscan` and
#'   `Matrix` are required.
#' @section Scaling:
#' The calculation builds dense `N` by `N` distances and shared-neighbor
#' similarities internally, even for global scores. A single double matrix at
#' 5,000 rows uses about 191 MiB. Several coexist, and `ret_extra = TRUE` adds
#' square matrices of partner losses and visit counts. For large datasets,
#' evaluate the same representative subset in both spaces and report its size.
#' Reducing `n_iter` saves repeated clustering work. The dense graphs still
#' cover every supplied observation.
#' @return A named vector of `steadiness` and `cohesiveness`. With
#'   `ret_extra = TRUE`, returns `scores` and a row-aligned data frame
#'   `pointwise` with both map values. Each map accumulates weighted losses
#'   by partner, averages repeated events with that partner, then rescales
#'   across observations. Map values depend on the sampled groups. A point
#'   with no recorded loss receives one. Map means do not recover the
#'   global scores. The other direction's global loss can amplify a map
#'   value's distance from one, so map values can fall below zero.
#' @references Jeon et al. (2022). Measuring and Explaining the Inter-Cluster
#'   Reliability of Multidimensional Projections.
#'   \doi{10.1109/TVCG.2021.3114833}
#'   The [authors' code](https://github.com/hj-n/steadiness-cohesiveness)
#'   accompanies the paper. For shared-neighbor similarities and pointwise
#'   map weights, Quadra follows [ZADU 0.5.4](https://github.com/hj-n/zadu).
#'   It uses R's random draws and the dbscan package's HDBSCAN partition.
#' @examplesIf requireNamespace("dbscan", quietly = TRUE) && requireNamespace("Matrix", quietly = TRUE)
#' x <- as.matrix(iris[, -5])
#' set.seed(42)
#' steadiness_cohesiveness(x, prcomp(x)$x[, 1:2], n_iter = 15)
#' @export
steadiness_cohesiveness <- function(
  Xin,
  Xout,
  k = NULL,
  n_iter = 150L,
  ret_extra = FALSE
) {
  if (
    !requireNamespace("dbscan", quietly = TRUE) ||
      !requireNamespace("Matrix", quietly = TRUE)
  ) {
    stop(
      "steadiness_cohesiveness requires 'dbscan' and 'Matrix'",
      call. = FALSE
    )
  }
  Xin <- prepare_input_matrix(Xin, name = "Xin", require_finite = TRUE)
  Xout <- prepare_input_matrix(Xout, name = "Xout", require_finite = TRUE)
  n <- nrow(Xin)
  if (n != nrow(Xout)) {
    stop(
      "Xin and Xout must have the same number of observations",
      call. = FALSE
    )
  }
  k <- validate_positive_integer(if (is.null(k)) floor(sqrt(n)) else k, "k")
  if (k >= n) {
    stop("k must be smaller than the number of observations", call. = FALSE)
  }
  n_iter <- validate_positive_integer(n_iter, "n_iter")
  ret_extra <- validate_scalar_logical(ret_extra, "ret_extra")
  reference <- snc_graph(Xin, k)
  embedding <- snc_graph(Xout, k)
  change <- 1 / (reference$similarity + 0.1) - 1 / (embedding$similarity + 0.1)
  steady <- snc_direction(embedding, reference, max(change), n_iter, ret_extra)
  cohesive <- snc_direction(
    reference,
    embedding,
    -min(change),
    n_iter,
    ret_extra
  )
  scores <- c(steadiness = steady$score, cohesiveness = cohesive$score)
  if (!ret_extra) {
    return(scores)
  }
  list(
    scores = scores,
    pointwise = data.frame(
      steadiness = 1 - steady$loss * max(1, 2 * (1 - cohesive$score)),
      cohesiveness = 1 - cohesive$loss * max(1, 2 * (1 - steady$score))
    )
  )
}

snc_graph <- function(x, k) {
  distances <- as.matrix(stats::dist(x))
  if (any(!is.finite(distances))) {
    stop("Nonfinite computed distance", call. = FALSE)
  }
  diag(distances) <- Inf
  n <- nrow(x)
  neighbors <- matrix(
    unlist(
      lapply(seq_len(n), function(i) {
        order(distances[i, ], method = "radix")[seq_len(k)]
      }),
      use.names = FALSE
    ),
    nrow = n,
    byrow = TRUE
  )
  weights <- Matrix::sparseMatrix(
    i = rep(seq_len(n), k),
    j = as.vector(neighbors),
    x = rep(seq(k + 1L, 2L), each = n),
    dims = c(n, n)
  )
  similarity <- as.matrix(Matrix::tcrossprod(weights))
  diag(similarity) <- 0
  peak <- max(similarity)
  if (peak > 0) {
    similarity <- similarity / peak
  }
  probability <- matrix(
    similarity[cbind(rep(seq_len(n), k), as.vector(neighbors))],
    nrow = n
  )
  list(
    neighbors = neighbors,
    similarity = similarity,
    probability = probability
  )
}

snc_walk <- function(graph, seeds, budget) {
  queue <- seeds[sample.int(length(seeds), 1L)]
  head <- 1L
  while (head <= length(queue) && length(queue) <= budget) {
    row <- queue[head]
    take <- stats::runif(ncol(graph$neighbors)) < graph$probability[row, ]
    # Condition the seed's first expansion on at least one accepted link.
    if (head == 1L && !any(take)) {
      next
    }
    queue <- c(queue, graph$neighbors[row, take])
    head <- head + 1L
  }
  unique(queue)
}

snc_partition <- function(similarity) {
  n <- nrow(similarity)
  if (n < 5L) {
    return(seq_len(n))
  }
  groups <- dbscan::hdbscan(
    stats::as.dist(1 / (similarity + 0.1)),
    minPts = 5L
  )$cluster
  noise <- which(groups == 0L)
  groups[noise] <- max(groups) + seq_along(noise)
  match(groups, unique(groups))
}

snc_group_distances <- function(similarity, groups, weights) {
  sums <- rowsum(t(rowsum(similarity, groups)), groups)
  1 / (sums / weights + 0.1)
}

snc_direction <- function(from, to, scale, n_iter, record) {
  n <- nrow(from$neighbors)
  loss <- numeric(n)
  seeds <- which(rowSums(from$probability) > 0)
  if (!length(seeds)) {
    stop("No walk can expand. Increase k.", call. = FALSE)
  }
  if (scale <= 0) {
    return(list(score = 1, loss = loss))
  }
  if (record) {
    partner_loss <- matrix(0, n, n)
    visits <- matrix(0L, n, n)
  }
  total <- weight <- 0
  for (iteration in seq_len(n_iter)) {
    ids <- snc_walk(from, seeds, max(1L, floor(0.3 * n)))
    opposite <- to$similarity[ids, ids, drop = FALSE]
    groups <- snc_partition(opposite)
    sizes <- tabulate(groups)
    weights <- outer(sizes, sizes)
    delta <- snc_group_distances(opposite, groups, weights) -
      snc_group_distances(
        from$similarity[ids, ids, drop = FALSE],
        groups,
        weights
      )
    diag(delta) <- 0
    events <- delta > 0
    weighted_loss <- pmax(delta, 0) * weights / scale
    keep <- lower.tri(delta) & events
    total <- total + sum(weighted_loss[keep])
    weight <- weight + sum(weights[keep])
    if (record) {
      partner_loss[ids, ids] <- partner_loss[ids, ids] +
        weighted_loss[groups, groups]
      visits[ids, ids] <- visits[ids, ids] + events[groups, groups]
    }
  }
  if (record) {
    loss <- rowSums(partner_loss / pmax(visits, 1L))
    if (max(loss) > 0) loss <- loss / max(loss)
  }
  list(score = if (weight > 0) 1 - total / weight else 1, loss = loss)
}
