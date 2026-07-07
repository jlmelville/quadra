#' Nearest Neighbor Preservation
#'
#' Calculates the overlap of two nearest neighbors graphs.
#'
#' As a measure of local structure preservation, the k-nearest neighbor graph of
#' the input data, and that of the output data is calculated. For each
#' observation the number of neighbors in common in each graph is calculated and
#' the mean value over the entire dataset is reported.
#'
#' There are two methods for calculating the nearest neighbor graph: a brute
#' force approach involving calculating all pairs distances or an approximate
#' nearest neighbors approach using the nearest neighbor descent method of Dong
#' and co-workers (2011). The brute force approach scales with the square of the
#' number of observations and is linear in the number of features. The nearest
#' neighbor descent method is likely to be more economical for high dimensional
#' input. Because the output data is often low-dimensional (e.g. 2D), brute
#' force neighbor search is often feasible for `Xout` even when it isn't for
#' `Xin`.
#'
#' For finer control of the nearest neighbor calculations, you can pass extra
#' arguments to those functions via the `nn_args_in` and `nn_args_out` lists.
#' See the documentation for the [rnndescent::brute_force_knn()] and
#' [rnndescent::nnd_knn()] functions for more details.
#'
#' Because the nearest neighbor search can be time-consuming, if you set
#' `ret_extra = TRUE`, the return value of this function is a list which
#' includes the nearest neighbor graph for `Xin` and for `Xout`. The graph can
#' be passed to `Xin` or `Xout` and re-used. This is useful if you are comparing
#' a fixed `Xin` with multiple `Xout`, e.g. where different dimensionality
#' reduction methods have been used, or parameters such as output dimensionality
#' or random number seed have been been modified.
#'
#' @section Nearest Neighbor Graph Format:
#'
#' Rather than provide the observations to `Xin` and `Xout`, pre-computed
#' nearest neighbor graphs can be provided. The format of the graph must be a
#' list containing:
#'
#' * `idx` a matrix with as many rows as observations in the input data and
#' `k` columns. The `i`th row of this matrix contains the row indices of the
#' nearest non-self neighbors of observation `i` in non-decreasing distance
#' order. Graphs calculated by this function are self-excluded. Older cached
#' self-inclusive graphs are detected and stripped with a warning.
#' * `dist` (optional) a matrix with the same dimensions as `idx`, containing the
#' equivalent distances. This information is not actually used by the
#' preservation function but is included in the output of the nearest neighbor
#' calculation. So if you provide your own nearest neighbor data, this matrix
#' does not need to be present.
#'
#' @param Xin the input data (usually high-dimensional), a matrix or data frame
#'   with one observation per row, or if `is_transposed = TRUE`, one observation
#'   per column. Alternatively, it can be a pre-computed nearest neighbor
#'   graph. In the latter case, `nn_method_in`, `metric_in` and `nn_args_in` are
#'   ignored. If `Xin` is a data-frame, non-numeric columns are ignored.
#' @param Xout the output data (usually lower dimensional than `Xin`), a matrix
#'   or data frame with one observation per row, or if `is_transposed = TRUE`,
#'   one observation per column. Alternatively, it can be a pre-computed nearest
#'   neighbor graph. In the latter case, `nn_method_out`, `metric_out` and
#'   `nn_args_out` are ignored. If `Xout` is a data-frame, non-numeric columns
#'   are ignored.
#' @param k the number of nearest neighbors to find. Can be a numeric vector,
#' in which case the preservation is calculated for each value separately.
#' @param nn_method_in the nearest neighbor method to calculate the neighbors of
#'   `Xin`. Can be one of `"brute"` (brute force calculation) or `"nnd"`, the
#'   nearest neighbor descent method of Dong and co-workers (2011).
#' @param metric_in the distance calculation to apply to `Xin`. One of
#'   `"euclidean"`, `"sqeuclidean"` (squared Euclidean), `"cosine"`, `"manhattan"`,
#'   `"correlation"` (1 minus the Pearson correlation), or `"hamming"`.
#' @param nn_method_out the nearest neighbor method to calculate the neighbors
#'   of `Xout`. See `nn_method_in` for details.
#' @param metric_out the distance metric to apply to `Xout`. See `metric_in` for
#'   details.
#' @param is_transposed if `TRUE` then `Xin` and `Xout` are assumed to have been
#'   passed in transposed format, i.e. with one observation per column.
#'   Otherwise, `Xin` and `Xout` will be transposed. For large datasets,
#'   transposing can be slow, so if this function will be called multiple times
#'   with the same input data, it is more efficient to transpose the input data
#'   once outside of this function and set `is_transposed = TRUE`.
#' @param n_threads the maximum number of threads to use.
#' @param verbose if `TRUE`, log information about the calculation to the
#'   console.
#' @param ret_extra if `TRUE`, additionally return the nearest neighbor graphs
#'   for `Xin` and `Xout`.
#' @param nn_args_in list of extra arguments to pass to the nearest neighbor
#'   methods, [rnndescent::brute_force_knn()] or [rnndescent::nnd_knn()],
#'   depending on the value of `nn_method_in`.
#' @param nn_args_out list of extra arguments to pass to the nearest neighbor
#'   methods, [rnndescent::brute_force_knn()] or [rnndescent::nnd_knn()],
#'   depending on the value of `nn_method_out`.
#' @return the mean value of the intersection of the neighborhoods per
#'   observation scaled between `0` (no neighbors in common) to `1` (all
#'   neighbors in common). For unrelated output neighborhoods, the expected
#'   preservation is approximately `k / (n - 1)`, where n is the number of
#'   observations. If `k` is a vector, then the return value is a vector of the
#'   preservations for each `k` in the order they were passed. If
#'   `ret_extra = TRUE`, then a list is returned containing:
#'
#'   * `nnp`: the vector of nearest neighbor preservation values.
#'   * `nn_in`: the nearest neighbor graph for `Xin`. See the
#'   'Nearest Neighbor Graph Format' section for details.
#'   * `nn_out`: the nearest neighbor graph for `Xout`. See the
#'   'Nearest Neighbor Graph Format' section for details.
#'   * `nnpv`: a list of vectors where each vector contains the individual
#'   neighbor preservation per observation. Items are named `nnp<k>`, where
#'   `<k>` refers to the values provided in the `k` parameter.
#'
#' @seealso The [rnndescent](https://github.com/jlmelville/rnndescent) package
#'   and the [rnndescent::brute_force_knn()] and [rnndescent::nnd_knn()]
#'   functions.
#' @references
#' Dong, W., Moses, C., & Li, K. (2011, March).
#' Efficient k-nearest neighbor graph construction for generic similarity measures.
#' In *Proceedings of the 20th international conference on World Wide Web*
#' (pp. 577-586).
#' ACM.
#' <https://doi.org/10.1145/1963405.1963487>.
#' @examples
#' iris_pca <- stats::prcomp(iris[, -5], rank. = 2, scale = FALSE, retx = TRUE)$x
#' nn_preservation(iris, iris_pca)
#'
#' # Calculate for multiple values of k
#' nn_preservation(iris, iris_pca, k = c(15, 30))
#'
#' # Return the nearest neighbor graphs and per-observation preservation
#' res <- nn_preservation(iris, iris_pca, k = c(15, 30), ret_extra = TRUE)
#'
#' # Plot the 2D PCA, coloring each point by neighbor preservation
#' plot(iris_pca, col = rainbow(15)[round(res$nnpv$nnp15 * 15)], pch = 19)
#'
#' # Re-use the input neighbor graph for these calculations
#' nn_preservation(res$nn_in, iris_pca, k = c(2, 5, 10), ret_extra = TRUE)
#'
#' # For small datasets, brute force search is more efficient than nearest
#' # neighbor descent
#' nn_preservation(res$nn_in, iris_pca,
#'   k = c(2, 5, 10), ret_extra = TRUE,
#'   nn_method_in = "brute"
#' )
#' @export
nn_preservation <- function(
  Xin,
  Xout,
  k = 15,
  nn_method_in = "nnd",
  metric_in = "sqeuclidean",
  nn_method_out = "brute",
  metric_out = "sqeuclidean",
  is_transposed = FALSE,
  n_threads = 0,
  verbose = FALSE,
  ret_extra = FALSE,
  nn_args_in = list(),
  nn_args_out = list()
) {
  k <- validate_positive_integer_vector(k, "k")
  max_k <- max(k)

  tsmessage("Getting neighbor graph for Xin")
  nn_in <-
    do.call(
      get_nn_graph,
      list(
        X = Xin,
        k = max_k,
        nn_method = nn_method_in,
        metric = metric_in,
        is_transposed = is_transposed,
        n_threads = n_threads,
        verbose = verbose,
        nn_args = nn_args_in,
        name = "Xin"
      )
    )
  check_nn_graph(nn_in, "Xin")

  tsmessage("Getting neighbor graph for Xout")
  nn_out <-
    do.call(
      get_nn_graph,
      list(
        X = Xout,
        k = max_k,
        nn_method = nn_method_out,
        metric = metric_out,
        is_transposed = is_transposed,
        n_threads = n_threads,
        verbose = verbose,
        nn_args = nn_args_out,
        name = "Xout"
      )
    )
  check_nn_graph(nn_out, "Xout")

  if (graph_dim(nn_in)[1] != graph_dim(nn_out)[1]) {
    stop(
      "Xin and Xout neighbor graphs must have the same number of rows",
      call. = FALSE
    )
  }

  overlap_counts <- nn_overlap_counts(
    idx = nn_out,
    ref_idx = nn_in,
    k = k,
    n_threads = n_threads
  )
  nnp_by_row <- sweep(overlap_counts, 2L, k, `/`)
  nnps <- colMeans(nnp_by_row)
  names(nnps) <- paste0("nnp", k)

  if (ret_extra) {
    nnpvs <- list()
    for (i in seq_along(k)) {
      nnpvs[[paste0("nnp", as.character(k[i]))]] <- nnp_by_row[, i]
    }
    res <- list(
      nn_in = nn_in,
      nn_out = nn_out,
      nnp = nnps,
      nnpv = nnpvs
    )
  } else {
    res <- nnps
  }
  res
}

#' Local Radius Correlation
#'
#' Compares the local scale around each observation in the input data and output
#' embedding. The local scale is summarized from nearest-neighbor distances,
#' either as the distance to the `k`th nearest non-self neighbor or as the mean
#' distance to the first `k` nearest non-self neighbors.
#'
#' This is a local scale diagnostic, not a density estimator. A high value means
#' observations with small or large local radii in the input data tend to have
#' small or large local radii in the output embedding. It complements
#' [nn_preservation()], which measures whether neighbor identities are
#' preserved.
#'
#' `Xin` and `Xout` can be raw observation matrices or pre-computed nearest
#' neighbor graphs. Unlike [nn_preservation()], supplied nearest-neighbor graphs
#' must contain a `dist` matrix as well as an `idx` matrix because this metric
#' uses neighbor distances. Graphs calculated by this function are
#' self-excluded. Older cached self-inclusive graphs are detected and stripped
#' with a warning.
#'
#' If either local-radius vector is constant, the correlation is undefined and
#' the corresponding result is `NA_real_`. If `log = TRUE`, all selected local
#' radius values must be positive, so duplicate points that produce zero local
#' radii should be handled before requesting a log-scale comparison.
#'
#' @param Xin the input data (usually high-dimensional), a matrix or data frame
#'   with one observation per row, or if `is_transposed = TRUE`, one observation
#'   per column. Alternatively, it can be a pre-computed nearest neighbor graph
#'   with `idx` and `dist` matrix elements.
#' @param Xout the output data (usually lower dimensional than `Xin`), a matrix
#'   or data frame with one observation per row, or if `is_transposed = TRUE`,
#'   one observation per column. Alternatively, it can be a pre-computed nearest
#'   neighbor graph with `idx` and `dist` matrix elements.
#' @param k the number of nearest neighbors to use. Can be a numeric vector, in
#'   which case the local radius correlation is calculated for each value
#'   separately.
#' @param statistic the local scale statistic. `"radius"` uses the distance to
#'   the `k`th nearest neighbor. `"mean"` uses the mean distance to the first
#'   `k` nearest neighbors.
#' @param method correlation method, either `"spearman"` or `"pearson"`.
#' @param log if `TRUE`, calculate the correlation on log local radii.
#' @inheritParams nn_preservation
#' @return A named numeric vector of local radius correlations, one value for
#'   each `k`. Items are named `lrc<k>`, where `<k>` refers to the values
#'   provided in the `k` parameter. If `ret_extra = TRUE`, then a list is
#'   returned containing:
#'
#'   * `lrc`: the vector of local radius correlations.
#'   * `nn_in`: the nearest neighbor graph for `Xin`.
#'   * `nn_out`: the nearest neighbor graph for `Xout`.
#'   * `scale_in`: a matrix of input local scale values, one column per `k`.
#'   * `scale_out`: a matrix of output local scale values, one column per `k`.
#'
#' @seealso [nn_preservation()] for neighbor-identity preservation and
#'   [trustworthiness()] for exact rank-penalty neighborhood preservation on
#'   distance matrices.
#' @examples
#' iris_pca <- stats::prcomp(iris[, -5], rank. = 2, scale = FALSE, retx = TRUE)$x
#' local_radius_correlation(
#'   iris[, -5],
#'   iris_pca,
#'   k = 15,
#'   nn_method_in = "brute"
#' )
#'
#' cached <- local_radius_correlation(
#'   iris[, -5],
#'   iris_pca,
#'   k = c(15, 30),
#'   nn_method_in = "brute",
#'   ret_extra = TRUE
#' )
#' local_radius_correlation(cached$nn_in, cached$nn_out, k = c(15, 30))
#' @export
local_radius_correlation <- function(
  Xin,
  Xout,
  k = 15,
  statistic = c("radius", "mean"),
  method = c("spearman", "pearson"),
  log = FALSE,
  nn_method_in = "nnd",
  metric_in = "sqeuclidean",
  nn_method_out = "brute",
  metric_out = "sqeuclidean",
  is_transposed = FALSE,
  n_threads = 0,
  verbose = FALSE,
  ret_extra = FALSE,
  nn_args_in = list(),
  nn_args_out = list()
) {
  k <- validate_positive_integer_vector(k, "k")
  statistic <- match.arg(statistic)
  method <- match.arg(method)
  if (!is.logical(log) || length(log) != 1L || is.na(log)) {
    stop("log must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(ret_extra) || length(ret_extra) != 1L || is.na(ret_extra)) {
    stop("ret_extra must be TRUE or FALSE", call. = FALSE)
  }
  max_k <- max(k)

  tsmessage("Getting neighbor graph for Xin")
  nn_in <-
    do.call(
      get_nn_graph,
      list(
        X = Xin,
        k = max_k,
        nn_method = nn_method_in,
        metric = metric_in,
        is_transposed = is_transposed,
        n_threads = n_threads,
        verbose = verbose,
        nn_args = nn_args_in,
        name = "Xin"
      )
    )
  check_nn_graph_dist(nn_in, "Xin")

  tsmessage("Getting neighbor graph for Xout")
  nn_out <-
    do.call(
      get_nn_graph,
      list(
        X = Xout,
        k = max_k,
        nn_method = nn_method_out,
        metric = metric_out,
        is_transposed = is_transposed,
        n_threads = n_threads,
        verbose = verbose,
        nn_args = nn_args_out,
        name = "Xout"
      )
    )
  check_nn_graph_dist(nn_out, "Xout")

  if (graph_dim(nn_in)[1] != graph_dim(nn_out)[1]) {
    stop(
      "Xin and Xout neighbor graphs must have the same number of rows",
      call. = FALSE
    )
  }

  scale_in <- local_radius_values(nn_in, k, statistic)
  scale_out <- local_radius_values(nn_out, k, statistic)
  colnames(scale_in) <- paste0("lrc", k)
  colnames(scale_out) <- paste0("lrc", k)

  if (log) {
    if (any(scale_in <= 0) || any(scale_out <= 0)) {
      stop(
        "local radius values must be positive when log = TRUE",
        call. = FALSE
      )
    }
    scale_in <- base::log(scale_in)
    scale_out <- base::log(scale_out)
  }

  lrc <- vapply(
    seq_along(k),
    function(i) local_scale_correlation(scale_in[, i], scale_out[, i], method),
    numeric(1)
  )
  names(lrc) <- paste0("lrc", k)

  if (ret_extra) {
    list(
      nn_in = nn_in,
      nn_out = nn_out,
      lrc = lrc,
      scale_in = scale_in,
      scale_out = scale_out
    )
  } else {
    lrc
  }
}

#' Mutual Neighbor Correlation
#'
#' Compares the mutual-neighbor count pattern in the input data and output
#' embedding. For each observation, its mutual-neighbor count is the number of
#' its first `k` nearest non-self neighbors that also include the observation
#' among their first `k` nearest non-self neighbors.
#'
#' This is a local graph diagnostic, not a replacement for nearest-neighbor
#' preservation. A high value means observations with many or few mutual-neighbor
#' relationships in the input graph tend to have many or few mutual-neighbor
#' relationships in the output graph. It complements [nn_preservation()], which
#' measures whether neighbor identities are preserved directly.
#'
#' `Xin` and `Xout` can be raw observation matrices or pre-computed nearest
#' neighbor graphs. Unlike [local_radius_correlation()], supplied
#' nearest-neighbor graphs only need an `idx` matrix because this metric uses
#' neighbor identities, not distances. Graphs calculated by this function are
#' self-excluded. Older cached self-inclusive graphs are detected and stripped
#' with a warning.
#'
#' If either mutual-neighbor count vector is constant, the correlation is
#' undefined and the corresponding result is `NA_real_`.
#'
#' @param k the number of nearest neighbors to use. Can be a numeric vector, in
#'   which case the mutual neighbor correlation is calculated for each value
#'   separately.
#' @param method correlation method, either `"pearson"` or `"spearman"`.
#' @inheritParams nn_preservation
#' @return A named numeric vector of mutual neighbor correlations, one value for
#'   each `k`. Items are named `mnc<k>`, where `<k>` refers to the values
#'   provided in the `k` parameter. If `ret_extra = TRUE`, then a list is
#'   returned containing:
#'
#'   * `mnc`: the vector of mutual neighbor correlations.
#'   * `nn_in`: the nearest neighbor graph for `Xin`.
#'   * `nn_out`: the nearest neighbor graph for `Xout`.
#'   * `mutual_neighbor_in`: a matrix of input mutual-neighbor counts, one
#'   column per `k`.
#'   * `mutual_neighbor_out`: a matrix of output mutual-neighbor counts, one
#'   column per `k`.
#'
#' @seealso [nn_preservation()] for neighbor-identity preservation and
#'   [local_radius_correlation()] for local scale preservation.
#' @examples
#' iris_pca <- stats::prcomp(iris[, -5], rank. = 2, scale = FALSE, retx = TRUE)$x
#' mutual_neighbor_correlation(
#'   iris[, -5],
#'   iris_pca,
#'   k = 15,
#'   nn_method_in = "brute"
#' )
#'
#' cached <- mutual_neighbor_correlation(
#'   iris[, -5],
#'   iris_pca,
#'   k = c(15, 30),
#'   nn_method_in = "brute",
#'   ret_extra = TRUE
#' )
#' mutual_neighbor_correlation(cached$nn_in, cached$nn_out, k = c(15, 30))
#' @export
mutual_neighbor_correlation <- function(
  Xin,
  Xout,
  k = 15,
  method = c("pearson", "spearman"),
  nn_method_in = "nnd",
  metric_in = "sqeuclidean",
  nn_method_out = "brute",
  metric_out = "sqeuclidean",
  is_transposed = FALSE,
  n_threads = 0,
  verbose = FALSE,
  ret_extra = FALSE,
  nn_args_in = list(),
  nn_args_out = list()
) {
  k <- validate_positive_integer_vector(k, "k")
  method <- match.arg(method)
  if (!is.logical(ret_extra) || length(ret_extra) != 1L || is.na(ret_extra)) {
    stop("ret_extra must be TRUE or FALSE", call. = FALSE)
  }
  max_k <- max(k)

  tsmessage("Getting neighbor graph for Xin")
  nn_in <-
    do.call(
      get_nn_graph,
      list(
        X = Xin,
        k = max_k,
        nn_method = nn_method_in,
        metric = metric_in,
        is_transposed = is_transposed,
        n_threads = n_threads,
        verbose = verbose,
        nn_args = nn_args_in,
        name = "Xin"
      )
    )
  check_nn_graph(nn_in, "Xin")

  tsmessage("Getting neighbor graph for Xout")
  nn_out <-
    do.call(
      get_nn_graph,
      list(
        X = Xout,
        k = max_k,
        nn_method = nn_method_out,
        metric = metric_out,
        is_transposed = is_transposed,
        n_threads = n_threads,
        verbose = verbose,
        nn_args = nn_args_out,
        name = "Xout"
      )
    )
  check_nn_graph(nn_out, "Xout")

  if (graph_dim(nn_in)[1] != graph_dim(nn_out)[1]) {
    stop(
      "Xin and Xout neighbor graphs must have the same number of rows",
      call. = FALSE
    )
  }

  mutual_neighbor_in <- mutual_neighbor_values(nn_in, k)
  mutual_neighbor_out <- mutual_neighbor_values(nn_out, k)
  colnames(mutual_neighbor_in) <- paste0("mnc", k)
  colnames(mutual_neighbor_out) <- paste0("mnc", k)

  mnc <- vapply(
    seq_along(k),
    function(i) {
      local_scale_correlation(
        mutual_neighbor_in[, i],
        mutual_neighbor_out[, i],
        method
      )
    },
    numeric(1)
  )
  names(mnc) <- paste0("mnc", k)

  if (ret_extra) {
    list(
      nn_in = nn_in,
      nn_out = nn_out,
      mnc = mnc,
      mutual_neighbor_in = mutual_neighbor_in,
      mutual_neighbor_out = mutual_neighbor_out
    )
  } else {
    mnc
  }
}

local_radius_values <- function(nn_graph, k, statistic) {
  dist <- nn_graph$dist
  vapply(
    k,
    function(ki) {
      if (statistic == "radius") {
        dist[, ki]
      } else {
        rowMeans(dist[, seq_len(ki), drop = FALSE])
      }
    },
    numeric(nrow(dist))
  )
}

mutual_neighbor_values <- function(nn_graph, k) {
  idx <- nn_graph$idx
  vapply(
    k,
    function(ki) mutual_neighbor_counts(idx, ki),
    integer(nrow(idx))
  )
}

mutual_neighbor_counts <- function(idx, k) {
  if (ncol(idx) < k) {
    stop("Not enough columns in idx for k = ", k, call. = FALSE)
  }
  idx_prefix <- idx[, seq_len(k), drop = FALSE]
  counts <- integer(nrow(idx_prefix))
  for (i in seq_len(nrow(idx_prefix))) {
    counts[i] <- sum(vapply(
      idx_prefix[i, ],
      function(j) any(idx_prefix[j, ] == i),
      logical(1)
    ))
  }
  counts
}

local_scale_correlation <- function(x, y, method) {
  if (length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(NA_real_)
  }
  unname(stats::cor(x = x, y = y, method = method))
}

calc_nn_graph <-
  function(
    X,
    k = 15,
    nn_method = "brute",
    metric = "sqeuclidean",
    n_threads = 0,
    verbose = FALSE,
    ...
  ) {
    varargs <- lmerge(
      list(
        data = X,
        k = k,
        verbose = verbose,
        n_threads = n_threads,
        metric = metric
      ),
      list(...)
    )

    nnfun <- switch(
      nn_method,
      brute = rnndescent::brute_force_knn,
      nnd = rnndescent::nnd_knn,
      stop("unknown method '", nn_method, "'")
    )
    do.call(nnfun, varargs)
  }

get_nn_graph <-
  function(
    X,
    k = 15,
    nn_method,
    metric,
    is_transposed,
    n_threads,
    verbose,
    nn_args,
    name = "Nearest-neighbor graph"
  ) {
    if (is_nn_graph(X)) {
      prepare_supplied_nn_graph(X, k = k, name = name)
    } else if (is.list(X) && !is.data.frame(X)) {
      check_nn_graph(X, "Nearest-neighbor graph")
    } else {
      X <- x2m(X)
      n_obs <- if (is_transposed) ncol(X) else nrow(X)
      check_k_for_n_obs(k, n_obs)
      if ("k" %in% names(nn_args)) {
        warning(
          "Ignoring 'k' in ",
          name,
          " nn_args; use the top-level k argument",
          call. = FALSE
        )
        nn_args$k <- NULL
      }
      nn_graph <- do.call(
        calc_nn_graph,
        lmerge(
          list(
            X = X,
            k = k + 1L,
            nn_method = nn_method,
            metric = metric,
            n_threads = n_threads,
            verbose = verbose,
            obs = ifelse(is_transposed, "C", "R")
          ),
          nn_args
        )
      )
      strip_self_neighbors(nn_graph, k = k, name = name)
    }
  }

graph_k <- function(nn_graph) {
  ncol(nn_graph$idx)
}

graph_dim <- function(nn_graph) {
  dim(nn_graph$idx)
}

is_nn_graph <- function(graph) {
  if (!is.list(graph) || is.null(graph$idx)) {
    return(FALSE)
  }
  idx <- graph$idx
  if (!is.matrix(idx)) {
    return(FALSE)
  }
  dist <- graph$dist
  if (is.null(dist)) {
    return(TRUE)
  }
  if (!is.matrix(dist)) {
    return(FALSE)
  }
  all(dim(idx) == dim(dist))
}

check_nn_graph <- function(graph, name = "graph") {
  if (!is_nn_graph(graph)) {
    stop(
      name,
      " must be a nearest-neighbor graph: a list with matrix element 'idx' ",
      "and optional matrix element 'dist' with matching dimensions",
      call. = FALSE
    )
  }
  idx <- graph$idx
  if (!is.numeric(idx)) {
    stop(name, " idx must be a numeric matrix", call. = FALSE)
  }
  if (anyNA(idx) || any(!is.finite(idx)) || any(idx != floor(idx))) {
    stop(name, " idx must contain finite integer indices", call. = FALSE)
  }
  if (any(idx < 1L) || any(idx > nrow(idx))) {
    stop(
      name,
      " idx values must be between 1 and the number of graph rows",
      call. = FALSE
    )
  }
  invisible(graph)
}

check_nn_graph_dist <- function(graph, name = "graph") {
  check_nn_graph(graph, name)
  dist <- graph$dist
  if (is.null(dist)) {
    stop(
      name,
      " must contain a 'dist' matrix for local radius correlation",
      call. = FALSE
    )
  }
  if (!is.numeric(dist)) {
    stop(name, " dist must be a numeric matrix", call. = FALSE)
  }
  if (anyNA(dist) || any(!is.finite(dist))) {
    stop(name, " dist must contain finite distances", call. = FALSE)
  }
  if (any(dist < 0)) {
    stop(name, " dist must contain non-negative distances", call. = FALSE)
  }
  invisible(graph)
}

check_k_for_n_obs <- function(k, n_obs) {
  if (k > n_obs - 1L) {
    stop(
      "k cannot be larger than the number of non-self observations",
      call. = FALSE
    )
  }
  invisible(k)
}

has_self_neighbors <- function(graph) {
  idx <- graph$idx
  any(idx == row(idx))
}

prepare_supplied_nn_graph <- function(
  graph,
  k,
  name = "Nearest-neighbor graph"
) {
  check_nn_graph(graph, name)
  check_k_for_n_obs(k, nrow(graph$idx))
  if (has_self_neighbors(graph)) {
    warning(
      name,
      " contains self-neighbors; stripping row self-indices. ",
      "Provide self-excluded graphs to avoid this warning.",
      call. = FALSE
    )
    return(strip_self_neighbors(graph, k = k, name = name))
  }
  if (graph_k(graph) < k) {
    stop(
      "Nearest-neighbor graph does not contain enough columns for requested k",
      call. = FALSE
    )
  }
  graph
}

strip_self_neighbors <- function(graph, k, name = "Nearest-neighbor graph") {
  check_nn_graph(graph, name)
  check_k_for_n_obs(k, nrow(graph$idx))

  idx <- graph$idx
  dist <- graph$dist
  n <- nrow(idx)
  idx_out <- matrix(idx[NA_integer_], nrow = n, ncol = k)
  dist_out <- if (!is.null(dist)) {
    matrix(dist[NA_integer_], nrow = n, ncol = k)
  }

  for (i in seq_len(n)) {
    non_self <- which(idx[i, ] != i)
    if (length(non_self) < k) {
      stop(
        name,
        " does not contain at least ",
        k,
        " non-self neighbors in every row",
        call. = FALSE
      )
    }
    keep <- non_self[seq_len(k)]
    idx_out[i, ] <- idx[i, keep]
    if (!is.null(dist)) {
      dist_out[i, ] <- dist[i, keep]
    }
  }

  graph$idx <- idx_out
  if (!is.null(dist)) {
    graph$dist <- dist_out
  }
  graph
}

nn_overlap_counts <- function(idx, ref_idx, k, n_threads = 0) {
  if (is.list(idx)) {
    idx <- idx$idx
  }
  if (!methods::is(idx, "matrix")) {
    stop("idx must be a matrix or nearest-neighbor graph", call. = FALSE)
  }

  if (is.list(ref_idx)) {
    ref_idx <- ref_idx$idx
  }
  if (!methods::is(ref_idx, "matrix")) {
    stop("ref_idx must be a matrix or nearest-neighbor graph", call. = FALSE)
  }

  k <- validate_positive_integer_vector(k, "k")
  max_k <- max(k)

  if (ncol(idx) < max_k) {
    stop("Not enough columns in idx for max(k) = ", max_k, call. = FALSE)
  }
  if (ncol(ref_idx) < max_k) {
    stop("Not enough columns in ref_idx for max(k) = ", max_k, call. = FALSE)
  }

  if (nrow(ref_idx) != nrow(idx)) {
    stop("idx and ref_idx must have the same number of rows", call. = FALSE)
  }

  neighbor_overlap_counts(idx, ref_idx, k, n_threads)
}
