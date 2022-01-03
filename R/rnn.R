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
#' * `idx` a matrix with as many as rows as observations in the input data and
#' `k` columns. The `i`th row of this matrix contains the row indices of the
#' nearest neighbors of observation `i` in non-decreasing distance order.
#' * `dist` (optional) a matrix with the same dimensions as `idx`, containing the
#' equivalent distances. This information is not actually used by the
#' preservation function but is included in the output of the nearest neighbor
#' calculation. So if you provide your own nearest neighbor data, this matrix
#' does not need to be present.
#'
#' @param Xin the input data (usually high-dimensional), a matrix or data frame
#'   with one observation per row. Alternatively, it can be a pre-computed
#'   nearest neighbor graph. In the latter case, `nn_method_in`, `metric_in` and
#'   `nn_args_in` are ignored. If `Xin` is a data-frame, non-numeric columns
#'   are ignored.
#' @param Xout the output data (usually lower dimensional than `Xin`), a matrix
#'   or data frame with one observation per row. Alternatively, it can be a
#'   pre-computed nearest neighbor graph. In the latter case, `nn_method_out`,
#'   `metric_out` and `nn_args_out` are ignored. If `Xout` is a data-frame,
#'   non-numeric columns are ignored.
#' @param k the number of nearest neighbors to find. Can be a numeric vector,
#' in which case the preservation is calculated for each value separately.
#' @param nn_method_in the nearest neighbor method to calculate the neighbors of
#'   `Xin`. Can be one of `"brute"` (brute force calculation) or `"nnd"`, the
#'   nearest neighbor descent method of Dong and co-workers (2011).
#' @param metric_in the distance calculation to apply to `Xin`. One of
#'   `"euclidean"`, `"l2sqr"` (squared Euclidean), `"cosine"`, `"manhattan"`,
#'   `"correlation"` (1 minus the Pearson correlation), or `"hamming"`.
#' @param nn_method_out the nearest neighbor method to calculate the neighbors
#'   of `Xout`. See `nn_method_in` for details.
#' @param metric_out the distance metric to apply to `Xout`. See `metric_in` for
#'   details.
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
#' @return the mean value of the intersection of the neighborhoods per observation
#'   scaled between `0` (no neighbors in common) to `1` (all neighbors in common).
#'   If `k` is a vector, then the return value is a vector of the preservations
#'   for each `k` in the order they were passed.
#'   If `ret_extra = TRUE`, then a list is returned containing:
#'
#'   * `nnp`: the vector of nearest neighbor preservation values.
#'   * `nn_in`: the nearest neighbor graph for `Xin`. See the
#'   'Nearest Neighbor Graph Format' section for details.
#'   * `nn_out`: the nearest neighbor graph for `Xout`. See the
#'   'Nearest Neighbor Graph Format' section for details.
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
#' @md
#' @examples
#' iris_pca <- stats::prcomp(iris[, -5], rank. = 2, scale = FALSE, retx = TRUE)$x
#' nn_preservation(iris, iris_pca)
#'
#' # Calculate for multiple values of k
#' nn_preservation(iris, iris_pca, k = c(15, 30))
#'
#' # Return the nearest neighbor graphs
#' res <- nn_preservation(iris, iris_pca, k = c(15, 30), ret_extra = TRUE)
#'
#' # Re-use the input neighbor graph for these calculations
#' nn_preservation(res$nn_in, iris_pca, k = c(2, 5, 10), ret_extra = TRUE)
#'
#' # For small datasets, brute force search is more efficient than nearest
#' # neighbor descent
#' nn_preservation(res$nn_in, iris_pca, k = c(2, 5, 10), ret_extra = TRUE,
#'                 nn_method_in = "brute")
#' @export
nn_preservation <- function(Xin,
                            Xout,
                            k = 15,
                            nn_method_in = "nnd",
                            metric_in = "l2sqr",
                            nn_method_out = "brute",
                            metric_out = "l2sqr",
                            n_threads = 0,
                            verbose = FALSE,
                            ret_extra = FALSE,
                            nn_args_in = list(),
                            nn_args_out = list()) {
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
        n_threads = n_threads,
        verbose = verbose,
        nn_args = nn_args_in
      )
    )
  stopifnot(is_nn_graph(nn_in))

  tsmessage("Getting neighbor graph for Xout")
  nn_out <-
    do.call(
      get_nn_graph,
      list(
        X = Xout,
        k = max_k,
        nn_method = nn_method_out,
        metric = metric_out,
        n_threads = n_threads,
        verbose = verbose,
        nn_args = nn_args_out
      )
    )
  stopifnot(is_nn_graph(nn_out))

  if (graph_dim(nn_in)[1] != graph_dim(nn_out)[1]) {
    stop("Number of observations in nn_in != nn_out")
  }

  nnps <-
    sapply(k, function(ki) {
      nn_accuracy(nn_out, ref_idx = nn_in, k = ki)
    })

  names(nnps) <- paste0("nnp", k)

  if (ret_extra) {
    res <- list(
      nn_in = nn_in,
      nn_out = nn_out,
      nnp = nnps
    )
  } else {
    res <- nnps
  }
  res
}

calc_nn_graph <-
  function(X,
           k = 15,
           nn_method = "brute",
           metric = "l2sqr",
           n_threads = 0,
           verbose = FALSE,
           ...) {
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

    nnfun <- switch(nn_method,
      brute = rnndescent::brute_force_knn,
      nnd = rnndescent::nnd_knn,
      stop("unknown method '", nn_method, "'")
    )
    do.call(nnfun, varargs)
  }

get_nn_graph <-
  function(X,
           k = 15,
           nn_method,
           metric,
           n_threads,
           verbose,
           nn_args) {
    if (is_nn_graph(X)) {
      if (graph_k(X) < k) {
        stop("Graph not large enough for requested k")
      }
      else {
        return(X)
      }
    } else {
      do.call(calc_nn_graph, lmerge(
        list(
          X = X,
          k = k,
          nn_method = nn_method,
          metric = metric,
          n_threads = n_threads,
          verbose = verbose
        ),
        nn_args
      ))
    }
  }

graph_k <- function(nn_graph) {
  ncol(nn_graph$idx)
}

graph_dim <- function(nn_graph) {
  dim(nn_graph$idx)
}

is_nn_graph <- function(graph) {
  if (!is.list(graph) || is.null(graph$idx) || is.null(graph$dist)) {
    return(FALSE)
  }
  idx <- graph$idx
  if (!is.matrix(idx)) {
    return(FALSE)
  }
  dist <- graph$dist
  if (!is.matrix(dist)) {
    return(FALSE)
  }
  all(dim(idx) == dim(dist))
}

check_graph <- function(idx, dist = NULL, k = NULL) {
  if (is.null(dist) && is.list(idx)) {
    dist <- idx$dist
  }
  if (is.list(idx)) {
    idx <- idx$idx
  }
  stopifnot(methods::is(idx, "matrix"))
  stopifnot(methods::is(dist, "matrix"))
  stopifnot(dim(idx) == dim(dist))
  if (is.null(k)) {
    k <- ncol(idx)
  }
  stopifnot(k > 0)
  list(idx = idx, dist = dist, k = k)
}

nn_intersect <-
  function(idx,
           ref_idx,
           k = NULL,
           include_self = TRUE,
           ref_k = NULL,
           verbose = TRUE) {
    if (is.list(idx)) {
      idx <- idx$idx
    }
    stopifnot(methods::is(idx, "matrix"))

    if (is.list(ref_idx)) {
      ref_idx <- ref_idx$idx
    }
    stopifnot(methods::is(ref_idx, "matrix"))

    if (is.null(k)) {
      k <- findk(idx, ref_idx)
    }
    if (is.null(ref_k)) {
      ref_k <- k
    }

    if (ncol(ref_idx) < k) {
      stop("Not enough columns in ref_idx for k = ", k)
    }

    n <- nrow(idx)
    if (nrow(ref_idx) != n) {
      stop("Not enough rows in ref_idx")
    }

    nbr_start <- 1
    nbr_end <- k

    ref_start <- nbr_start
    ref_end <- ref_k

    if (!include_self) {
      ref_start <- ref_start + 1
      if (nbr_end < ncol(ref_idx)) {
        ref_end <- ref_end + 1
      } else {
        nbr_end <- nbr_end - 1
      }
    }

    nbr_range <- nbr_start:nbr_end
    ref_range <- ref_start:ref_end

    total_intersect <- rep(0, )
    for (i in 1:n) {
      total_intersect[i] <-
        length(intersect(idx[i, nbr_range], ref_idx[i, ref_range]))
    }

    total_intersect
  }



findk <- function(idx, ref_idx) {
  if (is.list(idx)) {
    idx <- idx$idx
  }
  if (is.list(ref_idx)) {
    ref_idx <- ref_idx$idx
  }
  min(ncol(idx), ncol(ref_idx))
}

nn_accuracyv <-
  function(idx,
           ref_idx,
           k = NULL,
           include_self = TRUE,
           ref_k = NULL,
           verbose = TRUE) {
    nn_intersect(
      idx = idx,
      ref_idx = ref_idx,
      k = k,
      include_self = include_self,
      ref_k = ref_k,
      verbose = verbose
    ) /
      ifelse(is.null(k), findk(idx, ref_idx), k)
  }

nn_accuracy <-
  function(idx,
           ref_idx,
           k = NULL,
           include_self = TRUE,
           ref_k = NULL,
           verbose = TRUE) {
    vec <- nn_accuracyv(
      idx = idx,
      ref_idx = ref_idx,
      k = k,
      include_self = include_self,
      ref_k = ref_k,
      verbose = verbose
    )
    sum(vec) / length(vec)
  }
