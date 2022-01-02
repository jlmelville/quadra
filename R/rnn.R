calc_nn_graph <-
  function(X,
           k = 15,
           nn_method = "brute",
           metric = "euclidean",
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

get_nn_graph <- function(X, k = 15, nn_args) {
  if (is_nn_graph(X)) {
    nng(X, k = k)
  } else {
    do.call(calc_nn_graph, lmerge(list(X = X, k = k), nn_args))
  }
}

nn_preservation <- function(Xin,
                            Xout,
                            k = 15,
                            in_nn_args = list(),
                            out_nn_args = list(n_threads = 6, verbose = FALSE),
                            ret_extra = FALSE) {
  nn_in <-
    do.call(get_nn_graph, list(X = Xin, k = k, nn_args = in_nn_args))
  nn_out <-
    do.call(get_nn_graph, list(X = Xout, k = k, nn_args = out_nn_args))
  stopifnot(is_nn_graph(nn_in))
  stopifnot(is_nn_graph(nn_out))
  acc <- nn_accuracy(nn_out, ref_idx = nn_in, k = k)
  if (ret_extra) {
    res <- list(
      nn_in = nn_in,
      nn_out = nn_out,
      acc = acc
    )
  } else {
    res <- acc
  }
  res
}

nng <- function(graph, k) {
  list(idx = graph$idx[, 1:k], dist = graph$dist[, 1:k])
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
