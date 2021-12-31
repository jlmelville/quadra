get_dfun <- function(metric) {
  switch(metric,
    l2sqr = rnndescent::l2sqr_distance,
    euclidean = rnndescent::euclidean_distance,
    cosine = rnndescent::cosine_distance,
    manhattan = rnndescent::manhattan_distance,
    hamming = rnndescent::hamming_distance,
    correlation = rnndescent::correlation_distance,
    stop("Unknown metric '", metric, "'")
  )
}

rep_random_triplet_eval <- function(n = 5, seed = NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  varargs <- list(...)
  summary(replicate(n, {
    do.call(random_triplet_eval, varargs)
  }))
}

random_triplet_eval <-
  function(Xin,
           Xout,
           n_triplets = 5,
           in_metric = "euclidean",
           out_metric = "euclidean") {
    Xin <- x2m(Xin)
    Xout <- x2m(Xout)
    n <- nrow(Xin)
    if (n != nrow(Xout)) {
      stop("Xin and Xout must have the same number of rows")
    }

    # reading in Python triplets, assume 0-indexed
    # format is
    # 0, t1p1
    # 0, t1p2
    # 0, t2p1
    # etc
    if (methods::is(n_triplets, "matrix") ||
      methods::is(n_triplets, "data.frame")) {
      n_trips <- nrow(n_triplets)
      if (n_trips %% n != 0) {
        stop("Bad number of triplets")
      }
      n_trips <- n_trips / n
      n_triplets <- n_triplets + 1
      triplets <- matrix(t(as.matrix(n_triplets)), nrow = n_trips * 2)
      n_triplets <- n_trips
      avoid_self <- FALSE
    } else {
      # internal format is each column is t1p1, t1p2, t2p1... etc.
      # avoid redundant sampling as much as possible and never pick self as p1 or p2
      avoid_self <- TRUE
      if (n < (n_triplets * 2) + 1) {
        stop("n too small")
      }
      triplets <- replicate(n, {
        sample.int(n - 1, size = n_triplets * 2, replace = FALSE)
      })
    }

    in_dfun <- get_dfun(in_metric)
    out_dfun <- get_dfun(out_metric)
    acc <- 0
    for (i in 1:ncol(triplets)) {
      for (j in 1:n_triplets) {
        c1 <- j * 2 - 1
        p1 <- triplets[c1, i]
        p1 <- ifelse(avoid_self & p1 == i, p1 + 1, p1)

        p2 <- triplets[c1 + 1, i]
        p2 <- ifelse(avoid_self & p2 == i, p2 + 1, p2)

        if (i > nrow(Xin) || p1 > nrow(Xin) || p2 > nrow(Xin)) {
          browser()
        }
        dip1_in <- in_dfun(Xin[i, ], Xin[p1, ])
        dip2_in <- in_dfun(Xin[i, ], Xin[p2, ])

        dip1_out <- out_dfun(Xout[i, ], Xout[p1, ])
        dip2_out <- out_dfun(Xout[i, ], Xout[p2, ])

        if ((dip1_in < dip2_in) == (dip1_out < dip2_out)) {
          acc <- acc + 1
        }
      }
    }
    acc / (n * n_triplets)
  }

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
  if (is.list(X)) {
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
    if (is.list(ref_idx)) {
      ref_idx <- ref_idx$idx
    }

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
