#' Random Triplet Accuracy
#'
#' Evaluates the preservation of global structure of dimensionality reduction
#' results using the random triplet accuracy method of Wang and co-workers
#' (2020).
#'
#' The random triplet accuracy is calculated by randomly selecting three points
#' in the input data and calculating the distances for two sides of the
#' resulting triangle. This is repeated for the output data, and the relative
#' ordering of the distances are compared. The returned accuracy is the
#' proportion of triangles where the relative distances agree between the input
#' and output data.
#'
#' @param Xin the input data (usually high-dimensional), a matrix or data frame
#'   with one observation per row, or if `is_transposed = TRUE`, one observation
#'   per column.
#' @param Xout the output data (usually lower dimensional than `Xin`), a matrix
#'   or data frame with one observation per row, or if `is_transposed = TRUE`,
#'   one observation per column.
#' @param n_triplets the number of triplets per observation to generate.
#' @param metric_in the distance calculation to apply to `Xin`. One of
#'   `"euclidean"`, `"l2sqr"` (squared Euclidean), `"cosine"`, `"manhattan"`,
#'   `"correlation"` (1 minus the Pearson correlation), or `"hamming"`.
#' @param metric_out the distance metric to apply to `Xout`. See `metric_in` for
#'   details.
#' @param is_transposed if `TRUE` then `Xin` and `Xout` are assumed to have been
#'   passed in transposed format, i.e. with one observation per column.
#'   Otherwise, `Xin` and `Xout` will be transposed. For large datasets,
#'   transposing can be slow, so if this function will be called multiple times
#'   with the same input data, it is more efficient to transpose the input data
#'   once outside of this function and set `is_transposed = TRUE`.
#' @param n_threads the maximum number of threads to use.
#' @param grain_size the minimum number of observations to be processed per
#'   thread.
#' @return The triplet accuracy, ranging from 0 (no relative distances agree) to
#'   1 (all of them agree).
#' @references Wang, Y., Huang, H., Rudin, C., & Shaposhnik, Y. (2021).
#' Understanding how dimension reduction tools work: an empirical approach to
#' deciphering t-SNE, UMAP, TriMAP, and PaCMAP for data visualization.
#' *J Mach. Learn. Res*, *22*, 1-73. <https://jmlr.org/papers/v22/20-1061.html>.
#' @md
#' @examples
#' iris_pca2 <- stats::prcomp(iris[, -5], rank. = 2, scale = FALSE, retx = TRUE)$x
#' random_triplet_accuracy(iris, iris_pca2)
#'
#' # If you plan on comparing the results of multiple output methods, then
#' # pre-transposing the input data can save time
#' tiris <- t(iris[, -5])
#' iris_pca1 <- stats::prcomp(iris[, -5], rank. = 1, scale = FALSE, retx = TRUE)$x
#' iris_pca3 <- stats::prcomp(iris[, -5], rank. = 3, scale = FALSE, retx = TRUE)$x
#' random_triplet_accuracy(tiris, t(iris_pca1), is_transposed = TRUE)
#' random_triplet_accuracy(tiris, t(iris_pca2), is_transposed = TRUE)
#' random_triplet_accuracy(tiris, t(iris_pca3), is_transposed = TRUE)
#' @export
random_triplet_accuracy <-
  function(Xin,
           Xout,
           n_triplets = 5,
           metric_in = "l2sqr",
           metric_out = "l2sqr",
           is_transposed = FALSE,
           n_threads = 0,
           grain_size = 1) {
    metric_in <- validate_distance(metric_in)
    metric_out <- validate_distance(metric_out)

    Xin <- x2m(Xin)
    Xout <- x2m(Xout)
    if (!is_transposed) {
      Xin <- t(Xin)
      Xout <- t(Xout)
    }
    n_obs <- ncol(Xin)
    if (n_obs != ncol(Xout)) {
      stop("Xin and Xout must have the same number of observations")
    }

    triplets <-
      get_triplet_matrix(n_obs, n_triplets, zero_index = TRUE)

    triplet_sample(
      triplets,
      Xin,
      Xout,
      metric_in = metric_in,
      metric_out = metric_out,
      n_threads = n_threads
    )
  }


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

# n_triplets might be a pre-generated triplets matrix
get_triplet_matrix <- function(n_obs, n_triplets, zero_index) {
  if (is.matrix(n_triplets)) {
    triplets <- n_triplets
    if (ncol(triplets) != n_obs) {
      stop("Triplets matrix must have ", n_obs, " columns")
    }
    if (nrow(triplets) %% 2 != 0) {
      stop("Triplets matrix must have even number of rows")
    }
    if (min(triplets) < 0) {
      stop("Triplet matrix must have non-negative values")
    }
    max_trip_idx <- max(triplets)
    if (max_trip_idx > n_obs - 1) {
      stop("Triplet matrix elements must be in (0, ", n_obs - 1, ")")
    }
  } else {
    if (!is.numeric(n_triplets) || n_triplets < 1) {
      stop("n_triplets should be int > 0")
    }
    triplets <-
      create_triplet_matrix(n_obs, n_triplets, zero_index = TRUE)
  }
  triplets
}

# triplet matrix is n_obs x (2 x triplets)
# each column contains n_triplet pairs containing the indices of the other
# points of the triangle
# the other indices never contain i
create_triplet_matrix <-
  function(n, n_triplets, zero_index = FALSE) {
    res <- replicate(n * n_triplets, {
      sample.int(n - 1, size = 2, replace = FALSE)
    })
    dim(res) <- c(n_triplets * 2, n)
    res <- avoid_self_matrix(res)
    if (zero_index) {
      res <- res - 1
    }
    res
  }

# Ensure we never create a triangle with two identical indices
# Assumes that the triplet matrix contains indices in the range (1, nobs-1).
# for each column add 1 to each element which is >= the column index
avoid_self_matrix <- function(triplets) {
  # add the column index as a new row (1:ncol)
  tripi <- rbind(triplets, 1:ncol(triplets))

  # for each column, if the element >= the column index, add 1 to it
  tripi <-
    tripi + apply(tripi, 2, function(x) {
      x >= x[length(x)]
    })

  # undecorate the matrix
  tripi[-nrow(tripi), ]
}

random_triplet_eval_slow <-
  function(Xin,
           Xout,
           n_triplets,
           in_metric = "euclidean",
           out_metric = "euclidean",
           is_transposed = FALSE) {
    Xin <- x2m(Xin)
    Xout <- x2m(Xout)

    if (is_transposed) {
      Xin <- t(Xin)
      Xout <- t(Xout)
    }

    n <- nrow(Xin)
    if (n != nrow(Xout)) {
      stop("Xin and Xout must have the same number of rows")
    }

    triplets <- get_triplet_matrix(n, n_triplets, zero_index = TRUE)
    n_triplets <- nrow(triplets) / 2

    in_dfun <- get_dfun(in_metric)
    out_dfun <- get_dfun(out_metric)
    acc <- 0
    for (i in 1:ncol(triplets)) {
      for (j in 1:n_triplets) {
        c1 <- j * 2 - 1
        p1 <- triplets[c1, i]
        p1 <- ifelse(p1 == i, p1 + 1, p1)

        p2 <- triplets[c1 + 1, i]
        p2 <- ifelse(p2 == i, p2 + 1, p2)

        if (i > nrow(Xin) || p1 > nrow(Xin) || p2 > nrow(Xin)) {
          stop("Problem with indexing: ", paste(i, p1, p2, collapse = " "))
        }
        dip1_in <- in_dfun(Xin[i, ], Xin[p1, ])
        dip2_in <- in_dfun(Xin[i, ], Xin[p2, ])

        dip1_out <- out_dfun(Xout[i, ], Xout[p1, ])
        dip2_out <- out_dfun(Xout[i, ], Xout[p2, ])

        # if (i == 1) {
        # message(i, " ",
        #         p1, ": ", dip1_in, " ", p2, ": ", dip2_in,
        #         p1, ": ", dip1_out, " ", p2, ": ", dip2_out)
        # }


        if ((dip1_in < dip2_in) == (dip1_out < dip2_out)) {
          acc <- acc + 1
        }
      }
    }
    acc / (n * n_triplets)
  }
