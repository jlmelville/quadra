test_that("neighbor preservation", {
  named <- function(x, name) {
    names(x) <- name
    x
  }

  # fmt: skip
  n <- matrix(
    c(
      -0.3673, -1.595,
       1.2870,  0.5113,
       0.7557, -0.2293,
       1.5190, -0.4499,
       0.4918,  1.0850,
       1.7160, -2.0440,
      -0.1086,  2.2280,
      -0.04816, 0.4036,
       0.2425, -0.2668,
      -1.7340,  0.3950
    ),
    byrow = TRUE,
    nrow = 10
  )

  # fmt: skip
  m <- matrix(
    c(
       1.5060,  0.05989,  0.3542,
       0.8419, -1.14600, -0.2262,
       1.0310,  1.45700,  0.5978,
       0.4838, -0.005043, 1.7710,
      -0.3932, -0.01158,  1.4700,
      -0.5604, -0.31320, -0.5020,
       2.6240,  1.50900, -0.4661,
       0.1407, -1.11100,  0.1571,
       0.3512, -1.10300, -0.5345,
      -0.4729,  1.45300,  1.0870
    ),
    byrow = TRUE,
    nrow = 10
  )

  expect_equal(
    nn_preservation(
      m,
      n,
      k = 2,
      nn_method_in = "brute",
      nn_method_out = "brute"
    ),
    named(0.15, "nnp2")
  )
  expect_equal(
    nn_preservation(
      m,
      n,
      k = 2,
      nn_method_in = "brute",
      nn_method_out = "brute",
      n_threads = 2
    ),
    named(0.15, "nnp2")
  )

  cached_nn <-
    nn_preservation(
      m,
      n,
      k = 2,
      nn_method_in = "brute",
      nn_method_out = "brute",
      ret_extra = TRUE
    )
  expect_false(any(cached_nn$nn_in$idx == row(cached_nn$nn_in$idx)))
  expect_false(any(cached_nn$nn_out$idx == row(cached_nn$nn_out$idx)))
  expect_equal(ncol(cached_nn$nn_in$idx), 2)
  expect_equal(ncol(cached_nn$nn_out$idx), 2)
  expect_equal(
    nn_preservation(
      cached_nn$nn_in,
      cached_nn$nn_out,
      k = 2,
      nn_method_in = "brute",
      nn_method_out = "brute"
    ),
    named(0.15, "nnp2")
  )

  expect_equal(
    nn_preservation(
      m,
      n,
      k = c(2, 5),
      nn_method_in = "brute",
      nn_method_out = "brute"
    ),
    named(c(0.15, 0.5), c("nnp2", "nnp5"))
  )

  expect_equal(
    nn_preservation(
      t(m),
      t(n),
      k = c(2, 5),
      nn_method_in = "brute",
      nn_method_out = "brute",
      is_transposed = TRUE
    ),
    named(c(0.15, 0.5), c("nnp2", "nnp5"))
  )

  expect_equal(
    nn_preservation(
      as.data.frame(m),
      as.data.frame(n),
      k = 2,
      nn_method_in = "brute",
      nn_method_out = "brute"
    ),
    named(0.15, "nnp2")
  )
})

test_that("generated nearest-neighbor graphs exclude self-neighbors", {
  x <- matrix(c(0, 1, 3, 10), ncol = 1)

  brute_res <- nn_preservation(
    x,
    x,
    k = 2,
    nn_method_in = "brute",
    nn_method_out = "brute",
    ret_extra = TRUE
  )
  expect_false(any(brute_res$nn_in$idx == row(brute_res$nn_in$idx)))
  expect_false(any(brute_res$nn_out$idx == row(brute_res$nn_out$idx)))
  expect_equal(ncol(brute_res$nn_in$idx), 2)
  expect_equal(ncol(brute_res$nn_out$idx), 2)

  set.seed(1337)
  nnd_res <- nn_preservation(
    m,
    n,
    k = 2,
    nn_method_in = "nnd",
    nn_method_out = "nnd",
    ret_extra = TRUE
  )
  expect_false(any(nnd_res$nn_in$idx == row(nnd_res$nn_in$idx)))
  expect_false(any(nnd_res$nn_out$idx == row(nnd_res$nn_out$idx)))
  expect_equal(ncol(nnd_res$nn_in$idx), 2)
  expect_equal(ncol(nnd_res$nn_out$idx), 2)
})

test_that("cached self-excluded graphs reproduce raw-data values", {
  raw_res <- nn_preservation(
    m,
    n,
    k = c(2, 5),
    nn_method_in = "brute",
    nn_method_out = "brute",
    ret_extra = TRUE
  )

  expect_equal(
    nn_preservation(raw_res$nn_in, raw_res$nn_out, k = c(2, 5)),
    raw_res$nnp
  )
})

test_that("local radius correlation is one for preserved local scales", {
  x <- matrix(c(0, 1, 4, 10, 20), ncol = 1)

  expect_equal(
    local_radius_correlation(
      x,
      x,
      k = c(1, 2),
      nn_method_in = "brute",
      nn_method_out = "brute",
      method = "pearson"
    ),
    c(lrc1 = 1, lrc2 = 1)
  )

  cached <- local_radius_correlation(
    x,
    x * 2,
    k = 2,
    nn_method_in = "brute",
    nn_method_out = "brute",
    ret_extra = TRUE
  )

  expect_equal(cached$lrc, c(lrc2 = 1))
  expect_equal(
    local_radius_correlation(cached$nn_in, cached$nn_out, k = 2),
    cached$lrc
  )
  expect_equal(ncol(cached$scale_in), 1)
  expect_equal(ncol(cached$scale_out), 1)
})

test_that("local radius correlation supports graph distances and statistics", {
  # fmt: skip
  idx <- matrix(
    c(
      2, 3,
      1, 3,
      2, 4,
      3, 2
    ),
    nrow = 4,
    byrow = TRUE
  )
  # fmt: skip
  dist_in <- matrix(
    c(
      1, 4,
      1, 2,
      2, 5,
      3, 8
    ),
    nrow = 4,
    byrow = TRUE
  )
  # fmt: skip
  dist_out <- matrix(
    c(
      2, 2,
      1, 3,
      4, 4,
      2, 6
    ),
    nrow = 4,
    byrow = TRUE
  )
  graph_in <- list(idx = idx, dist = dist_in)
  graph_out <- list(idx = idx, dist = dist_out)

  expect_equal(
    local_radius_correlation(
      graph_in,
      graph_out,
      k = 1,
      method = "pearson"
    ),
    c(lrc1 = unname(stats::cor(dist_in[, 1], dist_out[, 1])))
  )
  expect_equal(
    local_radius_correlation(
      graph_in,
      graph_out,
      k = 2,
      statistic = "mean",
      method = "pearson"
    ),
    c(
      lrc2 = unname(
        stats::cor(rowMeans(dist_in), rowMeans(dist_out))
      )
    )
  )
})

test_that("local radius correlation validates graph distances", {
  # fmt: skip
  idx <- matrix(
    c(
      2, 3,
      1, 3,
      1, 2
    ),
    nrow = 3,
    byrow = TRUE
  )
  graph <- list(idx = idx, dist = matrix(1, nrow = 3, ncol = 2))

  expect_error(
    local_radius_correlation(list(idx = idx), graph, k = 1),
    "dist"
  )
  expect_error(
    local_radius_correlation(
      list(idx = idx, dist = matrix(NA_real_, nrow = 3, ncol = 2)),
      graph,
      k = 1
    ),
    "finite distances"
  )
})

test_that("local radius correlation handles constant and zero radii", {
  # fmt: skip
  idx <- matrix(
    c(
      2, 3,
      1, 3,
      1, 2
    ),
    nrow = 3,
    byrow = TRUE
  )
  constant_graph <- list(idx = idx, dist = matrix(1, nrow = 3, ncol = 2))

  expect_equal(
    local_radius_correlation(constant_graph, constant_graph, k = 1),
    c(lrc1 = NA_real_)
  )

  x <- matrix(c(0, 0, 2, 5), ncol = 1)
  expect_equal(
    local_radius_correlation(
      x,
      x,
      k = 1,
      nn_method_in = "brute",
      nn_method_out = "brute"
    ),
    c(lrc1 = 1)
  )
  expect_error(
    local_radius_correlation(
      x,
      x,
      k = 1,
      log = TRUE,
      nn_method_in = "brute",
      nn_method_out = "brute"
    ),
    "positive"
  )
})

test_that("old self-inclusive cached graphs are stripped with a warning", {
  old_graph <- rnndescent::brute_force_knn(m, k = 3, metric = "sqeuclidean")
  stripped_graph <- list(
    idx = old_graph$idx[, -1, drop = FALSE],
    dist = old_graph$dist[, -1, drop = FALSE]
  )

  expect_warning(
    res <- nn_preservation(old_graph, stripped_graph, k = 2),
    "contains self-neighbors"
  )
  expect_equal(res, c(nnp2 = 1))
})

test_that("idx-only nearest-neighbor graphs are accepted", {
  named <- function(x, name) {
    names(x) <- name
    x
  }

  # fmt: skip
  idx <- matrix(
    c(
      2, 3,
      1, 3,
      1, 2
    ),
    nrow = 3,
    byrow = TRUE
  )
  graph <- list(idx = idx)

  expect_equal(nn_preservation(graph, graph, k = 1), named(1, "nnp1"))
})

test_that("nearest-neighbor preservation returns per-row values", {
  # fmt: skip
  idx <- matrix(
    c(
      3, 2, 4,
      3, 4, 1,
      1, 2, 4,
      3, 1, 2
    ),
    nrow = 4,
    byrow = TRUE
  )
  # fmt: skip
  ref_idx <- matrix(
    c(
      2, 3, 4,
      1, 3, 4,
      4, 1, 2,
      3, 2, 1
    ),
    nrow = 4,
    byrow = TRUE
  )

  res <- nn_preservation(
    list(idx = ref_idx),
    list(idx = idx),
    k = c(2, 1, 2),
    ret_extra = TRUE
  )

  expect_equal(res$nnp, c(nnp2 = 0.625, nnp1 = 0.25, nnp2 = 0.625))
  expect_equal(res$nnpv$nnp2, c(1, 0.5, 0.5, 0.5))
  expect_equal(res$nnpv$nnp1, c(0, 0, 0, 1))
})

test_that("nearest-neighbor overlap counts preserve requested k order", {
  # fmt: skip
  idx <- matrix(
    c(
      3, 2, 4,
      3, 4, 1,
      1, 2, 4,
      3, 1, 2
    ),
    nrow = 4,
    byrow = TRUE
  )
  # fmt: skip
  ref_idx <- matrix(
    c(
      2, 3, 4,
      1, 3, 4,
      4, 1, 2,
      3, 2, 1
    ),
    nrow = 4,
    byrow = TRUE
  )

  expect_equal(
    nn_overlap_counts(idx, ref_idx, k = c(2, 1, 2)),
    matrix(
      c(
        2,
        1,
        1,
        1,
        0,
        0,
        0,
        1,
        2,
        1,
        1,
        1
      ),
      nrow = 4
    )
  )
})

test_that("nearest-neighbor overlap counts unique shared indices", {
  # fmt: skip
  idx <- matrix(
    c(
      2, 2, 3,
      1, 1, 3,
      1, 1, 2,
      1, 1, 2
    ),
    nrow = 4,
    byrow = TRUE
  )
  # fmt: skip
  ref_idx <- matrix(
    c(
      3, 4, 2,
      3, 4, 1,
      2, 4, 1,
      2, 3, 1
    ),
    nrow = 4,
    byrow = TRUE
  )

  expect_equal(
    nn_overlap_counts(idx, ref_idx, k = c(2, 3)),
    matrix(c(0, 0, 0, 0, 2, 2, 2, 2), nrow = 4)
  )
})

test_that("nearest-neighbor inputs are validated", {
  # fmt: skip
  idx <- matrix(
    c(
      2, 3,
      1, 3,
      1, 2
    ),
    nrow = 3,
    byrow = TRUE
  )
  graph <- list(idx = idx)

  expect_error(
    nn_preservation(graph, graph, k = 0),
    "k must contain positive integers"
  )
  expect_error(
    nn_preservation(graph, graph, k = 1.5),
    "k must contain positive integers"
  )
  expect_error(
    nn_preservation(graph, graph, k = 3),
    "non-self observations"
  )
  expect_error(
    nn_preservation(m[1:3, , drop = FALSE], n[1:3, , drop = FALSE], k = 3),
    "non-self observations"
  )
  expect_error(
    nn_preservation(list(idx = 1), graph, k = 1),
    "nearest-neighbor graph"
  )
  expect_error(
    nn_preservation(
      list(idx = idx, dist = idx[1, , drop = FALSE]),
      graph,
      k = 1
    ),
    "nearest-neighbor graph"
  )
})
