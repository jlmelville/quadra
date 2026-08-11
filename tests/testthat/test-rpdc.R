test_that("random pair distance correlation", {
  set.seed(1337)
  rdpc <- random_pair_distance_correlation(m, n)
  expect_lte(rdpc, 1.0)
  expect_gte(rdpc, -1.0)

  set.seed(1337)
  expect_equal(random_pair_distance_correlation(m, n), rdpc)

  set.seed(1337)
  expect_equal(
    random_pair_distance_correlation(t(m), t(n), is_transposed = TRUE),
    rdpc
  )

  expect_equal(random_pair_distance_correlation(m, m), 1.0)

  rdpct2 <- random_pair_distance_correlation(m, n, n_threads = 2)
  expect_lte(rdpct2, 1.0)
  expect_gte(rdpct2, -1.0)
})

test_that("random pair distance correlation supports correlation methods", {
  set.seed(987)
  randlist <- random_pair_distances(m, n, n_pairs = 50)

  set.seed(987)
  expect_equal(
    random_pair_distance_correlation(m, n, n_pairs = 50),
    unname(stats::cor(randlist$din, randlist$dout, method = "pearson"))
  )

  set.seed(987)
  expect_equal(
    random_pair_distance_correlation(
      m,
      n,
      n_pairs = 50,
      method = "spearman"
    ),
    unname(stats::cor(randlist$din, randlist$dout, method = "spearman"))
  )
})

test_that("random pair distance correlation supports documented distance metrics", {
  # fmt: skip
  x <- matrix(
    c(
      1, 2, 4,
      1, 2, 5,
      1, 3, 5,
      2, 3, 6,
      2, 4, 6,
      3, 4, 7
    ),
    nrow = 6,
    byrow = TRUE
  )

  for (metric in c(
    "sqeuclidean",
    "euclidean",
    "cosine",
    "hamming",
    "manhattan",
    "correlation"
  )) {
    set.seed(2468)
    expect_equal(
      random_pair_distance_correlation(
        x,
        x,
        n_pairs = 200,
        metric_in = metric,
        metric_out = metric
      ),
      1,
      tolerance = 1e-12,
      info = metric
    )
  }
})

test_that("unknown distance metrics error", {
  expect_error(
    random_pair_distance_correlation(m, n, metric_in = "not-a-metric"),
    "should be one of"
  )
})

test_that("unknown random pair correlation methods error", {
  expect_error(
    random_pair_distance_correlation(m, n, method = "kendall"),
    "should be one of"
  )
})

test_that("random pair sampling is reproducible for same thread count", {
  set.seed(2024)
  serial_1 <- random_pair_distance_correlation(
    m,
    n,
    n_pairs = 50,
    n_threads = 0
  )
  set.seed(2024)
  serial_2 <- random_pair_distance_correlation(
    m,
    n,
    n_pairs = 50,
    n_threads = 0
  )
  expect_equal(serial_2, serial_1)

  set.seed(2024)
  parallel_1 <- random_pair_distance_correlation(
    m,
    n,
    n_pairs = 50,
    n_threads = 2
  )
  set.seed(2024)
  parallel_2 <- random_pair_distance_correlation(
    m,
    n,
    n_pairs = 50,
    n_threads = 2
  )
  expect_equal(parallel_2, parallel_1)
})

test_that("random pair sampling excludes self-pairs", {
  x <- matrix(c(0, 1), ncol = 1)

  set.seed(20260811)
  sampled <- random_pair_distances(x, x, n_pairs = 20)

  expect_equal(sampled$din, rep(1, 20))
  expect_equal(sampled$dout, rep(1, 20))
})

test_that("random pair inputs are validated", {
  expect_error(
    random_pair_distance_correlation(
      m[1, , drop = FALSE],
      n[1, , drop = FALSE]
    ),
    "at least 2 observations"
  )
  expect_error(
    random_pair_distance_correlation(m, n[-1, ]),
    "same number of observations"
  )
  expect_error(
    random_pair_distance_correlation(m, n, n_pairs = 0),
    "n_pairs must be a positive integer"
  )
  expect_error(
    random_pair_distance_correlation(m, n, n_pairs = 1.5),
    "n_pairs must be a positive integer"
  )
  expect_error(
    random_pair_distance_correlation(
      data.frame(label = letters[1:3]),
      data.frame(label = letters[1:3])
    ),
    "at least one numeric column"
  )
  expect_error(
    random_pair_distance_correlation(m, n, is_transposed = NA),
    "is_transposed must be TRUE or FALSE"
  )
})

test_that("constant random pair distances use stats cor behavior", {
  x <- matrix(1, nrow = 4, ncol = 2)

  expect_warning(
    res <- random_pair_distance_correlation(x, x, n_pairs = 10),
    "standard deviation is zero"
  )
  expect_true(is.na(res))
})

test_that("dense-only sampled metrics reject sparse top-level inputs", {
  sampled_metrics <- list(
    correlation = function(Xin, Xout) {
      random_pair_distance_correlation(Xin, Xout, n_pairs = 10)
    },
    emd = function(Xin, Xout) {
      random_pair_distance_emd(Xin, Xout, n_pairs = 10)
    },
    stress = function(Xin, Xout) {
      random_pair_distance_stress(Xin, Xout, n_pairs = 10)
    },
    triplet = function(Xin, Xout) {
      random_triplet_accuracy(Xin, Xout, n_triplets = 2)
    }
  )

  for (metric_name in names(sampled_metrics)) {
    metric <- sampled_metrics[[metric_name]]
    for (input_name in c("Xin", "Xout")) {
      inputs <- list(Xin = m, Xout = n)
      inputs[[input_name]] <- sparse_m
      expect_error(
        do.call(metric, inputs),
        paste0(input_name, ".*sparse"),
        info = paste(metric_name, input_name)
      )
    }
  }
})

test_that("dense-only sampled metrics reject nonfinite top-level inputs", {
  explicit_triplets <- matrix(0, nrow = 2, ncol = nrow(m))
  cases <- list(
    list(
      metric = random_pair_distance_correlation,
      input = "Xin",
      value = NA_real_,
      args = list(n_pairs = 10)
    ),
    list(
      metric = random_pair_distance_correlation,
      input = "Xout",
      value = NaN,
      args = list(n_pairs = 10)
    ),
    list(
      metric = random_pair_distance_emd,
      input = "Xin",
      value = Inf,
      args = list(n_pairs = 10)
    ),
    list(
      metric = random_pair_distance_emd,
      input = "Xout",
      value = -Inf,
      args = list(n_pairs = 10)
    ),
    list(
      metric = random_pair_distance_stress,
      input = "Xin",
      value = NaN,
      args = list(n_pairs = 10)
    ),
    list(
      metric = random_pair_distance_stress,
      input = "Xout",
      value = Inf,
      args = list(n_pairs = 10)
    ),
    list(
      metric = random_triplet_accuracy,
      input = "Xin",
      value = -Inf,
      args = list(n_triplets = 2)
    ),
    list(
      metric = random_triplet_accuracy,
      input = "Xout",
      value = NA_real_,
      args = list(n_triplets = explicit_triplets)
    )
  )

  for (case in cases) {
    inputs <- list(Xin = m, Xout = n)
    inputs[[case$input]][1, 1] <- case$value
    expect_error(
      do.call(case$metric, c(inputs, case$args)),
      paste0(case$input, ".*finite values"),
      info = case$input
    )
  }
})

test_that("sampled input checks ignore nonnumeric data-frame columns", {
  xin <- data.frame(value = m[, 1], ignored = NA_character_)
  xout <- data.frame(value = n[, 1], ignored = NA_character_)

  set.seed(20260810)
  expected <- random_pair_distance_correlation(
    m[, 1, drop = FALSE],
    n[, 1, drop = FALSE],
    n_pairs = 20
  )
  set.seed(20260810)
  expect_equal(
    random_pair_distance_correlation(xin, xout, n_pairs = 20),
    expected
  )
})
