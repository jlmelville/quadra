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

test_that("random pair inputs are validated", {
  expect_error(
    random_pair_distance_correlation(
      m[1, , drop = FALSE],
      n[1, , drop = FALSE]
    ),
    "at least 2 observations"
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
})

test_that("constant random pair distances use stats cor behavior", {
  x <- matrix(1, nrow = 4, ncol = 2)

  expect_warning(
    res <- random_pair_distance_correlation(x, x, n_pairs = 10),
    "standard deviation is zero"
  )
  expect_true(is.na(res))
})
