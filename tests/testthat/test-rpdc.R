test_that("random pair distance correlation", {
  set.seed(1337)
  rdpc <- random_pair_distance_correlation(m, n)
  expect_lte(rdpc, 1.0)
  expect_gte(rdpc, -1.0)

  set.seed(1337)
  expect_equal(random_pair_distance_correlation(m, n), rdpc)

  set.seed(1337)
  expect_equal(random_pair_distance_correlation(t(m), t(n),
    is_transposed = TRUE
  ), rdpc)

  expect_equal(random_pair_distance_correlation(m, m), 1.0)

  rdpct2 <- random_pair_distance_correlation(m, n, n_threads = 2)
  expect_lte(rdpct2, 1.0)
  expect_gte(rdpct2, -1.0)
})

test_that("unknown C++ distance metrics error", {
  expect_error(
    random_distances(t(m), t(n), metric_in = "not-a-metric"),
    "Unknown distance metric"
  )
})

test_that("random pair inputs are validated", {
  expect_error(
    random_pair_distance_correlation(m[1, , drop = FALSE], n[1, , drop = FALSE]),
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
