test_that("random pair distance emd", {
  set.seed(42)
  expect_equal(
    random_pair_distance_emd(m, n, n_pairs = 100000),
    0.13,
    tolerance = 0.01
  )
  set.seed(42)
  expect_equal(
    random_pair_distance_emd(
      t(m),
      t(n),
      n_pairs = 100000,
      is_transposed = TRUE
    ),
    0.13,
    tolerance = 0.01
  )
  set.seed(42)
  expect_equal(
    random_pair_distance_emd(m, n, n_pairs = 100000, n_threads = 2),
    0.13,
    tolerance = 0.01
  )
  set.seed(42)
  expect_equal(
    random_pair_distance_emd(m, n, n_pairs = 100000, range_scale = FALSE),
    1.035,
    tolerance = 0.02
  )
})

test_that("constant random pair distances have zero EMD", {
  x <- matrix(1, nrow = 4, ncol = 2)

  expect_equal(random_pair_distance_emd(x, x, n_pairs = 10), 0)
})

test_that("random pair distance stress matches sampled distances", {
  set.seed(987)
  randlist <- random_pair_distances(m, n, n_pairs = 50)
  expected <- sqrt(
    mean((scale01(randlist$din) - scale01(randlist$dout))^2)
  )

  set.seed(987)
  expect_equal(
    random_pair_distance_stress(m, n, n_pairs = 50),
    expected
  )

  set.seed(987)
  expect_equal(
    random_pair_distance_stress(m, n, n_pairs = 50, range_scale = FALSE),
    sqrt(mean((randlist$din - randlist$dout)^2))
  )
})

test_that("random pair distance stress handles identity and scale", {
  x <- matrix(c(0, 1, 4, 10, 20), ncol = 1)

  expect_equal(random_pair_distance_stress(x, x, n_pairs = 20), 0)

  set.seed(42)
  expect_equal(random_pair_distance_stress(x, x * 2, n_pairs = 100), 0)

  set.seed(42)
  expect_gt(
    random_pair_distance_stress(
      x,
      x * 2,
      n_pairs = 100,
      range_scale = FALSE
    ),
    0
  )
})

test_that("random pair distance stress sampling is reproducible", {
  set.seed(2024)
  serial_1 <- random_pair_distance_stress(m, n, n_pairs = 50, n_threads = 0)
  set.seed(2024)
  serial_2 <- random_pair_distance_stress(m, n, n_pairs = 50, n_threads = 0)
  expect_equal(serial_2, serial_1)

  set.seed(2024)
  parallel_1 <- random_pair_distance_stress(m, n, n_pairs = 50, n_threads = 2)
  set.seed(2024)
  parallel_2 <- random_pair_distance_stress(m, n, n_pairs = 50, n_threads = 2)
  expect_equal(parallel_2, parallel_1)
})

test_that("random pair distance stress validates inputs", {
  expect_error(
    random_pair_distance_stress(
      m[1, , drop = FALSE],
      n[1, , drop = FALSE]
    ),
    "at least 2 observations"
  )
  expect_error(
    random_pair_distance_stress(m, n, n_pairs = 0),
    "n_pairs must be a positive integer"
  )
  expect_error(
    random_pair_distance_stress(m, n, metric_out = "not-a-metric"),
    "should be one of"
  )
})
