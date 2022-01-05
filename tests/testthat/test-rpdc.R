test_that("random pair distance correlation", {
  n <- matrix(c(
    -0.3673, -1.595,
    1.287, 0.5113,
    0.7557, -0.2293,
    1.519, -0.4499,
    0.4918, 1.085,
    1.716, -2.044,
    -0.1086, 2.228,
    -0.04816, 0.4036,
    0.2425, -0.2668,
    -1.734, 0.395
  ), byrow = TRUE, nrow = 10)

  m <- matrix(c(
    1.506, 0.05989, 0.3542,
    0.8419, -1.146, -0.2262,
    1.031, 1.457, 0.5978,
    0.4838, -0.005043, 1.771,
    -0.3932, -0.01158, 1.47,
    -0.5604, -0.3132, -0.502,
    2.624, 1.509, -0.4661,
    0.1407, -1.111, 0.1571,
    0.3512, -1.103, -0.5345,
    -0.4729, 1.453, 1.087
  ), byrow = TRUE, nrow = 10)

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
