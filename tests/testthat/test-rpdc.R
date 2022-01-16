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
