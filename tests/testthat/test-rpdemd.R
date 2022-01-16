test_that("random pair distance emd", {
  skip_if_not_installed("transport")

  set.seed(42)
  expect_equal(random_pair_distance_emd(m, n, n_pairs = 100000),
    0.13,
    tolerance = 0.01
  )
  set.seed(42)
  expect_equal(
    random_pair_distance_emd(t(m), t(n), n_pairs = 100000, is_transposed = TRUE),
    0.13,
    tolerance = 0.01
  )
  set.seed(42)
  expect_equal(random_pair_distance_emd(m, n, n_pairs = 100000, n_threads = 2),
    0.13,
    tolerance = 0.01
  )
})
