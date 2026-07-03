test_that("random pair distance emd", {
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
  set.seed(42)
  expect_equal(
    random_pair_distance_emd(m, n, n_pairs = 100000, range_scale = FALSE),
    1.06,
    tolerance = 0.02
  )
})

test_that("constant random pair distances have zero EMD", {
  x <- matrix(1, nrow = 4, ncol = 2)

  expect_equal(scale01(rep(3, 4)), rep(0, 4))
  expect_equal(random_pair_distance_emd(x, x, n_pairs = 10), 0)
})
