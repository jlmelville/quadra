test_that("random triplet accuracy", {
tm <- matrix(c(
  2, 5, 3, 4, 9, 4, 9, 3, 9, 6,
  4, 3, 4, 5, 8, 8, 8, 4, 4, 2,
  5, 0, 6, 2, 5, 3, 7, 1, 0, 7,
  6, 6, 8, 6, 6, 2, 5, 8, 3, 3
), byrow = TRUE, nrow = 4)

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

expect_equal(random_triplet_accuracy(m, n, tm), 0.45)
expect_equal(random_triplet_accuracy(t(m), t(n), tm, is_transposed = TRUE), 0.45)
expect_equal(random_triplet_accuracy(m, n, tm, n_threads = 2), 0.45)
})
