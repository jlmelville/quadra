test_that("neighbor preservation", {
  named <- function(x, name) {
    names(x) <- name
    x
  }

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

  expect_equal(
    nn_preservation(
      m,
      n,
      k = 2,
      nn_method_in = "brute",
      nn_method_out = "brute"
    ),
    named(0.55, "nnp2")
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
    named(0.55, "nnp2")
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
  expect_equal(
    nn_preservation(
      cached_nn$nn_in,
      cached_nn$nn_out,
      k = 2,
      nn_method_in = "brute",
      nn_method_out = "brute"
    ),
    named(0.55, "nnp2")
  )

  expect_equal(
    nn_preservation(
      m,
      n,
      k = c(2, 5),
      nn_method_in = "brute",
      nn_method_out = "brute"
    ),
    named(c(0.55, 0.48), c("nnp2", "nnp5"))
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
    named(c(0.55, 0.48), c("nnp2", "nnp5"))
  )
})
