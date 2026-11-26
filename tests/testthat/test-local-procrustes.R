test_that("local Procrustes follows reference neighbor identities", {
  x <- rbind(c(0, 0), c(1, 0), c(2, 1), c(0.4, 2))
  x <- rbind(x, sweep(x, 2, c(100, 20), "+"))
  rownames(x) <- letters[seq_len(nrow(x))]
  rigid <- sweep(cbind(-x[, 2], x[, 1]), 2, c(20, 17), "+")
  result <- local_procrustes(x, rigid, k = 3, ret_extra = TRUE)
  expect_equal(result$error, 0, tolerance = 1e-14)
  expect_identical(names(result$pointwise_error), rownames(x))
  expect_equal(result$n_defined, nrow(x))
  expect_equal(mean(result$pointwise_error), result$error)

  moved <- x
  moved[5:8, ] <- sweep(x[5:8, ], 2, c(40, 10), "+")
  expect_equal(
    local_procrustes(x, moved, k = 3, nn_in = result$nn_in),
    0,
    tolerance = 1e-14
  )

  swapped <- x
  swapped[c(2, 6), ] <- x[c(6, 2), ]
  expect_equal(sort(dist(x)), sort(dist(swapped)))
  expect_gt(local_procrustes(x, swapped, k = 3, nn_in = result$nn_in), 1)
  expect_equal(local_procrustes(x, 2 * x, k = 3, nn_in = result$nn_in), 1)
})

test_that("wide, unequal dimensions give the same residual as explicit padding", {
  set.seed(4)
  x <- matrix(rnorm(88), nrow = 8)
  y <- matrix(rnorm(56), nrow = 8)
  result <- local_procrustes(x, y, k = 3, ret_extra = TRUE)
  padded <- cbind(y, matrix(0, nrow(y), ncol(x) - ncol(y)))
  expected <- vapply(
    seq_len(nrow(x)),
    function(i) {
      ids <- c(i, result$nn_in$idx[i, ])
      a <- scale(x[ids, ], scale = FALSE)
      b <- scale(padded[ids, ], scale = FALSE)
      fit <- svd(crossprod(b, a))
      sum((a - b %*% fit$u %*% t(fit$v))^2) / sum(a^2)
    },
    numeric(1)
  )
  expect_equal(result$pointwise_error, expected, tolerance = 1e-12)
  expect_equal(
    local_procrustes(t(x), t(y), k = 3, is_transposed = TRUE),
    mean(expected)
  )
})

test_that("constant reference patches are undefined", {
  x <- cbind(c(0, 0, 1, 2))
  graph <- list(idx = cbind(c(2L, 1L, 4L, 3L)))
  result <- local_procrustes(x, x, k = 1, nn_in = graph, ret_extra = TRUE)
  expect_equal(result$pointwise_error, c(NA_real_, NA_real_, 0, 0))
  expect_identical(result$n_defined, 2L)
  expect_equal(result$error, 0)
  expect_identical(local_procrustes(x * 0, x, k = 1, nn_in = graph), NA_real_)
  expect_error(local_procrustes(x, x[-1, , drop = FALSE], k = 1), "same number")
  expect_error(
    local_procrustes(x, x, k = 1, nn_in = list(idx = rbind(2L, 1L))),
    "one row"
  )
})
