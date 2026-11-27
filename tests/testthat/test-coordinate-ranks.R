test_that("row-wise coordinate ranks reproduce dense profiles and their tie order", {
  # Repeated coordinates and tied distances exercise both kinds of ties.
  x <- cbind(c(0, 1, -1, 2, 0, 3, 4), c(0, 0, 0, 1, 0, 2, 1)) + 1e12
  rownames(x) <- letters[seq_len(nrow(x))]
  y <- cbind(as.integer(c(2, 0, 1, 3, 2, 1, 4)))
  ks <- c(3, 1)
  for (rows in list(seq_len(nrow(x)), c(5, 3, 1, 2, 4, 7, 6))) {
    input <- x[rows, , drop = FALSE]
    maps <- list(projected = y[rows, , drop = FALSE], identity = input)
    result <- exact_rank_profiles(input, maps, ks)
    expect_identical(names(result), names(maps))
    for (m in names(maps)) {
      din <- as.matrix(dist(input))
      dout <- as.matrix(dist(maps[[m]]))
      expected <- c(
        pointwise_rank_preservation(din, dout, ks),
        pointwise_rnx(din, dout, ks)
      )
      expect_equal(result[[m]], expected)
    }
  }
})

test_that("coordinate profiles use all dimensions and reject incompatible inputs", {
  x <- as.matrix(iris[1:12, 1:4])
  result <- exact_rank_profiles(x, list(iris[1:12, ]), 3)[[1]]
  expect_equal(as.numeric(result$auc), rep(1, nrow(x)))
  expect_error(exact_rank_profiles(x, list(), 1), "nonempty list")
  expect_error(exact_rank_profiles(x, list(x[-1, ]), 1), "same rows")
  expect_error(exact_rank_profiles(x, list(x), 6), "less than half")
  expect_error(exact_rank_profiles(x, list(x * Inf), 1), "finite")
  expect_error(exact_rank_profiles(x, list(x * 1e200), 1), "Nonfinite computed")
})

test_that("directional MRRE matches a direct calculation with first ties", {
  x <- cbind(c(0, 0, 1, 3, 4, 6, 9, 10, 8))
  rownames(x) <- letters[seq_len(nrow(x))]
  y <- cbind(c(0, 1, 0, 6, 3, 8, 9, 10, 4))
  ks <- c(4, 1, 3)
  maps <- list(moved = y, identity = x)
  profiles <- exact_rank_profiles(x, maps, ks, mrre = TRUE)
  baseline <- exact_rank_profiles(x, maps, ks)
  din <- as.matrix(dist(x))
  dout <- as.matrix(dist(y))

  direct <- function(a, b, k) {
    r <- seq_len(k)
    divisor <- sum((nrow(a) - 2 * r) / r)
    vapply(
      seq_len(nrow(a)),
      function(i) {
        ra <- rank(a[i, -i], ties.method = "first")
        rb <- rank(b[i, -i], ties.method = "first")
        sum(abs(rb[order(ra)[r]] - r) / r) / divisor
      },
      numeric(1)
    )
  }
  for (m in names(maps)) {
    expect_identical(profiles[[m]][names(baseline[[m]])], baseline[[m]])
    for (field in c("mrre_input", "mrre_output")) {
      expect_identical(
        dimnames(profiles[[m]][[field]]),
        dimnames(baseline[[m]]$continuity)
      )
    }
  }
  for (j in seq_along(ks)) {
    expect_equal(
      unname(profiles$moved$mrre_input[, j]),
      direct(din, dout, ks[j])
    )
    expect_equal(
      unname(profiles$moved$mrre_output[, j]),
      direct(dout, din, ks[j])
    )
  }
  expect_true(all(profiles$identity$mrre_input == 0))
  expect_true(all(profiles$identity$mrre_output == 0))
  expect_error(exact_rank_profiles(x, maps, ks, mrre = NA), "TRUE or FALSE")
})

test_that("MRRE uses reverse-order normalization without clipping", {
  x <- cbind(0:8)
  maps <- list(reversed = cbind(c(0, 8:1)), moved = cbind(c(0, 8:6, 1, 5:2)))
  p <- exact_rank_profiles(x, maps, 4, mrre = TRUE)
  expect_equal(p$reversed$mrre_input[[1]], 1)
  expect_equal(p$reversed$mrre_output[[1]], 1)
  expect_gt(p$moved$mrre_input[1, 1], 1)
})
