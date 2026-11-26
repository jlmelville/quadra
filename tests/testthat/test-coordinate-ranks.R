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
