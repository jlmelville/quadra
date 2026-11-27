test_that("Neighbor Dissimilarity uses both orders and separates gains from losses", {
  skip_if_not_installed("Matrix")
  a <- rbind(c(2, 3), c(1, 3), c(1, 4), c(2, 5), c(3, 4))
  b <- rbind(c(3, 4), c(3, 4), c(1, 2), c(1, 2), c(1, 2))
  shared <- function(idx) {
    n <- nrow(idx)
    outer(
      seq_len(n),
      seq_len(n),
      Vectorize(function(i, j) {
        length(intersect(idx[i, ], idx[j, ]))
      })
    )
  }
  delta <- (shared(b) - shared(a)) / 2
  gain <- sqrt(sum(pmax(delta, 0)^2))
  loss <- sqrt(sum(pmin(delta, 0)^2))
  expect_gt(gain, loss)
  f <- function(x, y, ...)
    neighbor_dissimilarity(list(idx = x), list(idx = y), k = 2, ...)
  out <- f(a, b, ret_extra = TRUE)
  expect_equal(out$nd, max(gain, loss))
  expect_equal(f(a, b), out$nd)
  expect_equal(out$gain, gain)
  expect_equal(out$loss, loss)
  expect_equal(f(b, a, ret_extra = TRUE)$loss, gain)
  ij <- as.matrix(out$pairs[, c("i", "j")])
  expect_true(all(ij[, 1] < ij[, 2]))
  expect_equal(out$pairs$shared_in, shared(a)[ij])
  expect_equal(out$pairs$shared_out, shared(b)[ij])
  expect_equal(nrow(unique(ij)), sum(upper.tri(delta) & delta != 0))
  expect_equal(f(a, a[, 2:1]), 0)
  expect_equal(nrow(f(a, a, ret_extra = TRUE)$pairs), 0)
  expect_warning(expect_equal(f(cbind(seq_len(nrow(a)), a), b), out$nd), "self")
  expect_error(f(a, rbind(b, c(1, 2))), "same number of rows")
})
