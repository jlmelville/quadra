test_that("label trustworthiness and continuity distinguish losses and gains", {
  x <- cbind(c(0, 2, 10, 12, 20, 22))
  y <- cbind(c(0, 2, 0, 2, 20, 22))
  labels <- factor(c("C", "C", "A", "A", "B", "B"))
  same <- label_trustworthiness_continuity(x, x, labels)
  expect_equal(same, c(trustworthiness = 1, continuity = 1))

  result <- label_trustworthiness_continuity(x, y, labels, ret_extra = TRUE)
  expect_identical(as.character(result$pairs$label_a), c("C", "C", "A"))
  expect_identical(as.character(result$pairs$label_b), c("A", "B", "B"))
  expect_equal(result$pairs$consistency_in, c(1, 1, 1))
  expect_equal(result$pairs$consistency_out, c(0.5, 1, 1))
  expect_equal(result$pairs$change, c(-0.5, 0, 0))
  expect_equal(result$scores, c(trustworthiness = 5 / 6, continuity = 1))
  expect_equal(
    label_trustworthiness_continuity(y, x, labels),
    c(trustworthiness = 1, continuity = 5 / 6)
  )
})

test_that("coincident centroids give half credit at unequal class sizes", {
  x <- cbind(c(0, 2, 12, 13, 14))
  y <- cbind(c(0, 2, 0, 1, 2))
  labels <- c(1, 1, 2, 2, 2)
  result <- label_trustworthiness_continuity(x, y, labels, ret_extra = TRUE)
  expect_equal(result$pairs$consistency_out, 0.5)
  shuffled <- c(5, 2, 3, 4, 1)
  expect_equal(
    label_trustworthiness_continuity(
      x[shuffled, , drop = FALSE],
      y[shuffled, , drop = FALSE],
      labels[shuffled]
    ),
    result$scores
  )
  expect_equal(
    label_trustworthiness_continuity(x, cbind(2 * x + 10, 3), labels),
    c(trustworthiness = 1, continuity = 1)
  )
})

test_that("label trustworthiness and continuity validate corresponding rows", {
  x <- cbind(1:3)
  expect_error(label_trustworthiness_continuity(x, x, 1:2), "number of rows")
  expect_error(label_trustworthiness_continuity(x, x, c(1, 2, NA)), "missing")
  expect_error(
    label_trustworthiness_continuity(x, x, rep(1, 3)),
    "two distinct"
  )
  expect_error(
    label_trustworthiness_continuity(x, cbind(1:2), 1:3),
    "same number"
  )
})
