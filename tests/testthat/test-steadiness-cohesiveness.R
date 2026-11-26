test_that("SNC computes reciprocal mean shared-neighbor similarities between groups", {
  sim <- rbind(c(0, .4, .3), c(.4, 0, .7), c(.3, .7, 0))
  sizes <- c(2, 1)
  result <- snc_group_distances(sim, c(1L, 1L, 2L), outer(sizes, sizes))
  expect_equal(unname(result[1, 2]), 1 / (.5 + .1))
  expect_equal(unname(result[2, 1]), 1 / (.5 + .1))
})

test_that("merging and splitting groups penalize opposite directions", {
  skip_if_not_installed("dbscan")
  skip_if_not_installed("Matrix")
  set.seed(450)
  clouds <- replicate(3, matrix(rnorm(160, sd = .45), 80, 2), simplify = FALSE)
  layout <- function(centers) {
    do.call(rbind, Map(function(z, v) sweep(z, 2, v, "+"), clouds, centers))
  }
  separate <- layout(list(c(-4, 0), c(4, 0), c(0, 7)))
  merged <- layout(list(c(-.1, 0), c(.1, 0), c(0, 7)))
  set.seed(451)
  merge <- steadiness_cohesiveness(separate, merged, ret_extra = TRUE)
  set.seed(451)
  expect_identical(steadiness_cohesiveness(separate, merged), merge$scores)
  expect_lt(merge$scores["steadiness"], merge$scores["cohesiveness"])
  set.seed(451)
  split <- steadiness_cohesiveness(merged, separate)
  expect_lt(split["cohesiveness"], split["steadiness"])
  expect_identical(nrow(merge$pointwise), nrow(separate))
  expect_true(all(is.finite(as.matrix(merge$pointwise))))
  expect_true(all(steadiness_cohesiveness(separate, separate) == 1))
  expect_error(
    steadiness_cohesiveness(separate, merged, k = 1),
    "increase k|Increase k"
  )
})
