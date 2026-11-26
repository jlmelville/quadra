test_that("class triplets use distinct indices and the right class directions", {
  labels <- factor(c("one", "two", "one", "single", "one"))
  set.seed(42)
  plan <- sample_class_triplets(labels, 150)
  expect_identical(dim(plan), c(150L, 3L))
  expect_true(all(plan[, 2L] != plan[, 3L]))
  expect_true(all(labels[plan[, 1L]] != labels[plan[, 2L]]))
  expect_true(all(labels[plan[, 2L]] == labels[plan[, 3L]]))
  expect_setequal(plan[, 1L], c(2L, 4L))
  set.seed(42)
  expect_identical(sample_class_triplets(labels, 150), plan)
  expect_error(sample_class_triplets(1:3), "two classes")
})

test_that("class triplets are weighted by the eligible observation counts", {
  labels <- rep(c("small", "large"), c(2, 4))
  set.seed(63)
  plan <- sample_class_triplets(labels, 4000)
  expect_lt(abs(mean(labels[plan[, 1L]] == "small") - 0.75), 0.03)
})

test_that("CADI reports directed class pairs and agrees with a cosine oracle", {
  x <- rbind(c(0, 0), c(-1, 0), c(1, 0), c(0, 1))
  y <- rbind(c(0, 0), c(-1, 0), c(1, 0), c(2, 0))
  labels <- c("A", "A", "B", "B")
  plan <- rbind(c(1L, 3L, 4L), c(3L, 1L, 2L), c(3L, 2L, 1L))
  set.seed(17)
  before <- .Random.seed
  result <- class_angular_distortion_index(
    x,
    y,
    labels,
    triplets = plan,
    ret_extra = TRUE
  )
  expect_identical(.Random.seed, before)
  expect_equal(result$cosine_in, c(0, 1, 1))
  expect_equal(result$cosine_out, c(1, 1, 1))
  expect_identical(result$pairs$anchor_class, c("A", "B"))
  expect_identical(result$pairs$endpoint_class, c("B", "A"))
  expect_identical(result$pairs$n_triplets, c(1L, 2L))
  expect_equal(result$pairs$cadi, c(1 / 4, 0))
  expect_equal(
    result$cadi,
    weighted.mean(result$pairs$cadi, result$pairs$n_triplets)
  )

  rigid <- sweep(3 * x[, 2:1], 2L, c(1, 2), "+")
  expect_equal(
    class_angular_distortion_index(x, rigid, labels, triplets = plan),
    0
  )
})

test_that("CADI validates the plan and uses zero cosine for a collapsed ray", {
  x <- rbind(c(0, 0), c(1, 0), c(2, 0))
  y <- x
  y[2L, ] <- y[1L, ]
  plan <- rbind(c(1L, 2L, 3L))
  result <- class_angular_distortion_index(
    x,
    y,
    c(1, 2, 2),
    triplets = plan,
    ret_extra = TRUE
  )
  expect_equal(result$cosine_in, 1)
  expect_equal(result$cosine_out, 0)
  expect_equal(result$cadi, 1 / 4)
  expect_error(
    class_angular_distortion_index(x, y, c(1, 2, 1), triplets = plan),
    "share a class"
  )
  expect_error(class_angular_distortion_index(x, y, c(1, 2)), "number of rows")
})

test_that("integer coordinates are promoted before calculating displacements", {
  x <- rbind(c(-2000000000L, 0L), c(2000000000L, 0L), c(-2000000000L, 1L))
  y <- rbind(c(0, 0), c(1, 0), c(2, 0))
  expect_equal(
    class_angular_distortion_index(x, y, c(1, 2, 2), triplets = rbind(1:3)),
    1 / 4
  )
})
