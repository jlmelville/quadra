test_that("Gap Index measures relative areas on the output triangulation", {
  skip_if_not_installed("geometry")
  x <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0.4, 0.6))
  y <- x
  y[5, ] <- c(0.2, 0.5)
  details <- gap_index(x, y, ret_extra = TRUE)
  expect_equal(gap_index(x, y), details$gap_index)
  expect_gt(details$gap_index, 0)
  expect_equal(sum(details$relative_area_in), 1)
  expect_equal(sum(details$relative_area_out), 1)

  # Compare with a direct planar area calculation for the same vertices.
  relative <- function(z) {
    t <- details$triangles
    a <- z[t[, 2L], , drop = FALSE] - z[t[, 1L], , drop = FALSE]
    b <- z[t[, 3L], , drop = FALSE] - z[t[, 1L], , drop = FALSE]
    area <- abs(a[, 1L] * b[, 2L] - a[, 2L] * b[, 1L])
    area / sum(area)
  }
  expect_equal(details$relative_area_in, relative(x))
  expect_equal(details$relative_area_out, relative(y))
  weight <- pmax(relative(x), relative(y))
  expect_equal(details$deformation, (relative(y) - relative(x)) / weight)
  expect_equal(
    details$gap_index,
    weighted.mean(abs(details$deformation), weight)
  )
  expect_true(any(details$deformation > 0) && any(details$deformation < 0))
})

test_that("global rescaling and planar shears preserve relative areas", {
  skip_if_not_installed("geometry")
  x <- rbind(c(0, 0), c(3, 0), c(1, 1), c(0, 2), c(0.4, 0.6))
  shear <- x %*% rbind(c(2, 0), c(0.8, 1))
  expect_equal(gap_index(x, x), 0)
  expect_equal(gap_index(x, shear), 0, tolerance = 1e-12)
  expect_equal(gap_index(cbind(x, 3), 5 * x), 0, tolerance = 1e-12)
})

test_that("Gap Index handles degenerate reference triangles and validates embeddings", {
  skip_if_not_installed("geometry")
  y <- rbind(c(0, 0), c(1, 0), c(0, 1))
  expect_equal(gap_index(y, y), 0)
  expect_error(gap_index(cbind(0:2, 0), y), "positive total triangle area")
  expect_error(gap_index(y, cbind(0:2, 0)), "non-collinear")
  expect_error(gap_index(y, y[c(1, 2, 1), ]), "distinct")
  expect_error(gap_index(y, y[, 1, drop = FALSE]), "two columns")
  expect_error(gap_index(y, y[-3, ]), "same number")
})
