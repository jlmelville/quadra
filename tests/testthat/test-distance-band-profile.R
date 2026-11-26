test_that("distance-band profiles use input radii and preserve supplied pair rows", {
  x <- cbind(c(0, 1, 2, 4))
  y <- cbind(c(0, 2, 1, 5))
  pairs <- t(combn(nrow(x), 2L))
  pairs <- rbind(pairs, pairs[1L, 2:1])

  set.seed(17)
  seed <- .Random.seed
  expect_warning(
    result <- random_pair_distance_band_profile(
      x,
      y,
      bands = c(0.4, 1),
      pairs = pairs,
      n_pairs = NA,
      ret_extra = TRUE
    ),
    "undefined"
  )
  expect_identical(.Random.seed, seed)
  expect_identical(unname(result$pairs), pairs)
  expect_equal(result$distance_in, c(1, 2, 4, 1, 3, 2, 1))
  expect_equal(result$distance_out, c(2, 1, 5, 1, 3, 4, 2))
  profile <- result$profile
  expect_equal(profile$radius, c(1, 4))
  expect_identical(profile$n_pairs, c(3L, 7L))
  expect_equal(profile$pair_fraction, c(3 / 7, 1))
  expect_true(is.na(profile$spearman[1L]))
  expect_equal(
    profile$spearman[2L],
    random_pair_distance_correlation(x, y, pairs = pairs, method = "spearman")
  )
  expect_equal(
    suppressWarnings(random_pair_distance_band_profile(
      t(x),
      t(y),
      bands = c(1, 4),
      band_type = "radius",
      pairs = pairs,
      is_transposed = TRUE,
      n_threads = 2
    ))[, -(1:2)],
    profile[, -(1:2)]
  )
})

test_that("fitted-scale stress uses a multiplier with no offset", {
  x <- cbind(0:3)
  y <- cbind(c(0, 2, 3, 4))
  pairs <- cbind(1L, 2:4)

  full <- random_pair_distance_band_profile(x, y, bands = 1, pairs = pairs)
  expect_equal(full$normalized_stress, sqrt(3 / 203))
  expect_equal(full$spearman, 1)
  expect_equal(
    random_pair_distance_stress(
      x,
      y,
      pairs = pairs,
      metric_in = "euclidean",
      metric_out = "euclidean"
    ),
    0,
    tolerance = 1e-15
  )
  scaled <- random_pair_distance_band_profile(
    x,
    10 * x,
    bands = 1,
    pairs = pairs
  )
  expect_equal(scaled$normalized_stress, 0, tolerance = 1e-15)
})

test_that("undefined bands report their counts and consolidate warnings", {
  x <- cbind(0:3)
  y <- x * 0
  warnings <- character()
  profile <- withCallingHandlers(
    random_pair_distance_band_profile(
      x,
      y,
      bands = c(0, 1, 3),
      band_type = "radius",
      pairs = cbind(1L, 2:4)
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings, 1L)
  expect_match(warnings, "undefined")
  expect_identical(profile$n_pairs, c(0L, 1L, 3L))
  expect_true(all(is.na(profile$spearman)))
  expect_true(all(is.na(profile$normalized_stress)))

  expect_error(
    random_pair_distance_band_profile(x, x, bands = c(1, 1)),
    "increasing"
  )
  expect_error(random_pair_distance_band_profile(x, x, bands = 0), "proportion")
  expect_error(
    random_pair_distance_band_profile(x, x, bands = -1, band_type = "radius"),
    "nonnegative"
  )
})
