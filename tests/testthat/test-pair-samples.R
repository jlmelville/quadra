pair_sqeuclidean <- function(x, pairs) {
  delta <-
    x[pairs[, "endpoint1"], , drop = FALSE] -
    x[pairs[, "endpoint2"], , drop = FALSE]
  rowSums(delta^2)
}

scale_to_unit_interval <- function(x) {
  x <- x - min(x)
  xmax <- max(x)
  if (xmax == 0) {
    return(rep(0, length(x)))
  }
  x / xmax
}

canonical_pairs <- function(pairs) {
  pairs <- matrix(
    as.integer(pairs),
    nrow = nrow(pairs),
    ncol = 2,
    dimnames = list(NULL, c("endpoint1", "endpoint2"))
  )
  pairs
}

test_that("sample_pairs returns reproducible one-based row pairs", {
  set.seed(20260823)
  pairs <- sample_pairs(n_obs = 7, n_pairs = 30)

  expect_type(pairs, "integer")
  expect_identical(dim(pairs), c(30L, 2L))
  expect_identical(colnames(pairs), c("endpoint1", "endpoint2"))
  expect_true(all(pairs >= 1L & pairs <= 7L))
  expect_true(all(pairs[, "endpoint1"] != pairs[, "endpoint2"]))

  set.seed(20260823)
  expect_identical(sample_pairs(7, 30), pairs)
})

test_that("sample_pairs validates its counts", {
  expect_error(sample_pairs(1, 10), "at least 2 observations")
  expect_error(sample_pairs(2.5, 10), "n_obs must be a positive integer")
  expect_error(sample_pairs(3, 0), "n_pairs must be a positive integer")
  expect_error(
    sample_pairs(.Machine$integer.max + 1, 10),
    "n_obs exceeds the supported integer range"
  )
})

test_that("explicit pair details match a direct R calculation", {
  # fmt: skip
  xin <- matrix(
    c(
      0, 0,
      1, 0,
      0, 2,
      3, 1
    ),
    nrow = 4,
    byrow = TRUE
  )
  # fmt: skip
  xout <- matrix(
    c(
      0, 0,
      2, 0,
      0, 1,
      3, 2
    ),
    nrow = 4,
    byrow = TRUE
  )
  # Repeated and reversed pairs are valid, and input row order is retained.
  # fmt: skip
  supplied <- matrix(
    c(
      4, 1,
      1, 2,
      3, 4,
      2, 3,
      4, 1,
      1, 4
    ),
    ncol = 2,
    byrow = TRUE
  )
  pairs <- canonical_pairs(supplied)
  distance_in <- pair_sqeuclidean(xin, pairs)
  distance_out <- pair_sqeuclidean(xout, pairs)

  pearson <- random_pair_distance_correlation(
    xin,
    xout,
    n_pairs = "ignored",
    pairs = supplied,
    ret_extra = TRUE
  )
  expect_named(
    pearson,
    c("correlation", "pairs", "distance_in", "distance_out")
  )
  expect_equal(pearson$correlation, stats::cor(distance_in, distance_out))
  expect_identical(pearson$pairs, pairs)
  expect_equal(pearson$distance_in, distance_in)
  expect_equal(pearson$distance_out, distance_out)
  expect_identical(
    random_pair_distance_correlation(
      t(xin),
      t(xout),
      is_transposed = TRUE,
      pairs = supplied,
      ret_extra = TRUE
    ),
    pearson
  )

  expect_equal(
    random_pair_distance_correlation(
      xin,
      xout,
      pairs = supplied,
      method = "spearman"
    ),
    stats::cor(distance_in, distance_out, method = "spearman")
  )

  scaled_in <- scale_to_unit_interval(distance_in)
  scaled_out <- scale_to_unit_interval(distance_out)
  expected_scaled_emd <- mean(abs(sort(scaled_in) - sort(scaled_out)))
  expected_raw_emd <- mean(abs(sort(distance_in) - sort(distance_out)))

  emd_scaled <- random_pair_distance_emd(
    xin,
    xout,
    pairs = supplied,
    ret_extra = TRUE
  )
  expect_named(emd_scaled, c("emd", "pairs", "distance_in", "distance_out"))
  expect_equal(emd_scaled$emd, expected_scaled_emd)
  expect_identical(emd_scaled$pairs, pairs)
  expect_equal(emd_scaled$distance_in, distance_in)
  expect_equal(emd_scaled$distance_out, distance_out)
  expect_equal(
    random_pair_distance_emd(
      xin,
      xout,
      pairs = supplied,
      range_scale = FALSE
    ),
    expected_raw_emd
  )

  expected_scaled_stress <- sqrt(mean((scaled_in - scaled_out)^2))
  expected_raw_stress <- sqrt(mean((distance_in - distance_out)^2))
  stress_scaled <- random_pair_distance_stress(
    xin,
    xout,
    pairs = supplied,
    ret_extra = TRUE
  )
  expect_named(
    stress_scaled,
    c("stress", "pairs", "distance_in", "distance_out")
  )
  expect_equal(stress_scaled$stress, expected_scaled_stress)
  expect_identical(stress_scaled$pairs, pairs)
  expect_equal(stress_scaled$distance_in, distance_in)
  expect_equal(stress_scaled$distance_out, distance_out)
  expect_equal(
    random_pair_distance_stress(
      xin,
      xout,
      pairs = supplied,
      range_scale = FALSE
    ),
    expected_raw_stress
  )

  reversed <- random_pair_distance_stress(
    xin,
    xout,
    pairs = supplied[rev(seq_len(nrow(supplied))), , drop = FALSE],
    ret_extra = TRUE
  )
  expect_equal(reversed$stress, stress_scaled$stress)
  expect_equal(reversed$distance_in, rev(stress_scaled$distance_in))
  expect_equal(reversed$distance_out, rev(stress_scaled$distance_out))
})

test_that("implicit pair details retain the evaluated sample", {
  metrics <- list(
    correlation = random_pair_distance_correlation,
    emd = random_pair_distance_emd,
    stress = random_pair_distance_stress
  )

  for (metric_name in names(metrics)) {
    metric <- metrics[[metric_name]]
    set.seed(4815)
    compact <- metric(m, n, n_pairs = 40, n_threads = 2)
    set.seed(4815)
    detailed <- metric(m, n, n_pairs = 40, n_threads = 2, ret_extra = TRUE)

    expect_equal(detailed[[metric_name]], compact, info = metric_name)
    expect_identical(dim(detailed$pairs), c(40L, 2L), info = metric_name)
    expect_identical(
      colnames(detailed$pairs),
      c("endpoint1", "endpoint2"),
      info = metric_name
    )
    expect_equal(
      detailed$distance_in,
      pair_sqeuclidean(m, detailed$pairs),
      info = metric_name
    )
    expect_equal(
      detailed$distance_out,
      pair_sqeuclidean(n, detailed$pairs),
      info = metric_name
    )
  }
})

test_that("supplied pairs bypass RNG and threads", {
  # fmt: skip
  pairs <- matrix(
    c(
      1, 2,
      3, 7,
      8, 4,
      5, 10
    ),
    ncol = 2,
    byrow = TRUE
  )
  metrics <- list(
    correlation = random_pair_distance_correlation,
    emd = random_pair_distance_emd,
    stress = random_pair_distance_stress
  )

  for (metric_name in names(metrics)) {
    metric <- metrics[[metric_name]]
    set.seed(13579)
    seed_before <- .Random.seed
    serial <- metric(m, n, pairs = pairs, n_threads = 0, ret_extra = TRUE)
    expect_identical(.Random.seed, seed_before, info = metric_name)

    expect_identical(
      metric(m, n, pairs = pairs, n_threads = 1, ret_extra = TRUE),
      serial,
      info = metric_name
    )
    expect_identical(
      metric(m, n, pairs = pairs, n_threads = 2, ret_extra = TRUE),
      serial,
      info = metric_name
    )
    expect_identical(.Random.seed, seed_before, info = metric_name)
  }
})

test_that("supplied pairs do not initialize the R RNG", {
  # fmt: skip
  pairs <- matrix(
    c(
      1, 2,
      3, 7,
      8, 4
    ),
    ncol = 2,
    byrow = TRUE
  )
  call_without_seed <- function() {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = ".Random.seed", envir = .GlobalEnv)
      }
    })

    if (had_seed) {
      rm(list = ".Random.seed", envir = .GlobalEnv)
    }
    result <- random_pair_distance_stress(
      m,
      n,
      pairs = pairs,
      n_threads = 2
    )
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    result
  }

  expect_type(call_without_seed(), "double")
})

test_that("explicit pair validation is shared by pair metrics", {
  # fmt: skip
  valid <- matrix(
    c(
      1, 2,
      3, 4
    ),
    ncol = 2,
    byrow = TRUE
  )
  invalid_pairs <- list(
    vector = list(value = c(1, 2), error = "nonempty numeric matrix"),
    data_frame = list(
      value = data.frame(endpoint1 = 1, endpoint2 = 2),
      error = "nonempty numeric matrix"
    ),
    character = list(
      value = matrix(c("1", "2"), nrow = 1),
      error = "nonempty numeric matrix"
    ),
    complex = list(
      value = matrix(c(1 + 0i, 2 + 0i), nrow = 1),
      error = "nonempty numeric matrix"
    ),
    empty = list(
      value = matrix(integer(), nrow = 0, ncol = 2),
      error = "nonempty numeric matrix"
    ),
    wrong_columns = list(
      value = matrix(1:3, nrow = 1),
      error = "2 columns"
    ),
    missing = list(value = matrix(c(1, NA), nrow = 1), error = "finite whole"),
    infinite = list(
      value = matrix(c(1, Inf), nrow = 1),
      error = "finite whole"
    ),
    fractional = list(
      value = matrix(c(1, 2.5), nrow = 1),
      error = "finite whole"
    ),
    too_small = list(value = matrix(c(0, 2), nrow = 1), error = "between 1"),
    too_large = list(value = matrix(c(1, 11), nrow = 1), error = "between 1"),
    self_pair = list(
      value = matrix(c(2, 2), nrow = 1),
      error = "distinct endpoints"
    )
  )
  metrics <- list(
    correlation = random_pair_distance_correlation,
    emd = random_pair_distance_emd,
    stress = random_pair_distance_stress
  )

  for (metric_name in names(metrics)) {
    metric <- metrics[[metric_name]]
    for (case_name in names(invalid_pairs)) {
      case <- invalid_pairs[[case_name]]
      expect_error(
        metric(m, n, pairs = case$value),
        case$error,
        info = paste(metric_name, case_name)
      )
    }
    expect_error(
      metric(m, n, pairs = valid, ret_extra = NA),
      "ret_extra must be TRUE or FALSE",
      info = metric_name
    )
  }
})

test_that("constant supplied distances retain metric conventions", {
  x <- matrix(1, nrow = 4, ncol = 2)
  # fmt: skip
  pairs <- matrix(
    c(
      1, 2,
      2, 3,
      3, 4
    ),
    ncol = 2,
    byrow = TRUE
  )

  expect_warning(
    correlation <- random_pair_distance_correlation(x, x, pairs = pairs),
    "standard deviation is zero"
  )
  expect_true(is.na(correlation))
  expect_equal(random_pair_distance_emd(x, x, pairs = pairs), 0)
  expect_equal(random_pair_distance_stress(x, x, pairs = pairs), 0)
})
