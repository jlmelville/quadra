test_that("nonfinite triplet distances fail before becoming ties or disagreements", {
  small <- matrix(c(0, 1, 2), ncol = 1)
  large <- small * 1e200
  explicit <- matrix(c(1L, 2L, 3L), nrow = 1)
  legacy <- rbind(c(1, 0, 0), c(2, 2, 1))
  for (threads in c(0, 2)) {
    for (input_large in c(FALSE, TRUE)) {
      xin <- if (input_large) large else small
      for (sample in list(1, legacy)) {
        expect_error(
          random_triplet_accuracy(
            xin,
            large,
            n_triplets = sample,
            n_threads = threads
          ),
          "nonfinite distance"
        )
      }
      expect_error(
        random_triplet_accuracy(
          xin,
          large,
          triplets = explicit,
          n_threads = threads
        ),
        "nonfinite distance"
      )
    }
  }
})

test_that("nonfinite pair distances fail before aggregation", {
  small <- matrix(c(0, 1, 2), ncol = 1)
  large <- small * 1e200
  pairs <- rbind(c(1L, 2L), c(1L, 3L))
  metrics <- list(
    random_pair_distance_correlation,
    random_pair_distance_emd,
    random_pair_distance_stress
  )
  for (fun in metrics) {
    for (threads in c(0, 2)) {
      for (sample in list(NULL, pairs)) {
        expect_error(
          fun(large, small, pairs = sample, n_pairs = 4, n_threads = threads),
          "nonfinite distance"
        )
        expect_error(
          fun(small, large, pairs = sample, n_pairs = 4, n_threads = threads),
          "nonfinite distance"
        )
      }
    }
  }
})

test_that("NaN angular distances also cause numerical-failure errors", {
  x <- rbind(c(1, 2), c(2, 1), c(3, 4))
  for (metric in c("cosine", "correlation")) {
    for (threads in c(0, 2)) {
      expect_error(
        random_pair_distance_stress(
          x,
          x * 1e200,
          pairs = rbind(c(1L, 2L), c(1L, 3L)),
          metric_in = metric,
          metric_out = metric,
          n_threads = threads
        ),
        "nonfinite distance"
      )
      expect_error(
        random_triplet_accuracy(
          x,
          x * 1e200,
          triplets = matrix(c(1L, 2L, 3L), nrow = 1),
          metric_in = metric,
          metric_out = metric,
          n_threads = threads
        ),
        "nonfinite distance"
      )
    }
  }
})

test_that("stress scales residuals before squaring", {
  pairs <- rbind(c(1L, 2L), c(1L, 3L))
  # Manhattan distances stay representable at both scales. Squared residuals do not.
  for (scale in c(1e200, 1e-200)) {
    x <- matrix(c(0, 1, 2) * scale, ncol = 1)
    result <- random_pair_distance_stress(
      x,
      x / 2,
      pairs = pairs,
      metric_in = "manhattan",
      metric_out = "manhattan",
      range_scale = FALSE,
      ret_extra = TRUE
    )
    expect_true(is.finite(result$stress))
    expect_equal(result$stress / scale, sqrt(5 / 8))
    expect_equal(result$distance_in / scale, c(1, 2))
    expect_equal(result$distance_out / scale, c(0.5, 1))
    expect_identical(
      random_pair_distance_stress(
        x,
        x,
        pairs = pairs,
        metric_in = "manhattan",
        metric_out = "manhattan",
        range_scale = FALSE
      ),
      0
    )
  }
})
