test_that("positive-integer controls reject unsupported ranges cleanly", {
  too_large <- .Machine$integer.max + 1
  public_calls <- list(
    scalar_k = function() nbr_pres(diag(3), diag(3), k = too_large),
    vector_k = function() nn_preservation(m, n, k = too_large),
    n_pairs = function() {
      random_pair_distance_correlation(m, n, n_pairs = too_large)
    },
    n_triplets = function() {
      random_triplet_accuracy(m, n, n_triplets = too_large)
    }
  )

  for (call_name in names(public_calls)) {
    call_control <- public_calls[[call_name]]
    expect_error(
      expect_warning(call_control(), NA),
      "supported integer range",
      info = call_name
    )
  }
})

test_that("positive-integer controls reject complex values cleanly", {
  public_calls <- list(
    scalar_k = function() nbr_pres(diag(3), diag(3), k = 1 + 0i),
    vector_k = function() nn_preservation(m, n, k = 1 + 0i),
    n_pairs = function() {
      random_pair_distance_correlation(m, n, n_pairs = 1 + 0i)
    },
    n_triplets = function() random_triplet_accuracy(m, n, n_triplets = 1 + 0i)
  )

  for (call_name in names(public_calls)) {
    call_control <- public_calls[[call_name]]
    expect_error(
      expect_warning(call_control(), NA),
      "positive integer",
      info = call_name
    )
  }
})
