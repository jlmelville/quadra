test_that("public thread controls require supported nonnegative integers", {
  # fmt: skip
  idx <- matrix(
    c(
      2, 3,
      1, 3,
      2, 1,
      3, 2
    ),
    nrow = 4,
    byrow = TRUE
  )
  graph <- list(idx = idx, dist = matrix(rep(c(1, 2), each = 4), nrow = 4))
  xin <- matrix(seq_len(12), nrow = 4)
  xout <- xin

  public_calls <- list(
    nbr_pres_knn = function(n_threads) {
      nbr_pres_knn(idx, idx, k = 1, n_threads = n_threads)
    },
    nn_preservation = function(n_threads) {
      nn_preservation(graph, graph, k = 1, n_threads = n_threads)
    },
    local_radius_correlation = function(n_threads) {
      local_radius_correlation(graph, graph, k = 1, n_threads = n_threads)
    },
    mutual_neighbor_correlation = function(n_threads) {
      mutual_neighbor_correlation(graph, graph, k = 1, n_threads = n_threads)
    },
    random_pair_distance_correlation = function(n_threads) {
      random_pair_distance_correlation(
        xin,
        xout,
        n_pairs = 5,
        n_threads = n_threads
      )
    },
    random_pair_distance_emd = function(n_threads) {
      random_pair_distance_emd(xin, xout, n_pairs = 5, n_threads = n_threads)
    },
    random_pair_distance_stress = function(n_threads) {
      random_pair_distance_stress(xin, xout, n_pairs = 5, n_threads = n_threads)
    },
    random_triplet_accuracy = function(n_threads) {
      random_triplet_accuracy(xin, xout, n_triplets = 2, n_threads = n_threads)
    }
  )
  invalid_counts <- list(
    -1,
    NA_real_,
    NaN,
    Inf,
    1.5,
    c(1, 2),
    "1",
    .Machine$integer.max + 1
  )

  for (call_name in names(public_calls)) {
    call_metric <- public_calls[[call_name]]
    for (invalid_count in invalid_counts) {
      expect_error(
        call_metric(invalid_count),
        "n_threads must be one finite, whole, nonnegative value",
        info = call_name
      )
    }
  }
})

test_that("zero and one thread retain serial execution", {
  xin <- matrix(seq_len(12), nrow = 4)

  set.seed(20260811)
  zero_threads <- random_pair_distances(xin, xin, n_pairs = 20, n_threads = 0)
  set.seed(20260811)
  one_thread <- random_pair_distances(xin, xin, n_pairs = 20, n_threads = 1)

  expect_equal(one_thread, zero_threads)
})
