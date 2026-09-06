test_that("mutual counts match an independent R scan for every observation", {
  reference_counts <- function(idx, k) {
    prefix <- idx[, seq_len(k), drop = FALSE]
    vapply(
      seq_len(nrow(idx)),
      function(i) {
        sum(vapply(prefix[i, ], function(j) i %in% prefix[j, ], logical(1)))
      },
      integer(1)
    )
  }
  set.seed(42)
  n <- 40L
  random <- t(vapply(
    seq_len(n),
    function(i) {
      sample(setdiff(seq_len(n), i), 8L)
    },
    integer(8)
  ))
  # Nearby points on a line have many reciprocal neighbors, including late hits.
  distances <- as.matrix(stats::dist(seq_len(n)))
  diag(distances) <- Inf
  reciprocal <- t(apply(distances, 1, order))[, seq_len(8)]
  for (idx in list(random, reciprocal)) {
    k <- c(8L, 1L, 4L)
    for (storage in c("integer", "double")) {
      storage.mode(idx) <- storage
      result <- mutual_neighbor_correlation(
        list(idx = idx),
        list(idx = random),
        k = k,
        ret_extra = TRUE
      )
      expected_in <- vapply(
        k,
        function(ki) reference_counts(idx, ki),
        integer(n)
      )
      expected_out <- vapply(
        k,
        function(ki) reference_counts(random, ki),
        integer(n)
      )
      expect_identical(unname(result$mutual_neighbor_in), expected_in)
      expect_identical(unname(result$mutual_neighbor_out), expected_out)
      expect_identical(colnames(result$mutual_neighbor_in), paste0("mnc", k))
      expected_cor <- vapply(
        seq_along(k),
        function(j) {
          if (
            length(unique(expected_in[, j])) < 2L ||
              length(unique(expected_out[, j])) < 2L
          ) {
            NA_real_
          } else {
            stats::cor(expected_in[, j], expected_out[, j], method = "pearson")
          }
        },
        numeric(1)
      )
      expect_equal(unname(result$mnc), expected_cor)
    }
  }
})
