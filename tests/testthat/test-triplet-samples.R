canonical_triplets <- function(triplets) {
  matrix(
    as.integer(triplets),
    nrow = nrow(triplets),
    ncol = 3,
    dimnames = list(NULL, c("anchor", "endpoint1", "endpoint2"))
  )
}

triplet_agreement <- function(xin, xout, triplets) {
  squared_distances <- function(x, endpoint) {
    delta <-
      x[triplets[, "anchor"], , drop = FALSE] -
      x[triplets[, endpoint], , drop = FALSE]
    rowSums(delta^2)
  }

  input_order <- sign(
    squared_distances(xin, "endpoint1") -
      squared_distances(xin, "endpoint2")
  )
  output_order <- sign(
    squared_distances(xout, "endpoint1") -
      squared_distances(xout, "endpoint2")
  )
  agreement <- input_order == output_order
  agreement[input_order == 0] <- NA
  agreement
}

triplet_accuracy <- function(agreement) {
  if (all(is.na(agreement))) {
    return(NA_real_)
  }
  mean(agreement, na.rm = TRUE)
}

test_that("sample_triplets returns reproducible anchor-major samples", {
  set.seed(20260823)
  triplets <- sample_triplets(n_obs = 7, n_triplets = 3)

  expect_type(triplets, "integer")
  expect_identical(dim(triplets), c(21L, 3L))
  expect_identical(colnames(triplets), c("anchor", "endpoint1", "endpoint2"))
  expect_identical(triplets[, "anchor"], rep(1:7, each = 3))
  expect_true(all(triplets >= 1L & triplets <= 7L))
  expect_true(all(triplets[, "anchor"] != triplets[, "endpoint1"]))
  expect_true(all(triplets[, "anchor"] != triplets[, "endpoint2"]))
  expect_true(all(triplets[, "endpoint1"] != triplets[, "endpoint2"]))

  set.seed(20260823)
  expect_identical(sample_triplets(7, 3), triplets)
})

test_that("sample_triplets validates its counts", {
  expect_error(sample_triplets(2, 5), "at least 3 observations")
  expect_error(sample_triplets(3.5, 5), "n_obs must be a positive integer")
  expect_error(sample_triplets(3, 0), "n_triplets must be a positive integer")
  expect_error(
    sample_triplets(.Machine$integer.max + 1, 5),
    "n_obs exceeds the supported integer range"
  )
  expect_error(
    sample_triplets(50000, 50000),
    "n_obs \\* n_triplets exceeds the supported integer range"
  )
})

test_that("explicit triplet details match a direct R calculation", {
  # fmt: skip
  xin <- matrix(
    c(
       0, 0,
       1, 0,
       0, 2,
       3, 0,
      -1, 0
    ),
    ncol = 2,
    byrow = TRUE
  )
  # fmt: skip
  xout <- matrix(
    c(
       0, 0,
       2, 0,
       0, 1,
       3, 0,
      -1, 0
    ),
    ncol = 2,
    byrow = TRUE
  )
  # Includes an input tie, an output tie, repeated and reversed endpoints,
  # unbalanced anchors, and a deliberately non-anchor-major row order.
  # fmt: skip
  supplied <- matrix(
    c(
      3, 1, 2,
      1, 2, 3,
      1, 2, 4,
      2, 1, 4,
      1, 2, 5,
      1, 3, 5,
      1, 2, 4,
      1, 4, 2
    ),
    ncol = 3,
    byrow = TRUE
  )
  triplets <- canonical_triplets(supplied)
  agreement <- triplet_agreement(xin, xout, triplets)

  detailed <- random_triplet_accuracy(
    xin,
    xout,
    n_triplets = "ignored",
    triplets = supplied,
    ret_extra = TRUE
  )
  expect_named(detailed, c("accuracy", "triplets", "agreement"))
  expect_equal(detailed$accuracy, triplet_accuracy(agreement))
  expect_identical(detailed$triplets, triplets)
  expect_identical(detailed$agreement, agreement)
  expect_identical(
    random_triplet_accuracy(
      t(xin),
      t(xout),
      is_transposed = TRUE,
      triplets = supplied,
      ret_extra = TRUE
    ),
    detailed
  )
  expect_equal(
    random_triplet_accuracy(xin, xout, triplets = supplied),
    detailed$accuracy
  )

  reversed <- random_triplet_accuracy(
    xin,
    xout,
    triplets = supplied[rev(seq_len(nrow(supplied))), , drop = FALSE],
    ret_extra = TRUE
  )
  expect_equal(reversed$accuracy, detailed$accuracy)
  expect_identical(reversed$agreement, rev(detailed$agreement))
})

test_that("legacy triplet matrices retain their layout and scalar result", {
  # fmt: skip
  legacy <- matrix(
    c(
      1, 0, 0, 0, 0,
      2, 2, 1, 1, 1,
      3, 3, 3, 2, 2,
      4, 4, 4, 4, 3
    ),
    nrow = 4,
    byrow = TRUE
  )
  # fmt: skip
  expected_triplets <- canonical_triplets(
    matrix(
      c(
        1, 2, 3,
        1, 4, 5,
        2, 1, 3,
        2, 4, 5,
        3, 1, 2,
        3, 4, 5,
        4, 1, 2,
        4, 3, 5,
        5, 1, 2,
        5, 3, 4
      ),
      ncol = 3,
      byrow = TRUE
    )
  )

  compact <- random_triplet_accuracy(m[1:5, ], n[1:5, ], legacy)
  detailed <- random_triplet_accuracy(
    m[1:5, ],
    n[1:5, ],
    legacy,
    ret_extra = TRUE
  )
  expect_identical(detailed$triplets, expected_triplets)
  expect_identical(
    detailed$agreement,
    triplet_agreement(m[1:5, ], n[1:5, ], expected_triplets)
  )
  expect_identical(detailed$accuracy, compact)

  expect_error(
    random_triplet_accuracy(
      m[1:5, ],
      n[1:5, ],
      n_triplets = legacy,
      triplets = expected_triplets
    ),
    "only one of n_triplets and triplets"
  )
})

test_that("implicit triplet details retain the evaluated sample", {
  set.seed(4815)
  compact <- random_triplet_accuracy(m, n, n_triplets = 4, n_threads = 2)
  seed_after_compact <- .Random.seed
  set.seed(4815)
  detailed <- random_triplet_accuracy(
    m,
    n,
    n_triplets = 4,
    n_threads = 2,
    ret_extra = TRUE
  )

  expect_identical(.Random.seed, seed_after_compact)
  expect_named(detailed, c("accuracy", "triplets", "agreement"))
  expect_identical(detailed$accuracy, compact)
  expect_identical(dim(detailed$triplets), c(40L, 3L))
  expect_identical(detailed$triplets[, "anchor"], rep(1:10, each = 4))
  expect_true(all(detailed$triplets >= 1L & detailed$triplets <= 10L))
  expect_true(
    all(detailed$triplets[, "anchor"] != detailed$triplets[, "endpoint1"])
  )
  expect_true(
    all(detailed$triplets[, "anchor"] != detailed$triplets[, "endpoint2"])
  )
  expect_true(
    all(detailed$triplets[, "endpoint1"] != detailed$triplets[, "endpoint2"])
  )
  expect_identical(
    detailed$agreement,
    triplet_agreement(m, n, detailed$triplets)
  )
})

test_that("supplied triplets bypass RNG and threads", {
  # fmt: skip
  triplets <- matrix(
    c(
      1, 2, 3,
      4, 1, 5,
      2, 3, 1,
      1, 2, 3
    ),
    ncol = 3,
    byrow = TRUE
  )

  set.seed(13579)
  seed_before <- .Random.seed
  serial <- random_triplet_accuracy(
    m[1:5, ],
    n[1:5, ],
    triplets = triplets,
    n_threads = 0,
    ret_extra = TRUE
  )
  expect_identical(.Random.seed, seed_before)
  expect_identical(
    random_triplet_accuracy(
      m[1:5, ],
      n[1:5, ],
      triplets = triplets,
      n_threads = 1,
      ret_extra = TRUE
    ),
    serial
  )
  expect_identical(
    random_triplet_accuracy(
      m[1:5, ],
      n[1:5, ],
      triplets = triplets,
      n_threads = 2,
      ret_extra = TRUE
    ),
    serial
  )
  expect_identical(.Random.seed, seed_before)

  # fmt: skip
  legacy <- matrix(
    c(
      1, 0, 0, 0, 0,
      2, 2, 1, 1, 1
    ),
    nrow = 2,
    byrow = TRUE
  )
  legacy_serial <- random_triplet_accuracy(
    m[1:5, ],
    n[1:5, ],
    n_triplets = legacy,
    n_threads = 0
  )
  expect_identical(
    random_triplet_accuracy(
      m[1:5, ],
      n[1:5, ],
      n_triplets = legacy,
      n_threads = 2
    ),
    legacy_serial
  )
  expect_identical(.Random.seed, seed_before)
})

test_that("supplied triplets do not initialize the R RNG", {
  # fmt: skip
  triplets <- matrix(
    c(
      1, 2, 3,
      4, 1, 5
    ),
    ncol = 3,
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
    result <- random_triplet_accuracy(
      m[1:5, ],
      n[1:5, ],
      triplets = triplets,
      n_threads = 2
    )
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    # fmt: skip
    legacy <- matrix(
      c(
        1, 0, 0, 0, 0,
        2, 2, 1, 1, 1
      ),
      nrow = 2,
      byrow = TRUE
    )
    random_triplet_accuracy(
      m[1:5, ],
      n[1:5, ],
      n_triplets = legacy,
      n_threads = 2
    )
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    result
  }

  expect_type(call_without_seed(), "double")
})

test_that("explicit triplet validation covers canonical samples", {
  # fmt: skip
  valid <- matrix(
    c(
      1, 2, 3,
      4, 1, 5
    ),
    ncol = 3,
    byrow = TRUE
  )
  invalid_triplets <- list(
    vector = list(value = c(1, 2, 3), error = "nonempty numeric matrix"),
    data_frame = list(
      value = data.frame(anchor = 1, endpoint1 = 2, endpoint2 = 3),
      error = "nonempty numeric matrix"
    ),
    character = list(
      value = matrix(c("1", "2", "3"), nrow = 1),
      error = "nonempty numeric matrix"
    ),
    complex = list(
      value = matrix(c(1 + 0i, 2 + 0i, 3 + 0i), nrow = 1),
      error = "nonempty numeric matrix"
    ),
    empty = list(
      value = matrix(integer(), nrow = 0, ncol = 3),
      error = "nonempty numeric matrix"
    ),
    wrong_columns = list(value = matrix(1:2, nrow = 1), error = "3 columns"),
    missing = list(
      value = matrix(c(1, 2, NA), nrow = 1),
      error = "finite whole"
    ),
    infinite = list(
      value = matrix(c(1, 2, Inf), nrow = 1),
      error = "finite whole"
    ),
    fractional = list(
      value = matrix(c(1, 2, 3.5), nrow = 1),
      error = "finite whole"
    ),
    too_small = list(value = matrix(c(0, 2, 3), nrow = 1), error = "between 1"),
    too_large = list(value = matrix(c(1, 2, 6), nrow = 1), error = "between 1"),
    anchor_endpoint = list(
      value = matrix(c(1, 1, 3), nrow = 1),
      error = "three distinct"
    ),
    duplicate_endpoints = list(
      value = matrix(c(1, 2, 2), nrow = 1),
      error = "three distinct"
    )
  )

  for (case_name in names(invalid_triplets)) {
    case <- invalid_triplets[[case_name]]
    expect_error(
      random_triplet_accuracy(m[1:5, ], n[1:5, ], triplets = case$value),
      case$error,
      info = case_name
    )
  }
  expect_error(
    random_triplet_accuracy(
      m[1:5, ],
      n[1:5, ],
      triplets = valid,
      ret_extra = NA
    ),
    "ret_extra must be TRUE or FALSE"
  )
})

test_that("all undefined triplets retain row-level outcomes", {
  # fmt: skip
  triplets <- matrix(
    c(
      1, 2, 3,
      2, 3, 4
    ),
    ncol = 3,
    byrow = TRUE
  )
  detailed <- random_triplet_accuracy(
    matrix(0, nrow = 4, ncol = 2),
    m[1:4, ],
    triplets = triplets,
    ret_extra = TRUE
  )

  expect_true(is.na(detailed$accuracy))
  expect_identical(detailed$agreement, c(NA, NA))
})
