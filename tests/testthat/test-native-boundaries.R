# Direct calls to generated internal wrappers are intentional here. These tests
# cover only malformed inputs that could otherwise reach unchecked indexing,
# narrowing, allocation, or thread setup. Public semantics are tested through
# exported R functions instead.

test_that("neighbor-overlap native boundary protects unchecked access", {
  idx <- matrix(c(2, 1), ncol = 1)

  expect_error(
    neighbor_overlap_counts(idx, idx[-1, , drop = FALSE], 1L),
    "same number of rows"
  )
  expect_error(
    neighbor_overlap_counts(idx, idx, integer()),
    "k must be nonempty"
  )
  expect_error(
    neighbor_overlap_counts(idx, idx, -1L),
    "k must contain nonnegative integers"
  )
  expect_error(
    neighbor_overlap_counts(idx, idx, 2L),
    "larger than the number of columns"
  )
  bad_idx <- idx
  bad_idx[1, 1] <- 0
  expect_error(
    neighbor_overlap_counts(bad_idx, idx, 1L),
    "finite integer indices"
  )
  expect_error(
    neighbor_overlap_counts(idx, idx, 1L, -1),
    "n_threads must be one finite, whole, nonnegative value"
  )
  expect_error(
    neighbor_overlap_counts(idx, idx, 1L, .Machine$integer.max + 1),
    "supported integer range"
  )
})

test_that("random-distance native boundary protects native setup", {
  xin <- matrix(c(0, 1), nrow = 1)

  expect_error(
    random_distances(
      matrix(numeric(), 0, 0),
      matrix(numeric(), 0, 0),
      n_pairs = 1
    ),
    "at least one feature and two observations"
  )
  expect_error(
    random_distances(xin, matrix(c(0, 1, 2), nrow = 1), n_pairs = 1),
    "same number of observations"
  )
  expect_error(
    random_distances(xin, xin, n_pairs = -1),
    "n_pairs must be finite, whole, and nonnegative"
  )
  expect_error(
    random_distances(xin, xin, n_pairs = 1, n_threads = -1),
    "n_threads must be one finite, whole, nonnegative value"
  )
  expect_error(
    random_distances(
      xin,
      xin,
      n_pairs = 1,
      n_threads = .Machine$integer.max + 1
    ),
    "supported integer range"
  )
})

test_that("exact-rank native boundaries protect matrix indexing", {
  empty <- matrix(numeric(), 0, 0)
  square <- diag(3)
  nonsquare <- matrix(1, nrow = 2, ncol = 3)

  expect_error(rnx_auc_direct(empty, empty), "at least three observations")
  expect_error(rnx_auc_direct(square, nonsquare), "same dimensions")
  expect_error(trustworthiness_exact(square, nonsquare, 1L), "same dimensions")
  expect_error(continuity_exact(nonsquare, nonsquare, 1L), "square distance")
})

test_that("explicit-triplet native boundary protects native indexing", {
  xin <- matrix(c(0, 1, 2), nrow = 1)
  # fmt: skip
  valid <- matrix(
    c(
      1, 0, 0,
      2, 2, 1
    ),
    nrow = 2,
    byrow = TRUE
  )

  expect_error(
    triplet_sample(valid, xin, matrix(c(0, 1, 2, 3), nrow = 1)),
    "same number of observations"
  )
  expect_error(
    triplet_sample(matrix(0, nrow = 2, ncol = 2), xin, xin),
    "one column per observation"
  )
  expect_error(
    triplet_sample(matrix(0, nrow = 3, ncol = 3), xin, xin),
    "even number of rows"
  )

  unsafe_triplets <- list(
    negative = replace(valid, 1, -1),
    out_of_range = replace(valid, 1, 3),
    nonfinite = replace(valid, 1, NA_real_),
    noninteger = replace(valid, 1, 1.5)
  )
  for (case_name in names(unsafe_triplets)) {
    expect_error(
      triplet_sample(unsafe_triplets[[case_name]], xin, xin),
      "finite integer indices in range",
      info = case_name
    )
  }

  expect_error(
    triplet_sample(valid, xin, xin, n_threads = -1),
    "n_threads must be one finite, whole, nonnegative value"
  )
  expect_error(
    triplet_sample(
      valid,
      xin,
      xin,
      n_threads = .Machine$integer.max + 1
    ),
    "supported integer range"
  )
})

test_that("random-triplet native boundary protects native setup", {
  xin <- matrix(c(0, 1, 2), nrow = 1)

  expect_error(
    random_triplet_sample(matrix(numeric(), 0, 0), matrix(numeric(), 0, 0)),
    "at least one feature and three observations"
  )
  expect_error(
    random_triplet_sample(xin, matrix(c(0, 1, 2, 3), nrow = 1)),
    "same number of observations"
  )
  expect_error(
    random_triplet_sample(xin, xin, n_triplets = -1),
    "n_triplets must be finite, whole, and nonnegative"
  )
  expect_error(
    random_triplet_sample(xin, xin, n_threads = -1),
    "n_threads must be one finite, whole, nonnegative value"
  )
  expect_error(
    random_triplet_sample(
      xin,
      xin,
      n_threads = .Machine$integer.max + 1
    ),
    "supported integer range"
  )
})
