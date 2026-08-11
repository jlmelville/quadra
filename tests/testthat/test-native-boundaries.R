# Direct calls to the generated internal wrappers are intentional here. Public
# wrappers mask these malformed inputs, while registered native routines must
# still reject them instead of performing unchecked indexing or unsigned casts.

test_that("neighbor-overlap native boundary validates shapes and threads", {
  idx <- matrix(c(2, 1), ncol = 1)

  expect_error(
    neighbor_overlap_counts(
      matrix(numeric(), 0, 0),
      matrix(numeric(), 0, 0),
      1L
    ),
    "at least two rows and one column"
  )
  expect_error(
    neighbor_overlap_counts(idx, idx[-1, , drop = FALSE], 1L),
    "same number of rows"
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

test_that("random-distance native boundary validates shapes and threads", {
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
    random_distances(xin, xin, n_pairs = 0),
    "n_pairs must be a positive integer"
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

test_that("exact-rank native boundaries validate distance-matrix shapes", {
  empty <- matrix(numeric(), 0, 0)
  square <- diag(3)
  nonsquare <- matrix(1, nrow = 2, ncol = 3)

  expect_error(rnx_auc_direct(empty, empty), "at least three observations")
  expect_error(rnx_auc_direct(square, nonsquare), "same dimensions")
  expect_error(trustworthiness_exact(square, nonsquare, 1L), "same dimensions")
  expect_error(continuity_exact(nonsquare, nonsquare, 1L), "square distance")
})

test_that("explicit-triplet native boundary validates shapes, indices, and threads", {
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
    triplet_sample(matrix(integer(), 0, 3), xin, xin),
    "nonempty"
  )
  expect_error(
    triplet_sample(valid, xin, matrix(c(0, 1, 2, 3), nrow = 1)),
    "same number of observations"
  )

  unsafe_triplets <- list(
    negative = replace(valid, 1, -1),
    out_of_range = replace(valid, 1, 3),
    anchor_endpoint = replace(valid, 1, 0),
    duplicate_endpoints = replace(valid, 2, 1)
  )
  for (case_name in names(unsafe_triplets)) {
    expect_error(
      triplet_sample(unsafe_triplets[[case_name]], xin, xin),
      "distinct.*in range",
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

test_that("random-triplet native boundary validates shapes and threads", {
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
    random_triplet_sample(xin, xin, n_triplets = 0),
    "n_triplets must be a positive integer"
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
