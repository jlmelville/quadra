test_that("random triplet accuracy", {
  # fmt: skip
  tm <- matrix(
    c(
      2, 5, 3, 4, 9, 4, 9, 3, 9, 6,
      4, 3, 4, 5, 8, 8, 8, 4, 4, 2,
      5, 0, 6, 2, 5, 3, 7, 1, 0, 7,
      6, 6, 8, 6, 6, 2, 5, 8, 3, 3
    ),
    byrow = TRUE,
    nrow = 4
  )

  # fmt: skip
  n <- matrix(
    c(
      -0.3673, -1.595,
       1.2870,  0.5113,
       0.7557, -0.2293,
       1.5190, -0.4499,
       0.4918,  1.0850,
       1.7160, -2.0440,
      -0.1086,  2.2280,
      -0.04816, 0.4036,
       0.2425, -0.2668,
      -1.7340,  0.3950
    ),
    byrow = TRUE,
    nrow = 10
  )

  # fmt: skip
  m <- matrix(
    c(
       1.5060,  0.05989,  0.3542,
       0.8419, -1.14600, -0.2262,
       1.0310,  1.45700,  0.5978,
       0.4838, -0.005043, 1.7710,
      -0.3932, -0.01158,  1.4700,
      -0.5604, -0.31320, -0.5020,
       2.6240,  1.50900, -0.4661,
       0.1407, -1.11100,  0.1571,
       0.3512, -1.10300, -0.5345,
      -0.4729,  1.45300,  1.0870
    ),
    byrow = TRUE,
    nrow = 10
  )

  expect_equal(random_triplet_accuracy(m, n, tm), 0.45)
  expect_equal(
    random_triplet_accuracy(t(m), t(n), tm, is_transposed = TRUE),
    0.45
  )
  expect_equal(random_triplet_accuracy(m, n, tm, n_threads = 1), 0.45)
  expect_equal(random_triplet_accuracy(m, n, tm, n_threads = 2), 0.45)
  expect_false("grain_size" %in% names(formals(random_triplet_accuracy)))
})

test_that("input-distance triplet ties are excluded", {
  # fmt: skip
  triplets <- matrix(
    c(
      1, 0, 0,
      2, 2, 1
    ),
    nrow = 2,
    byrow = TRUE
  )
  # fmt: skip
  xin <- matrix(
    c(
       0, 0,
       1, 0,
      -1, 0
    ),
    nrow = 3,
    byrow = TRUE
  )
  # fmt: skip
  xout <- matrix(
    c(
       0, 0,
       1, 0,
      -2, 0
    ),
    nrow = 3,
    byrow = TRUE
  )

  expect_equal(random_triplet_accuracy(xin, xout, triplets), 1)
  expect_true(is.na(random_triplet_accuracy(matrix(0, 3, 2), xout, triplets)))
})

test_that("random triplet inputs are validated", {
  expect_error(
    random_triplet_accuracy(m[1:2, ], n[1:2, ]),
    "at least 3 observations"
  )
  expect_error(
    random_triplet_accuracy(m, n, n_triplets = 0),
    "n_triplets must be a positive integer"
  )
  expect_error(
    random_triplet_accuracy(m, n, n_triplets = 1.5),
    "n_triplets must be a positive integer"
  )
  expect_error(
    random_triplet_accuracy(data.frame(label = letters[1:3]), n[1:3, ]),
    "at least one numeric column"
  )
})
