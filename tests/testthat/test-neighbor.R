test_that("distance-matrix neighborhood preservation excludes self", {
  din <- distance_matrix(matrix(c(0, 1, 3, 10), ncol = 1))
  dout <- distance_matrix(matrix(c(0, 10, 3, 1), ncol = 1))

  expect_equal(nbr_pres(din, din, k = 1), rep(1, 4))
  expect_equal(nbr_pres(din, dout, k = 1), rep(0, 4))
  expect_error(nbr_pres(din, dout, k = 4), "non-self observations")
})

test_that("nearest-neighbor matrix preservation uses the supplied neighbors", {
  # fmt: skip
  kin <- matrix(
    c(
      2, 3,
      1, 3,
      2, 1,
      3, 2
    ),
    nrow = 4,
    byrow = TRUE
  )
  # fmt: skip
  kout <- matrix(
    c(
      4, 3,
      3, 4,
      4, 1,
      1, 3
    ),
    nrow = 4,
    byrow = TRUE
  )

  expect_equal(nbr_pres_knn(kin, kin, k = 1), rep(1, 4))
  expect_equal(nbr_pres_knn(kin, kout, k = 1), rep(0, 4))
  expect_equal(nbr_pres_knn(kin, kout, k = 2), rep(0.5, 4))
})

test_that("co-ranking matrices exclude self-neighbors", {
  din <- distance_matrix(matrix(c(0, 1, 3, 10), ncol = 1))
  dout <- distance_matrix(matrix(c(0, 10, 3, 1), ncol = 1))

  expect_equal(coranking_matrix(din, din), diag(4, 3))
  # fmt: skip
  expected_crm <- matrix(
    c(
      0, 1, 3,
      1, 2, 1,
      3, 1, 0
    ),
    nrow = 3,
    byrow = TRUE
  )
  expect_equal(
    coranking_matrix(din, dout),
    expected_crm
  )
})

test_that("QNX and RNX use self-excluded co-ranking denominators", {
  din <- distance_matrix(matrix(c(0, 1, 3, 10), ncol = 1))
  dout <- distance_matrix(matrix(c(0, 10, 3, 1), ncol = 1))
  crm <- coranking_matrix(din, dout)

  expect_equal(qnx_crm(coranking_matrix(din, din), k = 1), 1)
  expect_equal(rnx_auc(din, din), 1)
  expect_equal(qnx_crm(crm, k = 1), 0)
  expect_equal(rnx_crm(crm, k = 1), -0.5)
  expect_equal(rnx_auc(din, dout), -0.5)
})

test_that("QNX and RNX match hand calculations on tiny co-ranking matrices", {
  # fmt: skip
  crm <- matrix(
    c(
      1, 1, 2,
      0, 2, 2,
      3, 1, 0
    ),
    nrow = 3,
    byrow = TRUE
  )

  expect_equal(qnx_crm(crm, k = 1), 1 / 4)
  expect_equal(qnx_crm(crm, k = 2), 1 / 2)
  expect_equal(rnx_crm(crm, k = 1), -1 / 8)
  expect_equal(rnx_crm(crm, k = 2), -1 / 2)
  expect_equal(rnx_auc_crm(crm), -1 / 4)
})
