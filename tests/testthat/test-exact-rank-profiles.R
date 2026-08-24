exact_rank_profile_fixture <- function() {
  # fmt: skip
  din <- matrix(
    c(
      0, 1, 1, 2, 2,
      1, 0, 3, 3, 2,
      1, 3, 0, 2, 2,
      2, 3, 2, 0, 1,
      2, 2, 2, 1, 0
    ),
    nrow = 5,
    byrow = TRUE
  )
  # fmt: skip
  dout <- matrix(
    c(
      0, 2, 1, 1, 2,
      2, 0, 2, 3, 3,
      1, 2, 0, 3, 3,
      1, 3, 3, 0, 2,
      2, 3, 3, 2, 0
    ),
    nrow = 5,
    byrow = TRUE
  )
  list(din = din, dout = dout)
}

rnx_profile_reference <- function(din, dout, k) {
  crm <- coranking_matrix(din, dout)
  stats::setNames(
    vapply(k, function(rank_k) rnx_crm(crm, rank_k), numeric(1)),
    paste0("rnx", k)
  )
}

test_that("rnx_curve matches the co-ranking reference", {
  fixture <- exact_rank_profile_fixture()
  full_k <- seq_len(nrow(fixture$din) - 2L)

  expect_equal(
    rnx_curve(fixture$din, fixture$dout),
    rnx_profile_reference(fixture$din, fixture$dout, full_k)
  )
  expect_equal(
    rnx_curve(fixture$din, fixture$dout, k = c(3, 1)),
    rnx_profile_reference(fixture$din, fixture$dout, c(3, 1))
  )
  expect_equal(
    rnx_curve(fixture$din, fixture$dout, k = 2),
    rnx_profile_reference(fixture$din, fixture$dout, 2)
  )
})

test_that("the full RNX curve reproduces RNX AUC", {
  fixture <- exact_rank_profile_fixture()
  curve <- rnx_curve(fixture$din, fixture$dout)
  k <- seq_along(curve)

  expect_equal(
    sum(unname(curve) / k) / sum(1 / k),
    rnx_auc(fixture$din, fixture$dout)
  )
  expect_equal(
    rnx_curve(fixture$din, fixture$din),
    c(rnx1 = 1, rnx2 = 1, rnx3 = 1)
  )
})

test_that("exact rank penalties accept unsorted k vectors", {
  fixture <- exact_rank_profile_fixture()
  k <- c(2, 1)
  expected_trustworthiness <- vapply(
    k,
    function(rank_k) trustworthiness(fixture$din, fixture$dout, rank_k),
    numeric(1)
  )
  expected_continuity <- vapply(
    k,
    function(rank_k) continuity(fixture$din, fixture$dout, rank_k),
    numeric(1)
  )

  expect_equal(
    trustworthiness(fixture$din, fixture$dout, k),
    stats::setNames(expected_trustworthiness, paste0("trustworthiness", k))
  )
  expect_equal(
    continuity(fixture$din, fixture$dout, k),
    stats::setNames(expected_continuity, paste0("continuity", k))
  )
  expect_equal(
    trustworthiness(fixture$din, fixture$din, k),
    c(trustworthiness2 = 1, trustworthiness1 = 1)
  )
  expect_equal(
    continuity(fixture$din, fixture$din, k),
    c(continuity2 = 1, continuity1 = 1)
  )
})

test_that("scalar exact rank penalty returns remain unchanged", {
  fixture <- exact_rank_profile_fixture()

  trust <- trustworthiness(fixture$din, fixture$dout, k = 1)
  cont <- continuity(fixture$din, fixture$dout, k = 1)

  expect_identical(trust, trustworthiness_exact(fixture$din, fixture$dout, 1L))
  expect_identical(cont, continuity_exact(fixture$din, fixture$dout, 1L))
  expect_null(names(trust))
  expect_null(names(cont))
})

test_that("exact rank profiles retain the first-tie policy after permutation", {
  fixture <- exact_rank_profile_fixture()
  permutation <- c(5, 2, 4, 1, 3)
  din <- fixture$din[permutation, permutation]
  dout <- fixture$dout[permutation, permutation]

  expect_equal(
    rnx_curve(din, dout),
    rnx_profile_reference(din, dout, seq_len(nrow(din) - 2L))
  )
  expect_equal(
    trustworthiness(din, dout, c(2, 1)),
    c(
      trustworthiness2 = trustworthiness(din, dout, 2),
      trustworthiness1 = trustworthiness(din, dout, 1)
    )
  )
  expect_equal(
    continuity(din, dout, c(2, 1)),
    c(
      continuity2 = continuity(din, dout, 2),
      continuity1 = continuity(din, dout, 1)
    )
  )
})

test_that("exact rank profiles can be worse than random", {
  din <- distance_matrix(matrix(c(0, 1, 3, 10), ncol = 1))
  dout <- distance_matrix(matrix(c(0, 10, 3, 1), ncol = 1))

  expect_lt(rnx_curve(din, dout, k = 1), 0)
})

test_that("exact rank profiles validate k", {
  fixture <- exact_rank_profile_fixture()
  din <- fixture$din
  dout <- fixture$dout

  expect_error(rnx_curve(din, dout, integer()), "positive integers")
  expect_error(rnx_curve(din, dout, c(1, 1)), "unique values")
  expect_error(rnx_curve(din, dout, k = 4), "less than")
  expect_error(rnx_curve(din, dout, k = 1.5), "positive integers")
  expect_error(
    rnx_curve(din, dout, k = c(1, .Machine$integer.max + 1)),
    "supported integer range"
  )
  expect_error(trustworthiness(din, dout, c(1, 1)), "unique values")
  expect_error(continuity(din, dout, c(1, 3)), "less than half")
  expect_error(
    trustworthiness(din, dout, c(1, .Machine$integer.max + 1)),
    "supported integer range"
  )
})

test_that("rnx_curve validates distance matrices", {
  fixture <- exact_rank_profile_fixture()
  nonfinite <- fixture$din
  nonfinite[1, 2] <- NA_real_

  expect_error(
    rnx_curve(nonfinite, fixture$dout),
    "`din` must contain only finite off-diagonal distances"
  )
  expect_error(
    rnx_curve(fixture$din, fixture$dout[-1, ]),
    "same dimensions"
  )
  expect_error(rnx_curve(diag(2), diag(2)), "at least three")
})

test_that("exact rank metrics do not touch the R RNG", {
  fixture <- exact_rank_profile_fixture()

  set.seed(13579)
  seed_before <- .Random.seed
  rnx_auc(fixture$din, fixture$dout)
  rnx_curve(fixture$din, fixture$dout)
  trustworthiness(fixture$din, fixture$dout, c(2, 1))
  continuity(fixture$din, fixture$dout, c(2, 1))
  expect_identical(.Random.seed, seed_before)

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
    rnx_auc(fixture$din, fixture$dout)
    rnx_curve(fixture$din, fixture$dout)
    trustworthiness(fixture$din, fixture$dout, c(2, 1))
    continuity(fixture$din, fixture$dout, c(2, 1))
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  }

  call_without_seed()
})
