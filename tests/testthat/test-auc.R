test_that("auc_mat averages valid rows by class", {
  dm <- matrix(0, nrow = 4, ncol = 4)
  labels <- c("a", "a", "b", "b")
  aucs <- c(0.8, 0.6, 0.5, 0.7)
  auc_row_fn <- function(dm, labels, i) aucs[[i]]

  res <- auc_mat(dm, labels, auc_row_fn)

  expect_equal(res$av_auc, 0.65)
  expect_equal(res$label_av$a, 0.7)
  expect_equal(res$label_av$b, 0.6)
})

test_that("auc_mat excludes undefined row AUCs", {
  dm <- matrix(0, nrow = 3, ncol = 3)
  labels <- c("a", "b", "b")
  aucs <- c(NaN, 0.4, 0.6)
  auc_row_fn <- function(dm, labels, i) aucs[[i]]

  res <- auc_mat(dm, labels, auc_row_fn)

  expect_equal(res$av_auc, 0.5)
  expect_true(is.na(res$label_av$a))
  expect_equal(res$label_av$b, 0.5)
})

test_that("auc_mat returns NA when all row AUCs are undefined", {
  dm <- matrix(0, nrow = 3, ncol = 3)
  labels <- rep("a", 3)
  auc_row_fn <- function(dm, labels, i) NaN

  res <- auc_mat(dm, labels, auc_row_fn)

  expect_true(is.na(res$av_auc))
  expect_true(is.na(res$label_av$a))
})
