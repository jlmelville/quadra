test_that("public AUC metrics average valid rows by class", {
  skip_if_not_installed("PRROC")

  # fmt: skip
  dm <- matrix(
    c(
      0, 1, 4, 5,
      1, 0, 5, 4,
      4, 5, 0, 1,
      5, 4, 1, 0
    ),
    nrow = 4,
    byrow = TRUE
  )
  labels <- c("a", "a", "b", "b")

  roc_res <- roc_auc(dm, labels)
  pr_res <- pr_auc(dm, labels)

  expect_equal(roc_res$av_auc, 1)
  expect_equal(roc_res$label_av$a, 1)
  expect_equal(roc_res$label_av$b, 1)

  expect_equal(pr_res$av_auc, 1)
  expect_equal(pr_res$label_av$a, 1)
  expect_equal(pr_res$label_av$b, 1)
})

test_that("public AUC metrics exclude undefined singleton-class rows", {
  skip_if_not_installed("PRROC")

  # fmt: skip
  dm <- matrix(
    c(
      0, 1, 4,
      1, 0, 5,
      4, 5, 0
    ),
    nrow = 3,
    byrow = TRUE
  )
  labels <- c("a", "a", "b")

  roc_res <- roc_auc(dm, labels)
  pr_res <- pr_auc(dm, labels)

  expect_equal(roc_res$av_auc, 1)
  expect_equal(roc_res$label_av$a, 1)
  expect_true(is.na(roc_res$label_av$b))

  expect_equal(pr_res$av_auc, 1)
  expect_equal(pr_res$label_av$a, 1)
  expect_true(is.na(pr_res$label_av$b))
})

test_that("public AUC metrics return NA when all rows are undefined", {
  skip_if_not_installed("PRROC")

  # fmt: skip
  dm <- matrix(
    c(
      0, 1, 2,
      1, 0, 3,
      2, 3, 0
    ),
    nrow = 3,
    byrow = TRUE
  )
  labels <- c("a", "b", "c")

  roc_res <- roc_auc(dm, labels)
  pr_res <- pr_auc(dm, labels)

  expect_true(is.na(roc_res$av_auc))
  expect_true(is.na(roc_res$label_av$a))
  expect_true(is.na(roc_res$label_av$b))
  expect_true(is.na(roc_res$label_av$c))

  expect_true(is.na(pr_res$av_auc))
  expect_true(is.na(pr_res$label_av$a))
  expect_true(is.na(pr_res$label_av$b))
  expect_true(is.na(pr_res$label_av$c))
})
