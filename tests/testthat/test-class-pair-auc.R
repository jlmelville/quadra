test_that("pair-restricted AUC matches single-metric calls and handles singletons", {
  skip_if_not_installed("PRROC")
  two <- class_pair_auc(as.matrix(dist(0:1)), c("a", "b"))
  expect_identical(two$n_defined_pr, 0L)
  expect_true(is.na(two$pr))
  labels <- factor(c("b", "a", "b", "c", "a", "d"))
  dm <- as.matrix(dist(c(0, 1, 2, 2, 4, 5)))
  result <- class_pair_auc(dm, labels)
  coords <- as.matrix(c(0, 1, 2, 2, 4, 5))
  centroids <- label_trustworthiness_continuity(coords, coords, labels, TRUE)
  expect_identical(result[1:2], centroids$pairs[1:2])
  for (i in seq_len(nrow(result))) {
    rows <- labels %in%
      c(as.character(result$label_a[i]), as.character(result$label_b[i]))
    if (sum(rows) >= 3L) {
      for (metric in c("roc", "pr")) {
        fn <- if (metric == "roc") roc_auc else pr_auc
        direct <- fn(dm[rows, rows], labels[rows], TRUE)
        expect_equal(result[[metric]][i], direct$av_auc)
        expect_identical(
          result[[paste0("n_defined_", metric)]][i],
          direct$n_defined
        )
      }
    } else {
      expect_true(is.na(result$pr[i]) && is.na(result$roc[i]))
      expect_identical(result$n_defined_roc[i], 0L)
    }
  }
  expect_identical(
    class_pair_auc(dm, labels, "pr"),
    result[c("label_a", "label_b", "pr", "n_defined_pr")]
  )
})
