retrieval_query_oracle <- function(dm, labels, i, curve) {
  other <- seq_len(nrow(dm))[-i]
  positive <- labels[other] == labels[i]
  if (!any(positive) || all(positive)) {
    return(NA_real_)
  }

  scores0 <- -dm[i, other[positive]]
  scores1 <- -dm[i, other[!positive]]
  if (curve == "roc") {
    return(as.numeric(PRROC::roc.curve(scores0, scores1)$auc))
  }
  as.numeric(PRROC::pr.curve(scores0, scores1)$auc.davis.goadrich)
}

test_that("public AUC metrics reject invalid distance matrices", {
  skip_if_not_installed("PRROC")

  for (auc_fn in list(roc_auc, pr_auc)) {
    expect_error(auc_fn(1:9, c("a", "a", "b")), "matrix")
    expect_error(
      auc_fn(matrix(seq_len(12), nrow = 3), c("a", "a", "b")),
      "square"
    )
    expect_error(
      auc_fn(matrix("0", nrow = 3, ncol = 3), c("a", "a", "b")),
      "numeric"
    )
    expect_error(
      auc_fn(matrix(numeric(), nrow = 0, ncol = 0), character()),
      "nonempty"
    )

    for (nonfinite in c(NA_real_, NaN, Inf, -Inf)) {
      dm <- diag(3)
      dm[1, 2] <- nonfinite
      expect_error(auc_fn(dm, c("a", "a", "b")), "finite")
    }

    expect_error(auc_fn(matrix(0, nrow = 1, ncol = 1), "a"), "at least 3")
    expect_error(
      auc_fn(matrix(0, nrow = 2, ncol = 2), c("a", "b")),
      "at least 3"
    )
  }
})

test_that("public AUC metrics reject invalid labels", {
  skip_if_not_installed("PRROC")

  dm <- diag(3)
  for (auc_fn in list(roc_auc, pr_auc)) {
    expect_error(auc_fn(dm, matrix(c("a", "a", "b"), ncol = 1)), "vector")
    expect_error(auc_fn(dm, list("a", "a", "b")), "vector")
    expect_error(auc_fn(dm, character()), "nonempty")
    expect_error(auc_fn(dm, c("a", "b")), "length")
    expect_error(auc_fn(dm, c("a", "a", "b", "b")), "length")
    expect_error(auc_fn(dm, c("a", NA_character_, "b")), "missing")
  }
})

test_that("public AUC metrics preserve perfect retrieval with factor labels", {
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
  labels <- factor(c("a", "a", "b", "b"), levels = c("b", "a", "unused"))

  for (auc_fn in list(roc_auc, pr_auc)) {
    result <- auc_fn(dm, labels)

    expect_identical(names(result), c("av_auc", "label_av"))
    expect_equal(result$av_auc, 1)
    expect_identical(names(result$label_av), c("a", "b"))
    expect_equal(result$label_av$a, 1)
    expect_equal(result$label_av$b, 1)
  }
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

  for (auc_fn in list(roc_auc, pr_auc)) {
    result <- auc_fn(dm, labels)

    expect_equal(result$av_auc, 1)
    expect_equal(result$label_av$a, 1)
    expect_true(is.na(result$label_av$b))
  }
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

  for (auc_fn in list(roc_auc, pr_auc)) {
    one_label <- auc_fn(dm, rep(TRUE, 3))
    expect_true(is.na(one_label$av_auc))
    expect_identical(one_label$label_av, setNames(list(NA_real_), "TRUE"))

    unique_labels <- auc_fn(dm, seq_len(3))
    expect_true(is.na(unique_labels$av_auc))
    expect_identical(names(unique_labels$label_av), c("1", "2", "3"))
    expect_true(all(vapply(unique_labels$label_av, is.na, logical(1))))
  }
})

test_that("public AUC metrics weight defined rows and preserve label order", {
  skip_if_not_installed("PRROC")

  # The first two classes have two and three defined rows, respectively. The
  # singleton final class is represented in the result but not a denominator.
  # fmt: skip
  dm <- matrix(
    c(
      0, 3, 4, 1, 5, 6,
      4, 0, 1, 5, 2, 6,
      4, 1, 0, 5, 2, 6,
      1, 3, 4, 0, 5, 6,
      1, 5, 6, 2, 0, 3,
      1, 2, 3, 4, 5, 0
    ),
    nrow = 6,
    byrow = TRUE
  )
  labels <- c("b", "a", "a", "b", "a", "c")

  for (auc_fn in list(roc_auc, pr_auc)) {
    result <- auc_fn(dm, labels)
    weighted_mean <- (2 * result$label_av$b + 3 * result$label_av$a) / 5
    unweighted_mean <- mean(c(result$label_av$b, result$label_av$a))

    expect_identical(names(result$label_av), c("b", "a", "c"))
    expect_true(is.na(result$label_av$c))
    expect_equal(result$av_auc, weighted_mean)
    expect_false(isTRUE(all.equal(result$av_auc, unweighted_mean)))
  }
})

test_that("public AUC metrics exclude diagonal distances", {
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
  changed_diagonal <- dm
  diag(changed_diagonal) <- c(-1e6, 1e6, -1e6, 1e6)
  labels <- c("a", "a", "b", "b")

  for (auc_fn in list(roc_auc, pr_auc)) {
    expect_equal(auc_fn(changed_diagonal, labels), auc_fn(dm, labels))
  }
})

test_that("public AUC details expose row-aligned query values", {
  skip_if_not_installed("PRROC")

  # Classes b and a have two and three defined queries; c is undefined.
  # fmt: skip
  dm <- matrix(
    c(
      0, 3, 4, 1, 5, 6,
      4, 0, 1, 5, 2, 6,
      4, 1, 0, 5, 2, 6,
      1, 3, 4, 0, 5, 6,
      1, 5, 6, 2, 0, 3,
      1, 2, 3, 4, 5, 0
    ),
    nrow = 6,
    byrow = TRUE
  )
  labels <- c("b", "a", "a", "b", "a", "c")
  cases <- list(
    list(metric = roc_auc, curve = "roc"),
    list(metric = pr_auc, curve = "pr")
  )

  for (case in cases) {
    compact <- case$metric(dm, labels)
    detailed <- case$metric(dm, labels, ret_extra = TRUE)
    expected_query <- vapply(
      seq_len(nrow(dm)),
      function(i) retrieval_query_oracle(dm, labels, i, case$curve),
      numeric(1)
    )

    expect_identical(case$metric(dm, labels, ret_extra = FALSE), compact)
    expect_identical(case$metric(dm, labels, TRUE), detailed)
    expect_identical(
      names(detailed),
      c(
        "av_auc",
        "label_av",
        "query_auc",
        "n_defined"
      )
    )
    expect_identical(detailed[1:2], compact)
    expect_equal(detailed$query_auc, expected_query)
    expect_identical(detailed$n_defined, 5L)

    defined <- is.finite(detailed$query_auc)
    total <- 0
    for (auc in detailed$query_auc[defined]) {
      total <- total + auc
    }
    expect_identical(detailed$av_auc, total / detailed$n_defined)
    expect_equal(
      unlist(detailed$label_av[c("b", "a")], use.names = FALSE),
      c(
        mean(detailed$query_auc[labels == "b"]),
        mean(detailed$query_auc[labels == "a"])
      )
    )
    expect_true(is.na(detailed$label_av$c))
  }
})

test_that("public AUC details preserve label and row identities", {
  skip_if_not_installed("PRROC")

  # fmt: skip
  dm <- matrix(
    c(
      0, 3, 4, 1, 5, 6,
      4, 0, 1, 5, 2, 6,
      4, 1, 0, 5, 2, 6,
      1, 3, 4, 0, 5, 6,
      1, 5, 6, 2, 0, 3,
      1, 2, 3, 4, 5, 0
    ),
    nrow = 6,
    byrow = TRUE
  )
  labels <- c("b", "a", "a", "b", "a", "c")
  factor_labels <- factor(labels, levels = c("c", "b", "a", "unused"))
  permutation <- c(6, 3, 1, 5, 2, 4)

  for (auc_fn in list(roc_auc, pr_auc)) {
    original <- auc_fn(dm, labels, ret_extra = TRUE)
    factor_result <- auc_fn(dm, factor_labels, ret_extra = TRUE)
    permuted <- auc_fn(
      dm[permutation, permutation],
      labels[permutation],
      ret_extra = TRUE
    )

    expect_identical(factor_result, original)
    expect_equal(permuted$query_auc, original$query_auc[permutation])
    expect_equal(permuted$av_auc, original$av_auc)
    expect_identical(
      names(permuted$label_av),
      unique(labels[permutation])
    )
    expect_equal(
      unlist(permuted$label_av[names(original$label_av)]),
      unlist(original$label_av)
    )
  }
})

test_that("public AUC details retain all undefined queries", {
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

  for (auc_fn in list(roc_auc, pr_auc)) {
    for (labels in list(rep(TRUE, 3), seq_len(3))) {
      compact <- auc_fn(dm, labels)
      detailed <- auc_fn(dm, labels, ret_extra = TRUE)

      expect_identical(detailed[1:2], compact)
      expect_identical(detailed$query_auc, rep(NA_real_, 3))
      expect_identical(detailed$n_defined, 0L)
    }
  }
})

test_that("public AUC metrics validate ret_extra", {
  skip_if_not_installed("PRROC")

  dm <- diag(3)
  labels <- c("a", "a", "b")
  for (auc_fn in list(roc_auc, pr_auc)) {
    expect_error(
      auc_fn(dm, labels, ret_extra = NA),
      "ret_extra must be TRUE or FALSE"
    )
    expect_error(
      auc_fn(dm, labels, ret_extra = 1),
      "ret_extra must be TRUE or FALSE"
    )
  }
})
