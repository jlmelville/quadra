#' Average ROC AUC for Label Retrieval
#'
#' Ranks other observations by distance from each query and treats matching
#' labels as relevant.
#'
#' `av_auc` is the mean ROC AUC over defined queries; `label_av` gives per-label
#' means in first-appearance order. Queries without both relevant and irrelevant
#' observations are excluded. An empty average is `NA_real_`.
#'
#' @note [PRROC::roc.curve()] computes the area for each query.
#' @references Bamber (1975) relates ROC area to ranking probability.
#' \doi{10.1016/0022-2496(75)90001-2}
#'
#' Venna et al. (2010) describe the retrieval view of
#' dimensionality reduction. <https://www.jmlr.org/papers/v11/venna10a.html>
#'
#' Grau, Grosse and Keilwagen (2015) describe the PRROC implementation.
#' \doi{10.1093/bioinformatics/btv153}
#'
#' @section Scaling:
#' For large datasets, pass `coords` and sample `query_rows`. Each query
#' computes one row of Euclidean distances and ranks all other observations.
#' The distance workspace grows linearly with the candidate pool. Scoring all
#' queries still takes quadratic distance work. If you already have `dm`,
#' you can sample the same queries from it. Subsetting the coordinates or
#' the matrix before the call also changes the candidates.
#' @param dm Finite square numeric distance matrix. Supply exactly one of `dm`
#'   and `coords`.
#' @param labels A label vector without missing values, with one label per row
#'   of the supplied data.
#' @param coords Dense Euclidean coordinates, one observation per row.
#'   Numeric columns are used for data frames.
#' @param query_rows Unique one-based row positions to score, in the order
#'   requested. The default scores every row. Only the query itself is
#'   excluded from its candidate pool.
#' @param ret_extra Whether to return query AUCs in the requested row order
#'   and the number of defined queries.
#' @examplesIf requireNamespace("PRROC", quietly = TRUE)
#' x <- as.matrix(iris[, 1:4])
#' roc_auc(coords = x, labels = iris$Species, query_rows = c(1, 51, 101))
#' @return A list whose first components are:
#'
#'   * `av_auc`: the mean AUC over defined queries.
#'   * `label_av`: a named list of per-label means in first-appearance order.
#'
#'   With `ret_extra = TRUE`, the list also contains `query_auc`, with
#'   `NA_real_` for undefined queries, and the integer `n_defined`. When
#'   supplied, `query_rows` identifies the returned queries. Means use the
#'   selected queries. Labels without a defined query have mean `NA_real_`.
#' @seealso [pr_auc()], [class_pair_auc()]
#' @export
roc_auc <- function(
  dm = NULL,
  labels,
  ret_extra = FALSE,
  coords = NULL,
  query_rows = NULL
) {
  retrieval_auc(dm, labels, ret_extra, coords, query_rows, roc_auc_scores)
}

#' Average Precision-Recall AUC for Label Retrieval
#'
#' Precision-recall counterpart to [roc_auc()], often more informative for
#' imbalanced labels.
#'
#' Aggregation and undefined-query handling are the same as for [roc_auc()].
#'
#' @note Quadra uses the Davis-Goadrich area returned by [PRROC::pr.curve()].
#'
#' @inheritParams roc_auc
#' @inheritSection roc_auc Scaling
#' @inherit roc_auc return
#' @references
#' Keilwagen, J., Grosse, I., & Grau, J. (2014).
#' Area under precision-recall curves for weighted and unweighted data.
#' *PLOS ONE*, *9*(3), e92209. \doi{10.1371/journal.pone.0092209}
#'
#' Davis, J., & Goadrich, M. (2006, June).
#' The relationship between Precision-Recall and ROC curves.
#' In *Proceedings of the 23rd international conference on Machine learning*
#' (pp. 233-240). ACM. \doi{10.1145/1143844.1143874}
#'
#' Grau, Grosse and Keilwagen (2015) describe the PRROC implementation.
#' \doi{10.1093/bioinformatics/btv153}
#' @export
pr_auc <- function(
  dm = NULL,
  labels,
  ret_extra = FALSE,
  coords = NULL,
  query_rows = NULL
) {
  retrieval_auc(dm, labels, ret_extra, coords, query_rows, pr_auc_scores)
}

retrieval_auc <- function(dm, labels, ret_extra, coords, query_rows, auc_fn) {
  if (!requireNamespace("PRROC", quietly = TRUE)) {
    stop("Label retrieval AUC requires the PRROC package", call. = FALSE)
  }
  ret_extra <- validate_scalar_logical(ret_extra, "ret_extra")
  if (is.null(dm) == is.null(coords)) {
    stop("Supply exactly one of dm and coords", call. = FALSE)
  }
  if (is.null(coords)) {
    validate_auc_inputs(dm, labels)
    distance_row <- function(i) {
      dm[i, ]
    }
  } else {
    coords <- prepare_input_matrix(
      coords,
      name = "coords",
      require_finite = TRUE
    )
    if (nrow(coords) < 3L) {
      stop("coords must have at least 3 rows", call. = FALSE)
    }
    validate_label_vector(labels, nrow(coords), "coords")
    if (is.integer(coords)) {
      storage.mode(coords) <- "double"
    }
    distance_row <- function(i) {
      euclidean_distances_from_row(coords, i)
    }
  }
  rows <- validate_query_rows(query_rows, length(labels))
  result <- summarize_retrieval_auc(
    distance_row,
    labels,
    auc_fn,
    ret_extra,
    rows
  )
  if (ret_extra && !is.null(query_rows)) {
    result$query_rows <- rows
  }
  result
}

validate_auc_inputs <- function(dm, labels, min_rows = 3L) {
  if (!is.matrix(dm)) {
    stop("`dm` must be a matrix", call. = FALSE)
  }
  if (!is.numeric(dm)) {
    stop("`dm` must be numeric", call. = FALSE)
  }
  if (length(dm) == 0L) {
    stop("`dm` must be nonempty", call. = FALSE)
  }
  if (nrow(dm) != ncol(dm)) {
    stop("`dm` must be square", call. = FALSE)
  }
  if (!all(is.finite(dm))) {
    stop("`dm` must contain only finite values", call. = FALSE)
  }
  if (nrow(dm) < min_rows) {
    stop(
      "`dm` must have at least ",
      min_rows,
      " rows and columns",
      call. = FALSE
    )
  }

  validate_label_vector(labels, nrow(dm), "dm")
}

pr_auc_scores <- function(positive, negative) {
  PRROC::pr.curve(positive, negative)$auc.davis.goadrich
}

roc_auc_scores <- function(positive, negative) {
  PRROC::roc.curve(positive, negative)$auc
}

summarize_retrieval_auc <- function(
  distance_row,
  labels,
  auc_fn,
  ret_extra,
  query_rows = seq_along(labels)
) {
  av_auc <- 0
  av_n <- 0
  query_auc <- if (ret_extra) rep(NA_real_, length(query_rows)) else NULL
  label_names <- unique(as.character(labels))
  ns <- stats::setNames(as.list(integer(length(label_names))), label_names)
  label_av <- stats::setNames(
    as.list(numeric(length(label_names))),
    label_names
  )
  for (j in seq_along(query_rows)) {
    i <- query_rows[j]
    label_index <- match(as.character(labels[i]), label_names)

    positive <- labels[-i] == labels[i]
    if (!any(positive) || all(positive)) {
      next
    }
    scores <- -distance_row(i)[-i]
    auc <- auc_fn(scores[positive], scores[!positive])
    if (is.numeric(auc) && length(auc) == 1L && isTRUE(is.finite(auc))) {
      av_auc <- av_auc + auc
      av_n <- av_n + 1
      label_av[[label_index]] <- label_av[[label_index]] + auc
      ns[[label_index]] <- ns[[label_index]] + 1L
      if (ret_extra) {
        query_auc[[j]] <- auc
      }
    }
  }
  for (label_index in seq_along(ns)) {
    if (ns[[label_index]] == 0L) {
      label_av[[label_index]] <- NA_real_
    } else {
      label_av[[label_index]] <- label_av[[label_index]] / ns[[label_index]]
    }
  }
  result <- list(
    av_auc = if (av_n == 0L) NA_real_ else av_auc / av_n,
    label_av = label_av
  )
  if (!ret_extra) {
    return(result)
  }
  c(
    result,
    list(query_auc = query_auc, n_defined = as.integer(av_n))
  )
}

#' Retrieval AUC Within Each Pair of Classes
#'
#' Restrict both queries and candidates to two classes at a time. Each
#' defined query receives equal weight within its class pair.
#'
#' @section Scaling:
#' The input is a dense `N` by `N` distance matrix. Each class pair also
#' requires a square submatrix for its members. For a particular pair, select
#' its two classes before computing distances. Keeping all members of those
#' classes gives the same pair score as extracting it from the full analysis.
#' @param dm Finite square numeric distance matrix. The query itself is excluded.
#' @param labels A label vector without missing values, with one label per row
#'   of `dm`.
#' @param metrics Which AUCs to report. Defaults to both `"roc"` and `"pr"`.
#' @return A data frame with one row per unordered class pair and columns
#'   `label_a` and `label_b` in the same first-appearance order as
#'   [label_trustworthiness_continuity()]. Each requested metric adds its
#'   mean (`roc` or `pr`) and defined-query count (`n_defined_roc` or
#'   `n_defined_pr`). Queries from singleton classes are undefined. A pair
#'   of singletons has count zero and mean `NA_real_`.
#' @note Quadra applies [roc_auc()] and [pr_auc()] to each two-class subset.
#' Both use PRROC. The PR score uses Davis and Goadrich's area.
#' @references Davis and Goadrich (2006). \doi{10.1145/1143844.1143874}
#'
#' Grau, Grosse and Keilwagen (2015). PRROC. \doi{10.1093/bioinformatics/btv153}
#' @examplesIf requireNamespace("PRROC", quietly = TRUE)
#' class_pair_auc(as.matrix(dist(iris[, -5])), iris$Species)
#' @export
class_pair_auc <- function(dm, labels, metrics = c("roc", "pr")) {
  metrics <- unique(match.arg(metrics, c("roc", "pr"), several.ok = TRUE))
  if (!requireNamespace("PRROC", quietly = TRUE)) {
    stop("class_pair_auc requires the 'PRROC' package", call. = FALSE)
  }
  validate_auc_inputs(dm, labels, min_rows = 2L)
  classes <- unique(labels)
  if (length(classes) < 2L) {
    stop("labels must contain at least two distinct labels", call. = FALSE)
  }
  groups <- split(seq_along(labels), match(labels, classes))
  pairs <- utils::combn(seq_along(classes), 2L)
  result <- data.frame(
    label_a = classes[pairs[1L, ]],
    label_b = classes[pairs[2L, ]]
  )
  row_functions <- list(roc = roc_auc_scores, pr = pr_auc_scores)
  for (metric in metrics) {
    result[[metric]] <- NA_real_
    result[[paste0("n_defined_", metric)]] <- 0L
  }
  for (i in seq_len(nrow(result))) {
    rows <- sort(c(groups[[pairs[1L, i]]], groups[[pairs[2L, i]]]))
    distances <- dm[rows, rows, drop = FALSE]
    for (metric in metrics) {
      summary <- summarize_retrieval_auc(
        function(i) {
          distances[i, ]
        },
        labels[rows],
        row_functions[[metric]],
        TRUE
      )
      result[[metric]][i] <- summary$av_auc
      result[[paste0("n_defined_", metric)]][i] <- summary$n_defined
    }
  }
  result
}
