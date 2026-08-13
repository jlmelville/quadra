#' Average Area Under the ROC Curve
#'
#' Embedding quality measure.
#'
#' The ROC curve plots the true positive rate vs false positive rate.
#' This function calculates the curve N times, where N is the number of the
#' observations. The label of the Nth observation is set as the positive class
#' and then the other observations are ranked according to their distance from
#' the Nth observation in the output coordinates (lower distances being better).
#' Observations with the same label as the Nth observation count as positive
#' observations. The final reported result is the average over all observations.
#' Rows with undefined AUC values, such as rows that cannot define both positive
#' and negative examples, are excluded from overall and per-label averages. If
#' no rows remain for an average, that average is `NA_real_`.
#'
#' `av_auc` weights each defined query equally; `label_av` reports per-label
#' means in first-appearance order.
#'
#' Perfect retrieval results in an AUC of 1. For random retrieval gives a value
#' of 0.5.
#'
#' @note Use of this function requires that the `PRROC` package be installed.
#'
#' @param dm A finite square numeric distance matrix. The diagonal is excluded
#'   from each retrieval problem.
#' @param labels A label vector without missing values, with one label per row
#'   of `dm`.
#' @return A list with components in this order:
#'
#'   * `av_auc`: the AUC averaged over defined observation rows.
#'   * `label_av`: a named list of per-label means in first-appearance order.
#' @export
roc_auc <- function(dm, labels) {
  if (!requireNamespace("PRROC", quietly = TRUE)) {
    stop("roc_auc function requires 'PRROC' package")
  }
  validate_auc_inputs(dm, labels)
  summarize_retrieval_auc(dm, labels, roc_auc_row)
}

#' Area Under the Precision-Recall Curve for an Embedding
#'
#' Embedding quality measure.
#'
#' The PR curve plots precision (also known as positive predictive value, PPV)
#' against recall (also known as the true positive rate). The area under the
#' curve provides similar information compared to the area under the ROC curve,
#' but may be more appropriate when classes are highly imbalanced.
#'
#' This function calculates the PR curve N times, where N is the number of the
#' observations. The label of the Nth observation is set as the positive class
#' and then the other observations are ranked according to their distance from
#' the Nth observation in the output coordinates (lower distances being better).
#' Observations with the same label as the Nth observation count as positive
#' observations. The final reported result is the average over all observations.
#' Rows with undefined AUC values, such as rows that cannot define both positive
#' and negative examples, are excluded from overall and per-label averages. If no
#' rows remain for an average, that average is `NA_real_`.
#'
#' `av_auc` weights each defined query equally; `label_av` reports per-label
#' means in first-appearance order.
#'
#' Perfect retrieval results in an AUC of 1. For random retrieval, the value
#' is the proportion of the positive class labels for that curve.
#'
#' @note Use of this function requires that the `PRROC` package be installed.
#'
#' @inheritParams roc_auc
#' @return A list with components in this order:
#'
#'   * `av_auc`: the PR AUC averaged over defined observation rows.
#'   * `label_av`: a named list of per-label means in first-appearance order.
#' @references
#' Keilwagen, J., Grosse, I., & Grau, J. (2014).
#' Area under precision-recall curves for weighted and unweighted data.
#' *PloS One*, *9*(3), e92209.
#'
#' Davis, J., & Goadrich, M. (2006, June).
#' The relationship between Precision-Recall and ROC curves.
#' In *Proceedings of the 23rd international conference on Machine learning*
#' (pp. 233-240). ACM.
#' @export
pr_auc <- function(dm, labels) {
  if (!requireNamespace("PRROC", quietly = TRUE)) {
    stop("pr_auc function requires 'PRROC' package")
  }
  validate_auc_inputs(dm, labels)
  summarize_retrieval_auc(dm, labels, pr_auc_row)
}

validate_auc_inputs <- function(dm, labels) {
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
  if (nrow(dm) < 3L) {
    stop("`dm` must have at least 3 rows and columns", call. = FALSE)
  }

  if (!is.atomic(labels) || !is.null(dim(labels))) {
    stop("`labels` must be a one-dimensional vector", call. = FALSE)
  }
  if (length(labels) == 0L) {
    stop("`labels` must be nonempty", call. = FALSE)
  }
  if (length(labels) != nrow(dm)) {
    stop("`labels` length must equal the number of rows in `dm`", call. = FALSE)
  }
  if (anyNA(labels)) {
    stop("`labels` must not contain missing values", call. = FALSE)
  }

  invisible(NULL)
}

# Area Under the PR Curve of an Observation
#
# Embedding quality measure.
#
# The PR curve plots precision (also known as positive predictive value, PPV)
# against recall (also known as the true positive rate). The area under the
# curve provides similar information compared to the area under the ROC curve,
# but may be more appropriate when classes are highly imbalanced.
#
# This function calculates the curve with the label of the specified
# observation set as the positive class. The other observations are then
# ranked according to their distance from the ith observation
# (lower distances being better). Observations with the same label as the
# specified observation count as the positive observations.
#
# Perfect retrieval results in an AUC of 1. Random retrieval gives a value
# of the proportion of positive class with respect to the entire data set
# (e.g. if there are 20 observations with the positive class label in a
# dataset of 100, then the random AUC is 0.2).
#
# @note Use of this function requires that the \code{PRROC} package be
# installed.
#
# @param dm Distance matrix.
# @param labels Vector of labels, of the same size as the number of rows
# (or columns) in the distance matrix.
# @param i The row of the distance matrix to use in the PR calculation.
# @return Area Under the curve.
# @references
#
# Keilwagen, J., Grosse, I., & Grau, J. (2014).
# Area under precision-recall curves for weighted and unweighted data.
# \emph{PloS One}, \emph{9}(3), e92209.
#
# Davis, J., & Goadrich, M. (2006, June).
# The relationship between Precision-Recall and ROC curves.
# In \emph{Proceedings of the 23rd international conference on Machine
# learning}
# (pp. 233-240). ACM.
pr_auc_row <- function(dm, labels, i) {
  other_ind <- seq_len(nrow(dm))[-i]
  is_positive <- labels[other_ind] == labels[i]
  if (!any(is_positive) || all(is_positive)) {
    return(NA_real_)
  }

  pos_dist <- dm[i, other_ind[is_positive]]
  neg_dist <- dm[i, other_ind[!is_positive]]

  as.numeric(
    PRROC::pr.curve(
      scores.class0 = -pos_dist,
      scores.class1 = -neg_dist
    )$auc.davis.goadrich
  )
}

# Area Under the ROC Curve of an Observation
#
# Embedding quality measure.
#
# The ROC curve plots the true positive rate vs false positive rate.
# This function calculates the curve with the label of the specified
# observation set as the positive class. The other observations are then
# ranked according to their distance from the ith observation
# (lower distances being better). Observations with the same label as the
# specified observation count as the positive observations.
#
# Perfect retrieval results in an AUC of 1. For random retrieval gives a value
# of 0.5.
#
# @note Use of this function requires that the \code{PRROC} package be
# installed.
#
# @param dm Distance matrix.
# @param labels Vector of labels, of the same size as the number of rows
# (or columns) in the distance matrix.
# @param i The row of the distance matrix to use in the ROC calculation.
# @return Area Under the curve.
roc_auc_row <- function(dm, labels, i) {
  other_ind <- seq_len(nrow(dm))[-i]
  is_positive <- labels[other_ind] == labels[i]
  if (!any(is_positive) || all(is_positive)) {
    return(NA_real_)
  }

  pos_dist <- dm[i, other_ind[is_positive]]
  neg_dist <- dm[i, other_ind[!is_positive]]

  as.numeric(
    PRROC::roc.curve(
      scores.class0 = -pos_dist,
      scores.class1 = -neg_dist
    )$auc
  )
}

# Average Area Under a Curve
#
# Embedding quality measure.
#
# This function calculates a curve using the specified function, repeating the
# procedure N times, where N is the number of the observations. Each time
# a different row of the distance matrix is used. The label of the Nth
# observation is set as the positive class and then the other observations are
# ranked according to their distance from the Nth observation in the output
# coordinates (lower distances being better). Observations with the same label
# as the Nth observation count as positive observations. The final reported
# result is the average over all observations.
#
# @note Use of this function requires that the \code{PRROC} package be
# installed.
#
# @param dm Distance matrix.
# @param labels Vector of labels, of the same size as the number of rows
# (or columns) in the distance matrix.
# @param auc_row_fn A function which can calculate the Area Under a Curve
# for a particular quality measure. Should have the signature
# \code{auc_row_fn(dm, labels, i)} where \code{i} is the ith row of the
# distance matrix, and should return a scalar value giving the area under the
# curve using the ith row of the distance matrix.
# @return A list containing:
# \item{av_auc}{Area Under the curve, averaged over each observation.}
# The list also contains the average AUC per class label, with each average
# being named after the class label.
summarize_retrieval_auc <- function(dm, labels, auc_row_fn) {
  av_auc <- 0
  av_n <- 0
  n <- nrow(dm)
  label_names <- unique(as.character(labels))
  ns <- stats::setNames(as.list(integer(length(label_names))), label_names)
  label_av <- stats::setNames(
    as.list(numeric(length(label_names))),
    label_names
  )
  for (i in seq_len(n)) {
    label_index <- match(as.character(labels[i]), label_names)

    auc <- auc_row_fn(dm, labels, i)
    if (is.numeric(auc) && length(auc) == 1L && isTRUE(is.finite(auc))) {
      av_auc <- av_auc + auc
      av_n <- av_n + 1
      label_av[[label_index]] <- label_av[[label_index]] + auc
      ns[[label_index]] <- ns[[label_index]] + 1L
    }
  }
  for (label_index in seq_along(ns)) {
    if (ns[[label_index]] == 0L) {
      label_av[[label_index]] <- NA_real_
    } else {
      label_av[[label_index]] <- label_av[[label_index]] / ns[[label_index]]
    }
  }
  list(
    av_auc = if (av_n == 0L) NA_real_ else av_auc / av_n,
    label_av = label_av
  )
}
