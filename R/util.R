stime <- function() {
  format(Sys.time(), "%T")
}

# message with a time stamp
# appears only if called from an environment where a logical verbose = TRUE
# OR force = TRUE
tsmessage <- function(
  ...,
  domain = NULL,
  appendLF = TRUE,
  force = FALSE,
  time_stamp = TRUE
) {
  verbose <- get0("verbose", envir = sys.parent())

  if (force || (!is.null(verbose) && verbose)) {
    msg <- ""
    if (time_stamp) {
      msg <- paste0(stime(), " ")
    }
    message(msg, ..., domain = domain, appendLF = appendLF)
    utils::flush.console()
  }
}

# convert data frame to matrix using numeric columns
x2m <- function(X) {
  if (is.data.frame(X)) {
    numeric_cols <- vapply(X, is.numeric, logical(1))
    if (!any(numeric_cols)) {
      stop(
        "Data frames must contain at least one numeric column",
        call. = FALSE
      )
    }
    m <- as.matrix(X[, numeric_cols, drop = FALSE])
  } else if (!methods::is(X, "matrix")) {
    m <- as.matrix(X)
  } else {
    m <- X
  }
  if (!is.numeric(m)) {
    stop("Input data must be numeric", call. = FALSE)
  }
  if (nrow(m) == 0 || ncol(m) == 0) {
    stop(
      "Input data must contain at least one row and one column",
      call. = FALSE
    )
  }
  m
}

validate_positive_integer <- function(x, name) {
  if (
    !is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x < 1 ||
      x != floor(x)
  ) {
    stop(name, " must be a positive integer", call. = FALSE)
  }
  as.integer(x)
}

validate_positive_integer_vector <- function(x, name) {
  if (
    !is.numeric(x) ||
      length(x) < 1L ||
      anyNA(x) ||
      any(!is.finite(x)) ||
      any(x < 1) ||
      any(x != floor(x))
  ) {
    stop(name, " must contain positive integers", call. = FALSE)
  }
  as.integer(x)
}

# Add the (named) values in l2 to l1.
# Use to override default values in l1 with user-supplied values in l2
lmerge <- function(l1, l2) {
  for (name in names(l2)) {
    l1[[name]] <- l2[[name]]
  }
  l1
}

supported_distances <- function() {
  c("sqeuclidean", "euclidean", "cosine", "hamming", "manhattan", "correlation")
}

validate_distance <- function(distance) {
  match.arg(distance, supported_distances())
}
