format_timestamp <- function() {
  format(Sys.time(), "%T")
}

# message with a time stamp
# appears only if verbose = TRUE or force = TRUE
report_progress <- function(
  ...,
  verbose,
  domain = NULL,
  appendLF = TRUE,
  force = FALSE,
  time_stamp = TRUE
) {
  if (force || (!is.null(verbose) && verbose)) {
    msg <- ""
    if (time_stamp) {
      msg <- paste0(format_timestamp(), " ")
    }
    message(msg, ..., domain = domain, appendLF = appendLF)
    utils::flush.console()
  }
}

# convert data frame to matrix using numeric columns
prepare_input_matrix <- function(
  X,
  name = "Input data",
  allow_sparse = FALSE,
  require_finite = FALSE
) {
  sparse_input <- methods::is(X, "sparseMatrix")
  if (sparse_input && !allow_sparse) {
    stop(
      name,
      " must be dense; sparse matrices are not supported",
      call. = FALSE
    )
  }

  if (is.data.frame(X)) {
    numeric_cols <- vapply(X, is.numeric, logical(1))
    if (!any(numeric_cols)) {
      stop(
        "Data frames must contain at least one numeric column",
        call. = FALSE
      )
    }
    m <- as.matrix(X[, numeric_cols, drop = FALSE])
  } else if (sparse_input) {
    m <- X
  } else if (!methods::is(X, "matrix")) {
    m <- as.matrix(X)
  } else {
    m <- X
  }
  if (!sparse_input && !is.numeric(m)) {
    stop("Input data must be numeric", call. = FALSE)
  }
  if (nrow(m) == 0 || ncol(m) == 0) {
    stop(
      "Input data must contain at least one row and one column",
      call. = FALSE
    )
  }
  if (require_finite && (anyNA(m) || any(!is.finite(m)))) {
    stop(name, " must contain only finite values", call. = FALSE)
  }
  m
}

validate_positive_integer <- function(x, name) {
  if (
    !is.numeric(x) ||
      is.complex(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x < 1 ||
      x != floor(x)
  ) {
    stop(name, " must be a positive integer", call. = FALSE)
  }
  if (x > .Machine$integer.max) {
    stop(name, " exceeds the supported integer range", call. = FALSE)
  }
  as.integer(x)
}

validate_positive_integer_vector <- function(x, name) {
  if (
    !is.numeric(x) ||
      is.complex(x) ||
      length(x) < 1L ||
      anyNA(x) ||
      any(!is.finite(x)) ||
      any(x < 1) ||
      any(x != floor(x))
  ) {
    stop(name, " must contain positive integers", call. = FALSE)
  }
  if (any(x > .Machine$integer.max)) {
    stop(
      name,
      " contains values that exceed the supported integer range",
      call. = FALSE
    )
  }
  if (anyDuplicated(x)) {
    stop(name, " must contain unique values", call. = FALSE)
  }
  as.integer(x)
}

validate_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE", call. = FALSE)
  }
  x
}

validate_n_threads <- function(n_threads) {
  if (
    !is.numeric(n_threads) ||
      is.complex(n_threads) ||
      length(n_threads) != 1L ||
      is.na(n_threads) ||
      !is.finite(n_threads) ||
      n_threads < 0 ||
      n_threads != floor(n_threads) ||
      n_threads > .Machine$integer.max
  ) {
    stop(
      paste(
        "n_threads must be one finite, whole, nonnegative value within",
        "the supported integer range"
      ),
      call. = FALSE
    )
  }
  as.integer(n_threads)
}

validate_nn_args <- function(x, name) {
  if (!is.list(x)) {
    stop(name, " must be a list", call. = FALSE)
  }
  if (length(x) == 0L) {
    return(x)
  }

  argument_names <- names(x)
  if (
    is.null(argument_names) ||
      anyNA(argument_names) ||
      any(!nzchar(argument_names))
  ) {
    stop(name, " must have nonempty, nonmissing names", call. = FALSE)
  }
  if (anyDuplicated(argument_names)) {
    stop(name, " must have unique names", call. = FALSE)
  }

  package_controls <- c(
    "X",
    "data",
    "k",
    "metric",
    "nn_method",
    "n_threads",
    "verbose",
    "obs"
  )
  collision <- intersect(argument_names, package_controls)
  if (length(collision) > 0L) {
    stop(
      name,
      " must not include package-owned control '",
      collision[1L],
      "'",
      call. = FALSE
    )
  }
  x
}

supported_distances <- function() {
  c("sqeuclidean", "euclidean", "cosine", "hamming", "manhattan", "correlation")
}

validate_distance <- function(distance) {
  match.arg(distance, supported_distances())
}
