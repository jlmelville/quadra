#ifndef QUADRA_NATIVE_VALIDATION_H
#define QUADRA_NATIVE_VALIDATION_H

#include <cmath>
#include <cstddef>
#include <limits>

#include <Rcpp.h>

namespace quadra {

// Public argument contracts belong in R. These guards protect independently
// callable registered routines before narrowing, allocation, or thread setup.
inline std::size_t checked_thread_count(double n_threads) {
  if (!std::isfinite(n_threads) || n_threads < 0 ||
      n_threads != std::floor(n_threads)) {
    Rcpp::stop(
        "n_threads must be one finite, whole, nonnegative value within the "
        "supported integer range");
  }
  if (n_threads > static_cast<double>((std::numeric_limits<int>::max)())) {
    Rcpp::stop("n_threads exceeds the supported integer range");
  }
  return static_cast<std::size_t>(n_threads);
}

inline std::size_t checked_nonnegative_count(double value, const char *name) {
  if (!std::isfinite(value) || value < 0 || value != std::floor(value)) {
    Rcpp::stop("%s must be finite, whole, and nonnegative", name);
  }
  if (value > static_cast<double>((std::numeric_limits<int>::max)())) {
    Rcpp::stop("%s exceeds the supported integer range", name);
  }
  return static_cast<std::size_t>(value);
}

} // namespace quadra

#endif // QUADRA_NATIVE_VALIDATION_H
