#include <algorithm>
#include <cstddef>
#include <vector>

#include <Rcpp.h>

#include "native-validation.h"
#include "neighbor-overlap-core.h"

using namespace Rcpp;

std::vector<quadra::detail::KQuery> prepare_k_queries(const IntegerVector &k,
                                                      std::size_t max_cols) {
  if (k.size() < 1) {
    stop("k must be nonempty");
  }

  std::vector<quadra::detail::KQuery> queries;
  queries.reserve(k.size());
  for (R_xlen_t i = 0; i < k.size(); ++i) {
    if (k[i] == NA_INTEGER || k[i] < 1) {
      stop("k must contain positive integers");
    }
    const auto value = static_cast<std::size_t>(k[i]);
    if (value > max_cols) {
      stop("k cannot be larger than the number of columns in idx or ref_idx");
    }
    queries.push_back({value, static_cast<std::size_t>(i)});
  }

  std::stable_sort(
      queries.begin(), queries.end(),
      [](const quadra::detail::KQuery &lhs, const quadra::detail::KQuery &rhs) {
        return lhs.value < rhs.value;
      });
  return queries;
}

std::vector<std::size_t> copy_neighbor_indices(const NumericMatrix &idx,
                                               std::size_t n_cols) {
  const std::size_t n_obs = idx.nrow();
  std::vector<std::size_t> copied(n_obs * n_cols);

  for (std::size_t col = 0; col < n_cols; ++col) {
    for (std::size_t row = 0; row < n_obs; ++row) {
      copied[quadra::detail::matrix_offset(row, col, n_obs)] =
          static_cast<std::size_t>(idx(row, col));
    }
  }
  return copied;
}

// [[Rcpp::export(rng = false)]]
IntegerMatrix neighbor_overlap_counts(const NumericMatrix &idx,
                                      const NumericMatrix &ref_idx,
                                      const IntegerVector &k,
                                      double n_threads = 0) {
  if (idx.nrow() != ref_idx.nrow()) {
    stop("idx and ref_idx must have the same number of rows");
  }

  const auto queries =
      prepare_k_queries(k, std::min(idx.ncol(), ref_idx.ncol()));
  const std::size_t n_obs = idx.nrow();
  const std::size_t max_k = queries.back().value;
  const std::size_t thread_count = quadra::checked_thread_count(n_threads);

  const auto idx_cpp = copy_neighbor_indices(idx, max_k);
  const auto ref_idx_cpp = copy_neighbor_indices(ref_idx, max_k);

  std::vector<int> counts(n_obs * static_cast<std::size_t>(k.size()), 0);
  quadra::detail::compute_overlap_counts(idx_cpp, ref_idx_cpp, n_obs, queries,
                                         counts, thread_count);

  IntegerMatrix result(n_obs, k.size());
  std::copy(counts.begin(), counts.end(), result.begin());
  return result;
}
