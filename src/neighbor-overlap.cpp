#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include <Rcpp.h>

#include "pforr.h"

using namespace Rcpp;

struct KQuery {
  std::size_t value;
  std::size_t output_col;
};

std::size_t matrix_offset(std::size_t row, std::size_t col, std::size_t nrow) {
  return row + (col * nrow);
}

std::vector<KQuery> prepare_k_queries(const IntegerVector& k,
                                      std::size_t max_cols) {
  if (k.size() < 1) {
    stop("k must contain positive integers");
  }

  std::vector<KQuery> queries;
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

  std::stable_sort(queries.begin(), queries.end(),
                   [](const KQuery& lhs, const KQuery& rhs) {
                     return lhs.value < rhs.value;
                   });
  return queries;
}

std::size_t checked_neighbor_index(double value, std::size_t n_obs,
                                   const char* name) {
  if (!R_finite(value) || value < 1.0 || value > static_cast<double>(n_obs) ||
      value != std::floor(value)) {
    stop("%s must contain finite integer indices between 1 and the number of "
         "rows",
         name);
  }
  return static_cast<std::size_t>(value);
}

std::vector<std::size_t> copy_neighbor_indices(const NumericMatrix& idx,
                                               std::size_t n_cols,
                                               const char* name) {
  const std::size_t n_obs = idx.nrow();
  std::vector<std::size_t> copied(n_obs * n_cols);

  for (std::size_t col = 0; col < n_cols; ++col) {
    for (std::size_t row = 0; row < n_obs; ++row) {
      copied[matrix_offset(row, col, n_obs)] =
          checked_neighbor_index(idx(row, col), n_obs, name);
    }
  }
  return copied;
}

void overlap_counts_inner(std::size_t begin, std::size_t end,
                          const std::vector<std::size_t>& idx,
                          const std::vector<std::size_t>& ref_idx,
                          std::size_t n_obs, const std::vector<KQuery>& queries,
                          std::vector<int>& counts) {
  std::vector<std::size_t> idx_seen(n_obs, 0);
  std::vector<std::size_t> ref_seen(n_obs, 0);
  std::size_t row_token = 1;

  for (std::size_t row = begin; row < end; ++row, ++row_token) {
    int overlap = 0;
    std::size_t query_pos = 0;

    for (std::size_t pos = 1; pos <= queries.back().value; ++pos) {
      const auto idx_value = idx[matrix_offset(row, pos - 1, n_obs)] - 1;
      if (idx_seen[idx_value] != row_token) {
        idx_seen[idx_value] = row_token;
        if (ref_seen[idx_value] == row_token) {
          ++overlap;
        }
      }

      const auto ref_value = ref_idx[matrix_offset(row, pos - 1, n_obs)] - 1;
      if (ref_seen[ref_value] != row_token) {
        ref_seen[ref_value] = row_token;
        if (idx_seen[ref_value] == row_token) {
          ++overlap;
        }
      }

      while (query_pos < queries.size() && queries[query_pos].value == pos) {
        counts[matrix_offset(row, queries[query_pos].output_col, n_obs)] =
            overlap;
        ++query_pos;
      }
    }
  }
}

// [[Rcpp::export]]
IntegerMatrix neighbor_overlap_counts(const NumericMatrix& idx,
                                      const NumericMatrix& ref_idx,
                                      const IntegerVector& k,
                                      std::size_t n_threads = 0) {
  if (idx.nrow() != ref_idx.nrow()) {
    stop("idx and ref_idx must have the same number of rows");
  }

  const auto queries =
      prepare_k_queries(k, std::min(idx.ncol(), ref_idx.ncol()));
  const std::size_t n_obs = idx.nrow();
  const std::size_t max_k = queries.back().value;

  const auto idx_cpp = copy_neighbor_indices(idx, max_k, "idx");
  const auto ref_idx_cpp = copy_neighbor_indices(ref_idx, max_k, "ref_idx");

  std::vector<int> counts(n_obs * static_cast<std::size_t>(k.size()), 0);
  auto worker = [&](std::size_t begin, std::size_t end) {
    overlap_counts_inner(begin, end, idx_cpp, ref_idx_cpp, n_obs, queries,
                         counts);
  };

  pforr::parallel_for(0, n_obs, worker, n_threads);

  IntegerMatrix result(n_obs, k.size());
  std::copy(counts.begin(), counts.end(), result.begin());
  return result;
}
