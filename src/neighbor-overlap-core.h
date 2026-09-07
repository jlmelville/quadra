#ifndef QUADRA_NEIGHBOR_OVERLAP_CORE_H
#define QUADRA_NEIGHBOR_OVERLAP_CORE_H

#include <cstddef>
#include <vector>

#include "pforr.h"

namespace quadra::detail {

struct KQuery {
  std::size_t value;
  std::size_t output_col;
};

inline std::size_t matrix_offset(std::size_t row, std::size_t col,
                                 std::size_t nrow) {
  return row + (col * nrow);
}

inline void overlap_counts_inner(std::size_t begin, std::size_t end,
                                 const std::vector<std::size_t> &idx,
                                 const std::vector<std::size_t> &ref_idx,
                                 std::size_t n_obs,
                                 const std::vector<KQuery> &queries,
                                 std::vector<int> &counts) {
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

inline void compute_overlap_counts(const std::vector<std::size_t> &idx,
                                   const std::vector<std::size_t> &ref_idx,
                                   std::size_t n_obs,
                                   const std::vector<KQuery> &queries,
                                   std::vector<int> &counts,
                                   std::size_t n_threads) {
  auto worker = [&](std::size_t begin, std::size_t end) {
    overlap_counts_inner(begin, end, idx, ref_idx, n_obs, queries, counts);
  };

  pforr::parallel_for(0, n_obs, worker, n_threads);
}

} // namespace quadra::detail

#endif // QUADRA_NEIGHBOR_OVERLAP_CORE_H
