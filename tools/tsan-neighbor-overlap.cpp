#include "neighbor-overlap-core.h"

#include <cstddef>
#include <iostream>
#include <vector>

namespace {

using quadra::detail::compute_overlap_counts;
using quadra::detail::KQuery;
using quadra::detail::matrix_offset;

bool overlap_counts_match(std::size_t n_threads,
                          const std::vector<std::size_t> &idx,
                          const std::vector<std::size_t> &ref_idx,
                          std::size_t n_obs, const std::vector<KQuery> &queries,
                          const std::vector<int> &expected) {
  std::vector<int> actual(n_obs * queries.size(), 0);
  compute_overlap_counts(idx, ref_idx, n_obs, queries, actual, n_threads);
  return actual == expected;
}

} // namespace

int main() {
  constexpr std::size_t n_obs = 257;
  constexpr std::size_t n_cols = 32;

  // Queries remain sorted by k while output columns deliberately do not.
  const std::vector<KQuery> queries{{8, 2}, {16, 0}, {32, 1}};
  std::vector<std::size_t> idx(n_obs * n_cols);
  std::vector<std::size_t> ref_idx(n_obs * n_cols);
  for (std::size_t col = 0; col < n_cols; ++col) {
    for (std::size_t row = 0; row < n_obs; ++row) {
      idx[matrix_offset(row, col, n_obs)] = ((row + col + 1) % n_obs) + 1;
      ref_idx[matrix_offset(row, col, n_obs)] = ((row + col + 5) % n_obs) + 1;
    }
  }

  std::vector<int> expected(n_obs * queries.size(), 0);
  compute_overlap_counts(idx, ref_idx, n_obs, queries, expected, 1);
  for (std::size_t row = 0; row < n_obs; ++row) {
    if (expected[matrix_offset(row, 0, n_obs)] != 12 ||
        expected[matrix_offset(row, 1, n_obs)] != 28 ||
        expected[matrix_offset(row, 2, n_obs)] != 4) {
      std::cerr << "serial overlap oracle has unexpected output mapping\n";
      return 1;
    }
  }

  const std::size_t thread_counts[] = {1, 2, 3, 4, 8};
  for (std::size_t repetition = 0; repetition < 200; ++repetition) {
    for (const std::size_t n_threads : thread_counts) {
      if (!overlap_counts_match(n_threads, idx, ref_idx, n_obs, queries,
                                expected)) {
        std::cerr << "parallel overlap mismatch at repetition " << repetition
                  << " with " << n_threads << " requested threads\n";
        return 1;
      }
    }
  }

  return 0;
}
