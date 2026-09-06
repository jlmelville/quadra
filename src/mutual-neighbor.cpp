#include <cstddef>
#include <vector>

#include <Rcpp.h>

// Indices have been validated and self-neighbors removed by resolve_nn_graph.
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector mutual_neighbor_counts_cpp(const Rcpp::NumericMatrix &idx,
                                               int k) {
  if (k == NA_INTEGER || k < 1 || k > idx.ncol()) {
    Rcpp::stop(
        "k must be positive and no larger than the number of columns in idx");
  }

  // The reciprocal scan needs whole rows; R matrices store columns
  // contiguously.
  const auto width = static_cast<std::size_t>(k);
  std::vector<int> neighbors(static_cast<std::size_t>(idx.nrow()) * width);
  for (int i = 0; i < idx.nrow(); ++i) {
    if (i % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }
    for (int neighbor = 0; neighbor < k; ++neighbor) {
      neighbors[static_cast<std::size_t>(i) * width + neighbor] =
          static_cast<int>(idx(i, neighbor)) - 1;
    }
  }

  Rcpp::IntegerVector counts(idx.nrow());
  for (int i = 0; i < idx.nrow(); ++i) {
    if (i % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }
    int count = 0;
    for (int neighbor = 0; neighbor < k; ++neighbor) {
      const int j = neighbors[static_cast<std::size_t>(i) * width + neighbor];
      const auto offset = static_cast<std::size_t>(j) * width;
      for (int reciprocal = 0; reciprocal < k; ++reciprocal) {
        if (neighbors[offset + reciprocal] == i) {
          ++count;
          break;
        }
      }
    }
    counts[i] = count;
  }
  return counts;
}
