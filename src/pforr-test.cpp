#include <stdexcept>

#include <Rcpp.h>

#include "pforr.h"

// [[Rcpp::export]]
void pforr_test_worker_exception(std::size_t n_threads = 2) {
  auto worker = [](std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      if (i == 1) {
        throw std::runtime_error("pforr worker exception");
      }
    }
  };

  pforr::parallel_for(static_cast<std::size_t>(0), static_cast<std::size_t>(4),
                      worker, n_threads, 1);
}
