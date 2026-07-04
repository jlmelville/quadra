#include <algorithm>
#include <cstddef>
#include <vector>

#include <Rcpp.h>

using namespace Rcpp;

bool rank_value_less(double lhs, double rhs) {
  const bool lhs_missing = ISNAN(lhs);
  const bool rhs_missing = ISNAN(rhs);
  if (lhs_missing || rhs_missing) {
    return !lhs_missing && rhs_missing;
  }
  return lhs < rhs;
}

std::vector<std::size_t> row_ranks_first_ties(const NumericMatrix& distances,
                                              std::size_t row) {
  const std::size_t n_obs = distances.nrow();

  std::vector<std::size_t> order;
  order.reserve(n_obs - 1);
  for (std::size_t col = 0; col < n_obs; ++col) {
    if (col != row) {
      order.push_back(col);
    }
  }

  std::stable_sort(
      order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
        return rank_value_less(distances(row, lhs), distances(row, rhs));
      });

  std::vector<std::size_t> ranks(n_obs, 0);
  for (std::size_t rank = 0; rank < order.size(); ++rank) {
    ranks[order[rank]] = rank + 1;
  }
  return ranks;
}

double rnx_auc_from_max_rank_histogram(const std::vector<double>& histogram,
                                       std::size_t n_obs) {
  const std::size_t n_ranks = n_obs - 1;
  if (n_ranks < 2) {
    return NA_REAL;
  }

  double top_left = 0;
  double num = 0;
  double den = 0;
  for (std::size_t k = 1; k < n_ranks; ++k) {
    top_left += histogram[k];
    const double qnx = top_left / (static_cast<double>(k) * n_obs);
    const double rnx = ((qnx * n_ranks) - k) / static_cast<double>(n_ranks - k);
    num += rnx / k;
    den += 1.0 / k;
  }
  return num / den;
}

double rank_penalty_score(const NumericMatrix& din, const NumericMatrix& dout,
                          std::size_t k, bool continuity) {
  const std::size_t n_obs = din.nrow();
  double penalty = 0;

  for (std::size_t i = 0; i < n_obs; ++i) {
    const auto ranks_in = row_ranks_first_ties(din, i);
    const auto ranks_out = row_ranks_first_ties(dout, i);

    for (std::size_t j = 0; j < n_obs; ++j) {
      if (j == i) {
        continue;
      }
      const std::size_t neighborhood_rank =
          continuity ? ranks_in[j] : ranks_out[j];
      const std::size_t penalty_rank = continuity ? ranks_out[j] : ranks_in[j];
      if (neighborhood_rank <= k && penalty_rank > k) {
        penalty += static_cast<double>(penalty_rank - k);
      }
    }
  }

  const double normalization = static_cast<double>(n_obs) *
                               static_cast<double>(k) *
                               static_cast<double>((2 * n_obs) - (3 * k) - 1);
  return 1.0 - ((2.0 * penalty) / normalization);
}

double exact_rank_penalty_metric(const NumericMatrix& din,
                                 const NumericMatrix& dout, int k,
                                 bool continuity) {
  if (din.nrow() != dout.nrow() || din.ncol() != dout.ncol() ||
      din.nrow() != din.ncol()) {
    stop("din and dout must be square distance matrices with the same "
         "dimensions");
  }
  if (k < 1) {
    stop("k must be a positive integer");
  }

  const std::size_t n_obs = din.nrow();
  const std::size_t rank_k = static_cast<std::size_t>(k);
  if (n_obs < 3) {
    stop("trustworthiness and continuity require at least three observations");
  }
  if ((2 * rank_k) >= n_obs) {
    stop("k must be less than half the number of observations");
  }

  return rank_penalty_score(din, dout, rank_k, continuity);
}

// [[Rcpp::export]]
double rnx_auc_direct(const NumericMatrix& din, const NumericMatrix& dout) {
  const std::size_t n_obs = din.nrow();
  const std::size_t n_ranks = n_obs - 1;

  std::vector<double> max_rank_histogram(n_ranks + 1, 0);
  for (std::size_t i = 0; i < n_obs; ++i) {
    const auto ranks_in = row_ranks_first_ties(din, i);
    const auto ranks_out = row_ranks_first_ties(dout, i);

    for (std::size_t j = 0; j < n_obs; ++j) {
      if (j == i) {
        continue;
      }
      const std::size_t max_rank = std::max(ranks_in[j], ranks_out[j]);
      ++max_rank_histogram[max_rank];
    }
  }

  return rnx_auc_from_max_rank_histogram(max_rank_histogram, n_obs);
}

// [[Rcpp::export]]
double trustworthiness_exact(const NumericMatrix& din,
                             const NumericMatrix& dout, int k) {
  return exact_rank_penalty_metric(din, dout, k, false);
}

// [[Rcpp::export]]
double continuity_exact(const NumericMatrix& din, const NumericMatrix& dout,
                        int k) {
  return exact_rank_penalty_metric(din, dout, k, true);
}
