#include <algorithm>
#include <cstddef>
#include <vector>

#include <Rcpp.h>

using namespace Rcpp;

std::size_t validate_exact_distance_matrices(const NumericMatrix &din,
                                             const NumericMatrix &dout) {
  if (din.nrow() != dout.nrow() || din.ncol() != dout.ncol()) {
    stop("din and dout must have the same dimensions");
  }
  if (din.nrow() != din.ncol()) {
    stop("din and dout must be square distance matrices");
  }
  if (din.nrow() < 3) {
    stop("exact rank metrics require at least three observations");
  }
  return static_cast<std::size_t>(din.nrow());
}

bool rank_value_less(double lhs, double rhs) {
  const bool lhs_missing = ISNAN(lhs);
  const bool rhs_missing = ISNAN(rhs);
  if (lhs_missing || rhs_missing) {
    return !lhs_missing && rhs_missing;
  }
  return lhs < rhs;
}

std::vector<std::size_t> row_ranks_first_ties(const NumericMatrix &distances,
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

std::vector<double> max_rank_histogram(const NumericMatrix &din,
                                       const NumericMatrix &dout,
                                       std::size_t n_obs) {
  const std::size_t n_ranks = n_obs - 1;
  std::vector<double> histogram(n_ranks + 1, 0);

  for (std::size_t i = 0; i < n_obs; ++i) {
    const auto ranks_in = row_ranks_first_ties(din, i);
    const auto ranks_out = row_ranks_first_ties(dout, i);

    for (std::size_t j = 0; j < n_obs; ++j) {
      if (j == i) {
        continue;
      }
      const std::size_t max_rank = std::max(ranks_in[j], ranks_out[j]);
      ++histogram[max_rank];
    }
  }
  return histogram;
}

double rnx_auc_from_max_rank_histogram(const std::vector<double> &histogram,
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

double rank_penalty_score(const NumericMatrix &din, const NumericMatrix &dout,
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

std::vector<std::size_t> validate_rank_penalty_ks(const IntegerVector &k,
                                                  std::size_t n_obs) {
  if (k.size() < 1) {
    stop("k must be nonempty");
  }

  std::vector<std::size_t> rank_ks;
  rank_ks.reserve(k.size());
  for (R_xlen_t i = 0; i < k.size(); ++i) {
    const int rank_k = k[i];
    if (rank_k == NA_INTEGER || rank_k < 1) {
      stop("k must contain positive integers");
    }
    const auto rank_k_size = static_cast<std::size_t>(rank_k);
    if ((2 * rank_k_size) >= n_obs) {
      stop("k must be less than half the number of observations");
    }
    rank_ks.push_back(rank_k_size);
  }
  return rank_ks;
}

NumericVector rank_penalty_scores(const NumericMatrix &din,
                                  const NumericMatrix &dout,
                                  const IntegerVector &k, bool continuity) {
  const std::size_t n_obs = validate_exact_distance_matrices(din, dout);
  const auto rank_ks = validate_rank_penalty_ks(k, n_obs);
  std::vector<double> penalties(rank_ks.size(), 0);

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
      for (std::size_t k_idx = 0; k_idx < rank_ks.size(); ++k_idx) {
        const std::size_t rank_k = rank_ks[k_idx];
        if (neighborhood_rank <= rank_k && penalty_rank > rank_k) {
          penalties[k_idx] += static_cast<double>(penalty_rank - rank_k);
        }
      }
    }
  }

  NumericVector scores(rank_ks.size());
  for (std::size_t k_idx = 0; k_idx < rank_ks.size(); ++k_idx) {
    const std::size_t rank_k = rank_ks[k_idx];
    const double normalization =
        static_cast<double>(n_obs) * static_cast<double>(rank_k) *
        static_cast<double>((2 * n_obs) - (3 * rank_k) - 1);
    scores[k_idx] = 1.0 - ((2.0 * penalties[k_idx]) / normalization);
  }
  return scores;
}

double exact_rank_penalty_metric(const NumericMatrix &din,
                                 const NumericMatrix &dout, int k,
                                 bool continuity) {
  const std::size_t n_obs = validate_exact_distance_matrices(din, dout);
  if (k < 1) {
    stop("k must be a positive integer");
  }

  const std::size_t rank_k = static_cast<std::size_t>(k);
  if ((2 * rank_k) >= n_obs) {
    stop("k must be less than half the number of observations");
  }

  return rank_penalty_score(din, dout, rank_k, continuity);
}

// [[Rcpp::export(rng = false)]]
double rnx_auc_direct(const NumericMatrix &din, const NumericMatrix &dout) {
  const std::size_t n_obs = validate_exact_distance_matrices(din, dout);
  const auto histogram = max_rank_histogram(din, dout, n_obs);
  return rnx_auc_from_max_rank_histogram(histogram, n_obs);
}

// [[Rcpp::export(rng = false)]]
NumericVector rnx_curve_direct(const NumericMatrix &din,
                               const NumericMatrix &dout,
                               const IntegerVector &k) {
  const std::size_t n_obs = validate_exact_distance_matrices(din, dout);
  const std::size_t n_ranks = n_obs - 1;
  if (k.size() < 1) {
    stop("k must be nonempty");
  }

  std::vector<std::size_t> rank_ks;
  rank_ks.reserve(k.size());
  for (R_xlen_t i = 0; i < k.size(); ++i) {
    const int rank_k = k[i];
    if (rank_k == NA_INTEGER || rank_k < 1) {
      stop("k must contain positive integers");
    }
    const auto rank_k_size = static_cast<std::size_t>(rank_k);
    if (rank_k_size >= n_ranks) {
      stop("k must be less than the number of non-self observations");
    }
    rank_ks.push_back(rank_k_size);
  }

  const auto histogram = max_rank_histogram(din, dout, n_obs);
  std::vector<double> full_curve(n_ranks, 0);
  double top_left = 0;
  for (std::size_t rank_k = 1; rank_k < n_ranks; ++rank_k) {
    top_left += histogram[rank_k];
    const double qnx = top_left / (static_cast<double>(rank_k) * n_obs);
    full_curve[rank_k] =
        ((qnx * n_ranks) - rank_k) / static_cast<double>(n_ranks - rank_k);
  }

  NumericVector curve(rank_ks.size());
  for (std::size_t i = 0; i < rank_ks.size(); ++i) {
    curve[i] = full_curve[rank_ks[i]];
  }
  return curve;
}

// [[Rcpp::export(rng = false)]]
double trustworthiness_exact(const NumericMatrix &din,
                             const NumericMatrix &dout, int k) {
  return exact_rank_penalty_metric(din, dout, k, false);
}

// [[Rcpp::export(rng = false)]]
double continuity_exact(const NumericMatrix &din, const NumericMatrix &dout,
                        int k) {
  return exact_rank_penalty_metric(din, dout, k, true);
}

// [[Rcpp::export(rng = false)]]
NumericVector trustworthiness_exact_multi(const NumericMatrix &din,
                                          const NumericMatrix &dout,
                                          const IntegerVector &k) {
  return rank_penalty_scores(din, dout, k, false);
}

// [[Rcpp::export(rng = false)]]
NumericVector continuity_exact_multi(const NumericMatrix &din,
                                     const NumericMatrix &dout,
                                     const IntegerVector &k) {
  return rank_penalty_scores(din, dout, k, true);
}
