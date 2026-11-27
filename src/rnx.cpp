#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
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

template <typename DistanceAt>
std::vector<std::size_t> first_tie_ranks(std::size_t n_obs, std::size_t row,
                                         DistanceAt distance_at) {
  std::vector<std::size_t> order;
  order.reserve(n_obs - 1);
  for (std::size_t col = 0; col < n_obs; ++col) {
    if (col != row) {
      order.push_back(col);
    }
  }

  std::stable_sort(order.begin(), order.end(),
                   [&](std::size_t lhs, std::size_t rhs) {
                     return rank_value_less(distance_at(lhs), distance_at(rhs));
                   });

  std::vector<std::size_t> ranks(n_obs, 0);
  for (std::size_t rank = 0; rank < order.size(); ++rank) {
    ranks[order[rank]] = rank + 1;
  }
  return ranks;
}

std::vector<std::size_t> row_ranks_first_ties(const NumericMatrix &distances,
                                              std::size_t row) {
  return first_tie_ranks(distances.nrow(), row,
                         [&](std::size_t col) { return distances(row, col); });
}

void accumulate_max_ranks(const std::vector<std::size_t> &ranks_in,
                          const std::vector<std::size_t> &ranks_out,
                          std::size_t row, std::vector<double> &histogram) {
  for (std::size_t j = 0; j < ranks_in.size(); ++j) {
    if (j != row) {
      ++histogram[std::max(ranks_in[j], ranks_out[j])];
    }
  }
}

void add_max_ranks(const NumericMatrix &din, const NumericMatrix &dout,
                   std::size_t row, std::vector<double> &histogram) {
  accumulate_max_ranks(row_ranks_first_ties(din, row),
                       row_ranks_first_ties(dout, row), row, histogram);
}

std::vector<double> max_rank_histogram(const NumericMatrix &din,
                                       const NumericMatrix &dout,
                                       std::size_t n_obs) {
  const std::size_t n_ranks = n_obs - 1;
  std::vector<double> histogram(n_ranks + 1, 0);

  for (std::size_t i = 0; i < n_obs; ++i) {
    if (i % 16 == 0) {
      Rcpp::checkUserInterrupt();
    }
    add_max_ranks(din, dout, i, histogram);
  }
  return histogram;
}

double rnx_at_k(double overlap_count, std::size_t n_queries,
                std::size_t n_ranks, std::size_t k) {
  const double qnx = overlap_count / (static_cast<double>(k) * n_queries);
  return ((qnx * n_ranks) - k) / static_cast<double>(n_ranks - k);
}

double rnx_auc_from_max_rank_histogram(const std::vector<double> &histogram,
                                       std::size_t n_queries,
                                       std::vector<double> *counts = nullptr) {
  const std::size_t n_ranks = histogram.size() - 1;
  double top_left = 0;
  double num = 0;
  double den = 0;
  for (std::size_t k = 1; k < n_ranks; ++k) {
    top_left += histogram[k];
    const double rnx = rnx_at_k(top_left, n_queries, n_ranks, k);
    if (counts != nullptr) {
      (*counts)[k] = top_left;
    }
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
    if (i % 16 == 0) {
      Rcpp::checkUserInterrupt();
    }
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
    if (i % 16 == 0) {
      Rcpp::checkUserInterrupt();
    }
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

std::vector<std::size_t> validate_rnx_curve_ks(const IntegerVector &k,
                                               std::size_t n_ranks) {
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
  return rank_ks;
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
  const auto rank_ks = validate_rnx_curve_ks(k, n_ranks);

  const auto histogram = max_rank_histogram(din, dout, n_obs);
  std::vector<double> full_curve(n_ranks, 0);
  double top_left = 0;
  for (std::size_t rank_k = 1; rank_k < n_ranks; ++rank_k) {
    top_left += histogram[rank_k];
    full_curve[rank_k] = rnx_at_k(top_left, n_obs, n_ranks, rank_k);
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

std::vector<double>
penalty_normalizers(std::size_t n_obs,
                    const std::vector<std::size_t> &rank_ks) {
  std::vector<double> normalizers;
  normalizers.reserve(rank_ks.size());
  for (const auto rank_k : rank_ks) {
    normalizers.push_back(static_cast<double>(rank_k) *
                          static_cast<double>((2 * n_obs) - (3 * rank_k) - 1));
  }

  return normalizers;
}

void set_row_penalties(const std::vector<std::size_t> &ranks_in,
                       const std::vector<std::size_t> &ranks_out,
                       const std::vector<std::size_t> &rank_ks,
                       const std::vector<double> &normalizers,
                       std::size_t result_row, NumericMatrix &trustworthiness,
                       NumericMatrix &continuity) {
  for (std::size_t col = 0; col < rank_ks.size(); ++col) {
    const auto rank_k = rank_ks[col];
    double false_penalty = 0;
    double missing_penalty = 0;
    for (std::size_t j = 0; j < ranks_in.size(); ++j) {
      // Self has rank zero in both spaces and incurs neither penalty.
      const auto ri = ranks_in[j];
      const auto ro = ranks_out[j];
      if (ro <= rank_k && ri > rank_k) {
        false_penalty += static_cast<double>(ri - rank_k);
      } else if (ri <= rank_k && ro > rank_k) {
        missing_penalty += static_cast<double>(ro - rank_k);
      }
    }
    trustworthiness(result_row, col) =
        1.0 - 2.0 * false_penalty / normalizers[col];
    continuity(result_row, col) =
        1.0 - 2.0 * missing_penalty / normalizers[col];
  }
}

// [[Rcpp::export(rng = false)]]
List pointwise_rank_preservation_exact(const NumericMatrix &din,
                                       const NumericMatrix &dout,
                                       const IntegerVector &k) {
  const auto n_obs = validate_exact_distance_matrices(din, dout);
  const auto rank_ks = validate_rank_penalty_ks(k, n_obs);
  NumericMatrix trustworthiness(n_obs, rank_ks.size());
  NumericMatrix continuity(n_obs, rank_ks.size());
  const auto normalizers = penalty_normalizers(n_obs, rank_ks);

  for (std::size_t i = 0; i < n_obs; ++i) {
    if (i % 16 == 0) {
      Rcpp::checkUserInterrupt();
    }
    const auto ranks_in = row_ranks_first_ties(din, i);
    const auto ranks_out = row_ranks_first_ties(dout, i);
    set_row_penalties(ranks_in, ranks_out, rank_ks, normalizers, i,
                      trustworthiness, continuity);
  }

  return List::create(Named("trustworthiness") = trustworthiness,
                      Named("continuity") = continuity);
}

// [[Rcpp::export(rng = false)]]
List pointwise_rnx_exact(const NumericMatrix &din, const NumericMatrix &dout,
                         const IntegerVector &k) {
  const auto n_obs = validate_exact_distance_matrices(din, dout);
  const std::size_t n_ranks = n_obs - 1;
  const auto rank_ks = validate_rnx_curve_ks(k, n_ranks);
  NumericVector auc(n_obs);
  NumericMatrix curve(n_obs, rank_ks.size());
  std::vector<double> histogram(n_ranks + 1, 0);
  std::vector<double> row_counts(rank_ks.empty() ? 0 : n_ranks);

  for (std::size_t i = 0; i < n_obs; ++i) {
    if (i % 16 == 0) {
      Rcpp::checkUserInterrupt();
    }
    std::fill(histogram.begin(), histogram.end(), 0);
    add_max_ranks(din, dout, i, histogram);
    auc[i] = rnx_auc_from_max_rank_histogram(
        histogram, 1, rank_ks.empty() ? nullptr : &row_counts);
    for (std::size_t col = 0; col < rank_ks.size(); ++col) {
      curve(i, col) =
          rnx_at_k(row_counts[rank_ks[col]], 1, n_ranks, rank_ks[col]);
    }
  }
  return List::create(Named("auc") = auc, Named("curve") = curve);
}

// Form Euclidean distances in coordinate order, as in stats::dist(). Updating
// all candidates together keeps the inner loop contiguous in R's columns.
template <typename Vector>
void fill_euclidean_distances(const NumericMatrix &x, std::size_t row,
                              Vector &distances) {
  const auto n_obs = static_cast<std::size_t>(x.nrow());
  std::fill(distances.begin(), distances.end(), 0.0);
  for (int dim = 0; dim < x.ncol(); ++dim) {
    const double *column = x.begin() + static_cast<std::size_t>(dim) * n_obs;
    const double anchor = column[row];
    for (std::size_t j = 0; j < n_obs; ++j) {
      const double delta = anchor - column[j];
      distances[j] += delta * delta;
    }
  }
  for (auto &value : distances) {
    if (!R_FINITE(value)) {
      stop("Nonfinite computed Euclidean distance");
    }
    value = std::sqrt(value);
  }
}

// [[Rcpp::export(rng = false)]]
NumericVector euclidean_distances_from_row(const NumericMatrix &x, int row) {
  if (row < 1 || row > x.nrow()) {
    stop("row must index a coordinate row");
  }
  NumericVector distances(x.nrow());
  fill_euclidean_distances(x, row - 1, distances);
  return distances;
}

std::vector<std::size_t> euclidean_row_ranks(const NumericMatrix &x,
                                             std::size_t row,
                                             std::vector<double> &distances) {
  fill_euclidean_distances(x, row, distances);
  return first_tie_ranks(distances.size(), row,
                         [&](std::size_t j) { return distances[j]; });
}

// Each selected rank occurs exactly once, so every nonzero slot is
// overwritten before the prefix sum for this row.
void set_relative_rank_errors(const std::vector<std::size_t> &selected,
                              const std::vector<std::size_t> &other,
                              const std::vector<std::size_t> &sizes,
                              const std::vector<double> &normalizers,
                              std::vector<double> &errors, int row,
                              NumericMatrix &result) {
  for (std::size_t j = 0; j < selected.size(); ++j) {
    const auto rank = selected[j];
    if (rank > 0 && rank < errors.size()) {
      errors[rank] = std::abs(static_cast<double>(other[j]) - rank) / rank;
    }
  }
  std::partial_sum(errors.begin(), errors.end(), errors.begin());
  for (std::size_t col = 0; col < sizes.size(); ++col) {
    const auto k = sizes[col];
    result(row, col) = errors[k] / normalizers[k];
  }
}

struct CoordinateRankProfile {
  NumericMatrix trustworthiness;
  NumericMatrix continuity;
  NumericVector auc;
  NumericMatrix curve;
  IntegerMatrix retained;
  NumericMatrix mrre_input;
  NumericMatrix mrre_output;
  CoordinateRankProfile(int n, int n_sizes, bool mrre)
      : trustworthiness(n, n_sizes), continuity(n, n_sizes), auc(n),
        curve(n, n_sizes), retained(n, n_sizes),
        mrre_input(n, mrre ? n_sizes : 0), mrre_output(n, mrre ? n_sizes : 0) {}
};

// [[Rcpp::export(rng = false)]]
List coordinate_rank_profiles(const NumericMatrix &input,
                              const List &embeddings, const IntegerVector &k,
                              bool mrre, const IntegerVector &query_rows) {
  const int n = input.nrow();
  if (n < 3 || input.ncol() < 1 || embeddings.size() < 1) {
    stop("Expected reference coordinates and at least one embedding");
  }
  const auto sizes = validate_rank_penalty_ks(k, n);
  const auto normalizers = penalty_normalizers(n, sizes);
  const auto max_k = mrre ? *std::max_element(sizes.begin(), sizes.end()) : 0;
  std::vector<double> errors(max_k + 1), mrre_normalizers(max_k + 1);
  for (std::size_t r = 1; r <= max_k; ++r) {
    mrre_normalizers[r] = mrre_normalizers[r - 1] + (n - 2.0 * r) / r;
  }
  std::vector<NumericMatrix> maps;
  std::vector<CoordinateRankProfile> profiles;
  for (R_xlen_t m = 0; m < embeddings.size(); ++m) {
    NumericMatrix map = as<NumericMatrix>(embeddings[m]);
    if (map.nrow() != n || map.ncol() < 1) {
      stop("Each embedding must have the same rows as the reference");
    }
    maps.push_back(map);
    profiles.emplace_back(query_rows.size(), sizes.size(), mrre);
  }
  std::vector<double> distances(n), histogram(n), row_counts(n);
  for (R_xlen_t q = 0; q < query_rows.size(); ++q) {
    if (query_rows[q] < 1 || query_rows[q] > n) {
      stop("query_rows contains a position outside the reference");
    }
    const int anchor = query_rows[q] - 1;
    checkUserInterrupt();
    const auto reference = euclidean_row_ranks(input, anchor, distances);
    for (std::size_t m = 0; m < maps.size(); ++m) {
      const auto embedded = euclidean_row_ranks(maps[m], anchor, distances);
      auto &p = profiles[m];
      set_row_penalties(reference, embedded, sizes, normalizers, q,
                        p.trustworthiness, p.continuity);
      if (mrre) {
        set_relative_rank_errors(reference, embedded, sizes, mrre_normalizers,
                                 errors, q, p.mrre_input);
        set_relative_rank_errors(embedded, reference, sizes, mrre_normalizers,
                                 errors, q, p.mrre_output);
      }
      std::fill(histogram.begin(), histogram.end(), 0.0);
      accumulate_max_ranks(reference, embedded, anchor, histogram);
      p.auc[q] = rnx_auc_from_max_rank_histogram(histogram, 1, &row_counts);
      for (std::size_t col = 0; col < sizes.size(); ++col) {
        const auto size = sizes[col];
        p.retained(q, col) = static_cast<int>(row_counts[size]);
        p.curve(q, col) = rnx_at_k(row_counts[size], 1, n - 1, size);
      }
    }
  }
  List result(embeddings.size());
  for (std::size_t m = 0; m < profiles.size(); ++m) {
    const auto &p = profiles[m];
    List profile = List::create(Named("trustworthiness") = p.trustworthiness,
                                Named("continuity") = p.continuity,
                                Named("auc") = p.auc, Named("curve") = p.curve);
    if (mrre) {
      profile["mrre_input"] = p.mrre_input;
      profile["mrre_output"] = p.mrre_output;
    }
    profile["retained"] = p.retained;
    result[m] = profile;
  }
  return result;
}
