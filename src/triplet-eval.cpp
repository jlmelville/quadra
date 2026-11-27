#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include <Rcpp.h>

#include "rnndescent/random.h"

#include "distance.h"
#include "native-validation.h"
#include "pforr.h"

using namespace Rcpp;

int validate_observation_matrices(const NumericMatrix &xin,
                                  const NumericMatrix &xout) {
  if (xin.ncol() != xout.ncol()) {
    stop("xin and xout must have the same number of observations");
  }
  if (xin.nrow() < 1 || xout.nrow() < 1 || xin.ncol() < 3) {
    stop("xin and xout must each have at least one feature and three "
         "observations");
  }
  return xin.ncol();
}

void validate_triplet_layout(const NumericMatrix &triplets, int n_obs) {
  if (triplets.ncol() != n_obs) {
    stop("triplets must have one column per observation");
  }
  if ((triplets.nrow() % 2) != 0) {
    stop("triplets must have an even number of rows");
  }
}

struct TripletCounts {
  std::size_t agreements{0};
  std::size_t comparisons{0};
};

struct AnchorTripletCounts {
  TripletCounts defined;
  std::size_t sampled{0};
};

struct TripletDetails {
  explicit TripletDetails(std::size_t n)
      : anchor(n), endpoint1(n), endpoint2(n), agreement(n) {}
  std::vector<int> anchor, endpoint1, endpoint2, agreement;
};

int compare_distances(double lhs, double rhs) {
  return (lhs > rhs) - (lhs < rhs);
}

std::size_t n_parallel_chunks(std::size_t begin, std::size_t end,
                              std::size_t n_threads) {
  const auto ranges =
      pforr::split_input_range(pforr::IndexRange(begin, end), n_threads, 1);
  return std::max<std::size_t>(ranges.size(), 1);
}

int update_triplet_counts(TripletCounts &counts, std::size_t anchor,
                          std::size_t p1, std::size_t p2,
                          const DistanceSpace &input,
                          const DistanceSpace &output) {
  const int input_order =
      compare_distances(input.between(anchor, p1), input.between(anchor, p2));
  if (input_order == 0) {
    return NA_LOGICAL;
  }

  const int output_order =
      compare_distances(output.between(anchor, p1), output.between(anchor, p2));
  ++counts.comparisons;
  if (input_order == output_order) {
    ++counts.agreements;
    return TRUE;
  }
  return FALSE;
}

double summarize_triplet_counts(const std::vector<TripletCounts> &counts) {
  TripletCounts total_counts;
  for (const auto &chunk_counts : counts) {
    total_counts.agreements += chunk_counts.agreements;
    total_counts.comparisons += chunk_counts.comparisons;
  }
  if (total_counts.comparisons == 0) {
    return NA_REAL;
  }
  return static_cast<double>(total_counts.agreements) /
         static_cast<double>(total_counts.comparisons);
}

double triplet_plan_sample(const std::vector<std::size_t> &anchor,
                           const std::vector<std::size_t> &endpoint1,
                           const std::vector<std::size_t> &endpoint2,
                           const DistanceSpace &input,
                           const DistanceSpace &output, std::size_t n_threads,
                           std::vector<int> *agreement) {

  std::vector<TripletCounts> counts(
      n_parallel_chunks(0, anchor.size(), n_threads));

  auto worker = [&](std::size_t begin, std::size_t end, std::size_t chunk_id) {
    TripletCounts chunk;
    for (std::size_t i = begin; i < end; ++i) {
      const int outcome = update_triplet_counts(chunk, anchor[i], endpoint1[i],
                                                endpoint2[i], input, output);
      if (agreement != nullptr) {
        (*agreement)[i] = outcome;
      }
    }
    counts[chunk_id].agreements += chunk.agreements;
    counts[chunk_id].comparisons += chunk.comparisons;
  };

  pforr::parallel_for_indexed(anchor.size(), worker, n_threads);
  return summarize_triplet_counts(counts);
}

std::size_t avoid_anchor_index(std::size_t idx, std::size_t anchor) {
  return idx >= anchor ? idx + 1 : idx;
}

TripletCounts random_triplet_sample_inner(
    std::size_t begin, std::size_t end, std::size_t n_triplets_per_observation,
    int n_obs, tdoann::RandomIntGenerator<uint64_t> &int_sampler,
    const DistanceSpace &input, const DistanceSpace &output,
    TripletDetails *details, std::vector<AnchorTripletCounts> &anchor_counts) {

  TripletCounts counts;
  for (std::size_t i = begin; i < end; i++) {
    TripletCounts &row_counts =
        anchor_counts.empty() ? counts : anchor_counts[i].defined;
    for (std::size_t j = 0; j < n_triplets_per_observation; j++) {
      const auto idxs = int_sampler.sample(n_obs - 1, 2);
      const std::size_t p1 = avoid_anchor_index(idxs[0], i);
      const std::size_t p2 = avoid_anchor_index(idxs[1], i);

      const int outcome =
          update_triplet_counts(row_counts, i, p1, p2, input, output);
      if (details != nullptr) {
        const std::size_t row = i * n_triplets_per_observation + j;
        details->anchor[row] = static_cast<int>(i + 1);
        details->endpoint1[row] = static_cast<int>(p1 + 1);
        details->endpoint2[row] = static_cast<int>(p2 + 1);
        details->agreement[row] = outcome;
      }
    }
    if (!anchor_counts.empty()) {
      anchor_counts[i].sampled = n_triplets_per_observation;
      counts.agreements += row_counts.agreements;
      counts.comparisons += row_counts.comparisons;
    }
  }
  return counts;
}

double random_triplet_sample(std::size_t n_triplets_per_observation, int n_obs,
                             const DistanceSpace &input,
                             const DistanceSpace &output, std::size_t n_threads,
                             TripletDetails *details,
                             std::vector<AnchorTripletCounts> &anchor_counts) {

  rnndescent::ParallelIntRNGAdapter<uint64_t, rnndescent::DQIntSampler>
      sampler_provider;
  sampler_provider.initialize();

  std::vector<TripletCounts> counts(n_parallel_chunks(0, n_obs, n_threads));

  auto worker = [&](std::size_t begin, std::size_t end, std::size_t chunk_id) {
    auto thread_sampler = sampler_provider.get_parallel_instance(chunk_id);
    const auto chunk_counts = random_triplet_sample_inner(
        begin, end, n_triplets_per_observation, n_obs, *thread_sampler, input,
        output, details, anchor_counts);
    counts[chunk_id].agreements += chunk_counts.agreements;
    counts[chunk_id].comparisons += chunk_counts.comparisons;
  };

  pforr::parallel_for_indexed(0, n_obs, worker, n_threads);
  return summarize_triplet_counts(counts);
}

DataFrame
anchor_triplet_result(const std::vector<AnchorTripletCounts> &anchor_counts) {
  // One row per observation in the input R matrix.
  const auto n = static_cast<int>(anchor_counts.size());
  IntegerVector anchor(n), sampled(n), defined(n), agreed(n);
  NumericVector accuracy(n);
  for (int i = 0; i < n; ++i) {
    const auto &row = anchor_counts[i];
    anchor[i] = i + 1;
    sampled[i] = static_cast<int>(row.sampled);
    defined[i] = static_cast<int>(row.defined.comparisons);
    agreed[i] = static_cast<int>(row.defined.agreements);
    accuracy[i] = row.defined.comparisons == 0
                      ? NA_REAL
                      : static_cast<double>(row.defined.agreements) /
                            static_cast<double>(row.defined.comparisons);
  }
  return DataFrame::create(_["anchor"] = anchor, _["n_sampled"] = sampled,
                           _["n_defined"] = defined, _["n_agree"] = agreed,
                           _["accuracy"] = accuracy);
}

SEXP triplet_result(double accuracy, const TripletDetails *details,
                    const std::vector<AnchorTripletCounts> &anchor_counts) {
  if (details == nullptr && anchor_counts.empty()) {
    return Rcpp::wrap(accuracy);
  }

  List result = List::create(_["accuracy"] = accuracy);
  if (details != nullptr) {
    const int n_triplets = static_cast<int>(details->anchor.size());
    IntegerMatrix triplets(n_triplets, 3);
    LogicalVector agreement_result(n_triplets);
    for (int i = 0; i < n_triplets; i++) {
      triplets(i, 0) = details->anchor[i];
      triplets(i, 1) = details->endpoint1[i];
      triplets(i, 2) = details->endpoint2[i];
      agreement_result[i] = details->agreement[i];
    }
    colnames(triplets) =
        CharacterVector::create("anchor", "endpoint1", "endpoint2");
    result["triplets"] = triplets;
    result["agreement"] = agreement_result;
  }
  if (!anchor_counts.empty()) {
    result["by_anchor"] = anchor_triplet_result(anchor_counts);
  }
  return result;
}

// [[Rcpp::export(rng = false)]]
double triplet_sample(const NumericMatrix &triplets, const NumericMatrix &xin,
                      const NumericMatrix &xout,
                      const std::string &metric_in = "sqeuclidean",
                      const std::string &metric_out = "sqeuclidean",
                      double n_threads = 0) {

  const int n_obs = validate_observation_matrices(xin, xout);
  validate_triplet_layout(triplets, n_obs);
  const std::size_t thread_count = quadra::checked_thread_count(n_threads);
  const DistanceSpace input(xin, metric_in), output(xout, metric_out);
  std::vector<std::size_t> endpoints;
  endpoints.reserve(triplets.size());
  for (const double endpoint : triplets) {
    if (!(endpoint >= 0 && endpoint < n_obs &&
          endpoint == std::floor(endpoint))) {
      stop("triplets must contain whole, zero-based observation indices");
    }
    endpoints.push_back(static_cast<std::size_t>(endpoint));
  }
  const std::size_t endpoints_per_observation = triplets.nrow();
  std::vector<TripletCounts> counts(n_parallel_chunks(0, n_obs, thread_count));

  auto worker = [&](std::size_t begin, std::size_t end, std::size_t chunk_id) {
    TripletCounts chunk;
    for (std::size_t i = begin; i < end; ++i) {
      const std::size_t offset = i * endpoints_per_observation;
      for (std::size_t j = 0; j < endpoints_per_observation; j += 2) {
        update_triplet_counts(chunk, i, endpoints[offset + j],
                              endpoints[offset + j + 1], input, output);
      }
    }
    counts[chunk_id].agreements += chunk.agreements;
    counts[chunk_id].comparisons += chunk.comparisons;
  };

  pforr::parallel_for_indexed(0, n_obs, worker, thread_count);
  return summarize_triplet_counts(counts);
}

// [[Rcpp::export(rng = false)]]
SEXP triplet_plan_sample(const IntegerMatrix &triplets,
                         const NumericMatrix &xin, const NumericMatrix &xout,
                         const std::string &metric_in = "sqeuclidean",
                         const std::string &metric_out = "sqeuclidean",
                         double n_threads = 0, bool ret_triplets = false,
                         bool ret_anchors = false) {

  const int n_obs = validate_observation_matrices(xin, xout);
  if (triplets.nrow() < 1 || triplets.ncol() != 3) {
    stop("triplets must have at least one row and exactly three columns");
  }
  const std::size_t thread_count = quadra::checked_thread_count(n_threads);
  const std::size_t triplet_count = static_cast<std::size_t>(triplets.nrow());

  std::vector<std::size_t> anchor(triplet_count);
  std::vector<std::size_t> endpoint1(triplet_count);
  std::vector<std::size_t> endpoint2(triplet_count);
  TripletDetails details(ret_triplets ? triplet_count : 0);
  if (ret_anchors && !ret_triplets) {
    details.agreement.resize(triplet_count);
  }

  for (std::size_t i = 0; i < triplet_count; i++) {
    const int anchor_i = triplets(i, 0);
    const int endpoint1_i = triplets(i, 1);
    const int endpoint2_i = triplets(i, 2);
    if (anchor_i == NA_INTEGER || endpoint1_i == NA_INTEGER ||
        endpoint2_i == NA_INTEGER || anchor_i < 1 || endpoint1_i < 1 ||
        endpoint2_i < 1 || anchor_i > xin.ncol() || endpoint1_i > xin.ncol() ||
        endpoint2_i > xin.ncol()) {
      stop("triplets must contain indices between 1 and the number of "
           "observations");
    }
    if (anchor_i == endpoint1_i || anchor_i == endpoint2_i ||
        endpoint1_i == endpoint2_i) {
      stop("triplets must contain three distinct indices per row");
    }
    anchor[i] = static_cast<std::size_t>(anchor_i - 1);
    endpoint1[i] = static_cast<std::size_t>(endpoint1_i - 1);
    endpoint2[i] = static_cast<std::size_t>(endpoint2_i - 1);
    if (ret_triplets) {
      details.anchor[i] = anchor_i;
      details.endpoint1[i] = endpoint1_i;
      details.endpoint2[i] = endpoint2_i;
    }
  }

  const DistanceSpace input(xin, metric_in), output(xout, metric_out);

  const double accuracy = triplet_plan_sample(
      anchor, endpoint1, endpoint2, input, output, thread_count,
      details.agreement.empty() ? nullptr : &details.agreement);

  std::vector<AnchorTripletCounts> anchor_counts(ret_anchors ? n_obs : 0);
  if (ret_anchors) {
    for (std::size_t i = 0; i < triplet_count; ++i) {
      auto &row = anchor_counts[anchor[i]];
      ++row.sampled;
      row.defined.comparisons += details.agreement[i] != NA_LOGICAL;
      row.defined.agreements += details.agreement[i] == TRUE;
    }
  }

  return triplet_result(accuracy, ret_triplets ? &details : nullptr,
                        anchor_counts);
}

// [[Rcpp::export(rng = false)]]
SEXP random_triplet_sample(const NumericMatrix &xin, const NumericMatrix &xout,
                           double n_triplets = 5,
                           const std::string &metric_in = "sqeuclidean",
                           const std::string &metric_out = "sqeuclidean",
                           double n_threads = 0, bool ret_triplets = false,
                           bool ret_anchors = false) {

  const int n_obs = validate_observation_matrices(xin, xout);
  const std::size_t triplet_count =
      quadra::checked_nonnegative_count(n_triplets, "n_triplets");
  const std::size_t thread_count = quadra::checked_thread_count(n_threads);

  if (ret_triplets &&
      triplet_count >
          static_cast<std::size_t>((std::numeric_limits<int>::max)()) / n_obs) {
    stop("detailed triplet sample exceeds the supported integer range");
  }
  const std::size_t comparison_count = n_obs * triplet_count;
  TripletDetails details(ret_triplets ? comparison_count : 0);

  const DistanceSpace input(xin, metric_in), output(xout, metric_out);

  std::vector<AnchorTripletCounts> anchor_counts(ret_anchors ? n_obs : 0);
  double accuracy;
  {
    Rcpp::RNGScope rng_scope;
    accuracy =
        random_triplet_sample(triplet_count, n_obs, input, output, thread_count,
                              ret_triplets ? &details : nullptr, anchor_counts);
  }

  return triplet_result(accuracy, ret_triplets ? &details : nullptr,
                        anchor_counts);
}
