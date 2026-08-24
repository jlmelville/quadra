#include <algorithm>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <Rcpp.h>

#include "rnndescent/random.h"
#include "tdoann/distance.h"

#include "distance.h"
#include "native-validation.h"
#include "pforr.h"

using namespace Rcpp;
using It = typename std::vector<double>::const_iterator;
using Dfun = double(It, It, It);
using TripIt = typename std::vector<std::size_t>::const_iterator;

std::size_t validate_observation_matrices(const NumericMatrix &xin,
                                          const NumericMatrix &xout) {
  if (xin.ncol() != xout.ncol()) {
    stop("xin and xout must have the same number of observations");
  }
  if (xin.nrow() < 1 || xout.nrow() < 1 || xin.ncol() < 3) {
    stop("xin and xout must each have at least one feature and three "
         "observations");
  }
  return static_cast<std::size_t>(xin.ncol());
}

void validate_triplet_layout(const NumericMatrix &triplets, std::size_t n_obs) {
  if (triplets.ncol() != static_cast<int>(n_obs)) {
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

int compare_distances(double lhs, double rhs) {
  return (lhs > rhs) - (lhs < rhs);
}

std::size_t n_parallel_chunks(std::size_t begin, std::size_t end,
                              std::size_t n_threads) {
  const auto ranges =
      pforr::split_input_range(pforr::IndexRange(begin, end), n_threads, 1);
  return std::max<std::size_t>(ranges.size(), 1);
}

int update_triplet_counts(TripletCounts &counts, std::size_t p1, std::size_t p2,
                          const It xin_i_begin, const It xin_i_end,
                          const It xin_begin, std::size_t xin_ncol,
                          const std::function<Dfun> &dfunin,
                          const It xout_i_begin, const It xout_i_end,
                          const It xout_begin, std::size_t xout_ncol,
                          const std::function<Dfun> &dfunout) {

  const auto input_distance_1 =
      dfunin(xin_i_begin, xin_i_end, xin_begin + p1 * xin_ncol);
  const auto input_distance_2 =
      dfunin(xin_i_begin, xin_i_end, xin_begin + p2 * xin_ncol);

  const int input_order = compare_distances(input_distance_1, input_distance_2);
  if (input_order == 0) {
    return NA_LOGICAL;
  }

  const auto output_distance_1 =
      dfunout(xout_i_begin, xout_i_end, xout_begin + p1 * xout_ncol);
  const auto output_distance_2 =
      dfunout(xout_i_begin, xout_i_end, xout_begin + p2 * xout_ncol);

  ++counts.comparisons;
  const int output_order =
      compare_distances(output_distance_1, output_distance_2);
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
  return total_counts.agreements /
         static_cast<double>(total_counts.comparisons);
}

TripletCounts legacy_triplet_sample_inner(
    std::size_t begin, std::size_t end, std::size_t n_triplets_per_observation,
    const TripIt triplets_begin, const It xin_begin, std::size_t xin_ncol,
    const std::function<Dfun> &dfunin, const It xout_begin,
    std::size_t xout_ncol, const std::function<Dfun> &dfunout) {

  TripletCounts counts;
  const std::size_t endpoints_per_observation = n_triplets_per_observation * 2;
  for (std::size_t i = begin; i < end; i++) {
    const TripIt observation_triplets =
        triplets_begin + i * endpoints_per_observation;
    const It xin_i_begin = xin_begin + i * xin_ncol;
    const It xin_i_end = xin_i_begin + xin_ncol;
    const It xout_i_begin = xout_begin + i * xout_ncol;
    const It xout_i_end = xout_i_begin + xout_ncol;

    for (std::size_t j = 0; j < n_triplets_per_observation; j++) {
      const TripIt endpoints = observation_triplets + j * 2;
      update_triplet_counts(counts, *endpoints, *(endpoints + 1), xin_i_begin,
                            xin_i_end, xin_begin, xin_ncol, dfunin,
                            xout_i_begin, xout_i_end, xout_begin, xout_ncol,
                            dfunout);
    }
  }
  return counts;
}

double legacy_triplet_sample(TripIt triplets_begin, TripIt triplets_end,
                             std::size_t n_obs, It xin_begin, It xin_end,
                             It xout_begin, It xout_end,
                             const std::function<Dfun> &dfunin,
                             const std::function<Dfun> &dfunout,
                             std::size_t n_threads) {

  const std::size_t n_triplets_per_observation =
      (triplets_end - triplets_begin) / n_obs / 2;
  const std::size_t xin_nfeat = (xin_end - xin_begin) / n_obs;
  const std::size_t xout_nfeat = (xout_end - xout_begin) / n_obs;
  std::vector<TripletCounts> counts(n_parallel_chunks(0, n_obs, n_threads));

  auto worker = [&](std::size_t begin, std::size_t end, std::size_t chunk_id) {
    const auto chunk_counts = legacy_triplet_sample_inner(
        begin, end, n_triplets_per_observation, triplets_begin, xin_begin,
        xin_nfeat, dfunin, xout_begin, xout_nfeat, dfunout);
    counts[chunk_id].agreements += chunk_counts.agreements;
    counts[chunk_id].comparisons += chunk_counts.comparisons;
  };

  pforr::parallel_for_indexed(0, n_obs, worker, n_threads);
  return summarize_triplet_counts(counts);
}

TripletCounts triplet_plan_sample_inner(
    std::size_t begin, std::size_t end, const std::vector<std::size_t> &anchor,
    const std::vector<std::size_t> &endpoint1,
    const std::vector<std::size_t> &endpoint2, const It xin_begin,
    std::size_t xin_ncol, const std::function<Dfun> &dfunin,
    const It xout_begin, std::size_t xout_ncol,
    const std::function<Dfun> &dfunout, std::vector<int> *agreement) {

  TripletCounts counts;
  for (std::size_t i = begin; i < end; i++) {
    const It xin_i_begin = xin_begin + anchor[i] * xin_ncol;
    const It xin_i_end = xin_i_begin + xin_ncol;
    const It xout_i_begin = xout_begin + anchor[i] * xout_ncol;
    const It xout_i_end = xout_i_begin + xout_ncol;

    const int outcome = update_triplet_counts(
        counts, endpoint1[i], endpoint2[i], xin_i_begin, xin_i_end, xin_begin,
        xin_ncol, dfunin, xout_i_begin, xout_i_end, xout_begin, xout_ncol,
        dfunout);
    if (agreement != nullptr) {
      (*agreement)[i] = outcome;
    }
  }
  return counts;
}

double triplet_plan_sample(const std::vector<std::size_t> &anchor,
                           const std::vector<std::size_t> &endpoint1,
                           const std::vector<std::size_t> &endpoint2,
                           std::size_t n_obs, It xin_begin, It xin_end,
                           It xout_begin, It xout_end,
                           const std::function<Dfun> &dfunin,
                           const std::function<Dfun> &dfunout,
                           std::size_t n_threads, std::vector<int> *agreement) {

  const std::size_t xin_nfeat = (xin_end - xin_begin) / n_obs;
  const std::size_t xout_nfeat = (xout_end - xout_begin) / n_obs;
  std::vector<TripletCounts> counts(
      n_parallel_chunks(0, anchor.size(), n_threads));

  auto worker = [&](std::size_t begin, std::size_t end, std::size_t chunk_id) {
    const auto chunk_counts = triplet_plan_sample_inner(
        begin, end, anchor, endpoint1, endpoint2, xin_begin, xin_nfeat, dfunin,
        xout_begin, xout_nfeat, dfunout, agreement);
    counts[chunk_id].agreements += chunk_counts.agreements;
    counts[chunk_id].comparisons += chunk_counts.comparisons;
  };

  pforr::parallel_for_indexed(anchor.size(), worker, n_threads);
  return summarize_triplet_counts(counts);
}

std::size_t avoid_anchor_index(std::size_t idx, std::size_t anchor) {
  return idx >= anchor ? idx + 1 : idx;
}

TripletCounts random_triplet_sample_inner(
    std::size_t begin, std::size_t end, std::size_t n_triplets_per_observation,
    std::size_t n_obs, tdoann::RandomIntGenerator<uint64_t> &int_sampler,
    const It xin_begin, std::size_t xin_ncol, const std::function<Dfun> &dfunin,
    const It xout_begin, std::size_t xout_ncol,
    const std::function<Dfun> &dfunout, std::vector<int> *returned_anchor,
    std::vector<int> *returned_endpoint1, std::vector<int> *returned_endpoint2,
    std::vector<int> *agreement) {

  TripletCounts counts;
  for (std::size_t i = begin; i < end; i++) {
    const It xin_i_begin = xin_begin + i * xin_ncol;
    const It xin_i_end = xin_i_begin + xin_ncol;
    const It xout_i_begin = xout_begin + i * xout_ncol;
    const It xout_i_end = xout_i_begin + xout_ncol;

    for (std::size_t j = 0; j < n_triplets_per_observation; j++) {
      const auto idxs = int_sampler.sample(n_obs - 1, 2);
      const std::size_t p1 = avoid_anchor_index(idxs[0], i);
      const std::size_t p2 = avoid_anchor_index(idxs[1], i);

      const int outcome = update_triplet_counts(
          counts, p1, p2, xin_i_begin, xin_i_end, xin_begin, xin_ncol, dfunin,
          xout_i_begin, xout_i_end, xout_begin, xout_ncol, dfunout);
      if (agreement != nullptr) {
        const std::size_t row = i * n_triplets_per_observation + j;
        (*returned_anchor)[row] = static_cast<int>(i + 1);
        (*returned_endpoint1)[row] = static_cast<int>(p1 + 1);
        (*returned_endpoint2)[row] = static_cast<int>(p2 + 1);
        (*agreement)[row] = outcome;
      }
    }
  }
  return counts;
}

double random_triplet_sample(
    std::size_t n_triplets_per_observation, std::size_t n_obs, It xin_begin,
    It xin_end, It xout_begin, It xout_end, const std::function<Dfun> &dfunin,
    const std::function<Dfun> &dfunout, std::size_t n_threads,
    std::vector<int> *returned_anchor, std::vector<int> *returned_endpoint1,
    std::vector<int> *returned_endpoint2, std::vector<int> *agreement) {
  const std::size_t xin_nfeat = (xin_end - xin_begin) / n_obs;
  const std::size_t xout_nfeat = (xout_end - xout_begin) / n_obs;

  rnndescent::ParallelIntRNGAdapter<uint64_t, rnndescent::DQIntSampler>
      sampler_provider;
  sampler_provider.initialize();

  std::vector<TripletCounts> counts(n_parallel_chunks(0, n_obs, n_threads));

  auto worker = [&](std::size_t begin, std::size_t end, std::size_t chunk_id) {
    auto thread_sampler = sampler_provider.get_parallel_instance(chunk_id);
    const auto chunk_counts = random_triplet_sample_inner(
        begin, end, n_triplets_per_observation, n_obs, *thread_sampler,
        xin_begin, xin_nfeat, dfunin, xout_begin, xout_nfeat, dfunout,
        returned_anchor, returned_endpoint1, returned_endpoint2, agreement);
    counts[chunk_id].agreements += chunk_counts.agreements;
    counts[chunk_id].comparisons += chunk_counts.comparisons;
  };

  pforr::parallel_for_indexed(0, n_obs, worker, n_threads);
  return summarize_triplet_counts(counts);
}

SEXP triplet_result(double accuracy, bool ret_triplets,
                    const std::vector<int> &anchor,
                    const std::vector<int> &endpoint1,
                    const std::vector<int> &endpoint2,
                    const std::vector<int> &agreement) {
  if (!ret_triplets) {
    return Rcpp::wrap(accuracy);
  }

  const int n_triplets = static_cast<int>(anchor.size());
  IntegerMatrix triplets(n_triplets, 3);
  LogicalVector agreement_result(n_triplets);
  for (int i = 0; i < n_triplets; i++) {
    triplets(i, 0) = anchor[i];
    triplets(i, 1) = endpoint1[i];
    triplets(i, 2) = endpoint2[i];
    agreement_result[i] = agreement[i];
  }
  colnames(triplets) =
      CharacterVector::create("anchor", "endpoint1", "endpoint2");

  return List::create(_["accuracy"] = accuracy, _["triplets"] = triplets,
                      _["agreement"] = agreement_result);
}

// [[Rcpp::export(rng = false)]]
double triplet_sample(const NumericMatrix &triplets, const NumericMatrix &xin,
                      const NumericMatrix &xout,
                      const std::string &metric_in = "sqeuclidean",
                      const std::string &metric_out = "sqeuclidean",
                      double n_threads = 0) {

  const std::size_t n_obs = validate_observation_matrices(xin, xout);
  validate_triplet_layout(triplets, n_obs);
  const std::size_t thread_count = quadra::checked_thread_count(n_threads);
  const std::function<Dfun> dfunin = create_dfun(metric_in);
  const std::function<Dfun> dfunout = create_dfun(metric_out);
  const auto triplets_cpp = Rcpp::as<std::vector<std::size_t>>(triplets);
  const auto xin_cpp = Rcpp::as<std::vector<double>>(xin);
  const auto xout_cpp = Rcpp::as<std::vector<double>>(xout);

  return legacy_triplet_sample(triplets_cpp.begin(), triplets_cpp.end(), n_obs,
                               xin_cpp.begin(), xin_cpp.end(), xout_cpp.begin(),
                               xout_cpp.end(), dfunin, dfunout, thread_count);
}

// [[Rcpp::export(rng = false)]]
SEXP triplet_plan_sample(const IntegerMatrix &triplets,
                         const NumericMatrix &xin, const NumericMatrix &xout,
                         const std::string &metric_in = "sqeuclidean",
                         const std::string &metric_out = "sqeuclidean",
                         double n_threads = 0, bool ret_triplets = false) {

  const std::size_t n_obs = validate_observation_matrices(xin, xout);
  if (triplets.nrow() < 1 || triplets.ncol() != 3) {
    stop("triplets must have at least one row and exactly three columns");
  }
  const std::size_t thread_count = quadra::checked_thread_count(n_threads);
  const std::size_t triplet_count = static_cast<std::size_t>(triplets.nrow());

  std::vector<std::size_t> anchor(triplet_count);
  std::vector<std::size_t> endpoint1(triplet_count);
  std::vector<std::size_t> endpoint2(triplet_count);
  std::vector<int> returned_anchor;
  std::vector<int> returned_endpoint1;
  std::vector<int> returned_endpoint2;
  std::vector<int> agreement;
  if (ret_triplets) {
    returned_anchor.resize(triplet_count);
    returned_endpoint1.resize(triplet_count);
    returned_endpoint2.resize(triplet_count);
    agreement.resize(triplet_count);
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
      returned_anchor[i] = anchor_i;
      returned_endpoint1[i] = endpoint1_i;
      returned_endpoint2[i] = endpoint2_i;
    }
  }

  const std::function<Dfun> dfunin = create_dfun(metric_in);
  const std::function<Dfun> dfunout = create_dfun(metric_out);
  const auto xin_cpp = Rcpp::as<std::vector<double>>(xin);
  const auto xout_cpp = Rcpp::as<std::vector<double>>(xout);

  const double accuracy = triplet_plan_sample(
      anchor, endpoint1, endpoint2, n_obs, xin_cpp.begin(), xin_cpp.end(),
      xout_cpp.begin(), xout_cpp.end(), dfunin, dfunout, thread_count,
      ret_triplets ? &agreement : nullptr);

  return triplet_result(accuracy, ret_triplets, returned_anchor,
                        returned_endpoint1, returned_endpoint2, agreement);
}

// [[Rcpp::export(rng = false)]]
SEXP random_triplet_sample(const NumericMatrix &xin, const NumericMatrix &xout,
                           double n_triplets = 5,
                           const std::string &metric_in = "sqeuclidean",
                           const std::string &metric_out = "sqeuclidean",
                           double n_threads = 0, bool ret_triplets = false) {

  const std::size_t n_obs = validate_observation_matrices(xin, xout);
  const std::size_t triplet_count =
      quadra::checked_nonnegative_count(n_triplets, "n_triplets");
  const std::size_t thread_count = quadra::checked_thread_count(n_threads);

  if (ret_triplets &&
      triplet_count >
          static_cast<std::size_t>((std::numeric_limits<int>::max)()) / n_obs) {
    stop("detailed triplet sample exceeds the supported integer range");
  }
  const std::size_t comparison_count = n_obs * triplet_count;
  std::vector<int> returned_anchor;
  std::vector<int> returned_endpoint1;
  std::vector<int> returned_endpoint2;
  std::vector<int> agreement;
  if (ret_triplets) {
    returned_anchor.resize(comparison_count);
    returned_endpoint1.resize(comparison_count);
    returned_endpoint2.resize(comparison_count);
    agreement.resize(comparison_count);
  }

  const std::function<Dfun> dfunin = create_dfun(metric_in);
  const std::function<Dfun> dfunout = create_dfun(metric_out);
  const auto xin_cpp = Rcpp::as<std::vector<double>>(xin);
  const auto xout_cpp = Rcpp::as<std::vector<double>>(xout);

  Rcpp::RNGScope rng_scope;
  const double accuracy = random_triplet_sample(
      triplet_count, n_obs, xin_cpp.begin(), xin_cpp.end(), xout_cpp.begin(),
      xout_cpp.end(), dfunin, dfunout, thread_count,
      ret_triplets ? &returned_anchor : nullptr,
      ret_triplets ? &returned_endpoint1 : nullptr,
      ret_triplets ? &returned_endpoint2 : nullptr,
      ret_triplets ? &agreement : nullptr);

  return triplet_result(accuracy, ret_triplets, returned_anchor,
                        returned_endpoint1, returned_endpoint2, agreement);
}
