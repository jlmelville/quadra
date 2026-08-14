#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
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

std::vector<std::size_t> copy_checked_triplets(const NumericMatrix &triplets,
                                               std::size_t n_obs) {
  if (triplets.ncol() != static_cast<int>(n_obs)) {
    stop("triplets must have one column per observation");
  }
  if ((triplets.nrow() % 2) != 0) {
    stop("triplets must have an even number of rows");
  }

  std::vector<std::size_t> copied;
  copied.reserve(triplets.size());
  for (int col = 0; col < triplets.ncol(); ++col) {
    for (int row = 0; row < triplets.nrow(); row += 2) {
      const double first = triplets(row, col);
      const double second = triplets(row + 1, col);
      const bool invalid = !R_finite(first) || !R_finite(second) ||
                           first != std::floor(first) ||
                           second != std::floor(second) || first < 0 ||
                           second < 0 || first >= static_cast<double>(n_obs) ||
                           second >= static_cast<double>(n_obs);
      if (invalid) {
        stop("triplets must contain finite integer indices in range");
      }
      copied.push_back(static_cast<std::size_t>(first));
      copied.push_back(static_cast<std::size_t>(second));
    }
  }
  return copied;
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

void update_triplet_counts(TripletCounts &counts, std::size_t p1,
                           std::size_t p2, const It xin_i_begin,
                           const It xin_i_end, const It xin_begin,
                           std::size_t xin_ncol,
                           const std::function<Dfun> &dfunin,
                           const It xout_i_begin, const It xout_i_end,
                           const It xout_begin, std::size_t xout_ncol,
                           const std::function<Dfun> &dfunout) {

  auto input_distance_1 =
      dfunin(xin_i_begin, xin_i_end, xin_begin + p1 * xin_ncol);
  auto input_distance_2 =
      dfunin(xin_i_begin, xin_i_end, xin_begin + p2 * xin_ncol);

  auto output_distance_1 =
      dfunout(xout_i_begin, xout_i_end, xout_begin + p1 * xout_ncol);
  auto output_distance_2 =
      dfunout(xout_i_begin, xout_i_end, xout_begin + p2 * xout_ncol);

  const int input_order = compare_distances(input_distance_1, input_distance_2);
  if (input_order == 0) {
    return;
  }

  ++counts.comparisons;
  const int output_order =
      compare_distances(output_distance_1, output_distance_2);
  if (input_order == output_order) {
    ++counts.agreements;
  }
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

TripletCounts triplet_sample_inner(std::size_t begin, std::size_t end,
                                   std::size_t n_triplets_per_observation,
                                   std::size_t n_obs,
                                   const TripIt triplets_begin,
                                   const It xin_begin, std::size_t xin_ncol,
                                   const std::function<Dfun> &dfunin,
                                   const It xout_begin, std::size_t xout_ncol,
                                   const std::function<Dfun> &dfunout) {

  TripletCounts counts;
  const std::size_t nt2 = n_triplets_per_observation * 2;
  for (std::size_t i = begin; i < end; i++) {
    const TripIt trip_obs_begin = triplets_begin + i * nt2;
    const It xin_i_begin = xin_begin + i * xin_ncol;
    const It xin_i_end = xin_i_begin + xin_ncol;

    const It xout_i_begin = xout_begin + i * xout_ncol;
    const It xout_i_end = xout_i_begin + xout_ncol;

    for (std::size_t j = 0; j < n_triplets_per_observation; j++) {
      const TripIt tripi = trip_obs_begin + (j * 2);
      const std::size_t p1 = *tripi;
      const std::size_t p2 = *(tripi + 1);

      update_triplet_counts(counts, p1, p2, xin_i_begin, xin_i_end, xin_begin,
                            xin_ncol, dfunin, xout_i_begin, xout_i_end,
                            xout_begin, xout_ncol, dfunout);
    }
  }
  return counts;
}

std::size_t avoid_anchor_index(std::size_t idx, std::size_t anchor) {
  return idx >= anchor ? idx + 1 : idx;
}

TripletCounts random_triplet_sample_inner(
    std::size_t begin, std::size_t end, std::size_t n_triplets_per_observation,
    std::size_t n_obs, tdoann::RandomIntGenerator<uint64_t> &int_sampler,
    const It xin_begin, std::size_t xin_ncol, const std::function<Dfun> &dfunin,
    const It xout_begin, std::size_t xout_ncol,
    const std::function<Dfun> &dfunout) {

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

      update_triplet_counts(counts, p1, p2, xin_i_begin, xin_i_end, xin_begin,
                            xin_ncol, dfunin, xout_i_begin, xout_i_end,
                            xout_begin, xout_ncol, dfunout);
    }
  }
  return counts;
}

double triplet_sample(TripIt triplets_begin, TripIt triplets_end,
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

  auto worker = [&](std::size_t begin, std::size_t end, std::size_t thread_id) {
    const auto chunk_counts = triplet_sample_inner(
        begin, end, n_triplets_per_observation, n_obs, triplets_begin,
        xin_begin, xin_nfeat, dfunin, xout_begin, xout_nfeat, dfunout);
    counts[thread_id].agreements += chunk_counts.agreements;
    counts[thread_id].comparisons += chunk_counts.comparisons;
  };

  pforr::parallel_for_indexed(0, n_obs, worker, n_threads);

  return summarize_triplet_counts(counts);
}

double random_triplet_sample(std::size_t n_triplets_per_observation,
                             std::size_t n_obs, It xin_begin, It xin_end,
                             It xout_begin, It xout_end,
                             const std::function<Dfun> &dfunin,
                             const std::function<Dfun> &dfunout,
                             std::size_t n_threads) {
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
        xin_begin, xin_nfeat, dfunin, xout_begin, xout_nfeat, dfunout);
    counts[chunk_id].agreements += chunk_counts.agreements;
    counts[chunk_id].comparisons += chunk_counts.comparisons;
  };

  pforr::parallel_for_indexed(0, n_obs, worker, n_threads);

  return summarize_triplet_counts(counts);
}

// [[Rcpp::export]]
double triplet_sample(const NumericMatrix &triplets, const NumericMatrix &xin,
                      const NumericMatrix &xout,
                      const std::string &metric_in = "sqeuclidean",
                      const std::string &metric_out = "sqeuclidean",
                      double n_threads = 0) {

  const std::size_t n_obs = validate_observation_matrices(xin, xout);
  const auto triplets_cpp = copy_checked_triplets(triplets, n_obs);
  const std::size_t thread_count = quadra::checked_thread_count(n_threads);

  std::function<Dfun> dfunin = create_dfun(metric_in);
  std::function<Dfun> dfunout = create_dfun(metric_out);

  auto xin_cpp = Rcpp::as<std::vector<double>>(xin);
  auto xout_cpp = Rcpp::as<std::vector<double>>(xout);

  return triplet_sample(triplets_cpp.begin(), triplets_cpp.end(), n_obs,
                        xin_cpp.begin(), xin_cpp.end(), xout_cpp.begin(),
                        xout_cpp.end(), dfunin, dfunout, thread_count);
}

// [[Rcpp::export]]
double random_triplet_sample(const NumericMatrix &xin,
                             const NumericMatrix &xout, double n_triplets = 5,
                             const std::string &metric_in = "sqeuclidean",
                             const std::string &metric_out = "sqeuclidean",
                             double n_threads = 0) {

  const std::size_t n_obs = validate_observation_matrices(xin, xout);
  const std::size_t triplet_count =
      quadra::checked_nonnegative_count(n_triplets, "n_triplets");
  const std::size_t thread_count = quadra::checked_thread_count(n_threads);

  std::function<Dfun> dfunin = create_dfun(metric_in);
  std::function<Dfun> dfunout = create_dfun(metric_out);

  auto xin_cpp = Rcpp::as<std::vector<double>>(xin);
  auto xout_cpp = Rcpp::as<std::vector<double>>(xout);

  return random_triplet_sample(triplet_count, n_obs, xin_cpp.begin(),
                               xin_cpp.end(), xout_cpp.begin(), xout_cpp.end(),
                               dfunin, dfunout, thread_count);
}
