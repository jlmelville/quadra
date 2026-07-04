#include <algorithm>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <Rcpp.h>

#include "rnndescent/random.h"
#include "tdoann/distance.h"

#include "distance.h"
#include "pforr.h"

using namespace Rcpp;
using It = typename std::vector<double>::const_iterator;
using Dfun = double(It, It, It);
using TripIt = typename std::vector<std::size_t>::const_iterator;

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

void update_triplet_counts(TripletCounts& counts, std::size_t p1,
                           std::size_t p2, const It xin_i_begin,
                           const It xin_i_end, const It xin_begin,
                           std::size_t xin_ncol,
                           const std::function<Dfun>& dfunin,
                           const It xout_i_begin, const It xout_i_end,
                           const It xout_begin, std::size_t xout_ncol,
                           const std::function<Dfun>& dfunout) {

  auto din_ip1 = dfunin(xin_i_begin, xin_i_end, xin_begin + p1 * xin_ncol);
  auto din_ip2 = dfunin(xin_i_begin, xin_i_end, xin_begin + p2 * xin_ncol);

  auto dout_ip1 =
      dfunout(xout_i_begin, xout_i_end, xout_begin + p1 * xout_ncol);
  auto dout_ip2 =
      dfunout(xout_i_begin, xout_i_end, xout_begin + p2 * xout_ncol);

  const int input_order = compare_distances(din_ip1, din_ip2);
  if (input_order == 0) {
    return;
  }

  ++counts.comparisons;
  const int output_order = compare_distances(dout_ip1, dout_ip2);
  if (input_order == output_order) {
    ++counts.agreements;
  }
}

double summarize_triplet_counts(const std::vector<TripletCounts>& counts) {
  TripletCounts total_counts;
  for (const auto& chunk_counts : counts) {
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
                                   std::size_t ntriplets_per_obs,
                                   std::size_t nobs,
                                   const TripIt triplets_begin,
                                   const It xin_begin, std::size_t xin_ncol,
                                   const std::function<Dfun>& dfunin,
                                   const It xout_begin, std::size_t xout_ncol,
                                   const std::function<Dfun>& dfunout) {

  TripletCounts counts;
  const std::size_t nt2 = ntriplets_per_obs * 2;
  for (std::size_t i = begin; i < end; i++) {
    const TripIt trip_obs_begin = triplets_begin + i * nt2;
    const It xin_i_begin = xin_begin + i * xin_ncol;
    const It xin_i_end = xin_i_begin + xin_ncol;

    const It xout_i_begin = xout_begin + i * xout_ncol;
    const It xout_i_end = xout_i_begin + xout_ncol;

    for (std::size_t j = 0; j < ntriplets_per_obs; j++) {
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
    std::size_t begin, std::size_t end, std::size_t ntriplets_per_obs,
    std::size_t nobs, tdoann::RandomIntGenerator<uint64_t>& int_sampler,
    const It xin_begin, std::size_t xin_ncol, const std::function<Dfun>& dfunin,
    const It xout_begin, std::size_t xout_ncol,
    const std::function<Dfun>& dfunout) {

  TripletCounts counts;
  for (std::size_t i = begin; i < end; i++) {
    const It xin_i_begin = xin_begin + i * xin_ncol;
    const It xin_i_end = xin_i_begin + xin_ncol;

    const It xout_i_begin = xout_begin + i * xout_ncol;
    const It xout_i_end = xout_i_begin + xout_ncol;

    for (std::size_t j = 0; j < ntriplets_per_obs; j++) {
      const auto idxs = int_sampler.sample(nobs - 1, 2);
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
                      std::size_t nobs, It xin_begin, It xin_end, It xout_begin,
                      It xout_end, const std::function<Dfun>& dfunin,
                      const std::function<Dfun>& dfunout,
                      std::size_t n_threads) {

  const std::size_t ntriplets_per_obs =
      (triplets_end - triplets_begin) / nobs / 2;
  const std::size_t xin_nfeat = (xin_end - xin_begin) / nobs;
  const std::size_t xout_nfeat = (xout_end - xout_begin) / nobs;

  std::vector<TripletCounts> counts(n_parallel_chunks(0, nobs, n_threads));

  auto worker = [&](std::size_t begin, std::size_t end, std::size_t thread_id) {
    const auto chunk_counts = triplet_sample_inner(
        begin, end, ntriplets_per_obs, nobs, triplets_begin, xin_begin,
        xin_nfeat, dfunin, xout_begin, xout_nfeat, dfunout);
    counts[thread_id].agreements += chunk_counts.agreements;
    counts[thread_id].comparisons += chunk_counts.comparisons;
  };

  pforr::parallel_for_indexed(0, nobs, worker, n_threads);

  return summarize_triplet_counts(counts);
}

double random_triplet_sample(std::size_t ntriplets_per_obs, std::size_t nobs,
                             It xin_begin, It xin_end, It xout_begin,
                             It xout_end, const std::function<Dfun>& dfunin,
                             const std::function<Dfun>& dfunout,
                             std::size_t n_threads) {
  const std::size_t xin_nfeat = (xin_end - xin_begin) / nobs;
  const std::size_t xout_nfeat = (xout_end - xout_begin) / nobs;

  rnndescent::ParallelIntRNGAdapter<uint64_t, rnndescent::DQIntSampler>
      sampler_provider;
  sampler_provider.initialize();

  std::vector<TripletCounts> counts(n_parallel_chunks(0, nobs, n_threads));

  auto worker = [&](std::size_t begin, std::size_t end, std::size_t chunk_id) {
    auto thread_sampler = sampler_provider.get_parallel_instance(chunk_id);
    const auto chunk_counts = random_triplet_sample_inner(
        begin, end, ntriplets_per_obs, nobs, *thread_sampler, xin_begin,
        xin_nfeat, dfunin, xout_begin, xout_nfeat, dfunout);
    counts[chunk_id].agreements += chunk_counts.agreements;
    counts[chunk_id].comparisons += chunk_counts.comparisons;
  };

  pforr::parallel_for_indexed(0, nobs, worker, n_threads);

  return summarize_triplet_counts(counts);
}

// [[Rcpp::export]]
double triplet_sample(const IntegerMatrix& triplets, const NumericMatrix& xin,
                      const NumericMatrix& xout,
                      const std::string& metric_in = "sqeuclidean",
                      const std::string& metric_out = "sqeuclidean",
                      std::size_t n_threads = 0) {

  std::function<Dfun> dfunin = create_dfun(metric_in);
  std::function<Dfun> dfunout = create_dfun(metric_out);

  auto triplets_cpp = Rcpp::as<std::vector<std::size_t>>(triplets);
  auto xin_cpp = Rcpp::as<std::vector<double>>(xin);
  auto xout_cpp = Rcpp::as<std::vector<double>>(xout);

  return triplet_sample(triplets_cpp.begin(), triplets_cpp.end(),
                        triplets.ncol(), xin_cpp.begin(), xin_cpp.end(),
                        xout_cpp.begin(), xout_cpp.end(), dfunin, dfunout,
                        n_threads);
}

// [[Rcpp::export]]
double random_triplet_sample(const NumericMatrix& xin,
                             const NumericMatrix& xout,
                             std::size_t n_triplets = 5,
                             const std::string& metric_in = "sqeuclidean",
                             const std::string& metric_out = "sqeuclidean",
                             std::size_t n_threads = 0) {

  std::function<Dfun> dfunin = create_dfun(metric_in);
  std::function<Dfun> dfunout = create_dfun(metric_out);

  auto xin_cpp = Rcpp::as<std::vector<double>>(xin);
  auto xout_cpp = Rcpp::as<std::vector<double>>(xout);

  return random_triplet_sample(n_triplets, xin.ncol(), xin_cpp.begin(),
                               xin_cpp.end(), xout_cpp.begin(), xout_cpp.end(),
                               dfunin, dfunout, n_threads);
}
