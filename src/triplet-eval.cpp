#include <functional>
#include <iostream>
#include <memory>
#include <string>

#include <Rcpp.h>
using namespace Rcpp;

#include "RcppPerpendicular.h"
#include "tdoann/distance.h"

template <typename TripIt, typename XIt, typename Xout>
std::size_t triplet_sample_inner(
    std::size_t begin, std::size_t end, std::size_t ntriplets_per_obs,
    std::size_t nobs, const TripIt triplets_begin, const XIt xin_begin,
    std::size_t xin_ncol,
    const std::function<Xout(const XIt, const XIt, const XIt)> &dfunin,
    const XIt xout_begin, std::size_t xout_ncol,
    const std::function<Xout(const XIt, const XIt, const XIt)> &dfunout) {

  std::size_t acc{0};
  const std::size_t nt2 = ntriplets_per_obs * 2;
  for (std::size_t i = begin; i < end; i++) {
    const TripIt trip_obs_begin = triplets_begin + i * nt2;
    const XIt xin_i_begin = xin_begin + i * xin_ncol;
    const XIt xin_i_end = xin_i_begin + xin_ncol;

    const XIt xout_i_begin = xout_begin + i * xout_ncol;
    const XIt xout_i_end = xout_i_begin + xout_ncol;

    for (std::size_t j = 0; j < ntriplets_per_obs; j++) {
      const TripIt tripi = trip_obs_begin + (j * 2);
      const std::size_t p1 = *tripi;
      const std::size_t p2 = *(tripi + 1);

      auto din_ip1 = dfunin(xin_i_begin, xin_i_end, xin_begin + p1 * xin_ncol);
      auto din_ip2 = dfunin(xin_i_begin, xin_i_end, xin_begin + p2 * xin_ncol);

      auto dout_ip1 =
          dfunout(xout_i_begin, xout_i_end, xout_begin + p1 * xout_ncol);
      auto dout_ip2 =
          dfunout(xout_i_begin, xout_i_end, xout_begin + p2 * xout_ncol);

      if ((din_ip1 < din_ip2) == (dout_ip1 < dout_ip2)) {
        ++acc;
      }
    }
  }
  return acc;
}

double triplet_sample(
    const int *triplets_begin, const int *triplets_end, std::size_t nobs,
    const double *xin_begin, const double *xin_end, const double *xout_begin,
    const double *xout_end,
    const std::function<double(const double *, const double *, const double *)>
        &dfunin,
    const std::function<double(const double *, const double *, const double *)>
        &dfunout,
    std::size_t n_threads, std::size_t grain_size) {

  const std::size_t ntriplets_per_obs =
      (triplets_end - triplets_begin) / nobs / 2;
  const std::size_t xin_nfeat = (xin_end - xin_begin) / nobs;
  const std::size_t xout_nfeat = (xout_end - xout_begin) / nobs;

  std::vector<std::size_t> accs(std::max(n_threads, std::size_t{1}));

  auto worker = [&](std::size_t begin, std::size_t end, std::size_t thread_id) {
    accs[thread_id] = triplet_sample_inner(
        begin, end, ntriplets_per_obs, nobs, triplets_begin, xin_begin,
        xin_nfeat, dfunin, xout_begin, xout_nfeat, dfunout);
  };

  RcppPerpendicular::pfor(nobs, worker, n_threads, grain_size);

  std::size_t acc{0};
  for (auto a : accs) {
    acc += a;
  }
  return acc / static_cast<double>(nobs * ntriplets_per_obs);
}

std::function<double(const double *, const double *, const double *)>
create_dfun(const std::string &metric) {
  if (metric == "euclidean") {
    return tdoann::euclidean<double, const double *>;
  } else if (metric == "l2sqr") {
    return tdoann::l2sqr<double, const double *>;
  } else if (metric == "cosine") {
    return tdoann::cosine<double, const double *>;
  } else if (metric == "hamming") {
    return tdoann::hamming<double, const double *>;
  } else if (metric == "manhattan") {
    return tdoann::manhattan<double, const double *>;
  } else if (metric == "correlation") {
    return tdoann::correlation<double, const double *>;
  } else {
    return tdoann::euclidean<double, const double *>;
  }
}

// [[Rcpp::export]]
double triplet_sample(const IntegerMatrix &triplets, const NumericMatrix &xin,
                      const NumericMatrix &xout,
                      const std::string &metric_in = "l2sqr",
                      const std::string &metric_out = "l2sqr",
                      std::size_t n_threads = 0, std::size_t grain_size = 1) {

  using Dfun = double(const double *, const double *, const double *);
  std::function<Dfun> dfunin = create_dfun(metric_in);
  std::function<Dfun> dfunout = create_dfun(metric_out);

  return triplet_sample(triplets.begin(), triplets.end(), triplets.ncol(),
                        xin.begin(), xin.end(), xout.begin(), xout.end(),
                        dfunin, dfunout, n_threads, grain_size);
}
