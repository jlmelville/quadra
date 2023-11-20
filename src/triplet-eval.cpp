#include <functional>
#include <memory>
#include <string>

#include <Rcpp.h>

#include "tdoann/distance.h"

#include "distance.h"
#include "quadrapforr.h"

using namespace Rcpp;
using It = typename std::vector<double>::const_iterator;
using Dfun = double(It, It, It);
using TripIt = typename std::vector<std::size_t>::const_iterator;

std::size_t triplet_sample_inner(std::size_t begin, std::size_t end,
                                 std::size_t ntriplets_per_obs,
                                 std::size_t nobs, const TripIt triplets_begin,
                                 const It xin_begin, std::size_t xin_ncol,
                                 const std::function<Dfun> &dfunin,
                                 const It xout_begin, std::size_t xout_ncol,
                                 const std::function<Dfun> &dfunout) {

  std::size_t acc{0};
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

double triplet_sample(TripIt triplets_begin, TripIt triplets_end,
                      std::size_t nobs, It xin_begin, It xin_end, It xout_begin,
                      It xout_end, const std::function<Dfun> &dfunin,
                      const std::function<Dfun> &dfunout,
                      std::size_t n_threads) {

  const std::size_t ntriplets_per_obs =
      (triplets_end - triplets_begin) / nobs / 2;
  const std::size_t xin_nfeat = (xin_end - xin_begin) / nobs;
  const std::size_t xout_nfeat = (xout_end - xout_begin) / nobs;

  std::vector<std::size_t> accs(std::max(n_threads, std::size_t{1}));

  auto worker = [&](std::size_t begin, std::size_t end, std::size_t thread_id) {
    accs[thread_id] += triplet_sample_inner(
        begin, end, ntriplets_per_obs, nobs, triplets_begin, xin_begin,
        xin_nfeat, dfunin, xout_begin, xout_nfeat, dfunout);
  };

  pforr::pfor(0, nobs, worker, n_threads);

  std::size_t acc{0};
  for (auto a : accs) {
    acc += a;
  }
  return acc / static_cast<double>(nobs * ntriplets_per_obs);
}

// [[Rcpp::export]]
double triplet_sample(const IntegerMatrix &triplets, const NumericMatrix &xin,
                      const NumericMatrix &xout,
                      const std::string &metric_in = "sqeuclidean",
                      const std::string &metric_out = "sqeuclidean",
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
