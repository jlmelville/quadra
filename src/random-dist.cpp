#include <iostream>
#include <vector>

#include <Rcpp.h>

#include "RcppPerpendicular.h"

#include "rnndescent/random.h"

#include "distance.h"

using namespace Rcpp;

template <typename Sampler, typename XIt, typename Xout>
void distance_sample_inner(
    std::size_t begin, std::size_t end, uint64_t seed, std::size_t n_obs,
    const XIt xin_begin, std::size_t xin_ncol,
    const std::function<Xout(const XIt, const XIt, const XIt)> &dfunin,
    const XIt xout_begin, std::size_t xout_ncol,
    const std::function<Xout(const XIt, const XIt, const XIt)> &dfunout,
    std::vector<Xout> &din, std::vector<Xout> &dout) {

  Sampler int_sampler(seed, end);
  for (std::size_t i = begin; i < end; i++) {
    auto idxs = int_sampler.sample(n_obs, 2);

    const XIt xin_i_begin = xin_begin + idxs[0] * xin_ncol;
    din[i] = dfunin(xin_i_begin, xin_i_begin + xin_ncol,
                    xin_begin + idxs[1] * xin_ncol);

    const XIt xout_i_begin = xout_begin + idxs[0] * xout_ncol;
    dout[i] = dfunout(xout_i_begin, xout_i_begin + xout_ncol,
                      xout_begin + idxs[1] * xout_ncol);
  }
}

template <typename Sampler>
void random_distances(
    std::size_t n_pairs, std::size_t n_obs, const double *xin_begin,
    const double *xin_end, const double *xout_begin, const double *xout_end,
    const std::function<double(const double *, const double *, const double *)>
        &dfunin,
    const std::function<double(const double *, const double *, const double *)>
        &dfunout,
    std::vector<double> &din, std::vector<double> &dout,
    std::size_t n_threads) {
  const std::size_t xin_nfeat = (xin_end - xin_begin) / n_obs;
  const std::size_t xout_nfeat = (xout_end - xout_begin) / n_obs;

  uint64_t seed = Sampler::get_seed();

  auto worker = [&](std::size_t begin, std::size_t end) {
    distance_sample_inner<Sampler>(begin, end, seed, n_obs, xin_begin,
                                   xin_nfeat, dfunin, xout_begin, xout_nfeat,
                                   dfunout, din, dout);
  };

  RcppPerpendicular::parallel_for(n_pairs, worker, n_threads);
}

// [[Rcpp::export]]
List random_distances(NumericMatrix xin, NumericMatrix xout,
                      const std::string &metric_in = "euclidean",
                      const std::string &metric_out = "euclidean",
                      std::size_t n_pairs = 10000, std::size_t n_threads = 0,
                      bool verbose = false) {
  using Dfun = double(const double *, const double *, const double *);
  std::function<Dfun> dfunin = create_dfun(metric_in);
  std::function<Dfun> dfunout = create_dfun(metric_out);

  std::vector<double> din(n_pairs);
  std::vector<double> dout(n_pairs);

  random_distances<rnndescent::DQIntSampler<uint32_t>>(
      n_pairs, xin.ncol(), xin.begin(), xin.end(), xout.begin(), xout.end(),
      dfunin, dfunout, din, dout, n_threads);

  NumericVector res_in(din.begin(), din.end());
  NumericVector res_out(dout.begin(), dout.end());

  return List::create(_("din") = res_in, _("dout") = res_out);
}
