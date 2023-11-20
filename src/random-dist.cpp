#include <iostream>
#include <vector>

#include <Rcpp.h>

#include "quadrapforr.h"

#include "rnndescent/random.h"

#include "distance.h"

using namespace Rcpp;

using It = typename std::vector<double>::const_iterator;
using Dfun = double(It, It, It);

void distance_sample_inner(std::size_t begin, std::size_t end,
                           tdoann::RandomIntGenerator<uint64_t> &int_sampler,
                           std::size_t n_obs, const It xin_begin,
                           std::size_t xin_ncol,
                           const std::function<Dfun> &dfunin,
                           const It xout_begin, std::size_t xout_ncol,
                           const std::function<Dfun> &dfunout,
                           std::vector<double> &din,
                           std::vector<double> &dout) {

  for (std::size_t i = begin; i < end; i++) {
    auto idxs = int_sampler.sample(n_obs, 2);

    const It xin_i_begin = xin_begin + idxs[0] * xin_ncol;
    din[i] = dfunin(xin_i_begin, xin_i_begin + xin_ncol,
                    xin_begin + idxs[1] * xin_ncol);

    const It xout_i_begin = xout_begin + idxs[0] * xout_ncol;
    dout[i] = dfunout(xout_i_begin, xout_i_begin + xout_ncol,
                      xout_begin + idxs[1] * xout_ncol);
  }
}

void random_distances(std::size_t n_pairs, std::size_t n_obs, It xin_begin,
                      It xin_end, It xout_begin, It xout_end,
                      const std::function<Dfun> &dfunin,
                      const std::function<Dfun> &dfunout,
                      std::vector<double> &din, std::vector<double> &dout,
                      std::size_t n_threads) {
  const std::size_t xin_nfeat = (xin_end - xin_begin) / n_obs;
  const std::size_t xout_nfeat = (xout_end - xout_begin) / n_obs;

  rnndescent::ParallelIntRNGAdapter<uint64_t, rnndescent::DQIntSampler>
      sampler_provider;

  sampler_provider.initialize();

  auto worker = [&](std::size_t begin, std::size_t end) {
    auto thread_sampler = sampler_provider.get_parallel_instance(end);
    distance_sample_inner(begin, end, *thread_sampler, n_obs, xin_begin,
                          xin_nfeat, dfunin, xout_begin, xout_nfeat, dfunout,
                          din, dout);
  };

  pforr::parallel_for(n_pairs, worker, n_threads);
}

// [[Rcpp::export]]
List random_distances(NumericMatrix xin, NumericMatrix xout,
                      const std::string &metric_in = "euclidean",
                      const std::string &metric_out = "euclidean",
                      std::size_t n_pairs = 10000, std::size_t n_threads = 0,
                      bool verbose = false) {

  std::function<Dfun> dfunin = create_dfun(metric_in);
  std::function<Dfun> dfunout = create_dfun(metric_out);

  std::vector<double> din(n_pairs);
  std::vector<double> dout(n_pairs);

  auto xin_cpp = Rcpp::as<std::vector<double>>(xin);
  auto xout_cpp = Rcpp::as<std::vector<double>>(xout);

  random_distances(n_pairs, xin.ncol(), xin_cpp.begin(), xin_cpp.end(),
                   xout_cpp.begin(), xout_cpp.end(), dfunin, dfunout, din, dout,
                   n_threads);

  NumericVector res_in(din.begin(), din.end());
  NumericVector res_out(dout.begin(), dout.end());

  return List::create(_("din") = res_in, _("dout") = res_out);
}
