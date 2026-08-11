#include <vector>

#include <Rcpp.h>

#include "pforr.h"

#include "rnndescent/random.h"

#include "distance.h"
#include "native-validation.h"

using namespace Rcpp;

using It = typename std::vector<double>::const_iterator;
using Dfun = double(It, It, It);

void distance_sample_inner(std::size_t begin, std::size_t end,
                           tdoann::RandomIntGenerator<uint64_t>& int_sampler,
                           std::size_t n_obs, const It xin_begin,
                           std::size_t xin_ncol,
                           const std::function<Dfun>& dfunin,
                           const It xout_begin, std::size_t xout_ncol,
                           const std::function<Dfun>& dfunout,
                           std::vector<double>& din,
                           std::vector<double>& dout) {

  for (std::size_t i = begin; i < end; i++) {
    // DQIntSampler::sample(n_obs, 2) draws distinct indices without
    // replacement.
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
                      const std::function<Dfun>& dfunin,
                      const std::function<Dfun>& dfunout,
                      std::vector<double>& din, std::vector<double>& dout,
                      std::size_t n_threads) {
  const std::size_t xin_nfeat = (xin_end - xin_begin) / n_obs;
  const std::size_t xout_nfeat = (xout_end - xout_begin) / n_obs;

  rnndescent::ParallelIntRNGAdapter<uint64_t, rnndescent::DQIntSampler>
      sampler_provider;

  sampler_provider.initialize();

  auto worker = [&](std::size_t begin, std::size_t end, std::size_t chunk_id) {
    auto thread_sampler = sampler_provider.get_parallel_instance(chunk_id);
    distance_sample_inner(begin, end, *thread_sampler, n_obs, xin_begin,
                          xin_nfeat, dfunin, xout_begin, xout_nfeat, dfunout,
                          din, dout);
  };

  pforr::parallel_for_indexed(n_pairs, worker, n_threads);
}

// [[Rcpp::export]]
List random_distances(NumericMatrix xin, NumericMatrix xout,
                      const std::string& metric_in = "euclidean",
                      const std::string& metric_out = "euclidean",
                      double n_pairs = 10000, double n_threads = 0) {
  if (xin.ncol() != xout.ncol()) {
    stop("xin and xout must have the same number of observations");
  }
  if (xin.nrow() < 1 || xout.nrow() < 1 || xin.ncol() < 2) {
    stop("xin and xout must each have at least one feature and two "
         "observations");
  }

  const std::size_t pair_count =
      quadra::validate_positive_count(n_pairs, "n_pairs");
  const std::size_t thread_count = quadra::validate_n_threads(n_threads);

  std::function<Dfun> dfunin = create_dfun(metric_in);
  std::function<Dfun> dfunout = create_dfun(metric_out);

  std::vector<double> din(pair_count);
  std::vector<double> dout(pair_count);

  auto xin_cpp = Rcpp::as<std::vector<double>>(xin);
  auto xout_cpp = Rcpp::as<std::vector<double>>(xout);

  random_distances(pair_count, xin.ncol(), xin_cpp.begin(), xin_cpp.end(),
                   xout_cpp.begin(), xout_cpp.end(), dfunin, dfunout, din, dout,
                   thread_count);

  NumericVector res_in(din.begin(), din.end());
  NumericVector res_out(dout.begin(), dout.end());

  return List::create(_("din") = res_in, _("dout") = res_out);
}
