#include <vector>

#include <Rcpp.h>

#include "pforr.h"

#include "rnndescent/random.h"

#include "distance.h"
#include "native-validation.h"

using namespace Rcpp;

using It = typename std::vector<double>::const_iterator;
using Dfun = double(It, It, It);

void distance_sample_inner(
    std::size_t begin, std::size_t end,
    tdoann::RandomIntGenerator<uint64_t> &int_sampler, std::size_t n_obs,
    const It xin_begin, std::size_t xin_ncol, const std::function<Dfun> &dfunin,
    const It xout_begin, std::size_t xout_ncol,
    const std::function<Dfun> &dfunout, std::vector<double> &input_distances,
    std::vector<double> &output_distances, std::vector<int> *endpoint1,
    std::vector<int> *endpoint2) {

  for (std::size_t i = begin; i < end; i++) {
    // DQIntSampler::sample(n_obs, 2) draws distinct indices without
    // replacement.
    auto idxs = int_sampler.sample(n_obs, 2);

    if (endpoint1 != nullptr) {
      (*endpoint1)[i] = static_cast<int>(idxs[0] + 1);
      (*endpoint2)[i] = static_cast<int>(idxs[1] + 1);
    }

    const It xin_i_begin = xin_begin + idxs[0] * xin_ncol;
    input_distances[i] = dfunin(xin_i_begin, xin_i_begin + xin_ncol,
                                xin_begin + idxs[1] * xin_ncol);

    const It xout_i_begin = xout_begin + idxs[0] * xout_ncol;
    output_distances[i] = dfunout(xout_i_begin, xout_i_begin + xout_ncol,
                                  xout_begin + idxs[1] * xout_ncol);
  }
}

void random_distances(std::size_t n_pairs, std::size_t n_obs, It xin_begin,
                      It xin_end, It xout_begin, It xout_end,
                      const std::function<Dfun> &dfunin,
                      const std::function<Dfun> &dfunout,
                      std::vector<double> &input_distances,
                      std::vector<double> &output_distances,
                      std::size_t n_threads, std::vector<int> *endpoint1,
                      std::vector<int> *endpoint2) {
  const std::size_t xin_nfeat = (xin_end - xin_begin) / n_obs;
  const std::size_t xout_nfeat = (xout_end - xout_begin) / n_obs;

  rnndescent::ParallelIntRNGAdapter<uint64_t, rnndescent::DQIntSampler>
      sampler_provider;

  sampler_provider.initialize();

  auto worker = [&](std::size_t begin, std::size_t end, std::size_t chunk_id) {
    auto thread_sampler = sampler_provider.get_parallel_instance(chunk_id);
    distance_sample_inner(begin, end, *thread_sampler, n_obs, xin_begin,
                          xin_nfeat, dfunin, xout_begin, xout_nfeat, dfunout,
                          input_distances, output_distances, endpoint1,
                          endpoint2);
  };

  pforr::parallel_for_indexed(n_pairs, worker, n_threads);
}

void pair_distances_inner(std::size_t begin, std::size_t end,
                          const std::vector<std::size_t> &endpoint1,
                          const std::vector<std::size_t> &endpoint2,
                          const It xin_begin, std::size_t xin_ncol,
                          const std::function<Dfun> &dfunin,
                          const It xout_begin, std::size_t xout_ncol,
                          const std::function<Dfun> &dfunout,
                          std::vector<double> &input_distances,
                          std::vector<double> &output_distances) {
  for (std::size_t i = begin; i < end; i++) {
    const It xin_i_begin = xin_begin + endpoint1[i] * xin_ncol;
    input_distances[i] = dfunin(xin_i_begin, xin_i_begin + xin_ncol,
                                xin_begin + endpoint2[i] * xin_ncol);

    const It xout_i_begin = xout_begin + endpoint1[i] * xout_ncol;
    output_distances[i] = dfunout(xout_i_begin, xout_i_begin + xout_ncol,
                                  xout_begin + endpoint2[i] * xout_ncol);
  }
}

void pair_distances(std::size_t n_pairs, std::size_t n_obs, It xin_begin,
                    It xin_end, It xout_begin, It xout_end,
                    const std::function<Dfun> &dfunin,
                    const std::function<Dfun> &dfunout,
                    const std::vector<std::size_t> &endpoint1,
                    const std::vector<std::size_t> &endpoint2,
                    std::vector<double> &input_distances,
                    std::vector<double> &output_distances,
                    std::size_t n_threads) {
  const std::size_t xin_nfeat = (xin_end - xin_begin) / n_obs;
  const std::size_t xout_nfeat = (xout_end - xout_begin) / n_obs;

  auto worker = [&](std::size_t begin, std::size_t end, std::size_t) {
    pair_distances_inner(begin, end, endpoint1, endpoint2, xin_begin, xin_nfeat,
                         dfunin, xout_begin, xout_nfeat, dfunout,
                         input_distances, output_distances);
  };

  pforr::parallel_for_indexed(n_pairs, worker, n_threads);
}

// [[Rcpp::export(rng = false)]]
List random_distances(NumericMatrix xin, NumericMatrix xout,
                      const std::string &metric_in = "euclidean",
                      const std::string &metric_out = "euclidean",
                      double n_pairs = 10000, double n_threads = 0,
                      Nullable<IntegerMatrix> pairs = R_NilValue,
                      bool ret_pairs = false) {
  if (xin.ncol() != xout.ncol()) {
    stop("xin and xout must have the same number of observations");
  }
  if (xin.nrow() < 1 || xout.nrow() < 1 || xin.ncol() < 2) {
    stop("xin and xout must each have at least one feature and two "
         "observations");
  }

  const std::size_t thread_count = quadra::checked_thread_count(n_threads);

  const bool has_pairs = pairs.isNotNull();
  std::size_t pair_count = 0;
  std::vector<std::size_t> pair_endpoint1;
  std::vector<std::size_t> pair_endpoint2;
  std::vector<int> returned_endpoint1;
  std::vector<int> returned_endpoint2;

  if (has_pairs) {
    IntegerMatrix pair_matrix(pairs);
    if (pair_matrix.nrow() < 1 || pair_matrix.ncol() != 2) {
      stop("pairs must have at least one row and exactly two columns");
    }
    pair_count = static_cast<std::size_t>(pair_matrix.nrow());
    pair_endpoint1.resize(pair_count);
    pair_endpoint2.resize(pair_count);
    if (ret_pairs) {
      returned_endpoint1.resize(pair_count);
      returned_endpoint2.resize(pair_count);
    }

    for (std::size_t i = 0; i < pair_count; i++) {
      const int endpoint1 = pair_matrix(i, 0);
      const int endpoint2 = pair_matrix(i, 1);
      if (endpoint1 == NA_INTEGER || endpoint2 == NA_INTEGER || endpoint1 < 1 ||
          endpoint2 < 1 || endpoint1 > xin.ncol() || endpoint2 > xin.ncol()) {
        stop("pairs must contain indices between 1 and the number of "
             "observations");
      }
      pair_endpoint1[i] = static_cast<std::size_t>(endpoint1 - 1);
      pair_endpoint2[i] = static_cast<std::size_t>(endpoint2 - 1);
      if (ret_pairs) {
        returned_endpoint1[i] = endpoint1;
        returned_endpoint2[i] = endpoint2;
      }
    }
  } else {
    pair_count = quadra::checked_nonnegative_count(n_pairs, "n_pairs");
    if (ret_pairs) {
      returned_endpoint1.resize(pair_count);
      returned_endpoint2.resize(pair_count);
    }
  }

  std::function<Dfun> dfunin = create_dfun(metric_in);
  std::function<Dfun> dfunout = create_dfun(metric_out);

  std::vector<double> input_distances(pair_count);
  std::vector<double> output_distances(pair_count);

  auto xin_cpp = Rcpp::as<std::vector<double>>(xin);
  auto xout_cpp = Rcpp::as<std::vector<double>>(xout);

  if (has_pairs) {
    pair_distances(pair_count, xin.ncol(), xin_cpp.begin(), xin_cpp.end(),
                   xout_cpp.begin(), xout_cpp.end(), dfunin, dfunout,
                   pair_endpoint1, pair_endpoint2, input_distances,
                   output_distances, thread_count);
  } else {
    Rcpp::RNGScope rng_scope;
    random_distances(pair_count, xin.ncol(), xin_cpp.begin(), xin_cpp.end(),
                     xout_cpp.begin(), xout_cpp.end(), dfunin, dfunout,
                     input_distances, output_distances, thread_count,
                     ret_pairs ? &returned_endpoint1 : nullptr,
                     ret_pairs ? &returned_endpoint2 : nullptr);
  }

  NumericVector res_in(input_distances.begin(), input_distances.end());
  NumericVector res_out(output_distances.begin(), output_distances.end());

  List result = List::create(_("din") = res_in, _("dout") = res_out);
  if (ret_pairs) {
    IntegerMatrix returned_pairs(pair_count, 2);
    for (std::size_t i = 0; i < pair_count; i++) {
      returned_pairs(i, 0) = returned_endpoint1[i];
      returned_pairs(i, 1) = returned_endpoint2[i];
    }
    colnames(returned_pairs) =
        CharacterVector::create("endpoint1", "endpoint2");
    result["pairs"] = returned_pairs;
  }

  return result;
}
