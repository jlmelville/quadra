#include <vector>

#include <Rcpp.h>

#include "pforr.h"

#include "rnndescent/random.h"

#include "distance.h"
#include "native-validation.h"

using namespace Rcpp;

void random_distances(int n_obs, const DistanceSpace &input,
                      const DistanceSpace &output,
                      std::vector<double> &input_distances,
                      std::vector<double> &output_distances,
                      std::size_t n_threads, std::vector<int> *endpoint1,
                      std::vector<int> *endpoint2) {
  rnndescent::ParallelIntRNGAdapter<uint64_t, rnndescent::DQIntSampler>
      sampler_provider;
  sampler_provider.initialize();

  auto worker = [&](std::size_t begin, std::size_t end, std::size_t chunk_id) {
    auto sampler = sampler_provider.get_parallel_instance(chunk_id);
    for (std::size_t i = begin; i < end; ++i) {
      // DQIntSampler draws distinct indices without replacement.
      const auto idxs = sampler->sample(n_obs, 2);
      if (endpoint1 != nullptr) {
        (*endpoint1)[i] = static_cast<int>(idxs[0] + 1);
        (*endpoint2)[i] = static_cast<int>(idxs[1] + 1);
      }
      input_distances[i] = input.between(idxs[0], idxs[1]);
      output_distances[i] = output.between(idxs[0], idxs[1]);
    }
  };
  pforr::parallel_for_indexed(input_distances.size(), worker, n_threads);
}

void pair_distances(const std::vector<std::size_t> &endpoint1,
                    const std::vector<std::size_t> &endpoint2,
                    const DistanceSpace &input, const DistanceSpace &output,
                    std::vector<double> &input_distances,
                    std::vector<double> &output_distances,
                    std::size_t n_threads) {
  auto worker = [&](std::size_t begin, std::size_t end, std::size_t) {
    for (std::size_t i = begin; i < end; ++i) {
      input_distances[i] = input.between(endpoint1[i], endpoint2[i]);
      output_distances[i] = output.between(endpoint1[i], endpoint2[i]);
    }
  };
  pforr::parallel_for_indexed(endpoint1.size(), worker, n_threads);
}

// [[Rcpp::export(rng = false)]]
List random_distances(const NumericMatrix &xin, const NumericMatrix &xout,
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
  int pair_count = 0;
  std::vector<std::size_t> pair_endpoint1;
  std::vector<std::size_t> pair_endpoint2;
  std::vector<int> returned_endpoint1;
  std::vector<int> returned_endpoint2;

  if (has_pairs) {
    IntegerMatrix pair_matrix(pairs);
    if (pair_matrix.nrow() < 1 || pair_matrix.ncol() != 2) {
      stop("pairs must have at least one row and exactly two columns");
    }
    pair_count = pair_matrix.nrow();
    pair_endpoint1.resize(pair_count);
    pair_endpoint2.resize(pair_count);
    if (ret_pairs) {
      returned_endpoint1.resize(pair_count);
      returned_endpoint2.resize(pair_count);
    }

    for (int i = 0; i < pair_count; i++) {
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

  const DistanceSpace input(xin, metric_in), output(xout, metric_out);

  std::vector<double> input_distances(pair_count);
  std::vector<double> output_distances(pair_count);

  if (has_pairs) {
    pair_distances(pair_endpoint1, pair_endpoint2, input, output,
                   input_distances, output_distances, thread_count);
  } else {
    Rcpp::RNGScope rng_scope;
    random_distances(xin.ncol(), input, output, input_distances,
                     output_distances, thread_count,
                     ret_pairs ? &returned_endpoint1 : nullptr,
                     ret_pairs ? &returned_endpoint2 : nullptr);
  }

  NumericVector res_in(input_distances.begin(), input_distances.end());
  NumericVector res_out(output_distances.begin(), output_distances.end());

  List result = List::create(_("din") = res_in, _("dout") = res_out);
  if (ret_pairs) {
    IntegerMatrix returned_pairs(pair_count, 2);
    for (int i = 0; i < pair_count; i++) {
      returned_pairs(i, 0) = returned_endpoint1[i];
      returned_pairs(i, 1) = returned_endpoint2[i];
    }
    colnames(returned_pairs) =
        CharacterVector::create("endpoint1", "endpoint2");
    result["pairs"] = returned_pairs;
  }

  return result;
}
