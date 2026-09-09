#include <Rcpp.h>
#include <omp.h>
#include <random>
#include "gengamma_orig.h"

using namespace Rcpp;

// [[Rcpp::export]]
List simulate_renewal_multi_omp(
    std::vector<double> time_vec,
    NumericVector modulant_mat_flat,
    int n_ind,
    double shape,
    double k
) {
  // Create an Rcpp List to store the results
  List all_results(n_ind);

  double t_diff = time_vec[1] - time_vec[0];
  int n_time = time_vec.size();

  // Fix 1: Pull the raw C++ pointer out to handle multi-threaded indexing
  const double* raw_modulant = modulant_mat_flat.begin();

#pragma omp parallel for schedule(static)
  for (int i = 0; i < n_ind; i++) {

    // Pinpoint the direct starting memory location for this row
    const double* ind_modulant = &raw_modulant[i * n_time];

    std::vector<double> event_times;
    event_times.reserve(1000);

    double last_event = 0;
    event_times.push_back(last_event);

    // Fix 2: Modern C++ thread-safe Random Number Generator
    // This instantiates a local, isolated engine for each individual loop (thread)
    std::mt19937 rng(i + 1); // Seed it with a unique index per individual
    std::uniform_real_distribution<double> dist(0.0, 1.0); // Generates uniform numbers in [0,1)

    // Draw the initial exponential threshold using the Inverse CDF method: -log(1 - u)
    // Which is mathematically identical to -log(u) for u in (0, 1]
    double u = std::max(dist(rng), 1e-15); //
    double H_threshold = -log(u);

    double H_level = 0;

    gengamma_orig::density gen_pdf;
    gengamma_orig::cdf gen_cdf;

    for (int ti = 0; ti < n_time; ti++) {
      double dt = time_vec[ti] - last_event;
      if (dt < 0) {
        continue;
      }

      double current_scale = ind_modulant[ti];

      double pdf = gen_pdf(dt, current_scale, shape, k);
      double cdf = gen_cdf(dt, current_scale, shape, k);
      double surv = std::max(1.0 - cdf, 1e-15);

      double h_t = t_diff * pdf / surv;
      H_level += h_t;

      if (H_level >= H_threshold) {
        event_times.push_back(time_vec[ti]);
        last_event = time_vec[ti];
        H_level = 0;

        // Draw next thread-safe random threshold
        double u_next = std::max(dist(rng), 1e-15); //
        H_threshold = -log(u_next);
      }
    }

    // Protect R's global object assignment context from simultaneous cross-thread writes
#pragma omp critical
{
  all_results[i] = wrap(event_times);
}
  }

  return all_results;
}
