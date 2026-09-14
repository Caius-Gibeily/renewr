#include <Rcpp.h>
#include <RcppParallel.h>
#include <random>
#include <cmath>
#include <vector>
#include "gengamma_orig.h"

using namespace Rcpp;
using namespace RcppParallel;

struct ThreadData {
  std::vector<int> sample;
  std::vector<int> ind;
  std::vector<int> group;
  std::vector<double> event_times;
  std::vector<double> dt;
};

struct FlexibleRenewalSimulator : public Worker {
  const std::vector<double>& time_vec;
  const double* flat_modulant_ptr;
  const IntegerVector& groups_vec;
  const NumericVector& shape_vec;
  const NumericVector& k_vec;
  const double t_diff;
  const int n_time;
  const int n_ind;
  const double max_x;
  const bool use_samples;
  const double start_time;

  std::vector<ThreadData>& thread_storage;

  FlexibleRenewalSimulator(const std::vector<double>& time_vec,
                           const double* flat_modulant_ptr,
                           const IntegerVector& groups_vec,
                           const NumericVector& shape_vec,
                           const NumericVector& k_vec,
                           int n_time, int n_ind, double t_diff, double max_x,
                           bool use_samples, std::vector<ThreadData>& thread_storage,
                           double start_time)
    : time_vec(time_vec), flat_modulant_ptr(flat_modulant_ptr), groups_vec(groups_vec),
      shape_vec(shape_vec), k_vec(k_vec), t_diff(t_diff), n_time(n_time),
      n_ind(n_ind), max_x(max_x), use_samples(use_samples), thread_storage(thread_storage),
      start_time(start_time) {}

  void operator()(std::size_t begin, std::size_t end) {
    gengamma_orig::density gen_pdf;
    gengamma_orig::cdf gen_cdf;

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (std::size_t task_id = begin; task_id < end; task_id++) {
      int s_idx = use_samples ? (task_id / n_ind) : 0;
      int i_idx = use_samples ? (task_id % n_ind) : task_id;

      ThreadData& local = thread_storage[task_id];
      local.ind.reserve(200);
      local.group.reserve(200);
      local.event_times.reserve(200);
      local.dt.reserve(200);
      if (use_samples) { local.sample.reserve(200); }

      double current_shape = shape_vec[s_idx];
      double current_k = k_vec[s_idx];
      int current_group = groups_vec[i_idx];

      double last_event = start_time;

      double u = std::max(dist(rng), 1e-15);
      double H_threshold = -std::log(u);
      double H_level = 0;
      double current_scale;
      double dt_val = 0;

      std::size_t task_index = use_samples ? (s_idx * n_ind + i_idx) : i_idx;
      std::size_t time_series_offset = task_index * n_time;

      int warmup = 1;

      int ti_start = static_cast<int>(std::round(start_time / t_diff));

      for (int ti = ti_start; ti < n_time; ti++) {

        double t_current;   // real elapsed time this iteration represents

        if (ti < 0) {
          warmup = 1;
          t_current = ti * t_diff;
        } else {
          warmup = 0;
          t_current = time_vec[ti];
        }

        dt_val = t_current - last_event;

        if (dt_val < 0) { continue; }

        if (warmup == 1) {
          current_scale = flat_modulant_ptr[time_series_offset];
        } else {
          current_scale = flat_modulant_ptr[time_series_offset + ti];
        }

        double pdf = gen_pdf(dt_val, current_scale, current_shape, current_k);
        double cdf = gen_cdf(dt_val, current_scale, current_shape, current_k);

        double surv = 1.0 - cdf;
        if (surv < 1e-12) surv = 1e-12;

        double h_t = t_diff * pdf / surv;
        H_level += h_t;

        if (H_level >= H_threshold) {
          double ev_time = t_current;

          if (ev_time > 0) {
            if (use_samples) { local.sample.push_back(s_idx + 1); }
            local.ind.push_back(i_idx + 1);
            local.group.push_back(current_group);
            local.event_times.push_back(ev_time);
            local.dt.push_back(ev_time - last_event);
          }

          last_event = ev_time;
          H_level = 0;

          double u_next = std::max(dist(rng), 1e-15);
          H_threshold = -std::log(u_next);
        }
      }

      if (use_samples) { local.sample.push_back(s_idx + 1); }
      local.ind.push_back(i_idx + 1);
      local.group.push_back(current_group);
      local.event_times.push_back(max_x);
      local.dt.push_back(max_x - last_event);
    }
  }
};

// [[Rcpp::export]]
DataFrame simulate_renewal_flexible(
    std::vector<double> time_vec,
    NumericVector modulant_mat_flat,
    IntegerVector groups_vec,
    NumericVector shape_vec,
    NumericVector k_vec,
    int n_ind,
    int n_samples,
    double max_x,
    bool use_samples,
    double start_time
) {
  double t_diff = time_vec[1] - time_vec[0];
  int n_time = time_vec.size();

  // Total parallel tasks changes dynamically depending on whether samples are active
  int total_tasks = use_samples ? (n_ind * n_samples) : n_ind;

  const double* flat_ptr = modulant_mat_flat.begin();
  std::vector<ThreadData> thread_storage(total_tasks);

  FlexibleRenewalSimulator simulator(time_vec, flat_ptr, groups_vec, shape_vec, k_vec,
                                     n_time, n_ind, t_diff, max_x, use_samples, thread_storage, start_time);

  parallelFor(0, total_tasks, simulator);

  std::size_t total_rows = 0;

  for (int i = 0; i < total_tasks; i++) {
    total_rows += thread_storage[i].event_times.size();
  }

  IntegerVector final_sample(use_samples ? total_rows : 0);
  IntegerVector final_ind(total_rows);
  IntegerVector final_group(total_rows);
  NumericVector final_event_times(total_rows);
  NumericVector final_dt(total_rows);

  std::size_t write_offset = 0;
  for (int i = 0; i < total_tasks; i++) {
    std::size_t sz = thread_storage[i].event_times.size();
    if (sz == 0) continue;

    if (use_samples) {
      std::copy(thread_storage[i].sample.begin(), thread_storage[i].sample.end(), final_sample.begin() + write_offset);
    }
    std::copy(thread_storage[i].ind.begin(), thread_storage[i].ind.end(), final_ind.begin() + write_offset);
    std::copy(thread_storage[i].group.begin(), thread_storage[i].group.end(), final_group.begin() + write_offset);
    std::copy(thread_storage[i].event_times.begin(), thread_storage[i].event_times.end(), final_event_times.begin() + write_offset);
    std::copy(thread_storage[i].dt.begin(), thread_storage[i].dt.end(), final_dt.begin() + write_offset);

    write_offset += sz;
  }

  // Conditionally include the "sample" column in the output table
  if (use_samples) {
    return DataFrame::create(
      _["sample"]      = final_sample,
      _["ind"]         = final_ind,
      _["group"]       = final_group,
      _["event_times"] = final_event_times,
      _["dt"]          = final_dt
    );
  } else {
    return DataFrame::create(
      _["ind"]         = final_ind,
      _["group"]       = final_group,
      _["event_times"] = final_event_times,
      _["dt"]          = final_dt
    );
  }
}
