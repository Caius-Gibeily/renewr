// spectral basis set, phi
matrix phi(int N, int M, real L, vector x) {
  matrix[N, M] res;

  for (m in 1:M) {
    res[, m] = inv_sqrt(L) * sin(pi() * m * (x + L) / (2 * L));
  }

  return res;
}

matrix phi_periodic(int N, int M, real w0, vector x) {

  row_vector[M / 2] k = linspaced_row_vector(M / 2, 1, M / 2);

  matrix[N, M / 2] w0xk = (w0 * x) * k;
  return append_col(cos(w0xk), sin(w0xk));
}

// Basis indices
vector basis_indices(int M,real L) {
  vector[M] indices = linspaced_vector(M, 1, M);

  return square(pi() / (2 * L) * indices);
}

// Matern 1/2 - sourced from
// https://epiforecasts.io/EpiNow2/stan/gaussian__process_8stan_source.html#l00090

vector diagSPD_Matern12(real alpha, real rho, int M, real L) {
  vector[M] denom = 1 / rho + rho * basis_indices(M, L);
  return alpha * sqrt(2 ./ denom);
}
// Matern 3/2
vector diagSPD_Matern32(real alpha, real rho, int M, real L) {

  real factor = 2 * alpha * (sqrt(3) / rho)^1.5;
  vector[M] denom = 3 / square(rho) + basis_indices(M, L);

  return factor ./ denom;
}
// Matern 5/2
vector diagSPD_Matern52(real alpha, real rho, int M, real L) {
  real factor = 16 * pow(sqrt(5) / rho, 5);
  vector[M] denom = 3 * pow(5 / square(rho) + basis_indices(M, L), 3);
  return alpha * sqrt(factor ./ denom);
}

// Periodic
vector diagSPD_Periodic(real alpha, real rho, int M) {

  real a = inv_square(rho);
  vector[M / 2] indices = linspaced_vector(M / 2, 1, M / 2);
  vector[M / 2] q = exp(
    log(alpha) + 0.5 *
    (log(2) - a + to_vector(log_modified_bessel_first_kind(indices, a)))
    );
    return append_row(q,q);
}

vector diagSPD_EQ(real alpha, real rho, int M, real L) {
  vector[M] indices = linspaced_vector(M, 1, M);
  real factor = alpha * sqrt(sqrt(2 * pi()) * rho);
  real exponent = -0.25 * (rho * pi() / 2 / L)^2;
  return factor * exp(exponent * square(indices));
}

real gengamma_lpdf(real x, real shape, real k, real scale) {

  if (x <= 0.0) {
    return negative_infinity();
  }
  real logdens = log(shape) - lgamma(k) + (shape * k - 1.0) * log(x) - (shape * k) * log(scale) - pow(x/scale,shape);
  return logdens;

}

vector get_diagSPD(real alpha, real rho, int M, real L, int kernel) {

  vector[M] diag_S;
  if (kernel == 1) {
    diag_S = diagSPD_EQ(alpha,rho,M,L);
  } else if (kernel == 2) {
    diag_S = diagSPD_Matern12(alpha,rho,M,L);
  } else if (kernel == 3) {
    diag_S = diagSPD_Matern32(alpha,rho,M,L);
  } else if (kernel == 4) {
    diag_S = diagSPD_Matern52(alpha,rho,M,L);
  } else if (kernel == 5) {
    diag_S = diagSPD_Periodic(alpha,rho,M);
  }
  return diag_S;
}

matrix sum_to_zero_groups(int G, int I, int M, array[] int g_membership, matrix z_ind_raw, array[] int I_per_group) {
  matrix[I, M] z_ind;

  int pos = 1;
  for (g in 1:G) {

    int n = I_per_group[g];

    int idx = 0;
    int last_i = 0;


    for (i in 1:I) {
      if (g_membership[i] == g) {
        idx += 1;
        if (idx < n) {
          z_ind[i, ] = z_ind_raw[pos, ];
          pos += 1;
        } else {
          last_i = i;
        }
      }
    }

    // Constrained row
    for (m in 1:M) {
      real s = 0;

      for (i in 1:I) {
        if (g_membership[i] == g && i != last_i)
        s += z_ind[i, m];
      }

      z_ind[last_i, m] = -s;
    }
  }

  return z_ind;
}

matrix constrain_groups(int G, int I, int M, array[] int g_membership, matrix z_ind_raw, array[] int I_per_group) {
  matrix[I, M] z_ind;
  int pos = 1;

  for (g in 1:G) {
    int n = I_per_group[g];

    matrix[n, M] group_raw;
    for (row_idx in 1:n) {
      group_raw[row_idx, ] = z_ind_raw[pos, ];
      pos += 1;
    }

    real scale_factor = sqrt(n / (n - 1.0));
    matrix[n, M] group_centered;

    for (m in 1:M) {
      real col_mean = mean(group_raw[, m]);
      group_centered[, m] = (group_raw[, m] - col_mean) * scale_factor;
    }

    int idx = 1;
    for (i in 1:I) {
      if (g_membership[i] == g) {
        z_ind[i, ] = group_centered[idx, ];
        idx += 1;
      }
    }
  }

  return z_ind;
}


vector sum_to_zero_mu(int l0, int l1, array[] int l1_membership, vector mu_raw, array[] int l0_per_l1) {
  vector[l0] mu_constrained;
  int pos = 1;

  for (l in 1:l1) {
    int n = l0_per_l1[l];
    int idx = 0;
    int last_m = 0;
    real group_sum = 0.0;

    for (m in 1:l0) {
      if (l1_membership[m] == l) {
        idx += 1;
        if (idx < n) {
          mu_constrained[m] = mu_raw[pos];
          group_sum += mu_raw[pos];
          pos += 1;
        } else {
          last_m = m;
        }
      }
    }

    mu_constrained[last_m] = -group_sum;
  }

  return mu_constrained;
}

matrix qr_decomp(int N) {
  matrix[N, N - 1] Q_R;
  for (i in 1:(N - 1)) {
    for (j in 1:N) {
      if (j <= i) {
        Q_R[j, i] = 1.0 / sqrt(i * (i + 1.0));
      } else if (j == i + 1) {
        Q_R[j, i] = -i / sqrt(i * (i + 1.0));
      } else {
        Q_R[j, i] = 0.0;
      }
    }
  }
  return Q_R;
}

vector constrain_mu(int l0, int l1, array[] int l1_membership, vector mu_raw, array[] int l0_per_l1) {
  vector[l0] mu_constrained;
  int pos = 1;

  for (l in 1:l1) {
    int n = l0_per_l1[l];

    vector[n] group_raw;
    for (i in 1:n) {
      group_raw[i] = mu_raw[pos];
      pos += 1;
    }

    real group_mean = mean(group_raw);

    real scale_factor = sqrt(n / (n - 1.0));
    vector[n] group_centered = (group_raw - group_mean) * scale_factor;

    int idx = 1;
    for (m in 1:l0) {
      if (l1_membership[m] == l) {
        mu_constrained[m] = group_centered[idx];
        idx += 1;
      }
    }
  }

  return mu_constrained;
}


void apply_prior_lp(real param, real dist, real arg1, real arg2) {
  if (dist == 1) {
    target += normal_lpdf(param | arg1, arg2);
  } else if (dist == 2) {
    target += lognormal_lpdf(param | arg1, arg2);
  } else if (dist == 3) {
    target += cauchy_lpdf(param | arg1, arg2);
  } else if (dist == 4) {
    target += exponential_lpdf(param | arg1) - exponential_lccdf(1.0 | arg1);
  }
}

real apply_prior_rng(real dist, real arg1, real arg2, int apply_floor) {
  real param = -1;

  if (apply_floor == 1 || apply_floor == 0) {
    real floor_value = (apply_floor == 1) ? 1.0 : 0.0;
    while (param < floor_value) {
      if (dist == 1) {
        param = normal_rng(arg1, arg2);
      } else if (dist == 2) {
        param = lognormal_rng(arg1, arg2);
      } else if (dist == 3) {
        param = cauchy_rng(arg1, arg2);
      } else if (dist == 4) {
        param = exponential_rng(arg1);
      }
    }
  } else {
    if (dist == 1) {
      param = normal_rng(arg1, arg2);
    } else if (dist == 2) {
      param = lognormal_rng(arg1, arg2);
    } else if (dist == 3) {
      param = cauchy_rng(arg1, arg2);
    } else if (dist == 4) {
      param = exponential_rng(arg1);
    }
  }

  return param;
}

matrix exponential_likelihood(int N_total, vector log_qw, array[] vector eta_quad, vector dt) {

  matrix[N_total, 3] log_kernel;

  for (j in 1:3) {
    vector[N_total] scale = exp(eta_quad[j]);
    for (n in 1:N_total) {
      log_kernel[n, j] = log_qw[j] + exponential_lpdf(dt[n] | inv(scale[n]));
    }
  }
  return log_kernel;
}

matrix gamma_likelihood(int N_total, vector log_qw, array[] vector eta_quad,
  vector dt, real k) {

  matrix[N_total, 3] log_kernel;

  for (j in 1:3) {
    vector[N_total] rate = exp(-eta_quad[j]);
    for (n in 1:N_total) {
        log_kernel[n, j] = log_qw[j] + gamma_lpdf(dt[n] | k, rate[n]);
    }
  }
  return log_kernel;
}

matrix weibull_likelihood(int N_total, vector log_qw, array[] vector eta_quad, vector dt, real shape) {

  matrix[N_total, 3] log_kernel;

  for (j in 1:3) {
    vector[N_total] scale = exp(eta_quad[j]);
    for (n in 1:N_total) {

        log_kernel[n, j] = log_qw[j] + weibull_lpdf(dt[n] | shape, scale[n]);

    }
  }
  return log_kernel;
}

matrix lognormal_likelihood(int N_total, vector log_qw, array[] vector eta_quad, vector dt, real sigma_lognormal) {

  matrix[N_total, 3] log_kernel;

  for (j in 1:3) {
    vector[N_total] mu_lognormal = eta_quad[j];
    for (n in 1:N_total) {

        log_kernel[n, j] = log_qw[j] + lognormal_lpdf(dt[n] | mu_lognormal[n], sigma_lognormal);

    }
  }
  return log_kernel;
}

matrix gengamma_likelihood(int N_total, vector log_qw, array[] vector eta_quad, vector dt, real k, real shape) {

  matrix[N_total, 3] log_kernel;

  for (j in 1:3) {
    vector[N_total] scale = exp(eta_quad[j]);
    for (n in 1:N_total) {
      log_kernel[n, j] = log_qw[j] + gengamma_lpdf(dt[n] | shape, k, scale[n]);
    }
  }
  return log_kernel;
}
real get_rate_t(int ind, vector mu_ind, vector beta_ind_i,
  vector beta_group, vector t, int M, real L, real w0, int kernel) {

  matrix[1, M] PHI;
  if (kernel != 5) {
    PHI = phi(1, M, L, t);
  } else {
    PHI = phi_periodic(1, M, w0, t);
  }

  vector[1] f_group = PHI * beta_group;
  vector[1] f_ind = PHI * beta_ind_i;
  vector[1] eta_t = mu_ind[ind] + f_group + f_ind;

  return exp(-eta_t[1]); // inv(exp(x)) is exp(-x)
}

// Numerically stable hazard functions using exp(lpdf - lccdf)
real exponential_h(real dt, real rate) {
  return exp(exponential_lpdf(dt | rate) - exponential_lccdf(dt | rate));
}

real gamma_h(real dt, real shape, real rate) {
  return exp(gamma_lpdf(dt | shape, rate) - gamma_lccdf(dt | shape, rate));
}

real weibull_h(real dt, real shape, real scale) {
  return exp(weibull_lpdf(dt | shape, scale) - weibull_lccdf(dt | shape, scale));
}

real lognormal_h(real dt, real mu_lognormal, real sigma_lognormal) {
  return exp(lognormal_lpdf(dt | mu_lognormal, sigma_lognormal) - lognormal_lccdf(dt | mu_lognormal, sigma_lognormal));
}

