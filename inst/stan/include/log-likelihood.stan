vector[N_total] log_kernel;

if (family == 1) log_kernel = exponential_likelihood(N_total, n_quad, qw, eta_quad,
  dt_quad, censored);
else if (family == 2) log_kernel = gamma_likelihood(N_total, n_quad, qw, eta_quad,
  dt_quad, k[1], censored);
else if (family == 3) log_kernel = weibull_likelihood(N_total, n_quad, qw, eta_quad,
  dt_quad, shape[1], censored);
else if (family == 4) log_kernel = lognormal_likelihood(N_total, n_quad, qw, eta_quad,
  dt_quad, sigma_lognormal[1], censored);
//else if (family == 5) log_kernel = gengamma_likelihood(N_total, qw, eta_quad, dt, k[1], shape[1]);

