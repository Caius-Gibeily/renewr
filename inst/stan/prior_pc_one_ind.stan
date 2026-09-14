functions {
#include include/functions.stan
}

data {
  int<lower=1> M;
  int<lower=1> I;
  int<lower=1> N_params;
  array[N_params, 2] real params;
  vector[N_params] distributions;
  int<lower=1, upper=4> kernel;
  int<lower=1, upper=5> family;

  real<lower=0> L_factor;
  real<lower=1> duration;
  real<lower=0> w0;

  int<lower=0, upper=N_params> include_k;
  int<lower=0, upper=N_params> include_shape;
  int<lower=0, upper=N_params> include_sigma_lognormal;
}

transformed data {
  real L = L_factor * duration / 2;
  vector[100] t_grid = linspaced_vector(100, 0, duration);
  matrix[100, M] PHI;
  if (kernel != 5) {
    PHI = phi(100, M, L, duration, t_grid);
  } else {
    PHI = phi(100, M, w0, duration, t_grid);
  }
}

generated quantities {

  vector[M] z_ind;

  real<lower=0> mu_ind;
  real<lower=0> sigma_ind;
  real<lower=0> mu_raw_ind;

  real<lower=0> alpha_ind;
  real<lower=0> rho_ind;

  vector<lower=1>[include_k ? 1 : 0] k;
  vector<lower=1>[include_shape ? 1 : 0] shape;
  vector<lower=0>[include_sigma_lognormal ? 1 : 0] sigma_lognormal;

  vector[M] diag_S_ind;

  for (m in 1:M) {
    z_ind[m] = std_normal_rng();
  }

  mu_ind = apply_prior_rng(distributions[1], params[1,1], params[1,2], -1);

  alpha_ind = apply_prior_rng(distributions[2], params[2, 1], params[2, 2], 0);
  rho_ind   = apply_prior_rng(distributions[3], params[3, 1], params[3, 2], 0);

  if (include_k != 0) {
    k[1] = apply_prior_rng(distributions[include_k], params[include_k, 1], params[include_k, 2], 1);
  }
  if (include_shape != 0) {
    shape[1] = apply_prior_rng(distributions[include_shape], params[include_shape, 1], params[include_shape, 2], 1);
  }
  if (include_sigma_lognormal != 0) {
    sigma_lognormal[1] = apply_prior_rng(distributions[include_sigma_lognormal], params[include_sigma_lognormal, 1], params[include_sigma_lognormal, 2], 0);
  }

  diag_S_ind = get_diagSPD(alpha_ind, rho_ind, M, L, kernel);


  vector[M] beta_ind = diag_S_ind .* z_ind;
}
