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

  matrix[I,I-1] Q_R = qr_decomp(I);
}

generated quantities {

  vector[M] z_group;
  matrix[I-1, M] z_ind_raw;
  matrix[I, M] z_ind;

  real mu_group;
  real<lower=0> sigma_ind;
  vector[I-1] mu_raw_ind;

  real<lower=0> alpha_group;
  real<lower=0> alpha_ind;
  real<lower=0> rho_group;
  real<lower=0> rho_ind;

  vector<lower=1>[include_k ? 1 : 0] k;
  vector<lower=1>[include_shape ? 1 : 0] shape;
  vector<lower=0>[include_sigma_lognormal ? 1 : 0] sigma_lognormal;

  vector[M] diag_S_group;
  vector[M] diag_S_ind;


  mu_group = apply_prior_rng(distributions[1], params[1,1], params[1,2], -1);
  sigma_ind = apply_prior_rng(1, 0.5, 0.2, 0);

  for (i in 1:I-1) {
    mu_raw_ind[i] = std_normal_rng();
  }

  vector[I] mu_raw_ind_std = Q_R * mu_raw_ind; //constrain_mu(I, 1, rep_array(1, I), mu_raw_ind, {I});
  vector[I] mu_ind = mu_group + mu_raw_ind_std * sigma_ind;

  alpha_group = apply_prior_rng(distributions[2], params[2, 1], params[2, 2],0);
  alpha_ind = apply_prior_rng(distributions[3], params[3, 1], params[3, 2], 0);
  rho_group = apply_prior_rng(distributions[4], params[4, 1], params[4, 2], 0);
  rho_ind   = apply_prior_rng(distributions[5], params[5, 1], params[5, 2], 0);

  if (include_k != 0) {
    k[1] = apply_prior_rng(distributions[include_k], params[include_k, 1], params[include_k, 2], 1);
  }
  if (include_shape != 0) {
    shape[1] = apply_prior_rng(distributions[include_shape], params[include_shape, 1], params[include_shape, 2], 1);
  }
  if (include_sigma_lognormal != 0) {
    sigma_lognormal[1] = apply_prior_rng(distributions[include_sigma_lognormal], params[include_sigma_lognormal, 1], params[include_sigma_lognormal, 2], 0);
  }

  diag_S_group = get_diagSPD(alpha_group, rho_group, M, L, kernel);
  diag_S_ind   = get_diagSPD(alpha_ind, rho_ind, M, L, kernel);

  for (m in 1:M) {
    z_group[m] = std_normal_rng();
  }

  for (i in 1:I-1) {
    for (m in 1:M) {
      z_ind_raw[i, m] = std_normal_rng();
    }
  }

  z_ind = Q_R * z_ind_raw; //constrain_groups(G, I, M, rep_array(1, I), z_ind_raw, {I});
  //z_ind[1:(I - 1), ] = z_ind_raw;
  //for (m in 1:M) {
  //  z_ind[I, m] = -sum(z_ind_raw[, m]);
  //}

  vector[M] beta_group = diag_S_group .* z_group;
  matrix[I, M] beta_ind = diag_post_multiply(z_ind, diag_S_ind);
}
