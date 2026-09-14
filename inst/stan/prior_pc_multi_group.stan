functions {
#include include/functions.stan
}

data {

  // Individual/group data
  int<lower=1> I;

  int<lower=1> G;

  array[I] int<lower=1,upper=G> g_membership;
  array[G] int I_per_group;

  int<lower=1> N_params;

  // prior metadata
  array[N_params, 2] real params;
  vector[N_params] distributions;
  int<lower=1, upper=4> kernel;
  int<lower=1, upper=5> family;

  // Hilbert space data
  int<lower=1> M;
  real<lower=0> L_factor;
  real<lower=1> duration;
  real<lower=0> w0;

  // Booleans for survival params
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
  matrix[G,G-1] Q_R = qr_decomp(G);
}

generated quantities {
  // Initialisation
  vector[M] z_global;
  matrix[G-1, M] z_group_raw;
  matrix[I, M] z_ind_raw;
  matrix[G, M] z_group;
  matrix[I, M] z_ind;
  real mu_global;
  real<lower=0> sigma_ind;
  real<lower=0> sigma_group;
  vector[G-1] mu_raw_group;
  vector[I] mu_raw_ind;
  real<lower=0> alpha_global;
  real<lower=0> alpha_group;
  real<lower=0> alpha_ind;
  real<lower=0> rho_global;
  real<lower=0> rho_group;
  real<lower=0> rho_ind;
  vector<lower=1>[include_k ? 1 : 0] k;
  vector<lower=1>[include_shape ? 1 : 0] shape;
  vector<lower=0>[include_sigma_lognormal ? 1 : 0] sigma_lognormal;
  vector[M] diag_S_global;
  vector[M] diag_S_group;
  vector[M] diag_S_ind;

  z_global = to_vector(normal_rng(rep_vector(0.0, M), 1.0));

  for (g in 1:(G - 1)) {
    z_group_raw[g, ] = to_row_vector(normal_rng(rep_vector(0.0, M), 1.0));
  }

  for (i in 1:I) {
    z_ind_raw[i, ] = to_row_vector(normal_rng(rep_vector(0.0, M), 1.0));
  }

  z_group = Q_R * z_group_raw;

  z_ind = constrain_groups(G, I, M, g_membership, z_ind_raw, I_per_group);


  mu_global = apply_prior_rng(distributions[1], params[1,1], params[1,2], -1);

  alpha_global = apply_prior_rng(distributions[2], params[2, 1], params[2, 2], 0);
  alpha_group = apply_prior_rng(distributions[3], params[3, 1], params[3, 2], 0);
  alpha_ind = apply_prior_rng(distributions[4], params[4, 1], params[4, 2], 0);
  rho_global = apply_prior_rng(distributions[5], params[5, 1], params[5, 2], 0);
  rho_group = apply_prior_rng(distributions[6], params[6, 1], params[6, 2], 0);
  rho_ind = apply_prior_rng(distributions[7], params[7, 1], params[7, 2], 0);

  if (include_k != 0) {
    k[1] = apply_prior_rng(distributions[include_k], params[include_k, 1], params[include_k, 2], 1);
  }
  if (include_shape != 0) {
    shape[1] = apply_prior_rng(distributions[include_shape], params[include_shape, 1], params[include_shape, 2], 1);
  }
  if (include_sigma_lognormal != 0) {
    sigma_lognormal[1] = apply_prior_rng(distributions[include_sigma_lognormal], params[include_sigma_lognormal, 1], params[include_sigma_lognormal, 2], 0);
  }

  // Hierarchical mu
  sigma_group = exponential_rng(1);
  sigma_ind = exponential_rng(1);


  mu_raw_group = to_vector(normal_rng(rep_vector(0.0, G - 1), 1.0));
  vector[G] mu_raw_group_std = Q_R * mu_raw_group;
  vector[G] mu_group = mu_global + sigma_group * mu_raw_group_std;


  mu_raw_ind = to_vector(normal_rng(rep_vector(0.0, I), 1.0));
  vector[I] mu_raw_ind_std = constrain_mu(I, G, g_membership, mu_raw_ind, I_per_group);
  vector[I] mu_ind = mu_group[g_membership] + sigma_ind * mu_raw_ind_std;

  diag_S_global = get_diagSPD(alpha_global, rho_global, M, L, kernel);
  diag_S_group  = get_diagSPD(alpha_group, rho_group, M, L, kernel);
  diag_S_ind    = get_diagSPD(alpha_ind, rho_ind, M, L, kernel);

  vector[M] beta_global = diag_S_global .* z_global;
  matrix[G, M] beta_group = diag_post_multiply(z_group, diag_S_group);
  matrix[I, M] beta_ind = diag_post_multiply(z_ind, diag_S_ind);
}
