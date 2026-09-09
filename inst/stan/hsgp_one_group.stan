functions {
#include include/functions.stan
}


data {
#include include/data.stan
  int<lower=1> I;
  array[N_total] int<lower=1,upper=I> ind_id;

}

transformed data {
#include include/transformed_data.stan
  matrix[I,I-1] Q_R = qr_decomp(I);

}



parameters {

  // population GP
  vector[M] z_group;
  matrix[I-1,M] z_ind_raw;

  real<lower=0> rho_group;
  //real<lower=0> alpha_group;

  // Individual-level
  real<lower=0> rho_ind;
  //real<lower=0> alpha_ind;
  real<lower=0> alpha;
  real<lower=0,upper=1>omega;

  // indect intercept
  //vector[I] mu;
  real mu_group;
  real<lower=0> sigma_ind;
  vector[I-1] mu_raw_ind;
  // renewal shape
  //real log_k_minus1;
  //vector<lower=1>[S] k;

  vector<lower=1>[include_k ? 1 : 0] k;
  vector<lower=1>[include_shape ? 1 : 0] shape;
  vector<lower=0>[include_sigma_lognormal ? 1 : 0] sigma_lognormal;


}

transformed parameters {

  vector[I] mu_raw_ind_std = Q_R * mu_raw_ind; //mu_raw_ind_std = sum_to_zero_mu(I, 1, rep_array(1,I), mu_raw_ind, {I}); //mu_raw_ind_std = Q_R * mu_raw_ind;
  vector[I] mu_ind = mu_group + mu_raw_ind_std * sigma_ind;

  matrix[I,M] z_ind;

  real alpha_group = alpha * sqrt(omega);
  real alpha_ind = alpha * sqrt(1 - omega);
  //real k = 1.0 + exp(log_k_minus1);
  vector[M] diag_S_group;
  vector[M] diag_S_ind;

  diag_S_group = get_diagSPD(alpha_group,rho_group,M,L,kernel);
  diag_S_ind = get_diagSPD(alpha_ind,rho_ind,M,L,kernel);


  //z_ind[1:(I - 1), ] = z_ind_raw;
  z_ind = Q_R * z_ind_raw;
  //for (m in 1:M) {
  //  z_ind[I, m] = -sum(z_ind_raw[, m]);
  //}

  vector[M] beta_group = diag_S_group .* z_group;
  matrix[I, M] beta_ind = diag_post_multiply(z_ind,diag_S_ind);

}

model {

  // Priors
  z_group ~ std_normal();
  to_vector(z_ind_raw) ~ std_normal();

  mu_raw_ind ~ std_normal();
  sigma_ind ~ normal(0.5,0.2);

  omega ~ beta(6,4);

  apply_prior_lp(mu_group, distributions[1], params[1, 1], params[1, 2]);

  apply_prior_lp(alpha, distributions[2], params[2, 1], params[2, 2]);
  //apply_prior_lp(alpha_ind, distributions[3], params[3, 1], params[3, 2]);

  apply_prior_lp(rho_group, distributions[4], params[4, 1], params[4, 2]);
  apply_prior_lp(rho_ind, distributions[5], params[5, 1], params[5, 2]);


  if (include_k != 0) apply_prior_lp(k[1], distributions[include_k],
  params[include_k, 1], params[include_k, 2]);
  if (include_shape != 0) apply_prior_lp(shape[1], distributions[include_shape],
  params[include_shape, 1], params[include_shape, 2]);
  if (include_sigma_lognormal != 0) apply_prior_lp(sigma_lognormal[1], distributions[include_sigma_lognormal],
  params[include_sigma_lognormal, 1], params[include_sigma_lognormal, 2]);

  // GP instantiation
  //array[3] vector[N_total] eta_quad;
  matrix[n_quad+1,N_total] eta_quad;
  for (j in 1:n_quad+1) {
    vector[N_total] f_group = PHI_quad[j] * beta_group;
    vector[N_total] f_ind = rows_dot_product(PHI_quad[j], beta_ind[ind_id, ]);
    vector[N_total] f_aggregated = mu_ind[ind_id] + f_group + f_ind;
    eta_quad[j] = f_aggregated';
  }

  // Likelihood
#include include/log-likelihood.stan
  target += sum(log_kernel) + 1e-8;
}

generated quantities {

  // Log likelihood computation
  vector[N_total] log_lik;

  matrix[n_quad+1,N_total] eta_quad;

  for (j in 1:n_quad+1) {
    vector[N_total] f_group = PHI_quad[j] * beta_group;
    vector[N_total] f_ind = rows_dot_product(PHI_quad[j], beta_ind[ind_id, ]);

    vector[N_total] f_aggregated = mu_ind[ind_id] + f_group + f_ind;
    eta_quad[j] = f_aggregated';
  }

#include include/log-likelihood.stan
  log_lik = log_kernel;

}
