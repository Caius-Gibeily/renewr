functions {
#include include/functions.stan
}


data {
  // model settings
#include include/data.stan
}

transformed data {
#include include/transformed_data.stan
}

parameters {

  // population GP
  vector[M] z_ind;
  real<lower=0> rho_ind;
  real<lower=0> alpha_ind;

  // subject intercept
  real mu_ind;

  // renewal shape
  //real log_k_minus1;
  //vector<lower=1>[S] k;

  vector<lower=1>[include_k ? 1 : 0] k;
  vector<lower=1>[include_shape ? 1 : 0] shape;
  vector<lower=0>[include_sigma_lognormal ? 1 : 0] sigma_lognormal;


}

transformed parameters {


  //real k = 1.0 + exp(log_k_minus1);
  vector[M] diag_S_ind;

  diag_S_ind = get_diagSPD(alpha_ind,rho_ind,M,L,kernel);

  vector[M] beta_ind = diag_S_ind .* z_ind;

}

model {

  // Priors
  z_ind ~ std_normal();

  apply_prior_lp(mu_ind, distributions[1], params[1, 1], params[1, 2]);
  apply_prior_lp(alpha_ind, distributions[2], params[2, 1], params[2, 2]);
  apply_prior_lp(rho_ind, distributions[3], params[3, 1], params[3, 2]);

  if (include_k != 0) apply_prior_lp(k[1], distributions[include_k],
    params[include_k, 1], params[include_k, 2]);
  if (include_shape != 0) apply_prior_lp(shape[1], distributions[include_shape],
    params[include_shape, 1], params[include_shape, 2]);
  if (include_sigma_lognormal != 0) apply_prior_lp(sigma_lognormal[1], distributions[include_sigma_lognormal],
    params[include_sigma_lognormal, 1], params[include_sigma_lognormal, 2]);

  matrix[n_quad+1,N_total] eta_quad;

  for (j in 1:n_quad+1) {

    vector[N_total] f_ind = PHI_quad[j] * beta_ind;

    eta_quad[j] = mu_ind + f_ind';
  }

#include include/log-likelihood.stan

  target += sum(log_kernel);
}

