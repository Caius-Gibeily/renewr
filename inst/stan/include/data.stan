
int<lower=1> M;
int<lower=1> N_params;
int<lower=3,upper=5>n_quad;

array[N_params,2] real params;
vector[N_params] distributions;
int<lower=1, upper=5> kernel;
int<lower=1, upper=5> family;

real<lower=0> L_factor;
real<lower=0> w0;

int<lower=0, upper=N_params> include_k;
int<lower=0, upper=N_params> include_shape;
int<lower=0, upper=N_params> include_sigma_lognormal;
// data
int<lower=1> N_total;
real<lower=0> duration;

vector[N_total] t_ev;
vector[N_total] dt;
array[N_total] int censored;
