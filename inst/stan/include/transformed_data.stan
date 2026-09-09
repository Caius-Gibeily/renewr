
real L = L_factor * duration;

// n-point Gauss-Legendre quadrature
vector[n_quad+1] qx = get_qx(n_quad);
vector[n_quad] qw = get_qw(n_quad);

//vector[3] log_qw = log(qw);

// precompute basis matrices
array[n_quad+1] matrix[N_total, M] PHI_quad;
matrix[N_total, M] PHI_end;
// dt / 2 * xi + (2*t_ev + dt) / 2
matrix[n_quad+1,N_total] dt_quad;
vector[N_total] t_quad;

if (kernel != 5) {
  for (j in 1:n_quad+1) {
    dt_quad[j] = dt' * qx[j]; //t_ev - (1.0 - qx[j]) .* dt;
    t_quad = t_ev - dt + dt_quad[j]';
    PHI_quad[j] = phi(N_total,M,L,t_quad);
  }

} else {
  for (j in 1:n_quad+1) {
    dt_quad[j] = dt' * qx[j]; //t_ev - (1.0 - qx[j]) .* dt;
    t_quad = t_ev - dt + dt_quad[j]';
    PHI_quad[j] = phi_periodic(N_total,M,w0,t_quad);
  }

}
