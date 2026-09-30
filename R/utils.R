#' @noRd
.phi2 <- function(x, M, L) {
  # major bug: L was being evaluated on [0, S] domain where S is the duration
  # the formula was developed for [-L, L] domain centred on 0.
  #L is specified relative to the half-duration L <- L/2
  outer(x, seq_len(M), function(x,m)
    sin(pi * m * (x + L - max(x)/2) / (2 * L)) / sqrt(L)
  )
}


#' @noRd
.get_M_bases <- function(duration, L_factor, rho_list,
                         kernel = "squared_exp") {
  idx <- match(kernel, c("squared_exp", "matern12", "matern32", "matern52"))

  if (idx == "matern12") {
    warning("No heuristic available for Mátern 1/2 kernel. Setting M = 100
            but it is recommended that you experiment with alternative numbers of M bases.
            Please check ?gp_fit for further guidance.")
    M <- 100
    return(M)
  } else {

    scale_factors <- c(1.75, NA, 3.42, 2.65)
    L_limits <- c(3.2, NA, 4.5, 4.1)
    S <- duration / 2
    if (L_factor >= 1.2 | L_factor >= L_limits[idx] * min(rho_list,na.rm = TRUE) / S) {
      M <- scale_factors[idx] * L_factor * S / min(rho_list,na.rm = TRUE)
    } else {
      stop("Please choose a higher L_factor for stability.")
    }
    return(ceiling(M))
  }
}
