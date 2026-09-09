simulate_gp_traces <- function(n_ind = 5, n_group = 1, duration = 120, alpha = 0.8,
                               mu_ind = 0, mu_group = 0, mu_global = 0, sigma_ind = 0.3, sigma_group = 0.5,
                               omega = 0.8, alpha_ind, alpha_group, alpha_global, rho_ind = 10,
                               rho_group = 5 / 4 * rho_ind, rho_global = 5 / 4 * rho_group,
                               kernel = "squared_exp", L_factor = 1.2, M, seed = NULL) {
  if (length(n_ind) > 1 & is.numeric(n_ind)) {
    n_group <- length(n_ind)
    type <- "multi_group"
  } else {
    n_ind <- rep(n_ind, n_group)
  }

  if (n_group == 1) {
    if (length(n_ind) == 1 & n_ind == 1) {
      type <- "one_ind"
    } else {
      type <- "one_group"
    }
  } else type <- "multi_group"

  if (missing(M)) {
    if (type == "one_ind") {
      rho_list <- rho_ind
    } else if (type == "one_group") {
      rho_list <- c(rho_ind, rho_group)
    } else {
      rho_list <- c(rho_ind, rho_group, rho_global)
    }
    M <- .get_M_bases(duration, L_factor, rho_list, kernel)
  }
  if (missing(omega)) {
    omega <- switch(type,
      "one_ind" = 1,
      "one_group" = c(1/4, 3/4),
      "multi_group" = c(1/6, 2/6, 3/6)
    )
  }

  if ((type == "one_ind" & missing(alpha_ind)) |
    (type == "one_group" & missing(alpha_ind) & missing(alpha_group)) |
    (type == "multi_group" & missing(alpha_ind) & missing(alpha_group) &
      missing(alpha_global))) {
    if (type == "one_ind") {
      alpha_ind <- alpha
    } else if (type == "one_group") {
      if (length(omega) > 2) {
        warning("Using only the first two elements to assign alpha_ind and alpha_group")
        omega <- omega[1:2] / sum(omega[1:2])
      } else if (length(omega) == 1) {
        if (omega > 1 | omega < 0) stop("Please ensure omega is >= 0 or <= 1")
        omega <- c(1 - omega, omega)
      }

      if (round(sum(omega),5) != 1) {
        warning("Sum of omega weights does not sum to 1. Normalising weights.")
        print(omega)
        omega <- omega / sum(omega)
      }
      alpha_group <- alpha * sqrt(omega[2])
      alpha_ind <- alpha * sqrt(omega[1])
    } else if (type == "multi_group") {
      if (length(omega) < 3) {
        stop("Please specify three weights for a multi-group simulation or
                    set alpha_ind, alpha_group and alpha_global separately")
      }
      if (sum(omega) != 1) {
        warning("Sum of omega weights does not sum to 1. Normalising weights.")
        omega <- omega / sum(omega)
      }
      alpha_global <- omega[3]
      alpha_group <- omega[2]
      alpha_ind <- omega[1]
    }
  }

  L <- duration * L_factor
  t_grid <- seq(0, duration, by = 1)
  phi <- .phi(t_grid, M, L)

  sim_data <- list()

  diagS_ind <- .get_diagSPD(
    alpha_ind, rho_ind,
    M, L, kernel
  )
  if (missing(alpha_ind)) stop("Please specify an alpha_ind or set alpha")
  if (type == "one_ind") {
    alpha_group <- alpha_global <- rho_group <- rho_global <- mu_group <- mu_global <- NULL

    z_ind <- rnorm(M)
    beta_ind <- z_ind * diagS_ind
    f_ind <- phi %*% beta_ind

    eta_ind <- mu_ind + beta_ind

    sim_data$traces <- list(
      ind_traces = .format_df(n_ind, t_grid, eta_ind,level="ind"),
      type = type
    )

  } else if (type == "one_group") {
    if (any(missing(alpha_ind),missing(alpha_group))) {
      stop("Please specify an alpha_ind and alpha_group or set alpha and omega to partition alpha between individual and group levels")
    }
    alpha_global <- rho_global <- mu_global <- NULL
    diagS_group <- .get_diagSPD(
      alpha_group, rho_group,
      M, L, kernel
    )
    z_group <- rnorm(M)
    beta_group <- z_group * diagS_group

    ## sum to zero requirement
    QR <- .qr_helmert(n_ind)
    mu_ind_raw <- rnorm(n_ind - 1)
    mu_ind <- mu_group + (QR %*% mu_ind_raw) * sigma_ind

    z_ind_raw <- matrix(rnorm(M * (n_ind - 1)),
      nrow = n_ind - 1, ncol = M
    )
    z_ind <- QR %*% z_ind_raw
    beta_ind <- t(t(z_ind) * diagS_ind)
    ## group-level random coefficients

    f_group <- phi %*% z_group
    f_ind <- phi %*% t(beta_ind)
    eta_ind <- matrix(nrow = n_ind, ncol = length(t))
    for (i in 1:n_ind)

    eta_ind <- f_ind + matrix(rep(mu_ind, each = nrow(f_ind)),nrow=nrow(f_ind)) + f_group[,rep(1,n_ind)]
    eta_group <- mu_group + f_group

    sim_data$traces <- list(
      ind_traces = .format_df(n_ind, t_grid, eta_ind,level="ind"),
      group_traces = .format_df(n_group, t_grid, eta_group,level="group"),
      type = type
    )


  } else if (type == "multi_group") {
    if (any(missing(alpha_ind),missing(alpha_group),missing(alpha_global))) {
      stop("Please specify an alpha_ind, alpha_group and alpha_global or set alpha and omega to partition alpha across the three levels.")
    }

    diagS_group <- .get_diagSPD(
      alpha_group, rho_group,
      M, L, kernel
    )
    diagS_global <- .get_diagSPD(
      alpha_global, rho_global,
      M, L, kernel
    )
    z_global <- rnorm(M)
    beta_global <- z_global * diagS_global

    QR <- .qr_helmert(n_group)
    ind_by_group <- rep(seq_along(n_ind), times = n_ind)
    ###########

    mu_group_raw <- rnorm(n_group - 1)
    mu_group <- mu_global + (QR %*% mu_group_raw) * sigma_group

    mu_ind_raw <- rnorm(sum(n_ind) - n_group)


    mu_ind <- mu_group[ind_by_group] + .constrain_by_hierarchy(n_ind, mu_ind_raw) * sigma_ind
    ###########

    z_group_raw <- matrix(rnorm(M * (n_group - 1)),
      nrow = n_group - 1, ncol = M
    )
    z_group <- QR %*% z_group_raw
    beta_group <- t(t(z_group) * diagS_group)
    ###########

    z_ind_raw <- matrix(rnorm(M * (sum(n_ind) - n_group)),
      nrow = sum(n_ind) - n_group, ncol = M
    )
    z_ind <- .constrain_by_hierarchy(n_ind, z_ind_raw)
    beta_ind <- t(t(z_ind) * diagS_ind)
    ###########

    f_global <- phi %*% beta_global
    f_group <- phi %*% t(beta_group)
    f_ind <- phi %*% t(beta_ind)


    eta_ind <- f_ind + matrix(rep(mu_ind, each = nrow(f_ind)),nrow=nrow(f_ind)) +
      f_group[,ind_by_group] +
      f_global[,rep(1,sum(n_ind))]

    eta_group <- f_group + matrix(rep(mu_group, each = nrow(f_group)),
                                  nrow=nrow(f_group)) + f_global[,rep(1,n_group)]

    eta_global <- mu_global + f_global

    sim_data$traces <- list(
      ind_traces = .format_df(n_ind, t_grid, eta_ind, level = "ind"),
      group_traces = .format_df(n_group, t_grid, eta_group, level = "group"),
      global_trace = .format_df(n_group, t_grid, eta_global, level = "global"),
      type = type
    )
  }



  hsgp_hyperparams <- list(rho_ind = rho_ind,
                      rho_group = rho_group,
                      rho_global = rho_global,
                      alpha_ind = alpha_ind,
                      alpha_group = alpha_group,
                      alpha_global = alpha_global,
                      M = M,
                      L = L,
                      kernel = kernel)

  Filter(Negate(is.null), hsgp_hyperparams)

  mu_ind_list <- as.list(setNames(mu_ind,
                              paste0("mu_ind[", seq.int(sum(n_ind)), "]")))

  if (type %in% c("one_group","multi_group")) {
    mu_group_list <- as.list(setNames(mu_group,
                                    paste0("mu_group[", seq.int(n_group), "]")))
    mu_params <- append(mu_ind_list, mu_group_list)
    if (type == "multi_group") {
      mu_params$mu_global <- mu_global
    }
  } else mu_params <- mu_ind_list


  class(hsgp_hyperparams) <- c("hsgp_hyperparams",class(hsgp_hyperparams))
  class(mu_params) <- c("mu_params",class(mu_params))
  sim_data$sim_params$hsgp_hyperparams <- hsgp_hyperparams
  sim_data$sim_params$mu_params <- mu_params

  sim_data$sim_parameters$seed_traces <- seed

  class(sim_data) <- append("sim_traces", class(sim_data))
  return(sim_data)
}

.format_df <- function(N, t_grid, eta,level) {


  len_t <- length(t_grid)

  if (level == "ind") {
    total_ind <- sum(N)

    ind_vec <- factor(rep.int(seq_len(total_ind), len_t))
    t_vec   <- rep(t_grid, each = total_ind)

    ind_by_group <- rep.int(seq_along(N), times = N)
    group_vec    <- rep.int(ind_by_group, len_t)

    formatted_df <- data.frame(
      ind = ind_vec,
      group = group_vec,
      t = t_vec,
      eta = as.vector(t(eta)),
      stringsAsFactors = FALSE
    )

  } else if (level == "group") {
    total_groups <- N

    group_vec <- rep.int(seq_len(total_groups),len_t)
    t_vec <- rep(t_grid, each = total_groups)

    formatted_df <- data.frame(
      group = group_vec,
      t = t_vec,
      eta = as.vector(t(eta)),
      stringsAsFactors = FALSE
    )

  } else if (level == "global") {

    formatted_df <- data.frame(
      t = t_grid,
      eta = as.vector(eta),
      stringsAsFactors = FALSE
    )
  }

  return(as.data.frame(formatted_df))
}

.constrain_by_hierarchy <- function(a_per_b, obj) {
  ids_raw <- rep(seq_along(a_per_b), times = a_per_b - 1)
  ids <- rep(seq_along(a_per_b), times = a_per_b)

  dims <- if (is.null(dim(obj)[[2]])) {
    c(sum(a_per_b), 1)
  } else {
    c(sum(a_per_b), dim(obj)[[2]])
  }

  constrained <- matrix(nrow = sum(a_per_b), ncol = dims[[2]])
  for (i in unique(ids)) {
    QR <- .qr_helmert(a_per_b[i])
    if (length(dim(obj)) > 1) {
      constrained[ids == i, ] <- QR %*% obj[ids_raw == i, ]
    } else {
      constrained[ids == i] <- QR %*% obj[ids_raw == i]
    }
  }
  return(constrained)
}

.get_diagSPD <- function(alpha, rho, M, L, kernel) {
  if (kernel == "squared_exp") {
    diagS <- .diagSPD_EQ(alpha, rho, M, L)
  } else if (kernel == "matern12") {
    diagS <- .diagSPD_Matern12(alpha, rho, M, L)
  } else if (kernel == "matern32") {
    diagS <- .diagSPD_Matern32(alpha, rho, M, L)
  } else if (kernel == "matern52") {
    diagS <- .diagSPD_Matern52(alpha, rho, M, L)
  }
  return(diagS)
}

.get_M_bases <- function(duration, L_factor, rho_list,
                         kernel = "squared_exp") {
  idx <- match(kernel, c("squared_exp", "matern12", "matern32", "matern52"))
  scale_factors <- c(1.75, 100, 3.42, 2.65)
  L_limits <- c(3.2, 8, 4.5, 4.1)
  S <- duration / 2
  if (L_factor >= 1.2 | L_factor >= L_limits[idx] * min(rho_list) / S) {
    M <- scale_factors[idx] * L_factor * S / min(rho_list)
  } else {
    stop("Please choose a higher L_factor for stability.")
  }
  return(ceiling(M))
}
#' @noRd
.qr_helmert <- function(N) {
  QR <- matrix(0, nrow = N, ncol = N - 1)

  for (i in 1:(N - 1)) {
    for (j in 1:N) {
      if (j <= i) {
        QR[j, i] <- 1 / sqrt(i * (i + 1))
      } else if (j == i + 1) {
        QR[j, i] <- -i / sqrt(i * (i + 1))
      }
    }
  }
  return(QR)
}

#' @noRd
.diagSPD_EQ <- function(alpha, rho, M, L) {
  indices <- seq_len(M)
  factor <- alpha * sqrt(sqrt(2 * pi) * rho)
  exponent <- -0.25 * (rho * pi / (2 * L))^2
  factor * exp(exponent * indices^2)
}

#' @noRd
.diagSPD_Matern52 <- function(alpha, rho, M, L) {
  factor <- 16 * (sqrt(5) / rho)^5
  indices <- (pi / (2 * L) * seq_len(M))^2

  denom <- 3 * ((5 / rho^2) + indices)^3
  return(alpha * sqrt(factor / denom))
}

#' @noRd
.diagSPD_Matern12 <- function(alpha, rho, M, L) {
  indices <- 1:M
  factor <- 2.0

  denom <- rho * ((1.0 / rho)^2 + (pi * indices / (2 * L))^2)
  return(alpha * sqrt(factor * (1 / denom)))
}

#' @noRd
.diagSPD_Matern32 <- function(alpha, rho, M, L) {
  M_series <- seq_len(M)
  indices <- ((pi / (2 * L)) * M_series)^2

  factor <- 2 * alpha * (sqrt(3) / rho)^1.5
  denom <- 3 / (rho^2) + indices
  return(factor / denom)
}

#' @noRd
.diagSPD_Periodic <- function(alpha, rho, M, ...) {
  a <- 1 / (rho^2)
  indices <- 1:M

  log_bessel <- log(besselI(x = a, nu = indices, expon.scaled = TRUE)) + a

  q <- exp(log(alpha) + 0.5 * (log(2) - a + log_bessel))

  # append_row(q, q) is equivalent to concatenating the vector with itself
  return(append_row(q, q))
}
