#' Simulate Latent Spline Traces
#' @description
#' A short description...
#'
#' @param duration Integer. The length of the full time series to simulate (in seconds). The default is 100 seconds
#' @param n_ind Integer or integer vector. Number of individuals per group. If 1 integer is supplied,
#' @param n_groups Integer. The number of groups to simulate. If the number of groups is greater than 1,
#' group splines are constrained to sum to zero. The default is 1 group.
#' @param df Integer. Degrees of freedom for the cubic spline basis set. The default is 5
#' @param sd_group Numeric. Standard deviation for group-level coefficients.
#' If the number of individuals is greater than 1, the s.d. of the group-level coefficients determines how  standard deviation.
#' @param sd_ind Numeric. Individual-level coefficient standard deviation.
#' @param mu_coeffs Numeric vector. Mean coefficients for the spline basis.
#' @param seed Numeric. Set a seed for reproducibility. If one is not specified, seed defaults to NULL
#'
#' @return An object of class `sim_traces`, containing a list of `ind_traces`,
#' `group_traces` (trace common to individuals in the same group),
#' `global_traces` (trace common to all individuals; NULL if simulation has only
#' 1 individual or 1 group) and `sim_parameters` (simulation parameters).
#' @seealso The methods [simulate_gp_traces()] and [simulate_boxcar_traces()]
#' @examples
#' spline_traces <- simulate_spline_traces(n_ind = 5, n_groups = 2)
#' plot(spline_traces)
#'
#' @export
simulate_spline_traces <- function(duration = 100, n_ind = 10, n_groups = 1,
                                         df = 5, sd_group = 0.5, sd_ind = 0.25,
                                         mu_coeffs = stats::rnorm(df, 0, 1), seed = NULL) {
  set.seed(seed)

  sim_struct <- .create_sim_struct(n_groups = n_groups, n_ind = n_ind)
  if (n_ind == 1) {
    alpha_global <- alpha_group <- 0
    type <- "one_ind"
  } else if (n_groups == 1) {
    alpha_global <- 0
    type <- "one_group"
  } else {
    type <- "multi_group"
  }

  if (n_ind == 1 || n_groups == 1) mu_coeffs <- 0

  group_traces <- .simulate_n(duration,
    n = n_groups, class = "spline",
    df = df, mu_coeffs = mu_coeffs, sd_coeffs = sd_group, level = "group"
  )

  ind_traces <- .simulate_n(duration,
    n = n_ind * n_groups, class = "spline",
    df = df, mu_coeffs = group_traces$coeffs[, sim_struct$group, drop = FALSE],
    sd_coeffs = sd_ind, g_membership = sim_struct |> pull(group), level = "ind"
  )

  ind_traces$traces <- ind_traces$traces |>
    dplyr::left_join(sim_struct, by = "ind") |>
    dplyr::relocate(group, ind)

  args_list <- mget(names(formals()), sys.frame(sys.nframe()))

  exclude <- c("ind_traces", "group_traces", "type", "mu_coeffs", "seed")
  args_list <- args_list[!names(args_list) %in% exclude]

  sim_traces <- list(
    traces = ind_traces$traces,
    group_traces = group_traces$traces,
    global_coeffs = mu_coeffs,
    group_coeffs = group_traces$coeffs,
    ind_coeffs = ind_traces$coeffs,
    type = type
  )
  sim_traces$sim_parameters <- args_list
  sim_traces$sim_parameters$seed_traces <- seed

  class(sim_traces) <- c("sim_traces", class(sim_traces))
  return(sim_traces)
}

#' Simulate Boxcar Trajectories
#' @inheritParams simulate_spline_traces
#' @param nodes Numeric vector. Node locations for step intervals.
#' @param steps Numeric vector. Step sizes for the boxcar transitions.
#' @param cycles Integer. Number of repeating cycles over the duration.
#' @param sd_group Numeric. Group-level step scaling standard deviation.
#' @param sd_ind Numeric. Individual-level step scaling standard deviation.
#' @inherit simulate_spline_traces return
#' @seealso The methods [simulate_spline_traces()] and
#' [simulate_gp_traces()]
#' @examples
#' boxcar_traces <- simulate_boxcar_traces(n_ind = 5, n_groups = 2)
#' plot(boxcar_traces)
#' @export
simulate_boxcar_traces <- function(duration = 100, n_ind = 10, n_groups = 1,
                                         nodes = seq(1, duration, length.out = 3), steps = c(-1, 1),
                                         cycles = 2, sd_group = 0.5, sd_ind = 0.25) {
  if (n_ind == 1) {
    alpha_global <- alpha_group <- 0
    type <- "one_ind"
  } else if (n_groups == 1) {
    alpha_global <- 0
    type <- "one_group"
  } else {
    type <- "multi_group"
  }

  sim_struct <- .create_sim_struct(n_groups = n_groups, n_ind = n_ind)

  scale_factors_group <- stats::rnorm(n_groups, 1, sd_group)
  if (n_groups == 1 || n_ind == 1) scale_factors_group <- 1

  scale_factors_ind <- scale_factors_group[sim_struct$group] + stats::rnorm(n_groups * n_ind, 0, sd_ind)

  group_traces <- .simulate_n(duration,
    n = n_groups, class = "boxcar",
    nodes = nodes, steps = steps, cycles = cycles,
    scale_factor = scale_factors_group, level = "group"
  )

  ind_traces <- .simulate_n(duration,
    n = n_ind * n_groups, class = "boxcar",
    nodes = nodes, steps = steps,
    cycles = cycles, scale_factor = scale_factors_ind, level = "ind"
  )

  ind_traces$traces <- ind_traces$traces |>
    dplyr::left_join(sim_struct, by = "ind") |>
    dplyr::relocate(ind, group)

  args_list <- mget(names(formals()), sys.frame(sys.nframe()))

  exclude <- c(
    "ind_traces", "group_traces",
    "scale_factors_ind", "scale_factors_group", "type"
  )
  args_list <- args_list[!names(args_list) %in% exclude]

  sim_traces <- list(
    traces = ind_traces$traces,
    steps = ind_traces$level_steps,
    group_traces = group_traces$traces,
    group_steps = group_traces$level_steps,
    type = type
  )
  sim_traces$sim_parameters <- args_list

  class(sim_traces) <- c("sim_traces", class(sim_traces))
  return(sim_traces)
}


#' Simulate Gaussian Process Trajectories
#'
#' @inheritParams simulate_spline_traces
#' @param kernel Character. Choice of covariance kernel, K, to simulate Gaussian Process
#' trajectories, comprising "squared_exp", "matern", "periodic".
#' The squared exponential (or radial basis function) is the default kernel.
#' If "matern" is specified, the default Mátern smoothness parameter, nu, is 3/2.
#' "matern12", "matern32" and "matern52" may also be specified, setting smoothness to
#' 1/2, 3/2 and 5/2, respectively. Note that currently only smoothness terms at these intervals
#' are supported in the package. Squared exponential, Mátern and periodic kernels are
#' parameterised by a combination of length-scale, \eqn{\rho} (or \eqn{\ell}) and amplitude-scale,
#' \eqn{\alpha}.
#' @param alpha_global,alpha_group,alpha_ind Numeric. Amplitude parameters of the
#' covariance kernel. `alpha_global`, `alpha_group` and `alpha_ind` respectively
#' set the amplitude terms at global, group and individual-level trajectories.
#' By default, alpha terms at each hierarchy are set to be 60% of the preceding level,
#' rounded to the nearest tenth, starting from a maximum alpha of 0.8.
#' If the simulation is multi-group, `alpha_global = 0.8`,`alpha_group = 0.5` and
#' `alpha_ind = 0.3`. If single-group and multi-individual, `alpha_group = 0.8` and
#' `alpha_ind = 0.5`. If only one alpha term is specified, the others are extrapolated
#' in the same geometric series. Setting `alpha_ind = 0.8`, for example, would set
#' `alpha_group = 1.3`.
#' @param rho_global,rho_group,rho_ind Numeric. Length-scale parameters of the covariance
#' kernel. `rho_global`, `rho_group` and `rho_ind` respectively
#' set the length-scale terms at global, group and individual-level trajectories. By default,
#' length-scale parameters are set at the same scale at 15. If only one rho term is specified,
#' the others are set to the same value.
#' @param nu Numeric. Mátern smoothness parameter. Typical smoothness terms are 1/2, 3/2 and 5/2
#' and these are the only intervals that can currently be modelled by the package.
#' The default is 3/2
#' @param P Numeric. Period duration parameter for periodic kernels.
#'
#' @inherit simulate_spline_traces return
#' @seealso The methods [simulate_spline_traces()] and
#' [simulate_boxcar_traces()]
#' @examples
#' gp_traces<- simulate_gp_traces(n_groups = 2, n_ind = 5, kernel = "matern32")
#' plot(gp_traces)
#'
#' @export
simulate_gp_traces <- function(duration = 100, n_groups = 1, n_ind = 10,
                                     kernel = "squared_exp",
                                     alpha_global = 0.5, rho_global = duration / 2,
                                     alpha_group = 0.5, rho_group = duration / 2,
                                     alpha_ind = 0.25, rho_ind = duration / 2,
                                     nu = 3 / 2, P = 1, seed = NULL) {
  set.seed(seed)

  sim_struct <- .create_sim_struct(n_groups = n_groups, n_ind = n_ind)
  if (n_ind == 1) {
    alpha_global <- alpha_group <- 0
    type <- "one_ind"
  } else if (n_groups == 1) {
    alpha_global <- 0
    type <- "one_group"
  } else {
    type <- "multi_group"
  }

  global_trace <- .simulate_n(duration,
    n = 1, class = "gp",
    kernel = kernel, alpha = alpha_global,
    rho = rho_global, nu = nu, P = P, level = "global"
  )

  group_traces <- .simulate_n(duration,
    n = n_groups, class = "gp",
    kernel = kernel, alpha = alpha_group,
    rho = rho_group, nu = nu, P = P, level = "group"
  )

  if (type == "multi_group") {
    group_traces <- group_traces %>%
      dplyr::group_by(x) |>
      dplyr::mutate(y = y - mean(y)) |>
      dplyr::ungroup() |>
      dplyr::mutate(y = y + global_trace$y)
  } else {
    group_traces <- group_traces |> dplyr::mutate(y = y + global_trace$y)
  }

  ind_traces <- .simulate_n(duration,
    n = n_groups * n_ind, class = "gp",
    kernel = kernel, alpha = alpha_ind,
    rho = rho_ind, nu = nu, P = P, level = "ind"
  ) |>
    dplyr::left_join(sim_struct, by = "ind") |>
    dplyr::left_join(group_traces,
      by = c("group", "x"),
      suffix = c("_ind", "_group")
    )
  if (type != "one_ind") {
    ind_traces <- ind_traces |>
      dplyr::group_by(x, group) |>
      dplyr::mutate(y_ind = y_ind - mean(y_ind)) |>
      dplyr::ungroup() |>
      dplyr::mutate(y = y_ind + y_group) |>
      dplyr::relocate(ind, .after = "group")
  } else {
    ind_traces <- ind_traces |>
      dplyr::mutate(y = y_ind + y_group) |>
      dplyr::relocate(ind, .after = "group")
  }

  args_list <- as.list(environment())
  exclude <- c(
    "ind_traces", "group_traces",
    "global_trace", "type"
  )
  args_list <- args_list[!names(args_list) %in% exclude]

  sim_traces <- list(
    traces = ind_traces,
    group_traces = group_traces,
    global_trace = global_trace,
    type = type
  )
  sim_traces$sim_parameters <- args_list
  sim_traces$sim_parameters$seed_traces <- seed

  class(sim_traces) <- append("sim_traces", class(sim_traces))
  return(sim_traces)
}

.create_sim_struct <- function(n_groups, n_ind) {
  tibble::rowid_to_column(
    tidyr::expand_grid(group = as.factor(seq_len(n_groups)), subgroup_id = as.factor(seq_len(n_ind))),
    var = "ind"
  ) |>
    mutate(ind = as.factor(ind)) |>
    select(-subgroup_id)
}
#' @noRd
.simulate_n <- function(duration, n, class, df = 5, mu_coeffs = 0, sd_coeffs = 0.5,
                        level = "ind", nodes, steps, cycles, scale_factor = 1,
                        g_membership = NULL, ...) {
  class <- match.arg(class, c("spline", "boxcar", "gp"))

  if (class == "spline") {
    if (length(sd_coeffs) != 1 && length(sd_coeffs) != df) {
      stop("Spline standard deviations length must equal 1 or df.")
    }
    if (length(mu_coeffs) != 1 && .get_first_dim(mu_coeffs) != df) {
      stop("Spline mean coefficients length must equal 1 or df.")
    }

    # normalise coefficients to scale to 0
    coeffs <- matrix(stats::rnorm(n * df, 0, 1), nrow = df, ncol = n)
    if (is.matrix(mu_coeffs)) {
      for (g in unique(g_membership)) {
        coeffs[, g_membership == g] <- mu_coeffs[, g_membership == g] +
          t(scale(t(coeffs[, g_membership == g]))) * sd_coeffs
      }
    }


    traces <- purrr::map(1:n, \(i) .sim_splines(duration, df, coeffs[, i, drop = FALSE], ...)) |>
      purrr::list_rbind(names_to = level) |>
      dplyr::mutate(across(all_of(level), as.factor))

    return(list(traces = traces, coeffs = coeffs))
  } else if (class == "boxcar") {
    level_steps <- as.matrix(steps) %*% t(scale_factor)
    traces <- purrr::map(1:n, \(i) .sim_boxcar(duration, nodes, level_steps[, i], cycles)) |>
      purrr::list_rbind(names_to = level) |>
      dplyr::mutate(across(all_of(level), as.factor))

    return(list(traces = traces, scales = level_steps))
  } else if (class == "gp") {
    traces <- purrr::map(1:n, \(i) .sim_GP(duration, ...)) |>
      purrr::list_rbind(names_to = level) |>
      dplyr::mutate(across(all_of(level), as.factor))

    return(traces)
  }
}

#' @noRd
.sim_splines <- function(duration, df, coeffs) {
  basis_matrix <- splines::bs(seq(1, duration), df = df)
  y <- as.vector(basis_matrix %*% coeffs)
  return(tibble::tibble(x = seq(1, duration), y = y))
}

#' @noRd
.sim_GP <- function(duration, kernel = "squared_exp", alpha = 0.5,
                    rho = duration / 5, nu = 3 / 2, P = duration / 2) {
  cov_mat <- .kernel(duration, kernel, alpha, rho, nu, P)
  y <- MASS::mvrnorm(n = 1, mu = rep(0, duration), Sigma = cov_mat)
  return(tibble::tibble(x = seq(1, duration), y = y))
}

#' @noRd
.sim_boxcar <- function(duration, nodes, steps, cycles) {
  if (nodes[length(nodes)] < duration) cycles <- 1
  if (any(diff(nodes) / cycles < 1)) warning("Number of cycles smaller than step grid resolution")

  interval_id <- findInterval(1:(nodes[length(nodes)] / cycles),
    nodes / cycles,
    rightmost.closed = TRUE
  ) |>
    rep_len(duration)

  return(tibble::tibble(x = 1:duration, y = steps[interval_id]))
}

#' @noRd
.kernel <- function(duration, kernel, alpha, rho, nu, P) {
  X <- seq(1, duration, by = 1)
  kernel <- match.arg(kernel, c("squared_exp", "matern", "matern12", "matern32", "matern52", "rational_quad", "periodic"))
  if (kernel == "squared_exp") {
    K <- rkriging::Gaussian.Kernel(rho)
    cov_mat <- alpha^2 * rkriging::Evaluate.Kernel(K, X)
  } else if (kernel == "matern") {
    K <- rkriging::Matern.Kernel(rho, nu = nu)
    cov_mat <- alpha^2 * rkriging::Evaluate.Kernel(K, X)
  } else if (kernel == "matern12") {
    K <- rkriging::Matern12.Kernel(rho)
    cov_mat <- alpha^2 * rkriging::Evaluate.Kernel(K, X)
  } else if (kernel == "matern32") {
    K <- rkriging::Matern32.Kernel(rho)
    cov_mat <- alpha^2 * rkriging::Evaluate.Kernel(K, X)
  } else if (kernel == "matern52") {
    K <- rkriging::Matern52.Kernel(rho)
    cov_mat <- alpha^2 * rkriging::Evaluate.Kernel(K, X)
  } else if (kernel == "rational_quad") {
    K <- rkriging::RQ.Kernel(rho, alpha = 1)
    cov_mat <- alpha^2 * rkriging::Evaluate.Kernel(K, X)
  } else if (kernel == "periodic") {
    cov_mat <- alpha^2 * .periodic_kernel(duration, rho, P)
  }
  return(cov_mat)
}

#' @noRd
.periodic_kernel <- function(duration, rho, P) {
  grid <- seq(1, duration, by = 1)
  dist_mat <- abs(outer(grid, grid, "-"))
  cov_mat <- exp(-(2 * sin(pi * dist_mat / P)^2) / rho^2)
  return(cov_mat)
}

#' @noRd
.get_first_dim <- function(x) {
  if (is.null(dim(x))) length(x) else dim(x)[1]
}
