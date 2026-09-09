#' Simulate inhomogeneous renewal processes
#'
#' @param sigma Numeric. Individual-level coefficient standard deviation.
#' @param Q Numeric. Mean coefficients for the spline basis.
#' @param baseline Numeric. Group-level coefficient standard deviation.
#' @param baseline_sd Numeric. Individual-level coefficient standard deviation.
#' @param lerp Numeric. Mean coefficients for the spline basis.
#'
#' @return A tibble containing simulated event times for each group.
#' @export

gp_fit <- function(events, ...) {
  UseMethod("gp_fit")
}

#' @export
#' @method gp_fit data.frame
gp_fit.data.frame <- function(events,...) {
  event_data <- list(events = events)
  gp_fit(event_data,...)
}

#' @export
#' @method gp_fit list
gp_fit.list <- function(event_data, duration, group_id = "group", ind_id = "ind",subs=NULL,
                       event_times = "event_times", dt = "dt",run = c("prior_pc","fit"),

                       M = 40, family = c("exponential","gamma","weibull","lognormal","gengamma"),
                       kernel = c("squared_exp", "matern12","matern32","matern52","periodic"),
                       priors = NULL, L_factor = 1.2, w0 = 0.5, n_quad = 3,
                       chains = 4, iter = 2000, warmup = 1000, adapt_delta = 0.95, max_treedepth = 12,old=FALSE,...) {
  kernel <- match.arg(kernel)
  family <- match.arg(family)

  if (kernel == "periodic") {
    warning("Periodic kernel uses half the number of basis functions (M). Doubling M")
    M = M * 2
  }

  events <- event_data$events

  ind_str <- tryCatch(rlang::as_name(rlang::ensym(ind_id)),
                      error = function(e) NULL)
  group_str <- tryCatch(rlang::as_name(rlang::ensym(group_id)),
                        error = function(e) NULL)
  event_times_str <- tryCatch(rlang::as_name(rlang::ensym(event_times)),
                     error = function(e) NULL)

  dt_str <- tryCatch(rlang::as_name(rlang::ensym(dt)),
                    error = function(e) NULL)

  if (missing(duration)) {
    duration <- events |>
      reframe(max_t = max(.data[[event_times_str]])) |>
      pull(max_t)
    warning("Setting duration to the latest event time. If this is incorrect, please set duration.")
  }

  if (is.null(ind_str) || !ind_str %in% names(events)) {
    warning("Individual ID not provided or not found in 'event_data'. Defaulting to 'one_ind' model")
    subclass <- "one_ind"
    I <- 1
    ind_id <- rep(1, nrow(events))
    g_id <- rep(1, nrow(events))
    g_membership <- 1
  } else {
    I <- events |>
      dplyr::distinct(.data[[ind_str]]) |>
      dplyr::pull() |>
      length()

    if (I == 1) {
      subclass <- "one_ind"
      ind_id <- rep(1, nrow(events))
      G <- 1
      g_id <- rep(1, nrow(events))
      g_membership <- 1
    } else {
      ind_id <- events |> dplyr::select(.data[[ind_str]]) |>
        dplyr::pull()
      if (is.null(group_str) || !group_str %in% names(events)) {
        warning("Group ID not provided or not found in 'events'. Defaulting to 'one_group' model.")
        subclass <- "one_group"
        G <- 1
        g_id <- rep(1, nrow(events))
        g_membership <- rep(1, I)

      } else {

        G <- events |>
          dplyr::distinct(.data[[group_str]]) |>
          dplyr::pull() |>
          length()
        g_id <- events |> dplyr::select(.data[[group_str]]) |>
          dplyr::pull()

        if (G == 1) {
          subclass <- "one_group"
          g_membership <- rep(1, I)
        } else {
          subclass <- "multi_group"
          g_membership <- events |>
            distinct(.data[[group_str]],.data[[ind_str]]) |>
            pull(group)

        }
      }
    }
  }

  if (is.null(event_times_str) || !event_times_str %in% names(events)) {
    stop("Please make sure that event times are present and that the column ID you specified has the same name.")
  }

  if (is.null(dt_str) || !dt_str %in% names(events)) {
    events <- events |>
      mutate(dt = c(0,diff(.data[[event_times_str]])))
  }

  model_data <- list(
    old = old,
    event_times = events |>
      pull(all_of(event_times_str)),
    traces = event_data$traces, # for simulations
    censored = events$censored,
    group_traces = event_data$group_traces,
    global_trace = event_data$global_trace,
    events = events |> rename(ind = .data[[ind_str]],
                              event_times = .data[[event_times_str]]),
    dt = events |> pull(all_of(dt_str)),
    ind_id = as.numeric(as.character(ind_id)),
    g_id = as.numeric(as.character(g_id)),
    g_membership = as.numeric(as.character(g_membership)),
    I = I,
    G = G,
    settings = list(duration = duration,
                    M = M,
                    L_factor = L_factor,
                    w0 = w0,
                    kernel = match.arg(kernel),
                    family = match.arg(family),
                    priors = priors,
                    n_quad = n_quad),
    stan_runtime = list(chains = chains,
                        iter = iter, warmup = warmup,
                        adapt_delta = adapt_delta,
                        max_treedepth = max_treedepth),
    sim_parameters = event_data$sim_parameters
  )

  if ("prior_pc" %in% run & "fit" %in% run) {
    class(model_data) <- c("gp_prior_pc",subclass)
    gp_fit(model_data) |> gp_fit()
  } else if ("prior_pc" %in% run) {
    class(model_data) <- c("gp_prior_pc",subclass)
    gp_fit(model_data)
  } else if ("fit" %in% run) {
    class(model_data) <- c("gp_model",subclass)
    gp_fit(model_data)
  }
}

#' @export
#' @method gp_fit gp_model
gp_fit.gp_model <- function(model_data, ...) {

  if ("one_ind" %in% class(model_data)) {
    message("Fitting a single-individual GP model")
    prior_frame <- parse_priors(model_data)
    stan_obj <- stanmodels$hsgp_one_ind
  } else if ("one_group" %in% class(model_data)) {
    message("Fitting a single-group hierarchical GP model")
    prior_frame <- parse_priors(model_data)
    if (model_data$old) {
      stan_obj <- stanmodels$hsgp_one_group2
    } else {
      stan_obj <- stanmodels$hsgp_one_group
    }
  } else if ("multi_group" %in% class(model_data)) {
    message("Fitting a multi-group hierarchical GP model")
    prior_frame <- parse_priors(model_data)
    stan_obj <- stanmodels$hsgp_multi_group
  }

  model_data$settings <- append(model_data$settings,
                                 list(variables=prior_frame |>
                                        pull(prior_variable)))
  flags <- .get_flags(model_data)

  stan_dat <- list(
    N_total = length(model_data$event_times),
    I = model_data$I,
    ind_id = model_data$ind_id,
    G = model_data$G,
    g_id = model_data$g_id,
    g_membership = model_data$g_membership,
    I_per_group = tabulate(model_data$g_membership),

    M = model_data$settings$M,
    L_factor = model_data$settings$L_factor,
    w0 = model_data$settings$w0,
    duration = model_data$settings$duration,
    n_quad = model_data$settings$n_quad,
    t_ev = model_data$event_times,
    dt = model_data$dt,
    censored = model_data$censored,

    distributions = as.double(prior_frame$distribution_id),
    N_params = nrow(prior_frame),
    params = as.matrix(prior_frame[c("param_1","param_2")]),
    kernel = match(model_data$settings$kernel,c("squared_exp",
                                                "matern12","matern32","matern52","periodic")),
    family = match(model_data$settings$family, c("exponential","gamma","weibull",
                                                 "lognormal","gengamma"))
  )
  stan_dat <- append(stan_dat, flags)

  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)

  fit <- rstan::sampling(
    object = stan_obj,
    data = stan_dat,
    init = .build_inits(stan_dat,
                        model_data$stan_runtime$chains),
    chains = model_data$stan_runtime$chains,
    iter = model_data$stan_runtime$iter,
    warmup = model_data$stan_runtime$warmup,
    control = list(adapt_delta = model_data$stan_runtime$adapt_delta,
                   max_treedepth = model_data$stan_runtime$max_treedepth))

  model_data$fit <- fit
  return(model_data)
}

#' @export
#' @method gp_fit gp_prior_pc
gp_fit.gp_prior_pc <- function(prior_pc_data, ...) {
  if ("one_ind" %in% class(prior_pc_data)) {
    message("Running prior predictive checks on model parameters")
    prior_frame <- parse_priors.one_ind(prior_pc_data)
    stan_obj <- stanmodels$prior_pc_one_ind
  } else if ("one_group" %in% class(prior_pc_data)) {
    message("Running prior predictive checks on model parameters")
    prior_frame <- parse_priors.one_group(prior_pc_data)
    stan_obj <- stanmodels$prior_pc_one_group
  } else if ("multi_group" %in% class(prior_pc_data)) {
    message("Running prior predictive checks on model parameters")
    prior_frame <- parse_priors.multi_group(prior_pc_data)
    stan_obj <- stanmodels$prior_pc_multi_group
  }

  prior_pc_data$settings <- append(prior_pc_data$settings,
                                   list(variables=prior_frame |>
                                          pull(prior_variable)))

  flags <- .get_flags(prior_pc_data)

  stan_dat <- list(
    I = prior_pc_data$I,
    ind_id = prior_pc_data$ind_id,
    G = prior_pc_data$G,
    g_id = prior_pc_data$g_id,
    g_membership = prior_pc_data$g_membership,
    I_per_group = tabulate(prior_pc_data$g_membership),

    M = prior_pc_data$settings$M,
    L_factor = prior_pc_data$settings$L_factor,
    w0 = prior_pc_data$settings$w0,
    duration = prior_pc_data$settings$duration,

    distributions = as.double(prior_frame$distribution_id),
    N_params = nrow(prior_frame),
    params = as.matrix(prior_frame[c("param_1","param_2")]),
    kernel = match(prior_pc_data$settings$kernel,c("squared_exp",
                                                   "matern12","matern32","matern52","periodic")),
    family = match(prior_pc_data$settings$family, c("exponential","gamma","weibull",
                                                    "lognormal","gengamma"))
  )

  stan_dat <- append(stan_dat, flags)

  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)

  fit <- rstan::sampling(
    object = stan_obj,
    data = stan_dat,
    chains = prior_pc_data$stan_runtime$chains,
    iter = prior_pc_data$stan_runtime$iter,
    algorithm = "Fixed_param")

  prior_pc_data$prior_pc <- fit
  class(prior_pc_data) <- c("full_config",class(prior_pc_data))

  return(prior_pc_data)
}

#' @export
#' @method gp_fit full_config
gp_fit.full_config <- function(model_data) {
  class(model_data) <- setdiff(class(model_data), c("gp_prior_pc","full_config"))
  class(model_data) <- c("gp_model",class(model_data))
  gp_fit(model_data)
}

#' @export
.filter_reindex <- function(tib, subs, ...) {
  cols <- c(...)
  tib <- tib |> dplyr::filter(if_any(all_of(cols),
                                    ~ .x %in% subs))
  for (col in cols) {
    tib[[col]] <- as.integer(forcats::fct_inorder(as.character(tib[[col]])))
  }
  return(tib)
}

#' @export
filter_data <- function(event_data, subs, by = c("ind", "group")) {
  by <- match.arg(by)

  if (by == "ind") {
    event_data$traces <- .filter_reindex(event_data$traces, subs, "ind")
    event_data$events <- .filter_reindex(event_data$events, subs, "ind")
  } else {
    event_data$traces <- .filter_reindex(event_data$traces, subs, "ind", "group")
    event_data$group_traces <- .filter_reindex(event_data$group_traces, subs, "group")
    event_data$events <- .filter_reindex(event_data$events, subs, "group", "ind")
  }

  event_data$n_ind <- dplyr::n_distinct(event_data$traces$ind)
  event_data$n_groups <- dplyr::n_distinct(event_data$traces$group)


  if (event_data$n_groups == 1) event_data$type <- "one_group"
  if (event_data$n_ind == 1)  event_data$type <- "one_ind"

  return(event_data)
}

#' @noRd
.get_flags <- function(model_data) {
  flags <- list()
  if (model_data$settings$family != "lognormal") {
    flags[["include_mu_lognormal"]] = 0
    flags[["include_sigma_lognormal"]] = 0
  }

  if (model_data$settings$family == "exponential") {
    flags[["include_k"]] = 0
    flags[["include_shape"]] = 0
  } else if (model_data$settings$family == "gamma") {
    flags[["include_k"]] = match("k",model_data$settings$variables)
    flags[["include_shape"]] = 0
  } else if (model_data$settings$family == "weibull") {
    flags[["include_k"]] = 0
    flags[["include_shape"]] = match("shape",model_data$settings$variables)
  } else if (model_data$settings$family == "gengamma") {
    flags[["include_k"]] = match("k",model_data$settings$variables)
    flags[["include_shape"]] = match("shape",model_data$settings$variables)
  }

  return(flags)

}

.build_inits <- function(stan_dat, n_chains) {
  mean_dt   <- mean(stan_dat$dt)
  median_dt <- median(stan_dat$dt)

  lapply(seq_len(n_chains), function(i) {
    list(
      rho_group = median_dt * exp(rnorm(1, 0, 0.1)),
      rho_ind = median_dt * exp(rnorm(1, 0, 0.1)),

      alpha_group = 0.3 * exp(rnorm(1, 0, 0.1)),
      alpha_ind = 0.3 * exp(rnorm(1, 0, 0.1)),

      alpha = 0.4 * exp(rnorm(1, 0, 0.1)),
      omega = 0.5,
      mu_group = -log(mean_dt) + rnorm(1, 0, 0.1),

      sigma_ind = 0.2,

      z_group = rnorm(stan_dat$M, 0, 0.1),
      z_ind_raw  = matrix(rnorm((stan_dat$I - 1) * stan_dat$M, 0, 0.1),
                          nrow = stan_dat$I - 1),
      mu_raw_ind = rnorm(stan_dat$I - 1, 0, 0.1),

      k = array(1.5 + abs(rnorm(1, 0, 0.3)))
      #shape = array(1.2 + abs(rnorm(1, 0, 0.3))),

      #sigma_lognormal = array(0.5 * exp(rnorm(1, 0, 0.1)))
    )
  })
}


