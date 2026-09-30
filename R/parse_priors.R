#' Parse priors
#' @noRd
.get_family_params <- function(family) {

  if (family == "gamma") return("k")
  else if (family == "weibull") return("shape")
  else if (family == "gengamma") return(c("k","shape"))
  else if (family == "lognormal") return(c("mu_lognormal","sigma_lognormal"))
}

#' @noRd
.process_priors <- function(model_settings, priors = NULL) {

  all_types <- c("one_ind","one_group","multi_group")
  dists <- c("normal", "lognormal", "cauchy",
             "invgamma", "gamma", "exp")


  # prepare prior frame
  duration <- model_settings$settings$duration
  type <- model_settings$settings$type


  survival_params <- .get_family_params(model_settings$settings$family)

  param_prefix <- c("mu","alpha","rho","sd",survival_params)
  params <- outer(c("mu_","alpha_","rho_","sd_"),
                  c("ind","group","global"),
                  paste0)

  valid_params <- c(as.vector(t(params[,1:match(type,all_types)])),survival_params)

  prior_frame <- dplyr::tibble(
    param = valid_params)

  prior_frame[paste0("hyperp",1:2)] <- NA
  prior_frame$dist_id <- NA
  #prior_frame <- readr::read_csv(path, show_col_types = FALSE)

  if (missing(priors) || is.null(priors)) {
    prior_frame <- .apply_default_priors(model_settings, prior_frame) |>
      dplyr::mutate(dist = dists[dist_id])

    return(prior_frame)
  }

  for (f in priors) {
    formula <- as.formula(f)
    prior_param <- as.character(formula[[2]])

    if (!(prior_param %in% valid_params) & !(prior_param %in% param_prefix)) {
      stop(paste0("Variable '", prior_param, "' is not valid for this model type. ",
                  "Allowed parameters: ", paste(valid_params, collapse = ", ")))
    }

    dist_data <- formula[[3]]
    dist_name <- as.character(dist_data[[1]])

    if (!(dist_name %in% dists)) {
      stop("Supported distributions: normal, lognormal, cauchy, invgamma, gamma, or exp.")
    }

    dist_args <- lapply(as.list(dist_data)[-1], function(arg) {
      val <- tryCatch(eval(arg), error = function(e) arg)
      if (!is.numeric(val) || length(val) != 1) {
        stop(sprintf(
          "Prior parameter '%s' for variable '%s' must be a single numeric constant (negative values are allowed).",
          deparse(arg), prior_param
        ))
      }
      val
    })
    #if (any(!sapply(dist_args, is.numeric))) {
    #  stop("Prior parameters must be numeric constants.")
    #}
    prior_index <- startsWith(prior_frame$param, prior_param)
    prior_frame[prior_index, "dist_id"] <- match(dist_name,dists)

    if (length(dist_args) == 1) {
      prior_frame[prior_index, "hyperp1"] <- dist_args[[1]]
      prior_frame[prior_index, "hyperp2"] <- NA
    } else {
      prior_frame[prior_index, c("hyperp1", "hyperp2")] <- dist_args
    }
  }
  prior_frame <- .apply_default_priors(model_settings,prior_frame) |>
    dplyr::mutate(dist = dists[dist_id])
  return(prior_frame)
}


#' @noRd
.apply_default_priors <- function(model_settings, prior_frame) {
  expected_evrange <- model_settings$settings$expected_evrange
  expected_lengthscale <- model_settings$settings$expected_lengthscale
  if (length(expected_evrange) != 2) {
    stop("Please ensure you specify a plausible range of event rates your data could
         a priori take. By default, this range is set to 2-50 events/min")
  }
  if (expected_evrange[1] > expected_evrange[2]) {
    stop("Please order the range of expected events from lowest to highest")
  }

  all_types <- c("one_ind","one_group","multi_group")
  duration <- model_settings$settings$duration
  type <- model_settings$settings$type

  ## survival

  ## mu - normal prior


  log_ers <- -log(expected_evrange/60)
  mu <- mean(log_ers)
  sd_mu <- (mu - log_ers[1]) / -2.326

  ## alphas - normal prior ~ N(0,dt_sd)
  dt_sd <- sd(model_settings$dt) / 3

  ## lengthscales / rho - inverse gamma rho ~ InvGamma(alpha,beta)
  invgamma_params <- .get_invgamma_hyperparams(model_settings$dt,
                                               duration = duration)
  a = 5
  b_ind = (a - 1) * expected_lengthscale

  b_group = (5/4 * a - 1) * 5/4 * expected_lengthscale
  b_global = (25/16 * a - 1) * 25/16 * expected_lengthscale

  hyperp1_prefix <- c(alpha = 0.5, mu = mu, sd = 0)
  hyperp1_exact  <- c(k = 0, shape = 0,
                      rho_ind = a,
                      rho_group = 5/4 * a,
                      rho_global = 25/16 * a)

  hyperp2_prefix <- c(alpha = 0.3, mu = sd_mu)
  hyperp2_exact  <- c(sd_ind = 0.5, sd_group = 0.1,
                   rho_ind =  b_ind,
                   rho_group = b_group,
                   rho_global = b_global,
                   k = 1, shape = 1)
  dist <- c(rho = 4, alpha = 1, mu = 1, sd = 1, k = 2, shape = 2)
  prior_frame <- prior_frame |>
    dplyr::mutate(
      prefix = sub("_.*", "", param),
      dist_id = dplyr::if_else(is.na(dist_id),
                               dist[prefix],dist_id),
      hyperp1 = dplyr::if_else(is.na(hyperp1),
                          dplyr::coalesce(hyperp1_prefix[prefix],
                        hyperp1_exact[param]),hyperp1),
      hyperp2 = dplyr::if_else(is.na(hyperp2),
                       dplyr::coalesce(hyperp2_prefix[prefix],
                                hyperp2_exact[param]),
                       hyperp2)
    ) |>
    select(-prefix)


  if ("mu_global" %in% prior_frame$param) {
    prior_frame <- prior_frame[!any(prior_frame$param %in%
                                        c("mu_ind","mu_group")) ,]
  } else if ("mu_group" %in% prior_frame$param) {
    prior_frame <- prior_frame[prior_frame$param != "mu_ind" ,]
  }

  if (model_settings$settings$family == "gamma") {
    prior_frame <- prior_frame[prior_frame$param != "shape" ,]
  } else if (model_settings$settings$family == "weibull") {
    prior_frame <- prior_frame[prior_frame$param != "k" ,]
  } else {
    prior_frame <- prior_frame[!(prior_frame$param %in% c("k","shape")) ,]
  }
  if (type != "multi_group") {
    prior_frame <- prior_frame[prior_frame$param != "sd_group" ,]
  }

  return(prior_frame)
}


############## inverse gamma priors on lengthscale parameters

# .get_invgamma_hyperparams <- function(dt, duration, l_factor = 3) {
#
#   # lower bound - set to 10% of dt mass
#   # upper bound - set to fraction of total duration
#   L <- quantile(dt,
#                 probs = 0.1,names=FALSE)
#   U <- duration / 2
#   print(c(L,U))
#   params <- .compute_invgamma(L,U)
#   return(params)
# }

# .compute_invgamma <- function(L, U, alpha = 6) {
#   beta_from_upper <- uniroot(function(b) invgamma::pinvgamma(U, alpha, b) - 0.99,
#                              lower = 1e-6, upper = 1e6)$root
#   list(alpha = alpha, beta = beta_from_upper,
#        check_lower = invgamma::pinvgamma(L, alpha, beta_from_upper))
# }
#' @noRd
.get_invgamma_hyperparams <- function(dt, duration,
                                               l_factor = 2.5, l_quantile = 0.2,
                                               group_mult = 5/4, global_mult = 5/4) {
  L_ind <- quantile(dt, probs = l_quantile, names = FALSE)
  U_ind <- duration / l_factor

  L_group <- L_ind * group_mult
  U_group <- U_ind * group_mult

  L_global <- L_group * global_mult
  U_global <- U_group * global_mult

  ind    <- .compute_invgamma(L_ind,    U_ind)
  group  <- .compute_invgamma(L_group,  U_group)
  global <- .compute_invgamma(L_global, U_global)

  list(ind = ind, group = group, global = global)
}

#' @noRd
.compute_invgamma <- function(L, U, alpha_floor = 6) {
  target <- function(pars) {
    alpha <- alpha_floor + exp(pars[1])
    beta  <- exp(pars[2])
    c(invgamma::pinvgamma(L, alpha, beta) - 0.01,
      invgamma::pinvgamma(U, alpha, beta) - 0.99)
  }

  start <- c(log(1), log((alpha_floor + 1) * sqrt(L * U)))
  sol <- nleqslv::nleqslv(start, target, control = list(xtol = 1e-10, ftol = 1e-10))

  alpha <- alpha_floor + exp(sol$x[1])
  beta  <- exp(sol$x[2])

  list(alpha = alpha, beta = beta,
       check_lower = invgamma::pinvgamma(L, alpha, beta),
       check_upper = invgamma::pinvgamma(U, alpha, beta),
       termcd = sol$termcd)
}
