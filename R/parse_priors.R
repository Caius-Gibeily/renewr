#' Parse priors
#' @return A tibble containing simulated event times for each group.
#' @export
parse_priors <- function(model_setup, ...) {
  UseMethod("parse_priors")
}

#' @export
#' @method parse_priors one_ind
parse_priors.one_ind <- function(model_setup) {
  valid_variables = c("mu_ind","alpha_ind","rho_ind") |>
    append(.get_family_priors(model_setup$settings$family))

  path <- system.file("extdata", "default_priors_one_ind.csv", package = "renewr")

  prior_frame <- .process_prior_frame(path = path,
                                      valid_variables = valid_variables,
                                      priors = model_setup$settings$priors)

  return(prior_frame)
}

#' @export
#' @method parse_priors one_group
parse_priors.one_group <- function(model_setup) {
  valid_variables = c("mu_group","alpha_group",
                      "alpha_ind", "rho_group","rho_ind") |>
    append(.get_family_priors(model_setup$settings$family))

  path <- system.file("extdata", "default_priors_one_group.csv", package = "renewr")

  prior_frame <- .process_prior_frame(path = path,
                                      valid_variables = valid_variables,
                                      priors = model_setup$settings$priors)
  return(prior_frame)
  }

#' @export
#' @method parse_priors multi_group
parse_priors.multi_group <- function(model_setup) {
  valid_variables = c("mu_global","alpha_global",
                      "alpha_group","alpha_ind", "rho_global","rho_group","rho_ind") |>
    append(.get_family_priors(model_setup$settings$family))

  path <- system.file("extdata", "default_priors_multi_group.csv", package = "renewr")

  prior_frame <- .process_prior_frame(path = path,
                                      valid_variables = valid_variables,
                                      priors = model_setup$settings$priors)
  return(prior_frame)

}

#############
# Helpers
#############
#' @noRd
.get_family_priors <- function(family) {

  if (family == "gamma") return("k")
  else if (family == "weibull") return("shape")
  else if (family == "gengamma") return(c("k","shape"))
  else if (family == "lognormal") return(c("mu_lognormal","sigma_lognormal"))
}

#' @noRd
.process_prior_frame <- function(path, valid_variables, priors = NULL) {
  # e.g.
  # prior_set <- alist(
  #   z_pop ~ normal(0,1),
  #   z_group ~ normal(0,1),
  #   z_ind ~ normal(0,1),
  #
  #   mu_group ~ normal(0,1),
  #   mu_ind ~ normal(0,1),
  #
  #   alpha_pop ~ normal(0,1),
  #   alpha_group ~ normal(0,1),
  #   alpha_ind ~ normal(0,1),
  #
  #   rho_pop ~ lognormal(0.3,0.2),
  #   rho_group ~ lognormal(0.3,0.2),
  #   rho_ind ~ lognormal(0.3,0.2),
  #
  #   k ~ exp(0.5),
  #   shape ~ exp(0.5)
  # )
  dists <- c("normal", "lognormal", "cauchy", "invgamma", "beta","exp")
  prior_frame <- readr::read_csv(path, show_col_types = FALSE)

  if (missing(priors) || is.null(priors)) {
    prior_frame <- prior_frame |>
      mutate(distribution_id = match(distribution,dists))
    return(prior_frame)
  }

  for (f in priors) {
    formula <- as.formula(f)
    prior_variable <- as.character(formula[[2]])

    if (!(prior_variable %in% valid_variables)) {
      stop(paste0("Variable '", prior_variable, "' is not valid for this model type. ",
                  "Allowed parameters: ", paste(valid_variables, collapse = ", ")))
    }

    dist_data <- formula[[3]]
    dist_name <- as.character(dist_data[[1]])

    if (!(dist_name %in% dists)) {
      stop("Supported distributions: normal, lognormal, cauchy, or exp.")
    }

    dist_args <- lapply(as.list(dist_data)[-1], function(arg) {
      val <- tryCatch(eval(arg), error = function(e) arg)
      if (!is.numeric(val) || length(val) != 1) {
        stop(sprintf(
          "Prior parameter '%s' for variable '%s' must be a single numeric constant (negative values are allowed).",
          deparse(arg), prior_variable
        ))
      }
      val
    })
    #if (any(!sapply(dist_args, is.numeric))) {
    #  stop("Prior parameters must be numeric constants.")
    #}

    prior_frame[prior_frame$prior_variable == prior_variable, "distribution"] <- dist_name

    if (length(dist_args) == 1) {
      prior_frame[prior_frame$prior_variable == prior_variable, "param_1"] <- dist_args[[1]]
      prior_frame[prior_frame$prior_variable == prior_variable, "param_2"] <- NA
    } else {
      prior_frame[prior_frame$prior_variable == prior_variable, c("param_1", "param_2")] <- dist_args
    }
  }
  prior_frame <- prior_frame |>
    filter(prior_variable %in% valid_variables) |>
    mutate(distribution_id = match(distribution,dists))
  return(prior_frame)
}
