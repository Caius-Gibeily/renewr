#' @export
model_loo <- function(model) {
  if (!inherits(model, "gp_model")) {
    stop("Please provide a fitted GP model.")
  }

  log_lik_array <- loo::extract_log_lik(
    model$fit,
    parameter_name = "log_lik",
    merge_chains = FALSE
  )

  loo_res <- loo::loo(log_lik_array,
                      r_eff = 1,
                      cores = parallel::detectCores())
  return(loo_res)
}

#' @export
model_compare <- function(...) {
  models <- list(...)

  loo_list <- lapply(models, model_loo)

  comp <- loo::loo_compare(loo_list)
  return(comp)
}
