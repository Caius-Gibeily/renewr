#' Compute negative log predictive density between a fitted GP model and ground truth
#' traces
#' @inheritParams plot.gp_model
#' @export
get_accuracy_metric <- function(model, level = "ind", resolution = 0.1,
                     prior = FALSE, dev_only = FALSE, rescale = TRUE) {
  if (!rescale) warning("NLPD should be computed between traces with the same scale.")

  summary_traces <- summarise_traces(model,
    level = level, prior = prior,
    dev_only = dev_only,
    resolution = resolution,
    rescale = rescale
  )


  traces <- switch(level,
    "ind" = model$traces,
    "group" = model$group_traces,
    "global" = model$global_trace
  )
  if (rescale) {
    scale_factor <- switch(model$sim_params$family,
      "exponential" = 1,
      "gamma" = model$sim_params$survival_params$k,
      "weibull" = gamma(1 + (1 / model$sim_params$survival_params$shape)),
      "gengamma" = gamma((model$sim_params$survival_params$shape + 1) / model$sim_params$survival_params$k) /
        gamma(model$sim_params$survival_params$shape / model$sim_params$survival_params$k)
    )
  } else {
    scale_factor <- 1
  }

  if (level == "global") {
    summary_traces$y_ground <- approx(t = traces$t,
      eta = traces$eta,
      tout = summary_traces$t
    )$y
  } else {
    summary_traces <- summary_traces |>
      dplyr::group_by(.data[[level]]) |>
      dplyr::group_modify(~ {
        current_id <- .y[[level]]
        traces_sub <- dplyr::filter(traces, .data[[level]] == current_id)

        .x$eta_ground <- approx(x = traces_sub$t,y = traces_sub$eta,
          xout = .x$t)$y

        .x
      }) |>
      dplyr::ungroup()
  }

  summary_traces <- summary_traces |>
    tidyr::drop_na(eta_ground, means, sds) |>
    dplyr::mutate(dnorms = dnorm(eta_ground + log(scale_factor),
      means, sds,
      log = TRUE
    )) |>
    dplyr::group_by(.data[[level]]) |>
    dplyr::summarise(
      mean_NLPD = -mean(dnorms),
      corr = cor(eta_ground + log(scale_factor), means),
      .groups = "drop"
    )
  return(summary_traces)
}
