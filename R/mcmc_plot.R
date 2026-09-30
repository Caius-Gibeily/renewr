

#' Plot estimates from MCMC draws
#' @description
#' A wrapper function on **bayesplot** to simplify visualisation of prior or posterior estimates from
#' MCMC draws.
#'
#' @param model gp_model. A fitted GP model or model with prior checks
#' @param pars Character or character vector. Set of parameters to plot, comprising
#' intercept parameters mu; GP parameters alpha and rho; and survival parameters k,
#' shape and sigma, if applicable. Characters may specify full parameter by name,
#' i.e. `"mu_group"` to plot all group intercepts or generically as `"mu"` to plot
#' all intercepts. By default, all parameters are shown.
#' @param prior Boolean. Plot prior estimates, if available, instead of posterior estimates.
#' The default is `FALSE`
#' @param type Character. Choice of output visualisation from `bayesplot::mcmc_` function set.
#' The type can be `"areas"` (see [bayesplot::mcmc_areas()]), `"dens"` (see [bayesplot::mcmc_dens()]),
#' `"dens_chains"` (see [bayesplot::mcmc_dens_chains()]), `"hist"` (see [bayesplot::mcmc_hist()]),
#' `"intervals"` (see [bayesplot::mcmc_intervals()]) and `"trace"` (see [bayesplot::mcmc_trace()]).
#' The default is `"intervals"`.
#' @param scheme Character. Colour scheme of plot. Option is passed to [bayesplot::color_scheme_set()].
#' The default is `"blue"`.
#' @param show_ground Boolean. If the model was fitted on simulated ground truth data, i.e.
#' [simulate_gp_traces()], setting `show_ground=TRUE` will show the respective parameters
#' underlying the simulated data as vertical lines. By default, this is true but will
#' not raise an error ground truth data are not available.
#' @param height Numeric. Height of the vertical line signposting the ground truth
#' parameter values if `show_ground=TRUE`.
#' @param linewidth Numeric. Width of the vertical line signposting the ground truth parameter
#' values if `show_ground=TRUE`.
#' @param color Character. The colour of the vertical line signposting the ground truth parameter
#' values if `show_ground=TRUE`.
#' @param ... Optional arguments passed to the selected mcmc_* plotting function.
#' @returns NULL
#'
#' @examples
#' # Simulate traces
#' gp_traces <- simulate_gp_traces(n_ind = 1, alpha_ind = 0.5,
#'     rho_ind = 10,
#'     seed = 123)
#' # Simulate events
#' events <- simulate_events(gp_traces, family = "weibull", baseline = 0.5)
#'
#' # Fit model
#'
#' fit <- gp_fit(events, family = "gamma",
#'     priors = alist(alpha_ind ~ normal(0.5, 0.3),
#'     rho_ind ~ lognormal(2.2, 0.4),
#'     mu_ind ~ normal(0.5, 0.2),
#'     shape ~ normal(1.5, 0.3)
#'     ))
#' # Plot mcmc intervals of the core model parameters
#' mcmc_plot(fit)
#' # Plot densities by chain to check chain convergence, setting red colour scheme
#' mcmc_plot(fit, type = "dens_chains", scheme = "red")
#' # Remove ground truth lines from plot
#' mcmc_plot(fit, show_ground = FALSE)

#' @seealso [bayesplot::mcmc_areas()], [bayesplot::mcmc_dens()],
#' [bayesplot::mcmc_dens_chains()], [bayesplot::mcmc_hist()],
#' [bayesplot::mcmc_intervals()]) and [bayesplot::mcmc_trace()]
#' @export
mcmc_plot <- function(model,pars=NULL,prior=FALSE,
                               type = c("areas", "dens", "dens_chains", "hist", "intervals", "trace"),
                      scheme="blue",show_ground = TRUE,height=0.3, linewidth=1,color="black",page="split",...) {
  if (length(type) > 1) type = "intervals" # default option

  type <- match.arg(type)

  if (!prior & !is.null(model$fit)) {
    model_dat <- model$fit
  } else if (prior & !is.null(model$prior_pc) | "gp_prior_pc" %in% class(model)) {
    model_dat <- model$prior_pc
  }

  all_pars <- names(model_dat)

  pars_filt <- grep(c("z_|beta|eta|log_|diag|mu_raw_|lp_"), all_pars, value = TRUE, invert = TRUE)

  pars_types <- list(mu = c("mu"), alpha = c("alpha"), rho = c("rho"), survival = c("k|shape"))
  if (page == "split") {
    if (is.null(pars)) {
      sub_pars <- lapply(pars_types, function(type) {
        grep(type, pars_filt, value = TRUE, invert = FALSE)
      })
    } else {
      matching_pars <- lapply(pars, function(type) {
        grep(type, pars_filt, value = TRUE, invert = FALSE)
      })
      sub_pars <- lapply(pars_types, function(type) {
        grep(type, purrr::flatten(matching_pars), value = TRUE, invert = FALSE)
      })
    }
  } else if (page == "single") {
    if (!is.null(pars)) {
      matching_pars <- lapply(pars, function(type) {
        grep(type, pars_filt, value = TRUE, invert = FALSE)
      }) |> list_c()
      sub_pars <- matching_pars
    } else {
      sub_pars <- pars_filt
    }

  }

  plot_fn <- switch(type,
                    "areas" = bayesplot::mcmc_areas,
                    "dens" = bayesplot::mcmc_dens,
                    "dens_chains" = bayesplot::mcmc_dens_chains,
                    "hist" = bayesplot::mcmc_hist,
                    "intervals" = bayesplot::mcmc_intervals,
                    "trace" = bayesplot::mcmc_trace,
                    stop("Unsupported plot type")
  )

  bayesplot::color_scheme_set(scheme = scheme)
  all_sim <- do.call(c,model$sim_params)

  if (page == "split") {
    plots <- lapply(sub_pars, function(p) {
      if (!rlang::is_empty(p)) {
        p_layer <- plot_fn(model_dat, pars = p, ...) +
          ggplot2::theme_minimal()
        pattern <- paste(gsub("[^A-Za-z_]","",p), collapse = "|")

        params_matched <- grep(pattern, names(all_sim))

        if (show_ground & !is.null(model$sim_params) & length(params_matched) != 0) {

          params <- tibble(
            parameter = p,
            xintercept = unlist(all_sim[params_matched]))

          y_levels <- ggplot2::ggplot_build(p_layer)$layout$panel_params[[1]]$y$get_labels()

          params <- params |>
            dplyr::filter(parameter %in% y_levels) |>
            dplyr::mutate(parameter = factor(parameter, levels = y_levels))


          p_layer +
            ggplot2::geom_linerange(
              data = params,
              aes(x = xintercept, y = parameter, ymin = as.numeric(parameter) - height,
                  ymax = as.numeric(parameter) + height),
              inherit.aes = FALSE,
              color = color, linewidth = linewidth
            )

        } else p_layer
      }

    }) |> purrr::discard(is.null)

    .ggplot_successive(plots)
  } else if (page == "single") {

      if (!rlang::is_empty(sub_pars)) {
        p_layer <- plot_fn(model_dat, pars = sub_pars, ...) +
          ggplot2::theme_minimal()

        clean_sub <- gsub("[^A-Za-z_]", "", sub_pars)
        clean_sim_names <- gsub("[^A-Za-z_]|^.*?\\.", "", names(all_sim))

        matched_idx <- match(clean_sub, clean_sim_names)

        valid <- !is.na(matched_idx)
        params_names <- sub_pars[valid]
        params_values <- matched_idx[valid]

        if (show_ground && !is.null(model$sim_params) && length(params_names) != 0) {

          y_levels <- ggplot2::ggplot_build(p_layer)$layout$panel_params[[1]]$y$get_labels()


          params <- tibble::tibble(
            parameter  = params_names,
            xintercept = unlist(all_sim[params_values])
          ) |>
            dplyr::filter(parameter %in% y_levels) |>
            dplyr::mutate(parameter = factor(parameter, levels = y_levels))

          p_layer <- p_layer +
            ggplot2::geom_linerange(
              data = params,
              ggplot2::aes(x = xintercept, y = parameter,
                           ymin = 1 - height, ymax = 1 + height),
              inherit.aes = FALSE,
              color = color, linewidth = linewidth
            )
        }
        ncol <- if (length(sub_pars) <= 2) {
          1} else length(sub_pars) %/% 2
        p_layer <- p_layer +
          ggplot2::scale_y_discrete(limits = NULL) +
          ggplot2::facet_wrap(~ parameter,ncol = ncol,
            scales = "free") + ggplot2::theme(
            axis.text.y  = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
          )

        } else p_layer

    return(p_layer)
    }

}
