#' Plot fitted GP model
#' @export
#' @param model gp_model. A fitted GP model or model with prior checks
#' @param level Character. The level (`"ind"`, `"group"`, `"global"`) at which
#' the posterior fitted Gaussian Processes are extracted. By default, the level is `"ind"`,
#' thereby showing the additive combination of group-level and individual GPs translated by
#' their respective baselines, `mu_ind[i]`, where i is the individual ID.
#' @param .width Numeric vector. Highest density posterior credible interval bands around
#' the mean GP trace to visualise. Strictly, widths show pointwise quantiles around
#' the mean. By default, `.width = c(0.5, 0.8, 0.99)`, marking 50%, 80% and 99% density
#' coverage (see [ggdist::interval_widths()]).
#' @param n_samples Integer. Number of samples to draw to compute prior/posterior
#' predictive checks. By default this is 1000 and can be no greater than the number of
#' post-warmup iterations times the number of chains (see [gp_fit()]). Using the default
#' parameters means this ceiling is 4000 samples.
#' @param prior Boolean. Determines whether to show GPs drawn from prior or posterior
#' parameters. By default, `prior=FALSE` (see [gp_fit()]).
#' @param show_cri Boolean. Controls whether credible intervals are plotted or not.
#' If `FALSE`, by default, `facet` will also be set to `FALSE`, overlaying all plots
#' on the same axes. By default, however, show_cri is `TRUE`.
#' @param facet Boolean. Controls whether plots are overlaid onto the same axes or
#' faceted by the level specified. By default, this is `TRUE` but is set to `FALSE`
#' if `show_cri` is `FALSE`. This behaviour can be overridden by setting `facet=TRUE`.
#' @param show_events Boolean. Determines whether event data, if available are shown
#' as subplots beneath each plot of GP traces. By default, if event data are available,
#' they are plotted (`show_events=TRUE`).
#' @param dev_only Boolean. Determines whether only the deviations relative to
#' the level above are to be plotted. For example, if `dev_only=TRUE` and `level="ind"`
#' in a single-group model, only \eqn{\mu_{ind[i]} + f_{ind[i]}} would be shown
#' for the ith individual, rather than the composite \eqn{\mu_{ind[i]} + f_{ind[i]} + f_{group[i]}}.
#' @param show_traces Boolean. Determines whether simulate ground truth traces are
#' to be plotted, if available.
#' @param palette Character. Colour scheme for plotting credible intervals. The choice
#' of palette is passed to [ggplot2::scale_fill_brewer()] (see [RColorBrewer::display.brewer.all()] for
#' the full set of colour palettes).
#' @inheritParams plot.sim_traces
#' @param sample_alpha Numeric. The alpha level for random MCMC samples drawn from the
#' prior or posterior model parameters. The default alpha is `0.6`
#' @param sample_colour Character. The line colour of the random MCMC samples drawn from the
#' prior or posterior model parameters. The default colour is `"darkgrey"`
#' @param sample_linewidth Numeric. The linewidth for random MCMC samples drawn from the
#' prior or posterior model parameters. The default linewidth is `0.9`.
#' @param rescale Boolean. The observed traces represent the time-modulated scale, \eqn{\theta(t)},
#' of the renewal process. However, for non-memoryless processes, namely, the exponential
#' the expected inter-event interval is also specified by shape parameters. `rescale`
#' applies a scaling constant, \eqn{c}, to the fitted GP traces (and simulated ground truth traces if applicable)
#' that corresponds with the survival family.
#' \deqn{\mathbb{E}[X] = \exp{\beta(t)} c \text{, where } \beta(t) = \ln{\theta(t)}}
#' By default, the traces are shown on the response scale and hence the inverse ln
#' is applied to \eqn{c}:
#' \deqn{\mathbb{E}[X] = \beta(t) + \ln{c}}
#' Where for the gamma family, \eqn{c = k}.
#' For Weibull, \eqn{c = \Gamma(1 + \frac{1}{\text{shape}})}.
#' For the generalised gamma, \eqn{c = \frac{\frac{\Gamma(k + 1)}{\Gamma(shape)}}{\frac{\Gamma(k)}{\Gamma(shape)}}}
#' Rescaling ensures that the expected inter-event interval (accounting for scale and shape)
#'is shown. For [get_nlpd()], rescale is by default `TRUE`.
#' @method plot gp_model
plot.gp_model <- function(model, level = c("ind", "group", "global"), .width = c(0.50, 0.8, 0.99), n_samples = 0, prior = FALSE,
                          show_ci = TRUE, facet = show_ci, show_events = TRUE, dev_only = FALSE,
                          show_traces = TRUE, palette = "Purples", width = 0.1, height = 0.2, size = 10,
                          sample_alpha = 0.6, sample_colour = "darkgrey", sample_linewidth = 0.9, rescale = FALSE) {
  level <- match.arg(level)
  model_classes <- class(model)

  if (any(c("one_ind", "one_group") %in% model_classes) && level == "global") {
    stop("Global trace not present for models with one individual or group.")
  }

  # tidy_model <- reconstruct_traces(model,level=level,.width=.width,dev_only=dev_only,from_prior=prior,...)
  tidy_model <- tidy_traces(model,
    level = level,
    .width = .width, dev_only = dev_only, prior = prior, rescale = rescale
  )

  if (rescale) {
    scale_factor <- switch(model$sim_parameters$family,
      "exponential" = 1,
      "gamma" = model$sim_parameters$`k[1]`,
      "weibull" = gamma(1 + (1 / model$sim_parameters$`shape[1]`)),
      "gengamma" = gamma((model$sim_parameters$`k[1]` + 1) / model$sim_parameters$`shape[1]`) /
        gamma(model$sim_parameters$`k[1]` / model$sim_parameters$`shape[1]`)
    )
  } else {
    scale_factor <- 1
  }


  if (n_samples > 0) {
    max_samples <- model$stan_runtime$chains *
      (model$stan_runtime$iter - model$stan_runtime$warmup)

    if (n_samples > max_samples) {
      stop("Please choose a number of samples less than the total number of non-warmup iterations across all chains")
    }

    traces_sampled <- draw_traces(model,
      level = level,
      n_samples = n_samples, prior = prior, rescale = rescale
    )
  }

  y_elements <- c(tidy_model$y)

  # if (show_ci && !is.null(model$traces)) y_elements <- c(tidy_model$ymin, tidy_model$ymax)
  # if (show_traces && !is.null(model$traces)) y_elements <- c(y_elements, model$traces$y_offset)
  # global_y <- quantile(y_elements + log(scale_factor),c(0.05,0.95))

  global_x <- c(0, model$settings$duration)

  facet_var <- level
  group_aes <- level

  build_base_plot <- function(fit_sub, trace_sub = NULL, ev_sub = NULL, sample_sub = NULL) {
    p <- ggplot2::ggplot() +
      # ggplot2::coord_cartesian(xlim = global_x,
      #                         ylim = global_y) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Time (t)", y = "Modulated intensity")

    if (show_ci && !is.null(fit_sub) && nrow(fit_sub) > 0) {
      p <- p + ggdist::geom_lineribbon(
        data = fit_sub,
        ggplot2::aes(
          x = x, y = y, ymin = ymin, ymax = ymax,
          fill = factor(.width), group = .data[[group_aes]]
        ),
        alpha = 0.6, linewidth = 0.5
      ) +
        ggplot2::scale_fill_brewer(palette = palette, direction = -1, name = "CrI Width")
    } else if (!is.null(fit_sub) && nrow(fit_sub) > 0) {
      p <- p + ggplot2::geom_line(
        data = fit_sub,
        ggplot2::aes(x = x, y = y, group = .data[[group_aes]]),
        linewidth = 0.8
      )
    }

    if (show_traces && !is.null(trace_sub) && nrow(trace_sub) > 0) {
      p <- p + ggplot2::geom_line(
        data = trace_sub,
        ggplot2::aes(x = x, y = y_offset + log(scale_factor), group = .data[[facet_var]]),
        inherit.aes = FALSE, linewidth = 0.7, color = "black", alpha = 0.5
      )
    }

    if (n_samples > 0 && !is.null(sample_sub) && nrow(sample_sub) > 0) {
      p <- p + ggplot2::geom_line(
        data = sample_sub,
        ggplot2::aes(
          x = x, y = y,
          group = interaction(sample, .data[[facet_var]])
        ),
        inherit.aes = FALSE, linewidth = sample_linewidth,
        color = sample_colour, alpha = sample_alpha
      )
    }

    if (show_events && !is.null(ev_sub) && nrow(ev_sub) > 0) {
      if (facet_var == "global") {
        sub_facet_var <- "group"
      } else if (facet_var == "group") {
        sub_facet_var <- "ind"
      } else {
        sub_facet_var <- "ind"
      }
      y_aes_events <- "ind"


      event_aes <- ggplot2::aes(
        x = event_times,
        y = factor(.data[[y_aes_events]]),
        group = .data[[sub_facet_var]],
        color = factor(.data[[y_aes_events]]),
        fill = factor(.data[[y_aes_events]])
      )
      y_label <- "Events"
      y_label <- tools::toTitleCase(y_aes_events)


      p_below <- ggplot2::ggplot(data = ev_sub, mapping = event_aes) +
        ggplot2::coord_cartesian(xlim = global_x) +
        ggplot2::geom_tile(width = width, height = height) +
        ggplot2::scale_color_viridis_d(guide = "none") +
        ggplot2::scale_fill_viridis_d(guide = "none") +
        ggplot2::labs(x = "Time (t)", y = y_label) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text = ggplot2::element_text(size = size))

      if (facet && facet_var == "ind") {
        p_below <- p_below + ggplot2::theme(axis.text.y = ggplot2::element_blank())
      }

      return((p / p_below) + patchwork::plot_layout(heights = c(3, 1)))
    } else {
      return(p)
    }
  }

  if (!is.null(model$traces) & level == "ind") {
    trace_data <- model$traces
  } else if (!is.null(model$traces) & level == "group") {
    trace_data <- model$group_traces
  } else if (!is.null(model$traces) & level == "global") {
    trace_data <- model$global_trace
  } else {
    trace_data <- NULL
  }

  event_data <- if (!is.null(model$events)) model$events else NULL
  sample_data <- if (n_samples > 0) traces_sampled$sampled_traces else NULL


  if (!facet || level == "global") {
    combined_plot <- build_base_plot(tidy_model, trace_data, event_data, sample_data)
    return(combined_plot)
  } else {
    plots_list <- list()
    unique_facets <- unique(tidy_model[[facet_var]])

    for (f_id in unique_facets) {
      fit_sub <- tidy_model |> dply_filter_equal(facet_var, f_id)
      trace_sub <- trace_data |> dply_filter_equal(facet_var, f_id)
      ev_sub <- event_data |> dply_filter_equal(facet_var, f_id)
      sample_sub <- sample_data |> dply_filter_equal(facet_var, f_id)

      p_f <- build_base_plot(fit_sub, trace_sub, ev_sub, sample_sub)
      p_f <- p_f + ggplot2::facet_wrap(ggplot2::vars(.data[[facet_var]]), scales = "fixed")

      plots_list[[paste0("Facet_", f_id)]] <- p_f
    }

    combined_plot <- patchwork::wrap_plots(plots_list) +
      patchwork::plot_layout(guides = "collect")

    return(combined_plot)
  }
}

model_plot <- getS3method("plot", "gp_model")
.S3method("plot", "gp_prior_pc", model_plot)


#' @noRd
dply_filter_equal <- function(df, col, value) {
  if (is.null(df) || nrow(df) == 0 || !(col %in% names(df))) {
    return(NULL)
  }
  df[df[[col]] == value, , drop = FALSE]
}
