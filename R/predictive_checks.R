#' Plot prior or posterior predictive checks on fitted model
#' @description
#' Plot a series of posterior predictive checks based on [ppc_get_eventrate()],
#' [ppc_get_interevent_dist()] and [ppc_get_inst_eventrate()].
#'
#' @inheritParams mcmc_plot
#' @inheritParams plot.gp_model
#' @param bw_kde Numeric or character. The smoothing bandwidth to use for the interevent interval
#' distribution. The default is `"nrd0"` (see [stats::density()])
#' @param bw_rate Numeric. The smoothing bandwidth to use for
#' computing the instantaneous event rate
#' @param scale_factor Numeric. The factor to scale event rate and instantaneous
#' event rate, both of which are by default in events/s.
#' @param resolution Numeric. The resolution of the time grid for computing instantaneous
#' event rate (see [ppc_get_inst_eventrate()]). By default, this is 0.01
#' @param ... Optional plotting arguments.
#' @returns A series of ggplots.
#'
#' @export
ppc_plot_all <- function(model, n_samples = 1000, palette = "Blues",prior=FALSE,
                              .width = c(0.5,0.8,0.99), bw_kde = 0.5, bw_rate = 5,
                         scale_factor = 1, resolution = 0.01,...) {

  ppc_events <- ppc_draw_events(model,n_samples = n_samples,prior=prior)

  if (!("gp_model") %in% class(model)) stop("Please use a fitted renewr model")

  # Overall blink rate
  p1 <- ppc_get_eventrate(model = model,
                          scale_factor=scale_factor,...)
  # Interevent distribution
  p2 <- ppc_get_interevent_dist(ppc_events = ppc_events,
                                    model = model,bw_kde = bw_kde, ...)
  # Instantaneous blink rate
  p3 <- ppc_get_inst_eventrate(ppc_events = ppc_events,
                                   model = model, bw_rate = bw_rate,
                               resolution = resolution,
                               scale_factor = scale_factor, ...)

  .ggplot_successive(p1,p2,p3)
}

#' Generate events from prior or posterior predictive distribution of a GP model
#' @description
#' Taking n random MCMC draws from the prior or posterior, events are generated,
#' using the baseline, GP and survival parameters of the draw at the individual level
#' across the full duration of the time domain.
#' @inheritParams plot.gp_model
#' @returns Generated events at the individual level given the prior or posterior
#' model parameters.
#' @seealso [draw_traces()], [ppc_get_eventrate()], [ppc_get_interevent_dist()],
#' [ppc_get_inst_eventrate().
#' @export
ppc_draw_events <- function(model = NULL, sampled_traces = NULL,
                            n_samples, prior = FALSE, resolution = 0.1) {

  if (is.null(sampled_traces) && is.null(model)) {
    stop("Please provide either posterior generated events or a fitted renewr model")
  } else if (is.null(sampled_traces)) {
    ppc_data <- draw_traces(model, n_samples = n_samples, prior = prior,
                            resolution = resolution)
  } else ppc_data <- sampled_traces


  samples <- ppc_data$sampled_traces
  params <- ppc_data$survival_params

  family <- if (!is.null(params$k) & !is.null(params$shape)) {
    "gengamma"
  } else if (!is.null(params$k)) {
    "gamma"
  } else if (!is.null(params$shape)) {
    "weibull"
  } else if (!is.null(params$sigma)) {
    "lognormal"
  } else {
    "exponential"
  }

  ppc_events <-  do.call(simulate_events,
                         c(list(trace_data = samples, family = family, resolution = 0.1),
                                params))


  return(ppc_events$events)
}


#' Plot overall event rate from prior or posterior predictive distribution
#' @inheritParams ppc_plot_all
#' @param ppc_events tibble. Optional tibble of pre-drawn prior/posterior predictive
#' events generated via [ppc_draw_events()]. If ppc_events is `NULL`, the default,
#' `model` cannot be `NULL` and events will be generated via the function on the passed model.
#' @inheritParams plot.gp_model
#' @param return_plot Boolean. Specify whether to return a ggplot or a list of data.
#' By default, this is TRUE
#' @returns Either a ggplot or a list of observed
#' overall blink rates per individual and blink rate distributions via PPC.
#' @seealso [ppc_get_interevent_dist()], [ppc_get_inst_eventrate()], [ppc_plot_all()].
#' @export
ppc_get_eventrate <- function(model = NULL,ppc_events=NULL,scale_factor=1,.width=c(0.5,0.8,0.99),
                              palette = "Purples",return_plot=TRUE, n_samples=1000, prior = FALSE,...) {
  if (is.null(ppc_events)) {
    ppc_events <- ppc_draw_events(model,n_samples = n_samples, prior = prior)

  }

  br_global <- model$events |>
    dplyr::group_by(ind) |>
    dplyr::reframe(ind = unique(ind),
                   group = unique(group),
                   br_global = n() / model$settings$duration * scale_factor)

  br_global_ppc <- ppc_events |>
    dplyr::group_by(ind,sample) |>
    dplyr::reframe(sample = unique(sample),
                   ind = unique(ind),
                   group = unique(group),
                   br_global = n() / model$settings$duration * scale_factor)

  if (return_plot) {
    p <- ggplot2::ggplot() +
      ggdist::stat_interval(data = br_global_ppc, ggplot2::aes(x = factor(ind), y = br_global),
                            .width = .width, size = 6) +
      ggplot2::geom_segment(data = br_global,
                            ggplot2::aes(x = as.numeric(factor(ind)) - 0.2,
                                         xend = as.numeric(factor(ind)) + 0.2,
                                         y = br_global,yend = br_global),
                            color = "#014636",size = 1.2) +
      ggplot2::scale_color_brewer(palette = palette, name = "CrI Width") +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Individual",
                    y = "Event rate (blinks/min)",
                    color = "CrI")
    return(p)
  } else {
  return(list(br_global=br_global,
              br_global_ppc = br_global_ppc))
  }
}

#' Plot inter-event distributions from prior posterior predictive distribution
#' @inheritParams ppc_plot_all
#' @inheritParams ppc_get_eventrate
#' @inheritParams plot.gp_model
#' @returns A ggplot or a list of data
#' @seealso [ppc_get_eventrate()], [ppc_get_inst_eventrate()], [ppc_plot_all()].
#' @export
ppc_get_interevent_dist <- function(model = NULL, ppc_events = NULL, n_samples = 1000, palette = "Purples",
                                    .width = c(0.5, 0.8, 0.99), bw_kde = 1.2, return_plot = TRUE, prior = FALSE) {
  if (is.null(ppc_events) && is.null(model)) {
    stop("Please provide either posterior generated events or a fitted Renewr model")
  } else if (is.null(ppc_events)) {
    ppc_events <- ppc_draw_events(model, n_samples = n_samples, prior = prior)
  }

  grid_bounds <- ppc_events |>
    dplyr::group_by(ind) |>
    dplyr::summarise(to_val = stats::quantile(dt, 0.99, na.rm = TRUE), .groups = "drop")

  ppc_events <- dplyr::left_join(ppc_events, grid_bounds, by = "ind")

  iei_dist <- ppc_events |>
    dplyr::group_by(ind, sample) |>
    dplyr::reframe({
      d <- stats::density(dt, n = 512, from = 0, to = to_val[1], bw = bw_kde)
      tibble::tibble(x = d$x, y = d$y)
    }) |>
    dplyr::group_by(ind, x) |>
    ggdist::mean_qi(y, .width = .width) |>
    dplyr::ungroup()

  if (return_plot) {
    fd_bw <- .fd_binwidth(model$events$dt)

    p <- ggplot2::ggplot() +
      ggplot2::geom_histogram(
        data = model$events,
        ggplot2::aes(x = dt, y = ggplot2::after_stat(density), group = ind),
        fill = "darkgrey", color = "white", alpha = 1,
        binwidth = fd_bw
      ) +
      ggdist::geom_lineribbon(data = iei_dist,ggplot2::aes(x = x, y = y,
                                                           ymin = .lower, ymax = .upper,
          fill = forcats::fct_rev(ordered(.width)),
          group = interaction(ind, .width)),
        alpha = 0.6) +
      ggplot2::facet_wrap(~ ind, scales = "free_x") +
      ggplot2::scale_fill_brewer(palette = palette, name = "CrI Width") +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Inter-event Interval (IBI (s))", y = "Density")

    return(p)
  } else {
    return(iei_dist)
  }
}


#' Plot binned count of blinks across time domain
#' @inheritParams ppc_plot_all
#' @inheritParams ppc_get_eventrate
#' @inheritParams plot.gp_model
#' @returns A ggplot or a list of data
#' @seealso [ppc_get_eventrate()], [ppc_get_inst_eventrate()], [ppc_plot_all()].
#' @export
ppc_get_binned_events <- function(model = NULL,ppc_events = NULL,
                                  bw_rate = 8,resolution = 0.1, n_samples = 1000,
                                  .width = c(0.5,0.8,0.99), return_plot = TRUE,
                                  palette = "Purples",show_events=TRUE, prior = FALSE, bin_width = 10,
                                  linewidth=2,ylabel = "Number of events") {

  if (is.null(ppc_events) & is.null(model)) {
    stop("Please provide either posterior generated events or a fitted Renewr model")
  } else if (is.null(ppc_events)) {
    ppc_events <- ppc_draw_events(model,n_samples = n_samples, prior = prior)
  }
  if (return_plot & is.null(model)) {
    stop("Please provide a fitted renewr model for plotting")
  }
  time_grid <- seq(0,max(ppc_events$event_times),
                   by = bin_width)

  all_samples <- unique(ppc_events$sample)

  ppc_binned <- ppc_events |>

    dplyr::mutate(event_window = cut(event_times, breaks = time_grid)) |>
    dplyr::filter(event_times < max(event_times)) |>
    dplyr::count(ind, event_window, sample) |>
    dplyr::group_by(ind, event_window) |>
    tidyr::complete(sample = all_samples,
      fill = list(n = 0)) |>
    ggdist::mean_qi(n, .width = .width)

  if (return_plot) {
    p_ppc_binned <- ggplot2::ggplot(ppc_binned,
                                    aes(x = as.factor(event_window),
                                                   y = n,group=ind)) +
      ggdist::geom_interval(aes(ymin=.lower,ymax=.upper)) +
      #ggplot2::geom_line(color = "black",linewidth=linewidth) +
      ggplot2::scale_color_brewer(palette = palette, name = "CrI Width") +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Time bin (s)", y = ylabel,
                    fill = "CrI")

    events_binned <- model$events |>
      dplyr::group_by(ind) |>
      dplyr::mutate(event_window = cut(event_times,time_grid)) |>
      dplyr::filter(event_times < max(event_times)) |>
      dplyr::count(event_window,ind) |>
      dplyr::ungroup() |>
      tidyr::complete(event_window, ind, fill = list(n = 0)) |>
      dplyr::group_by(event_window,ind) |>
      ggdist::mean_qi(n,.width = .width)

    p_ppc_binned <- p_ppc_binned +
      geom_segment(data = events_binned, ggplot2::aes(xend=as.numeric(factor(event_window))+0.4,
                                                      x=as.numeric(factor(event_window))-0.4,
                                                      y = n), linewidth = linewidth) +
      facet_wrap(~ ind) +
      theme(axis.text.x = element_text(angle = 90,
                                       hjust = 1))

    return(p_ppc_binned)
  } else {
    return(ppc_binned)
  }

}

#' Plot instantaneous event rates from prior or posterior predictive distributions
#' @inheritParams ppc_plot_all
#' @inheritParams ppc_get_eventrate
#' @inheritParams plot.gp_model
#' @returns A ggplot or a list of data
#' @seealso [ppc_get_interevent_dist()], [ppc_get_eventrate()], [ppc_plot_all()].
#' @export
ppc_get_inst_eventrate <- function(model = NULL,ppc_events = NULL,
                                       bw_rate = 8,resolution = 0.1, n_samples = 2000,
                                   .width = c(0.5,0.8,0.99), return_plot = TRUE,
                                   palette = "Purples",show_events=TRUE, prior = FALSE) {

  if (is.null(ppc_events) & is.null(model)) {
    stop("Please provide either posterior generated events or a fitted Renewr model")
  } else if (is.null(ppc_events)) {
    ppc_events <- ppc_draw_events(model,n_samples = n_samples, prior = prior,
                                  resolution = resolution) |>
      dplyr::filter(event_times < max(event_times))
  }
  if (return_plot & is.null(model)) {
    stop("Please provide a fitted renewr model for plotting")
  }

  time_grid <- seq(0,max(ppc_events$event_times),
                   by = 1)

  inst_rate_ppc <- ppc_events |>
    dplyr::group_by(ind, sample) |>
    dplyr::summarise(
      rate = list(sapply(time_grid, function(t)
        sum(dnorm(t - event_times, sd = bw_rate)))), .groups = "drop") |>
    dplyr::mutate(x = list(time_grid)) |>
    tidyr::unnest(c(x, rate)) |>
    dplyr::group_by(ind, x) |>
    ggdist::mean_qi(rate, .width = .width) |>
    dplyr::ungroup()

  inst_rate <- model$events |>
    dplyr::filter(event_times < max(event_times)) |>
    dplyr::group_by(ind) |>
    dplyr::summarise(rate = list(sapply(time_grid, function(t)
      sum(dnorm(t - event_times, sd = bw_rate)))), .groups = "drop") |>
    dplyr::mutate(x = list(time_grid)) |>
    tidyr::unnest(c(x,rate))

  if (return_plot) {
    p_inst_event_rate <- ggplot2::ggplot(inst_rate_ppc) +
      ggdist::geom_lineribbon(aes(x = x, y = rate, ymin = .lower, ymax = .upper,
                                  fill = forcats::fct_rev(ordered(.width)), group = interaction(ind, .width)),
                              alpha = 0.8, linewidth = 0.6) +

      ggplot2::scale_fill_brewer(palette = palette, name = "CrI Width") +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Time", y = "Instantaneous event rate (s)",
                    fill = "CrI") +
      ggplot2::geom_line(data=inst_rate,
                         ggplot2::aes(x = x, y = rate, group=ind),
                         linewidth=1,color = "black", alpha = 0.5) +
      ggplot2::facet_wrap(~ind, scales = "free") +
      ggplot2::theme_minimal() +
      ggplot2::coord_cartesian(xlim = c(0, model$settings$duration))
    if (!show_events) return(p_inst_event_rate)

      p_ob_events <- ggplot(model$events, ggplot2::aes(x = event_times, y = factor(1), group = ind)) +
      ggplot2::geom_tile(width = 0.3, height = 0.5) +
      ggplot2::labs(x = "Time", y = "Ind") +
      ggplot2::facet_wrap(~ind) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     strip.text = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank()) +
      ggplot2::coord_cartesian(xlim = c(0, model$settings$duration))

    p_ppc_events <- ggplot2::ggplot(ppc_events, ggplot2::aes(x = event_times, y = factor(sample),
                                                                     group = interaction(ind,sample))) +
      ggplot2::geom_tile(width = 0.5, height = 1) +
      ggplot2::labs(x = "Time", y = "Ind") +
      ggplot2::facet_wrap(~ind) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     strip.text = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank()) +
      ggplot2::coord_cartesian(xlim = c(0, model$settings$duration))

    p <- cowplot::plot_grid(p_inst_event_rate, p_ob_events,
                             p_ppc_events , align = "v",
                             ncol = 1, rel_heights = c(1.5,0.4,0.8),
                             rel_widths = c(1,0.5,1))
    return(p)
  } else {
  return(list(inst_rate_ppc = inst_rate_ppc, inst_rate = inst_rate))
  }
}

#' Plot empirical cumulative density functions from prior or posterior predictive distributions
#' @inheritParams ppc_plot_all
#' @inheritParams ppc_get_eventrate
#' @inheritParams plot.gp_model
#' @returns A ggplot or tibble of data
#' @seealso [ppc_get_interevent_dist()], [ppc_get_eventrate()], [ppc_plot_all()].
#' @export
ppc_get_ecdf <-  function(model = NULL, ppc_events = NULL, n_samples = 1000,
                          plot_density = TRUE, spaghetti = 0, palette = "Purples",.width = c(0.5, 0.8, 0.99),
                          bw_kde = 1.2, return_plot = TRUE, prior = FALSE,linewidth=1,color="blue",
                          resolution = 0.1) {
  if (is.null(ppc_events) & is.null(model)) {
    stop("Please provide either posterior generated events or a fitted Renewr model")
  } else if (is.null(ppc_events)) {
    ppc_events <- ppc_draw_events(model,n_samples = n_samples, prior = prior)
  }
  if (return_plot & is.null(model)) {
    stop("Please provide a fitted renewr model for plotting")
  }

  max_t <- max(ppc_events$event_times)
  ppc_events <- ppc_events[ppc_events$event_times != max(ppc_events$event_times),]

  t_grid <- seq(0,max_t, by = resolution)

  if (plot_density & spaghetti == 0) {
    ppc_ecdf <- ppc_events |>
      dplyr::group_by(ind, sample) |>
      dplyr::reframe(
        event_t_grid = t_grid,
        ecdf_each = ecdf(event_times)(t_grid)
      ) |>
      dplyr::group_by(ind, event_t_grid) |>
      ggdist::median_qi(ecdf_each, .width = .width)

    if (!return_plot) {
      return(ppc_ecdf)
    }

    p_ecdf <- ggplot(ppc_ecdf,aes(x = event_t_grid, y = ecdf_each)) +
      ggdist::geom_lineribbon(data = ppc_ecdf,ggplot2::aes(ymin = .lower,
                                                           ymax = .upper)) +
      ggplot2::scale_fill_brewer(palette = palette, name = "CrI Width") +
      ggplot2::theme_minimal() +
      facet_wrap(~ind)

  } else {
    if (plot_density & spaghetti == 0) {
      stop("Please set spaghetti, the number of posterior predictive ecdf draws to show
           to be greater than 0")
    }
    ppc_ecdf <- ppc_events |>
      dplyr::filter(sample %in% sample(1:n_samples,spaghetti)) |>
      dplyr::group_by(ind, sample) |>
      dplyr::reframe(
        event_t_grid = t_grid,
        ecdf_each = ecdf(event_times)(t_grid)
      )
    if (!return_plot) {
      return(ppc_ecdf)
    }
    p_ecdf <- ggplot(ppc_ecdf,aes(x = event_t_grid, y = ecdf_each,
                                  group = interaction(sample,ind),
                                  color = ind)) +
      ggplot2::geom_path(alpha=0.4,linewidth = max(linewidth / spaghetti * 5,0.1)) +
      ggplot2::theme_minimal() +
      facet_wrap(~ind)
  }
  p_ecdf <- p_ecdf + ggplot2::stat_ecdf(data = model$events,
                                        ggplot2::aes(x = event_times,group = ind),
                       geom = "step",
                       color = color,linewidth=linewidth,inherit.aes=FALSE) +
    ggplot2::facet_wrap(~ind) +
    ggplot2::labs(x = "Time", y = "Cumulative proportion of events",
                  fill = "CrI")
  return(p_ecdf)

}

#' @noRd
.fd_binwidth <- function(x) {
  2 * IQR(x, na.rm = TRUE) / (length(x)^(1/3))
}

#' @noRd
.ggplot_successive <- function(...) {
  input_list <- list(...)

  if (length(input_list) == 1 && is.list(input_list[[1]])) {
    plots <- input_list[[1]]
  } else {
    plots <- input_list
  }

  for (i in seq_along(plots)) {
    print(plots[[i]])

    if (i < length(plots)) {
      res <- readline(prompt = "Hit <Return> to see next plot (or 'q' to quit): ")
      if (tolower(trimws(res)) == "q") break
    }
  }
}


