#' Plot Traces of Simulated Data
#'
#' @param sim_data An object of class `sim_traces` or `sim_events`.
#' @param level Character vector specifying which levels to plot.
#'   Options include `"ind"` (individual), `"group"`, or `"global"`. Defaults to all three.
#' @param facet Logical. If `TRUE`, splits plots into separate panels by group or
#' individual.
#' @param width Numeric. Width of the event rasters. The default is 0.1. Passed to [ggplot2::geom_tile()]
#' @param height Numeric. Height of the event rasters. The default is 0.2. Passed to [ggplot2::geom_tile()]
#' @param size Numeric. Axis label size. Default is 10
#' @returns A ggplot showing plotted traces at the level(s) specified and event raster
#' subplots faceted by individual or group, if specified.
#' @examples
#' gp_traces <- simulate_gp_traces(n_ind = 5, alpha_group = 0.5, rho_group = 10,
#' seed = 123)
#' # Plot traces
#' plot(gp_traces)
#'
#' # Simulate a modulated point process drawn from a Weibull and plot the traces + events
#' events <- simulate_events(gp_traces, family = "weibull")
#' plot(events)
#'
#' @export
#' @method plot sim_traces
plot.sim_traces <- function(sim_data, level = c("ind", "group", "global"),
                            facet = TRUE, width = 0.1,
                            height = 0.2, size = 10) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for stacking subplots.")
  }

  plots_list <- list()
  level <- match.arg(level, several.ok = TRUE)

  ind_data <- sim_data$traces$ind_traces
  group_data <- sim_data$traces$group_traces
  global_data <- sim_data$traces$global_trace


  event_data <- NULL
  if ("sim_events" %in% class(sim_data) && !is.null(sim_data$events)) {
    event_data <- sim_data$events
  }

  n_groups <- sim_data$sim_parameters$n_groups
  purple_palette <- colorRampPalette(brewer.pal(9, "Purples"))(n_groups)

  build_base_plot <- function(ind_sub = NULL, grp_sub = NULL, ev_sub = NULL, title_suffix = "") {
    p <- ggplot() +
      theme_minimal() +
      labs(x = "Time (t)", y = "Eta(t)")

    if ("ind" %in% level && !is.null(ind_sub) && nrow(ind_sub) > 0) {
      i_x <- "t"
      i_y <- "eta"

      p <- p + geom_line(
        data = ind_sub,
        aes(
          x = .data[[i_x]], y = .data[[i_y]], group = ind,
          color = as.factor(ind)
        ), alpha = 0.6, linewidth = 0.5
      ) +
        scale_color_viridis_d(guide = "none")
    }

    if ("group" %in% level && !is.null(grp_sub) && nrow(grp_sub) > 0) {
      g_x <- "t"
      g_y <- "eta"

      p <- p + geom_line(
        data = grp_sub,
        aes(x = .data[[g_x]], y = .data[[g_y]], group = group),
        color = "black", linewidth = 1
      ) +
        labs(color = "Group")
    }

    if ("global" %in% level && !is.null(global_data) && nrow(global_data) > 0) {
      gl_x <- "t"
      gl_y <- "eta"

      p <- p + geom_line(
        data = global_data,
        aes(x = .data[[gl_x]], y = .data[[gl_y]]),
        color = "black", linewidth = 1.4, linetype = "dashed"
      )
    }

    if (title_suffix != "") {
      p <- p + ggtitle(title_suffix)
    }

    if (!is.null(ev_sub) && nrow(ev_sub) > 0) {
      p_below <- ggplot(
        data = ev_sub |> filter(censored != 1),
        aes(x = event_times,y = factor(ind),group = interaction(ind, group),
            color = interaction(ind, group),fill = interaction(ind, group))) +
        geom_tile(width = width, height = height) +
        scale_color_viridis_d(guide = "none") +
        scale_fill_viridis_d(guide = "none") +
        labs(x = "Time (t)", y = "Individual", ) +
        theme(axis.text = element_text(size = size)) +
        theme_minimal()

      # Use patchwork to stack them cleanly, allocating more space to the upper trace
      return(p / p_below + plot_layout(heights = c(2, 1)))
    }

    return(p)
  }


  unique_groups <- unique(ind_data$group)

  if (!facet) {
    p <- build_base_plot(ind_data, group_data,
                         event_data)

  } else {
    for (g in unique_groups) {
      i_sub <- ind_data[ind_data$group == g,]
      g_sub <- group_data[group_data$group == g,]

      ev_sub <- NULL
      if (!is.null(event_data)) {
        ev_sub <- event_data |> dplyr::filter(group == g)
      }

      p_g <- build_base_plot(i_sub, g_sub, ev_sub, title_suffix = paste("Group:", g))
      plots_list[[paste0("Group_", g)]] <- p_g
    }


    combined_plot <- patchwork::wrap_plots(plots_list, ) +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(title = "Trace Simulation")

    return(combined_plot)

  }
}
