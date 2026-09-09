#' Simulate modulated renewal processes from a collection of hierarchical traces
#' @description
#' A short description...
#'
#' @param trace_data A `sim_traces` object containing simulated traces, created
#' by any of the simulate_* functions (see [simulate_spline_traces()],
#' [simulate_gp_traces()], [simulate_boxcar_traces()]). An external
#' data.frame or tibble may also be passed.
#' @param group Character. The column name containing group ID. By default, this is
#' set to `"group"` and need not be modified if passing a `sim_traces` object.
#' @param ind Character. The column name containing group ID. By default, this is
#' set to `"ind"` and need not be modified if passing a `sim_traces` object.
#' @param family Character. Choice of survival family distribution to draw events from.
#' One of "exponential", "gamma", "weibull", "gengamma" (generalised gamma) and
#' "lognormal". The default is "gamma".
#' @param shape Numeric. The shape parameter of the Weibull and generalised gamma families.
#' The generalised gamma, in its 1962 Stacy parameterisation, takes both shape and k
#' parameters with the scale (equal to 1/rate)
#' @param k Numeric. The shape parameter of the gamma and generalised gamma family.
#' The generalised gamma, in its 1962 Stacy parameterisation, takes both shape and k
#' parameters with the scale (equal to 1/rate), which is modelled as a time-varying
#' smooth
#' @param sigma Numeric. The standard deviation of the log-normal survival family.
#' @param Q Numeric. In the modern generalised gamma parameterisation, Q is the shape
#' parameter, alongside location parameter, mu (modelled as a time-varying smooth) and scale,
#' sigma.
#' @param baseline Numeric or numeric vector. The vertical offset or \eqn{\mu} to apply to simulated
#' traces. If one value is supplied, the same baseline is applied to all individuals.
#' If the length of a passed vector baseline inputs equals the number of groups,
#' individual traces in each group will be adjusted by the respective value in the order
#' passed. Conversely, if the length equals the total number of individuals, each individual
#' trace will be translated by the respective value in the order of baselines passed and
#' order of individual IDs. # note: ensure this is correct
#' @param resolution Numeric. Time grid resolution for approximating the cumulative
#' hazard function. By default, linear interpolation is performed at a resolution of 0.005.
#' If the passed traces are already sufficiently time resolved, `resolution = NULL`
#' may be passed to skip linear interpolation.
#' @param seed Numeric. Set a seed for reproducibility. If one is not specified
#' but seed was specified during trace simulation via the simulate_*_traces functions,
#' seed is set to the same value. Otherwise, seed is NULL.
#' @returns An object of class `sim_events` containing all data from `sim_traces` object,
#' if passed, renewal process simulation parameters and simulated event data.
#' @examples
#' gp_traces <- simulate_gp_traces(n_ind = 5, alpha_group = 0.5, rho_group = 10,
#' seed = 123)
#' events <- simulate_events(gp_traces, family = "weibull")
#' plot(events)

#' @export
simulate_events <- function(trace_data, group = group, ind = ind, family = c(
                              "exponential", "gamma", "weibull",
                              "log-normal", "gengamma"
                            ), shape = 1, k = 2, sigma = 1,
                            Q = 0, baseline = 0, resolution = 0.01, seed = NULL) {
  # extract trace data from list if the object is already a list
  if (is.list(trace_data) & !is.data.frame(trace_data)) {
    traces <- trace_data$traces
  } # else, extract trace from dataframe and convert to a list container
  else if (is.data.frame(trace_data)) {
    traces <- trace_data
    trace_data <- list(traces = traces)
  }

  if (is.null(seed)) {
    if (!is.null(trace_data$sim_parameters$seed_traces)) {
      seed <- trace_data$sim_parameters$seed_traces
      set.seed(seed)
      warning("Setting seed to the one set in trace simulation. If this is unintended, please set seed in simulate_events().")
    } else {
      set.seed(seed)
    }
  } else {
    set.seed(seed)
  }

  trace_data$sim_parameters$seed_events <- seed
  trace_data$sim_parameters$family <- family

  sim_struct <- trace_data$traces <- traces |>
    group_by({{ group }}, {{ ind }}) |>
    group_keys()
  n_ind <- nrow(sim_struct)
  n_group <- unique(sim_struct$group) |> length()
  if (length(baseline) == 1) {
    # if a global mu offset is specified, translate all traces by this fixed value

    traces <- trace_data$traces <- traces |>
      mutate(
        y_offset = y + baseline,
        scale = 1 / exp(-(y_offset))
      )
    if (!is.null(trace_data$group_traces)) {
      # if the group trace is also present (multi-ind or multi-group models), add baseline to group traces
      trace_data$group_traces <- trace_data$group_traces |>
        mutate(y_offset = y + baseline)
    }
    if (!is.null(trace_data$global_traces)) {
      # if the global trace is also present (multi-group model), add baseline to group traces
      trace_data$global_trace <- trace_data$global_trace |>
        mutate(y_offset = y + baseline)

      trace_data$sim_parameters$mu_global <- baseline
    }

    mu_inds <- as.list(setNames(rep(baseline, n_ind), paste0("mu_ind[", sim_struct$ind, "]")))
    mu_groups <- as.list(setNames(
      rep(baseline, n_group),
      paste0(
        "mu_group[",
        unique(sim_struct$group), "]"
      )
    ))
  } else if (length(baseline) == n_ind) {
    # if the number of offsets supplied equals the number of individuals, add each
    # separately to each individual
    traces <- trace_data$traces <- traces |>
      group_by(ind) |>
      mutate(
        y_offset = y + baseline[ind],
        scale = exp(y_offset)
      ) |>
      dplyr::ungroup()

    # group baselines
    baseline_groups <- sim_struct |>
      dplyr::bind_cols(baseline = baseline) |>
      group_by(group) |>
      summarise(baseline = mean(baseline)) |>
      dplyr::pull(baseline)

    trace_data$group_traces <- trace_data$group_traces |>
      group_by(group) |>
      mutate(y_offset = y + baseline_groups[group]) |>
      dplyr::ungroup()

    mu_inds <- as.list(setNames(baseline, paste0("mu_ind[", sim_struct$ind, "]")))
    mu_groups <- as.list(setNames(
      baseline_groups,
      paste0(
        "mu_group[",
        unique(sim_struct$group), "]"
      )
    ))

    if (!is.null(trace_data$global_trace)) {
      baseline_global <- mean(baseline)

      trace_data$global_trace <- trace_data$global_trace |>
        mutate(y_offset = y + baseline_global)

      trace_data$sim_parameters$mu_global <- baseline_global
    }
  } else {
    stop("Please ensure a single baseline offset or a number of baseline offsets equal to the number of individuals")
  }

  trace_data$sim_parameters <- trace_data$sim_parameters |>
    append(c(mu_inds, mu_groups))


  if (family != "log-normal") {
    if (family == "exponential") {
      shape <- k <- 1
    } else if (family == "gamma") {
      shape <- 1
    } else if (family == "weibull") k <- 1

    event_times <- traces |>
      group_split(ind) |>
      map(\(x) .simulate_renewal(x,
        modulant = scale,
        shape = shape, k = k, resolution = resolution
      )) |>
      list_rbind(names_to = "ind") |>
      mutate(ind = as.factor(ind)) |>
      left_join(sim_struct, by = "ind") |>
      relocate(group, .after = ind)

    trace_data$sim_parameters <- trace_data$sim_parameters |>
      append(list(`shape[1]` = shape, `k[1]` = k))
  } else if (family == "log-normal") {
    if (is.null(mu)) {
      traces <- traces |>
        mutate(mu = log(scale) + log(k) / sqrt(shape))
    }

    event_times <- traces |>
      group_split(ind) |>
      map(\(x) .simulate_renewal(x,
        modulant = mu,
        sigma = sigma, Q = 0, resolution = resolution
      )) |>
      list_rbind(names_to = "ind") |>
      mutate(ind = as.factor(ind)) |>
      left_join(sim_struct, by = "ind") |>
      relocate(group, .after = ind)
  } else {
    stop(
      "Please choose a survival function from ",
      "exponential, gamma, weibull, log-normal, gengamma"
    )
  }
  if (is.list(trace_data)) {
    sim_data <- trace_data %>% append(
      list(events = event_times)
    )
    class(sim_data) <- c("sim_traces", class(sim_data))
  } else {
    sim_data <- list(
      traces = trace_data,
      events = event_times
    )
  }
  class(sim_data) <- c("sim_events", class(sim_data))
  return(sim_data)
}

#' @noRd
.simulate_renewal <- function(trace, modulant, shape, k, sigma, Q, resolution = 200) {
  if (!is.null(resolution)) {
    lerp <- 1 / resolution

    trace_parts <- trace |> pull({{ modulant }})
    trace_parts <- trace_parts |> approx(n = max(trace$x) * lerp)

    modulant <- trace_parts$y
    time <- trace_parts$x
  } else {
    modulant <- trace |> pull({{ modulant }})
    time <- trace$x
  }

  if (!missing(shape) && !missing(k)) {
    events <- simulate_renewal_multi_omp(time, modulant, 1, shape, k)
    print(events)
  } else if (!missing(sigma) && !missing(Q)) {
    events <- simulate_renewal(time, modulant, sigma, Q)
  }
  renewal_events <- tibble(
    event_times = events,
    dt = diff(c(0, events))
  ) |>
    dplyr::filter(event_times > 0)
  return(renewal_events)
}
