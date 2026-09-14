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
#' gp_traces <- simulate_gp_traces(
#'   n_ind = 5, alpha_group = 0.5, rho_group = 10,
#'   seed = 123
#' )
#' events <- simulate_events(gp_traces, family = "weibull")
#' plot(events)

#' @export
simulate_events <- function(
  trace_data, group = group, ind = ind, family = c("exponential", "gamma", "weibull",
                                                   "log-normal", "gengamma"),
  shape = 1, k = 2,sigma = 1, Q = 0,resolution = 0.01, seed_events = NULL) {
  if ("sim_traces" %in% class(trace_data)) {
    traces <- trace_data$traces$ind_traces
  } else if (is.data.frame(trace_data)) {
    traces <- trace_data
    trace_data <- list(traces = traces)
  }
  # for subsequent compatibility
  traces <- rename(traces,
    ind = ind,
    group = group
  )

  set.seed(seed_events)

  trace_data$sim_params$seed_events <- seed_events
  trace_data$sim_params$family <- family

  sim_struct <- unique(traces[c("ind", "group")])

  n_ind <- nrow(sim_struct)
  n_group <- max(unique(sim_struct$group))

  # acquire scale

  traces$scale <- exp(traces$eta)
  traces$rate <- 1 / traces$scale

  if (family != "log-normal") {
    if (family == "exponential") {
      shape <- k <- 1
    } else if (family == "gamma") {
      shape <- 1
    } else if (family == "weibull") {
      k <- 1
    }
  }

  trace_data$sim_params$survival_params <- list(
    shape = shape,
    k = k
  )

  ind_ids <- sim_struct$ind

  n_samples <- max(length(shape), length(k))
  use_samples <- n_samples > 1

  shape_vec <- rep_len(shape, n_samples)
  k_vec <- rep_len(k, n_samples)

  max_t <- if (!is.null(trace_data$sim_params$duration)) {
    trace_data$sim_params$duration
  } else {
    max(traces$t)
  }

  time_vec <- seq(0, max_t, by = resolution)
  n_time <- length(time_vec)

  traces_have_samples <- "sample" %in% names(traces)

  if (traces_have_samples) {
    sample_ids <- sort(unique(traces$sample))

    if (use_samples && length(sample_ids) != n_samples) {
      stop(
        "Number of unique `sample` values in traces (", length(sample_ids),
        ") does not match n_samples inferred from shape/k (", n_samples, ")."
      )
    }
    if (!use_samples && length(sample_ids) > 1) {
      stop(
        "traces contains multiple `sample` values but shape/k were supplied ",
        "as scalars -- please supply per-sample shape/k, or collapse traces ",
        "to a single sample before calling simulate_events()."
      )
    }


    traces_split <- split(traces,
      list(traces$sample, traces$ind),
      sep = "___"
    )
    ind_scale_by_sample <- lapply(sample_ids, function(s) {
      vapply(ind_ids, function(i) {
        key <- paste(s,i, sep = "___")
        tr <- traces_split[[key]]
        if (is.null(tr)) {
          stop(
            "No trace rows found for sample = ", s,
            ", ind = ", i, "."
          )
        }
        #tr <- tr[order(tr$t), ]
        stats::approx(tr$t, tr$scale, xout = time_vec, rule = 2)$y
      }, FUN.VALUE = numeric(n_time))
    })

    modulant_mat_flat <- as.vector(do.call(cbind, ind_scale_by_sample))

  } else {
    traces_split <- split(traces, traces$ind)

    ind_scale_mat <- vapply(ind_ids, function(i) {
      tr <- traces_split[[as.character(i)]]
      tr <- tr[order(tr$t), ]
      stats::approx(tr$t, tr$scale, xout = time_vec, rule = 2)$y
    }, FUN.VALUE = numeric(n_time))

    modulant_mat_flat <- if (use_samples) {
      rep(as.vector(ind_scale_mat), times = n_samples)
    } else {
      as.vector(ind_scale_mat)
    }
  }

  group_levels <- unique(sim_struct$group)
  groups_vec <- match(sim_struct$group, group_levels)


  # Parallel C++ loop
  events_raw <- simulate_renewal_flexible(
    time_vec = time_vec,
    modulant_mat_flat = modulant_mat_flat,
    groups_vec = as.integer(groups_vec),
    shape_vec = as.numeric(shape_vec),
    k_vec = as.numeric(k_vec),
    n_ind = n_ind,
    n_samples = n_samples,
    max_x = max_t,
    use_samples = use_samples,
    start_time = - 0.1 * max_t)

  events <- events_raw |>
    dplyr::mutate(
      ind   = ind_ids[ind],
      group = group_levels[group]
    )
  ## Determine censoring
  if (!use_samples) {

    events <- events |>
      dplyr::group_by(ind) |>
      dplyr::group_modify(\(ev_i, y) {
        current_group_val <- ev_i$group[1]
        ev_i |> dplyr::mutate(
          dt = diff(c(0, event_times)),
          censored = c(rep(0, length(event_times) - 1), 1)
        )
      }) |>
      ungroup()
  }


  if (use_samples) {
    events <- events |> dplyr::relocate(sample, .before = ind)
  }

  trace_data$events <- events
  trace_data


  if (is.list(trace_data)) {
    sim_data <- trace_data %>% append(
      list(events = events)
    )
    class(sim_data) <- c("sim_traces", class(sim_data))
  } else {
    sim_data <- list(
      traces = trace_data,
      events = events
    )
  }
  class(sim_data) <- c("sim_events", class(sim_data))
  return(sim_data)
}

.interp_points <- function(df, resolution, id, time, y) {
  t_old <- unique(df[[time]])
  new_times <- seq(min(t_old), max(t_old), length.out = max(t_old) * 1 / resolution)

  interp_list <- by(df, df[[id]], function(sub_df) {
    interp <- approx(
      x = sub_df[[time]],
      y = sub_df[[y]],
      xout = new_times
    )

    res_df <- data.frame(
      id = sub_df[[id]][1],
      t = interp$x,
      scale = interp$y
    )


    names(res_df) <- c(id, time, y)
    return(res_df)
  })

  interp_df <- do.call(rbind, interp_list)

  rownames(interp_df) <- NULL

  return(interp_df)
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
    events <- simulate_renewal_multi_parallel(time, modulant, 1, shape, k)
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
