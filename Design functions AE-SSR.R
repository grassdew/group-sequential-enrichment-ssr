# =============================================================================
# Adaptive enrichment + sample size re-estimation design functions
#
# Biomarker settings:
#   1) one binary biomarker, x in {0, 1}
#   2) one continuous biomarker, x in [0, 1]
#
# Designs:
#   GSD     : group sequential design, no enrichment, no SSR
#   GSE     : group sequential enrichment design, no SSR
#   GSD-SSR : group sequential design with sample size re-estimation
#   GSE-SSR : group sequential enrichment design with sample size re-estimation
#
# Main additions for the AE + SSR project:
#   - PPoS-based SSR at each non-final interim look
#   - candidate total sample size grid
#   - design-specific future accrual rule for PPoS
#   - additional output columns for SSR decisions and enrichment behavior
# =============================================================================

# ---- Basic utilities ---------------------------------------------------------

as_treatment_factor <- function(x) {
  if (is.factor(x)) {
    factor(x, levels = c("Control", "Treatment"))
  } else {
    factor(x, levels = c(0, 1), labels = c("Control", "Treatment"))
  }
}

assert_biomarker_type <- function(biomarker_type) {
  biomarker_type <- tolower(biomarker_type)
  if (!biomarker_type %in% c("binary", "continuous")) {
    stop("biomarker_type must be 'binary' or 'continuous'.")
  }
  biomarker_type
}

design_flags <- function(design_type) {
  design_type <- toupper(design_type)
  design_type <- gsub("_", "-", design_type)
  if (!design_type %in% c("GSD", "GSE", "GSD-SSR", "GSE-SSR")) {
    stop("design_type must be one of GSD, GSE, GSD-SSR, GSE-SSR.")
  }

  list(
    design_type = design_type,
    enrichment_enabled = design_type %in% c("GSE", "GSE-SSR"),
    ssr_enabled = design_type %in% c("GSD-SSR", "GSE-SSR")
  )
}

pad_numeric <- function(x, n, prefix) {
  out <- rep(NA_real_, n)
  if (length(x) > 0) out[seq_len(min(length(x), n))] <- x[seq_len(min(length(x), n))]
  names(out) <- paste0(prefix, seq_len(n))
  out
}

pad_character <- function(x, n, prefix) {
  out <- rep(NA_character_, n)
  if (length(x) > 0) out[seq_len(min(length(x), n))] <- x[seq_len(min(length(x), n))]
  names(out) <- paste0(prefix, seq_len(n))
  out
}

ceil_to_step <- function(x, step_size) {
  ceiling(x / step_size) * step_size
}

# ---- Biomarker summaries -----------------------------------------------------

count_biomarker_categories <- function(x_values, biomarker_type, n_bins = 20) {
  biomarker_type <- assert_biomarker_type(biomarker_type)

  if (length(x_values) == 0L) {
    if (biomarker_type == "binary") return(c(x0 = 0L, x1 = 0L))
    return(stats::setNames(integer(n_bins), paste0("bin", seq_len(n_bins))))
  }

  if (biomarker_type == "binary") {
    return(c(x0 = sum(x_values == 0, na.rm = TRUE),
             x1 = sum(x_values == 1, na.rm = TRUE)))
  }

  bin_breaks <- seq(0, 1, length.out = n_bins + 1L)
  bin_index <- as.numeric(cut(x_values, breaks = bin_breaks,
                              include.lowest = TRUE, labels = FALSE))
  stats::setNames(tabulate(bin_index, nbins = n_bins), paste0("bin", seq_len(n_bins)))
}

flatten_count_list <- function(count_list, prefix, max_looks,
                               biomarker_type, n_bins = 20) {
  biomarker_type <- assert_biomarker_type(biomarker_type)
  category_names <- if (biomarker_type == "binary") c("x0", "x1") else paste0("bin", seq_len(n_bins))

  out <- list()
  for (look in seq_len(max_looks)) {
    if (look <= length(count_list) && !is.null(count_list[[look]])) {
      counts <- count_list[[look]]
      for (nm in category_names) {
        out[[paste0(prefix, "_", nm, "_interim_", look)]] <- unname(counts[nm])
      }
    } else {
      for (nm in category_names) {
        out[[paste0(prefix, "_", nm, "_interim_", look)]] <- NA_integer_
      }
    }
  }
  out
}

# ---- Group sequential statistics --------------------------------------------

calculate_t_statistic <- function(interim_data) {
  treatment_y <- interim_data$y[interim_data$group == "Treatment"]
  control_y <- interim_data$y[interim_data$group == "Control"]

  n_treatment <- length(treatment_y)
  n_control <- length(control_y)
  if (n_treatment < 2 || n_control < 2) return(0)

  den <- sqrt(stats::var(treatment_y) / n_treatment +
                stats::var(control_y) / n_control)
  if (!is.finite(den) || den == 0) return(0)

  (mean(treatment_y) - mean(control_y)) / den
}

calculate_stratified_statistic <- function(test_statistics,
                                           incremental_sample_sizes,
                                           k = length(test_statistics),
                                           total_n = sum(incremental_sample_sizes[seq_len(k)])) {
  if (k <= 0 || length(test_statistics) == 0) return(0)
  idx <- seq_len(k)
  n_l <- incremental_sample_sizes[idx]
  n_l_sum <- total_n
  if (!is.finite(n_l_sum) || n_l_sum <= 0) return(0)
  sum(sqrt(n_l / n_l_sum) * test_statistics[idx])
}

assert_analysis_method <- function(analysis_method) {
  analysis_method <- tolower(analysis_method)
  if (!analysis_method %in% c("observed_stratified", "inverse_normal")) {
    stop("analysis_method must be 'observed_stratified' or 'inverse_normal'.")
  }
  analysis_method
}

calculate_inverse_normal_statistic <- function(test_statistics,
                                               information_rates,
                                               k = length(test_statistics)) {
  if (k <= 0 || length(test_statistics) == 0) return(0)

  idx <- seq_len(k)
  info_increments <- c(information_rates[1], diff(information_rates))
  weights <- sqrt(info_increments[idx])
  denom <- sqrt(sum(weights^2))
  if (!is.finite(denom) || denom <= 0) return(0)

  sum(weights * test_statistics[idx]) / denom
}

calculate_analysis_statistic <- function(test_statistics,
                                         incremental_sample_sizes,
                                         information_rates,
                                         analysis_method,
                                         k = length(test_statistics)) {
  analysis_method <- assert_analysis_method(analysis_method)

  if (analysis_method == "inverse_normal") {
    return(calculate_inverse_normal_statistic(
      test_statistics = test_statistics,
      information_rates = information_rates,
      k = k
    ))
  }

  calculate_stratified_statistic(
    test_statistics = test_statistics,
    incremental_sample_sizes = incremental_sample_sizes,
    k = k
  )
}

stopping_bounds <- function(alpha, beta, type_alpha, type_beta, informationRates) {
  if (!requireNamespace("rpact", quietly = TRUE)) {
    stop("Package 'rpact' is required for group sequential boundaries.")
  }

  gsd <- rpact::getDesignGroupSequential(
    typeOfDesign = type_alpha,
    alpha = alpha,
    typeBetaSpending = type_beta,
    bindingFutility = FALSE,
    informationRates = informationRates,
    sided = 1,
    beta = beta
  )

  list(
    efficacy_bound = gsd$criticalValues,
    futility_bound = c(gsd$futilityBounds,
                       gsd$criticalValues[length(gsd$criticalValues)])
  )
}

adaptive_threshold <- function(n_k, N, d, g) {
  d * (n_k / N)^g
}

# ---- Bayesian model for enrichment ------------------------------------------

update_model <- function(cumulative_data, biomarker_type,
                         prior_scale = 10,
                         chains = 4,
                         iter = 2000,
                         cores = 1) {
  biomarker_type <- assert_biomarker_type(biomarker_type)

  if (!requireNamespace("rstanarm", quietly = TRUE)) {
    stop("Package 'rstanarm' is required for enrichment modeling.")
  }

  cumulative_data$group <- as_treatment_factor(cumulative_data$group)
  cumulative_data$trt <- as.numeric(cumulative_data$group == "Treatment")

  if (biomarker_type == "binary") {
    return(rstanarm::stan_glm(
      y ~ group * x,
      data = cumulative_data,
      family = gaussian(),
      prior = rstanarm::normal(0, prior_scale),
      chains = chains,
      iter = iter,
      cores = cores,
      refresh = 0
    ))
  }

  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required for continuous-biomarker smooth terms.")
  }
  if (!"package:mgcv" %in% search()) {
    suppressPackageStartupMessages(library(mgcv))
  }

  rstanarm::stan_gamm4(
    formula = y ~ group + s(x) + s(x, by = trt),
    data = cumulative_data,
    family = gaussian(),
    prior = rstanarm::normal(0, prior_scale),
    chains = chains,
    iter = iter,
    cores = cores,
    refresh = 0
  )
}

make_prediction_data <- function(new_subject_data, treatment_group) {
  out <- new_subject_data
  out$group <- factor(treatment_group, levels = c("Control", "Treatment"))
  out$trt <- as.numeric(treatment_group == "Treatment")
  out
}

posterior_diff_draws <- function(model, new_subject_data, include_residual = TRUE) {
  if (is.null(model) || is.null(new_subject_data) || nrow(new_subject_data) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  nd_treat <- make_prediction_data(new_subject_data, "Treatment")
  nd_ctrl <- make_prediction_data(new_subject_data, "Control")

  if (include_residual) {
    y_treat <- rstanarm::posterior_predict(model, newdata = nd_treat)
    y_ctrl <- rstanarm::posterior_predict(model, newdata = nd_ctrl)
  } else {
    y_treat <- tryCatch(
      rstanarm::posterior_epred(model, newdata = nd_treat),
      error = function(e) rstanarm::posterior_predict(model, newdata = nd_treat)
    )
    y_ctrl <- tryCatch(
      rstanarm::posterior_epred(model, newdata = nd_ctrl),
      error = function(e) rstanarm::posterior_predict(model, newdata = nd_ctrl)
    )
  }

  y_treat - y_ctrl
}

calculate_pbi <- function(model, new_subject_data, epsilon) {
  diff_draws <- posterior_diff_draws(model, new_subject_data, include_residual = TRUE)
  if (length(diff_draws) == 0) return(numeric(0))
  colMeans(diff_draws > epsilon)
}

enroll_subjects <- function(model, subjects_to_screen, threshold,
                            enrolled_indices, enrolled_indices_current_stage,
                            max_enrollments, epsilon) {
  if (max_enrollments <= 0 || is.null(subjects_to_screen) || nrow(subjects_to_screen) == 0) {
    return(list(
      enrolled_subjects = NULL,
      enrolled_indices = enrolled_indices,
      enrolled_indices_current_stage = enrolled_indices_current_stage
    ))
  }

  pbis <- calculate_pbi(model, subjects_to_screen, epsilon)
  eligible_indices <- which(pbis > threshold)
  n_take <- min(length(eligible_indices), max_enrollments)

  if (n_take > 0) {
    take <- if (n_take < length(eligible_indices)) {
      sample(eligible_indices, n_take)
    } else {
      eligible_indices
    }
    enrolled_subjects <- subjects_to_screen[take, , drop = FALSE]
    enrolled_indices <- c(enrolled_indices, enrolled_subjects$subject_id)
    enrolled_indices_current_stage <- c(enrolled_indices_current_stage,
                                        enrolled_subjects$subject_id)
  } else {
    enrolled_subjects <- NULL
  }

  list(
    enrolled_subjects = enrolled_subjects,
    enrolled_indices = enrolled_indices,
    enrolled_indices_current_stage = enrolled_indices_current_stage
  )
}

# ---- PPoS helpers for SSR ----------------------------------------------------

estimate_effect_from_data <- function(data) {
  treatment_y <- data$y[data$group == "Treatment"]
  control_y <- data$y[data$group == "Control"]
  n_t <- length(treatment_y)
  n_c <- length(control_y)

  if (n_t < 2 || n_c < 2) {
    return(list(delta = 0, delta_sd = 1, sigma = stats::sd(data$y, na.rm = TRUE)))
  }

  delta <- mean(treatment_y) - mean(control_y)
  var_t <- stats::var(treatment_y)
  var_c <- stats::var(control_y)
  sigma <- sqrt(((n_t - 1) * var_t + (n_c - 1) * var_c) / (n_t + n_c - 2))
  se_delta <- sqrt(var_t / n_t + var_c / n_c)

  if (!is.finite(sigma) || sigma <= 0) sigma <- stats::sd(data$y, na.rm = TRUE)
  if (!is.finite(sigma) || sigma <= 0) sigma <- 1
  if (!is.finite(se_delta) || se_delta <= 0) se_delta <- sigma

  list(delta = delta, delta_sd = se_delta, sigma = sigma)
}

estimate_future_effect_for_pool <- function(model, future_pool, epsilon, threshold,
                                            enrichment_enabled,
                                            cumulative_data) {
  empirical <- estimate_effect_from_data(cumulative_data)

  if (!enrichment_enabled || is.null(model) || is.null(future_pool) || nrow(future_pool) == 0) {
    empirical$acceptance_rate <- ifelse(enrichment_enabled, NA_real_, 1)
    empirical$n_projection_pool <- ifelse(is.null(future_pool), 0L, nrow(future_pool))
    empirical$n_projection_accepted <- ifelse(is.null(future_pool), 0L, nrow(future_pool))
    return(empirical)
  }

  pbi <- calculate_pbi(model, future_pool, epsilon)
  accepted <- which(pbi > threshold)
  acceptance_rate <- mean(pbi > threshold)

  if (length(accepted) == 0 || !is.finite(acceptance_rate) || acceptance_rate <= 0) {
    empirical$acceptance_rate <- 0
    empirical$n_projection_pool <- nrow(future_pool)
    empirical$n_projection_accepted <- 0L
    return(empirical)
  }

  diff_draws <- posterior_diff_draws(model, future_pool[accepted, , drop = FALSE],
                                     include_residual = FALSE)
  draw_means <- rowMeans(diff_draws)
  delta <- mean(draw_means)
  delta_sd <- stats::sd(draw_means)

  if (!is.finite(delta_sd) || delta_sd <= 0) delta_sd <- empirical$delta_sd

  list(
    delta = delta,
    delta_sd = delta_sd,
    sigma = empirical$sigma,
    acceptance_rate = acceptance_rate,
    n_projection_pool = nrow(future_pool),
    n_projection_accepted = length(accepted)
  )
}

candidate_total_sample_sizes <- function(total_enrolled, look_index,
                                         current_target_n,
                                         info_rates,
                                         ssr_min_sample_size,
                                         ssr_max_sample_size,
                                         ssr_step_size,
                                         allow_sample_size_decrease = TRUE) {
  next_info_rate <- if (look_index < length(info_rates)) {
    info_rates[look_index + 1L]
  } else {
    1
  }

  # If N is reduced, the next scheduled milestone must still be ahead of the
  # enrollment already accrued at the current look.
  lower_by_next_info <- ceiling((total_enrolled + 1L) / next_info_rate)
  lower <- max(total_enrolled + 1L, ssr_min_sample_size, lower_by_next_info)

  if (!allow_sample_size_decrease) {
    lower <- max(lower, current_target_n)
  }

  upper <- max(ssr_max_sample_size, current_target_n, lower)

  grid <- seq(ceil_to_step(lower, ssr_step_size), upper, by = ssr_step_size)
  anchors <- c(current_target_n, ssr_max_sample_size)
  anchors <- anchors[anchors >= lower & anchors <= upper]

  sort(unique(as.integer(c(anchors, grid))))
}

compute_ppos_grid <- function(candidate_ns,
                              test_statistics,
                              incremental_sample_sizes,
                              info_rates,
                              current_look_index,
                              total_enrolled,
                              total_screened,
                              max_screen,
                              final_efficacy_bound,
                              effect_info,
                              enrichment_enabled,
                              analysis_method = "inverse_normal") {
  if (length(candidate_ns) == 0) {
    return(data.frame())
  }

  analysis_method <- assert_analysis_method(analysis_method)

  out <- data.frame(
    candidate_n = candidate_ns,
    remaining_enroll = NA_integer_,
    projected_total_screened = NA_real_,
    feasible = FALSE,
    ppos = NA_real_
  )

  acceptance_rate <- effect_info$acceptance_rate
  if (!is.finite(acceptance_rate) || is.na(acceptance_rate)) {
    acceptance_rate <- ifelse(enrichment_enabled, 0, 1)
  }

  sigma <- effect_info$sigma
  if (!is.finite(sigma) || sigma <= 0) sigma <- 1

  for (i in seq_len(nrow(out))) {
    n_total <- out$candidate_n[i]
    n_future <- max(0L, n_total - total_enrolled)

    if (n_future == 0) {
      projected_total_screened <- total_screened
      feasible <- TRUE
    } else if (acceptance_rate <= 0) {
      projected_total_screened <- Inf
      feasible <- FALSE
    } else {
      additional_screens <- if (enrichment_enabled) {
        ceiling(n_future / acceptance_rate)
      } else {
        n_future
      }
      projected_total_screened <- total_screened + additional_screens
      feasible <- projected_total_screened <= max_screen
    }

    out$remaining_enroll[i] <- n_future
    out$projected_total_screened[i] <- projected_total_screened
    out$feasible[i] <- feasible

    if (!feasible) next

    future_z_mean <- effect_info$delta * sqrt(n_future / (4 * sigma^2))
    future_z_var <- 1 + (n_future * effect_info$delta_sd^2 / (4 * sigma^2))
    if (!is.finite(future_z_var) || future_z_var <= 0) future_z_var <- 1

    if (analysis_method == "inverse_normal") {
      info_increments <- c(info_rates[1], diff(info_rates))
      k <- length(test_statistics)
      current_contribution <- if (k == 0) {
        0
      } else {
        sum(sqrt(info_increments[seq_len(k)]) * test_statistics)
      }

      future_weight <- sqrt(max(0, 1 - info_rates[current_look_index]))

      if (n_future == 0 || future_weight == 0) {
        out$ppos[i] <- as.numeric(current_contribution > final_efficacy_bound)
        next
      }

      final_mean <- current_contribution + future_weight * future_z_mean
      final_sd <- future_weight * sqrt(future_z_var)
    } else {
      current_contribution <- if (length(test_statistics) == 0) {
        0
      } else {
        sum(sqrt(incremental_sample_sizes / n_total) * test_statistics)
      }

      if (n_future == 0) {
        out$ppos[i] <- as.numeric(current_contribution > final_efficacy_bound)
        next
      }

      final_mean <- current_contribution + sqrt(n_future / n_total) * future_z_mean
      final_sd <- sqrt((n_future / n_total) * future_z_var)
    }

    if (!is.finite(final_sd) || final_sd <= 0) final_sd <- 1

    out$ppos[i] <- 1 - stats::pnorm((final_efficacy_bound - final_mean) / final_sd)
  }

  out
}

choose_ssr_target <- function(ppos_grid, ppos_target, ppos_futility) {
  if (nrow(ppos_grid) == 0 || all(!ppos_grid$feasible) ||
      all(is.na(ppos_grid$ppos[ppos_grid$feasible]))) {
    return(list(
      stop_for_futility = TRUE,
      selected_n = NA_integer_,
      selected_ppos = NA_real_,
      largest_n = NA_integer_,
      largest_ppos = NA_real_,
      decision = "SSR futility: no feasible candidate"
    ))
  }

  feasible_grid <- ppos_grid[ppos_grid$feasible & !is.na(ppos_grid$ppos), , drop = FALSE]
  target_hits <- feasible_grid[feasible_grid$ppos >= ppos_target, , drop = FALSE]

  largest_row <- feasible_grid[which.max(feasible_grid$candidate_n), , drop = FALSE]

  if (nrow(target_hits) > 0) {
    selected_row <- target_hits[which.min(target_hits$candidate_n), , drop = FALSE]
    return(list(
      stop_for_futility = FALSE,
      selected_n = selected_row$candidate_n,
      selected_ppos = selected_row$ppos,
      largest_n = largest_row$candidate_n,
      largest_ppos = largest_row$ppos,
      decision = "Target PPoS reached"
    ))
  }

  if (largest_row$ppos < ppos_futility) {
    return(list(
      stop_for_futility = TRUE,
      selected_n = largest_row$candidate_n,
      selected_ppos = largest_row$ppos,
      largest_n = largest_row$candidate_n,
      largest_ppos = largest_row$ppos,
      decision = "SSR futility: largest candidate below PPoS futility threshold"
    ))
  }

  selected_row <- feasible_grid[which.max(feasible_grid$ppos), , drop = FALSE]
  list(
    stop_for_futility = FALSE,
    selected_n = selected_row$candidate_n,
    selected_ppos = selected_row$ppos,
    largest_n = largest_row$candidate_n,
    largest_ppos = largest_row$ppos,
    decision = "Target PPoS not reached; selected maximum PPoS"
  )
}

# ---- Main design engine ------------------------------------------------------

run_one_biomarker_design <- function(trial_data,
                                     biomarker_type,
                                     design_type,
                                     max_screen,
                                     max_sample_size,
                                     info_rates,
                                     alpha,
                                     beta,
                                     type_alpha,
                                     type_beta,
                                     epsilon = 0.1,
                                     d = 0.5,
                                     g = 0.1,
                                     analysis_method = "inverse_normal",
                                     ssr_min_sample_size = ceiling(max_sample_size * 0.5),
                                     ssr_max_sample_size = ceiling(max_sample_size * 1.5),
                                     ssr_step_size = 50,
                                     allow_sample_size_decrease = TRUE,
                                     ssr_min_info_rate = 0,
                                     ppos_target = 0.8,
                                     ppos_futility = 0.1,
                                     projection_pool_size = 1000,
                                     model_prior_scale = 10,
                                     model_chains = 4,
                                     model_iter = 2000,
                                     model_cores = 1,
                                     n_bins = 20) {
  biomarker_type <- assert_biomarker_type(biomarker_type)
  flags <- design_flags(design_type)
  analysis_method <- assert_analysis_method(analysis_method)
  ssr_min_info_rate <- max(0, min(1, as.numeric(ssr_min_info_rate)))

  trial_data$group <- as_treatment_factor(trial_data$group)
  trial_data$trt <- as.numeric(trial_data$group == "Treatment")

  bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
  max_looks <- length(info_rates)
  original_target_n <- max_sample_size
  target_n <- max_sample_size

  total_screened <- 0L
  total_enrolled <- 0L
  start_idx <- 1L

  enrolled_indices <- integer(0)
  enrolled_indices_current_stage <- integer(0)
  screened_indices_current_stage <- integer(0)

  model <- NULL
  threshold <- 0
  test_statistics <- numeric(0)
  incremental_sample_sizes <- integer(0)
  interim_looks <- integer(0)

  early_stop <- FALSE
  conclusion <- NA_character_
  stop_reason <- NA_character_

  interim_enrollments <- numeric(0)
  interim_screenings <- numeric(0)
  cumulative_enrollments <- numeric(0)
  cumulative_screenings <- numeric(0)
  stratified_stats <- numeric(0)
  target_n_at_look <- numeric(0)
  threshold_at_look <- numeric(0)
  efficacy_bound_at_look <- numeric(0)
  futility_bound_at_look <- numeric(0)

  ppos_selected_n <- numeric(0)
  ppos_selected <- numeric(0)
  ppos_largest_n <- numeric(0)
  ppos_largest <- numeric(0)
  ppos_acceptance_rate <- numeric(0)
  ppos_projection_pool <- numeric(0)
  ppos_projection_accepted <- numeric(0)
  ssr_decision <- character(0)
  ssr_candidate_count <- numeric(0)
  ssr_feasible_candidate_count <- numeric(0)

  enrollment_count_list <- list()
  screening_count_list <- list()

  screen_batch_size <- max(1L, ceiling(original_target_n * info_rates[1]))

  record_counts <- function(look_index) {
    enrollment_count_list[[look_index]] <<- count_biomarker_categories(
      trial_data$x[enrolled_indices_current_stage],
      biomarker_type = biomarker_type,
      n_bins = n_bins
    )
    screening_count_list[[look_index]] <<- count_biomarker_categories(
      trial_data$x[screened_indices_current_stage],
      biomarker_type = biomarker_type,
      n_bins = n_bins
    )
  }

  while (total_enrolled < target_n &&
         total_screened < max_screen &&
         start_idx <= nrow(trial_data) &&
         length(interim_looks) < max_looks &&
         !early_stop) {

    look_index <- length(interim_looks) + 1L
    next_milestone <- min(ceiling(info_rates[look_index] * target_n), target_n)
    if (next_milestone <= total_enrolled) {
      next_milestone <- min(target_n, total_enrolled + 1L)
    }
    remaining_to_milestone <- next_milestone - total_enrolled
    if (remaining_to_milestone <= 0) break

    use_enrichment_rule <- flags$enrichment_enabled && !is.null(model)
    screen_count <- if (use_enrichment_rule) {
      min(max_screen - total_screened, screen_batch_size, nrow(trial_data) - start_idx + 1L)
    } else {
      min(max_screen - total_screened, remaining_to_milestone,
          nrow(trial_data) - start_idx + 1L)
    }
    if (screen_count <= 0) break

    end_idx <- start_idx + screen_count - 1L
    new_subjects <- trial_data[start_idx:end_idx, , drop = FALSE]

    total_screened <- total_screened + nrow(new_subjects)
    screened_indices_current_stage <- c(screened_indices_current_stage,
                                        new_subjects$subject_id)

    if (use_enrichment_rule) {
      enrollment_result <- enroll_subjects(
        model = model,
        subjects_to_screen = new_subjects,
        threshold = threshold,
        enrolled_indices = enrolled_indices,
        enrolled_indices_current_stage = enrolled_indices_current_stage,
        max_enrollments = remaining_to_milestone,
        epsilon = epsilon
      )
      enrolled_indices <- enrollment_result$enrolled_indices
      enrolled_indices_current_stage <- enrollment_result$enrolled_indices_current_stage
    } else {
      n_take <- min(nrow(new_subjects), remaining_to_milestone)
      if (n_take > 0) {
        subjects_to_enroll <- new_subjects[seq_len(n_take), , drop = FALSE]
        enrolled_indices <- c(enrolled_indices, subjects_to_enroll$subject_id)
        enrolled_indices_current_stage <- c(enrolled_indices_current_stage,
                                            subjects_to_enroll$subject_id)
      }
    }

    total_enrolled <- length(enrolled_indices)

    if (total_enrolled >= next_milestone) {
      look_index <- length(interim_looks) + 1L
      interim_looks <- c(interim_looks, total_enrolled)

      cumulative_data <- trial_data[enrolled_indices, , drop = FALSE]

      if (flags$enrichment_enabled) {
        model <- update_model(
          cumulative_data = cumulative_data,
          biomarker_type = biomarker_type,
          prior_scale = model_prior_scale,
          chains = model_chains,
          iter = model_iter,
          cores = model_cores
        )
      }

      interim_data <- trial_data[enrolled_indices_current_stage, , drop = FALSE]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))

      stratified_stat <- calculate_analysis_statistic(
        test_statistics = test_statistics,
        incremental_sample_sizes = incremental_sample_sizes,
        information_rates = info_rates,
        analysis_method = analysis_method,
        k = length(test_statistics)
      )

      interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
      interim_screenings <- c(interim_screenings, length(screened_indices_current_stage))
      cumulative_enrollments <- c(cumulative_enrollments, total_enrolled)
      cumulative_screenings <- c(cumulative_screenings, total_screened)
      stratified_stats <- c(stratified_stats, stratified_stat)
      target_n_at_look <- c(target_n_at_look, target_n)
      efficacy_bound_at_look <- c(efficacy_bound_at_look,
                                  bounds$efficacy_bound[look_index])
      futility_bound_at_look <- c(futility_bound_at_look,
                                  bounds$futility_bound[look_index])

      record_counts(look_index)

      if (total_enrolled < target_n) {
        if (stratified_stat > bounds$efficacy_bound[look_index]) {
          early_stop <- TRUE
          conclusion <- "Efficacy"
          stop_reason <- "Group sequential efficacy"
        } else if (stratified_stat < bounds$futility_bound[look_index]) {
          early_stop <- TRUE
          conclusion <- "Futility"
          stop_reason <- "Group sequential futility"
        }
      } else {
        conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[look_index],
                             "Efficacy", "Futility")
        stop_reason <- "Final analysis"
      }

      current_info_rate <- info_rates[look_index]
      ssr_timing_allowed <- current_info_rate >= ssr_min_info_rate

      if (!early_stop && flags$ssr_enabled && ssr_timing_allowed &&
          total_enrolled < target_n &&
          look_index < max_looks) {
        future_pool <- if (start_idx <= nrow(trial_data)) {
          trial_data[start_idx:nrow(trial_data), , drop = FALSE]
        } else {
          trial_data[0, , drop = FALSE]
        }
        if (nrow(future_pool) > projection_pool_size) {
          future_pool <- future_pool[seq_len(projection_pool_size), , drop = FALSE]
        }

        threshold_for_projection <- if (flags$enrichment_enabled) {
          adaptive_threshold(total_enrolled, target_n, d, g)
        } else {
          0
        }

        effect_info <- estimate_future_effect_for_pool(
          model = model,
          future_pool = future_pool,
          epsilon = epsilon,
          threshold = threshold_for_projection,
          enrichment_enabled = flags$enrichment_enabled,
          cumulative_data = cumulative_data
        )

        candidate_ns <- candidate_total_sample_sizes(
          total_enrolled = total_enrolled,
          look_index = look_index,
          current_target_n = target_n,
          info_rates = info_rates,
          ssr_min_sample_size = ssr_min_sample_size,
          ssr_max_sample_size = ssr_max_sample_size,
          ssr_step_size = ssr_step_size,
          allow_sample_size_decrease = allow_sample_size_decrease
        )

        ppos_grid <- compute_ppos_grid(
          candidate_ns = candidate_ns,
          test_statistics = test_statistics,
          incremental_sample_sizes = incremental_sample_sizes,
          info_rates = info_rates,
          current_look_index = look_index,
          total_enrolled = total_enrolled,
          total_screened = total_screened,
          max_screen = max_screen,
          final_efficacy_bound = bounds$efficacy_bound[max_looks],
          effect_info = effect_info,
          enrichment_enabled = flags$enrichment_enabled,
          analysis_method = analysis_method
        )

        ssr_choice <- choose_ssr_target(
          ppos_grid = ppos_grid,
          ppos_target = ppos_target,
          ppos_futility = ppos_futility
        )

        ppos_selected_n <- c(ppos_selected_n, ssr_choice$selected_n)
        ppos_selected <- c(ppos_selected, ssr_choice$selected_ppos)
        ppos_largest_n <- c(ppos_largest_n, ssr_choice$largest_n)
        ppos_largest <- c(ppos_largest, ssr_choice$largest_ppos)
        ppos_acceptance_rate <- c(ppos_acceptance_rate, effect_info$acceptance_rate)
        ppos_projection_pool <- c(ppos_projection_pool, effect_info$n_projection_pool)
        ppos_projection_accepted <- c(ppos_projection_accepted,
                                      effect_info$n_projection_accepted)
        ssr_decision <- c(ssr_decision, ssr_choice$decision)
        ssr_candidate_count <- c(ssr_candidate_count, nrow(ppos_grid))
        ssr_feasible_candidate_count <- c(ssr_feasible_candidate_count,
                                          sum(ppos_grid$feasible, na.rm = TRUE))

        if (ssr_choice$stop_for_futility) {
          early_stop <- TRUE
          conclusion <- "Futility"
          stop_reason <- ssr_choice$decision
        } else {
          target_n <- max(total_enrolled + 1L, as.integer(ssr_choice$selected_n))
        }
      } else {
        ppos_selected_n <- c(ppos_selected_n, NA_real_)
        ppos_selected <- c(ppos_selected, NA_real_)
        ppos_largest_n <- c(ppos_largest_n, NA_real_)
        ppos_largest <- c(ppos_largest, NA_real_)
        ppos_acceptance_rate <- c(ppos_acceptance_rate, NA_real_)
        ppos_projection_pool <- c(ppos_projection_pool, NA_real_)
        ppos_projection_accepted <- c(ppos_projection_accepted, NA_real_)
        ssr_decision <- c(ssr_decision, NA_character_)
        ssr_candidate_count <- c(ssr_candidate_count, NA_real_)
        ssr_feasible_candidate_count <- c(ssr_feasible_candidate_count, NA_real_)
      }

      if (flags$enrichment_enabled) {
        threshold <- adaptive_threshold(total_enrolled, target_n, d, g)
        threshold_at_look <- c(threshold_at_look, threshold)
      } else {
        threshold_at_look <- c(threshold_at_look, NA_real_)
      }

      screened_indices_current_stage <- integer(0)
      enrolled_indices_current_stage <- integer(0)
    }

    start_idx <- end_idx + 1L
  }

  if (is.na(conclusion) && total_enrolled > 0) {
    look_index <- min(length(test_statistics) + 1L, max_looks)

    if (length(enrolled_indices_current_stage) > 0 && length(test_statistics) < max_looks) {
      interim_looks <- c(interim_looks, total_enrolled)
      interim_data <- trial_data[enrolled_indices_current_stage, , drop = FALSE]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))

      stratified_stat <- calculate_analysis_statistic(
        test_statistics = test_statistics,
        incremental_sample_sizes = incremental_sample_sizes,
        information_rates = info_rates,
        analysis_method = analysis_method,
        k = length(test_statistics)
      )

      interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
      interim_screenings <- c(interim_screenings, length(screened_indices_current_stage))
      cumulative_enrollments <- c(cumulative_enrollments, total_enrolled)
      cumulative_screenings <- c(cumulative_screenings, total_screened)
      stratified_stats <- c(stratified_stats, stratified_stat)
      target_n_at_look <- c(target_n_at_look, target_n)
      threshold_at_look <- c(threshold_at_look, ifelse(flags$enrichment_enabled, threshold, NA_real_))
      efficacy_bound_at_look <- c(efficacy_bound_at_look, bounds$efficacy_bound[look_index])
      futility_bound_at_look <- c(futility_bound_at_look, bounds$futility_bound[look_index])
      ppos_selected_n <- c(ppos_selected_n, NA_real_)
      ppos_selected <- c(ppos_selected, NA_real_)
      ppos_largest_n <- c(ppos_largest_n, NA_real_)
      ppos_largest <- c(ppos_largest, NA_real_)
      ppos_acceptance_rate <- c(ppos_acceptance_rate, NA_real_)
      ppos_projection_pool <- c(ppos_projection_pool, NA_real_)
      ppos_projection_accepted <- c(ppos_projection_accepted, NA_real_)
      ssr_decision <- c(ssr_decision, NA_character_)
      ssr_candidate_count <- c(ssr_candidate_count, NA_real_)
      ssr_feasible_candidate_count <- c(ssr_feasible_candidate_count, NA_real_)
      record_counts(length(interim_looks))
    }

    k <- min(length(test_statistics), max_looks)
    final_stat <- calculate_analysis_statistic(
      test_statistics = test_statistics,
      incremental_sample_sizes = incremental_sample_sizes,
      information_rates = info_rates,
      analysis_method = analysis_method,
      k = k
    )
    conclusion <- ifelse(final_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
    stop_reason <- ifelse(total_enrolled >= target_n, "Final analysis",
                          "Analysis after screening limit")
  }

  if (is.na(conclusion)) {
    conclusion <- "Futility"
    stop_reason <- "No analyzable enrolled subjects"
  }

  status <- if (early_stop) {
    "Early stopping"
  } else if (total_enrolled >= target_n) {
    "Full enrollment"
  } else if (total_screened >= max_screen || start_idx > nrow(trial_data)) {
    "Screening limit reached"
  } else {
    "Stopped"
  }

  enrolled_x <- if (length(enrolled_indices) > 0) trial_data$x[enrolled_indices] else numeric(0)
  screened_x <- if (total_screened > 0) trial_data$x[seq_len(total_screened)] else numeric(0)

  overall_acceptance_rate <- ifelse(total_screened > 0, total_enrolled / total_screened, NA_real_)
  screen_exclusion_rate <- ifelse(is.na(overall_acceptance_rate), NA_real_,
                                  1 - overall_acceptance_rate)

  final_biomarker_summary <- if (biomarker_type == "binary") {
    list(
      enrolled_x0 = sum(enrolled_x == 0, na.rm = TRUE),
      enrolled_x1 = sum(enrolled_x == 1, na.rm = TRUE),
      screened_x0 = sum(screened_x == 0, na.rm = TRUE),
      screened_x1 = sum(screened_x == 1, na.rm = TRUE),
      enrolled_x1_rate = ifelse(length(enrolled_x) > 0, mean(enrolled_x == 1), NA_real_),
      screened_x1_rate = ifelse(length(screened_x) > 0, mean(screened_x == 1), NA_real_)
    )
  } else {
    list(
      enrolled_x_mean = ifelse(length(enrolled_x) > 0, mean(enrolled_x), NA_real_),
      enrolled_x_sd = ifelse(length(enrolled_x) > 1, stats::sd(enrolled_x), NA_real_),
      screened_x_mean = ifelse(length(screened_x) > 0, mean(screened_x), NA_real_),
      screened_x_sd = ifelse(length(screened_x) > 1, stats::sd(screened_x), NA_real_)
    )
  }

  ssr_evaluations <- sum(!is.na(ssr_decision))
  target_before_ssr <- rep(NA_real_, length(ppos_selected_n))
  n_target_compare <- min(length(target_before_ssr), length(target_n_at_look))
  if (n_target_compare > 0) {
    target_before_ssr[seq_len(n_target_compare)] <- target_n_at_look[seq_len(n_target_compare)]
  }
  target_updates <- sum(!is.na(ppos_selected_n) &
                          !is.na(target_before_ssr) &
                          ppos_selected_n != target_before_ssr,
                        na.rm = TRUE)
  target_decreases <- sum(!is.na(ppos_selected_n) &
                            !is.na(target_before_ssr) &
                            ppos_selected_n < target_before_ssr,
                          na.rm = TRUE)
  target_increases <- sum(!is.na(ppos_selected_n) &
                            !is.na(target_before_ssr) &
                            ppos_selected_n > target_before_ssr,
                          na.rm = TRUE)

  out <- c(
    list(
      total_screened = total_screened,
      total_enrolled = total_enrolled,
      status = status,
      conclusion = conclusion,
      early_stop = early_stop,
      stop_reason = stop_reason,
      enrichment_enabled = flags$enrichment_enabled,
      ssr_enabled = flags$ssr_enabled,
      analysis_method = analysis_method,
      sample_size_decrease_allowed = allow_sample_size_decrease,
      ssr_min_info_rate = ssr_min_info_rate,
      initial_max_sample_size = original_target_n,
      final_target_sample_size = target_n,
      target_sample_size_change = target_n - original_target_n,
      target_sample_size_ratio = target_n / original_target_n,
      sample_size_reduction = max(0, original_target_n - target_n),
      sample_size_reduction_rate = max(0, original_target_n - target_n) / original_target_n,
      ssr_evaluations = ssr_evaluations,
      ssr_target_updates = target_updates,
      ssr_target_decreases = target_decreases,
      ssr_target_increases = target_increases,
      ssr_any_decrease = target_decreases > 0,
      ssr_any_increase = target_increases > 0,
      overall_acceptance_rate = overall_acceptance_rate,
      screen_exclusion_rate = screen_exclusion_rate,
      n_looks_completed = length(test_statistics)
    ),
    final_biomarker_summary,
    as.list(pad_numeric(interim_enrollments, max_looks, "interim_enrollments_")),
    as.list(pad_numeric(interim_screenings, max_looks, "interim_screenings_")),
    as.list(pad_numeric(cumulative_enrollments, max_looks, "cumulative_enrolled_")),
    as.list(pad_numeric(cumulative_screenings, max_looks, "cumulative_screened_")),
    as.list(pad_numeric(stratified_stats, max_looks, "stratified_stat_")),
    as.list(pad_numeric(efficacy_bound_at_look, max_looks, "efficacy_bound_")),
    as.list(pad_numeric(futility_bound_at_look, max_looks, "futility_bound_")),
    as.list(pad_numeric(target_n_at_look, max_looks, "target_n_at_look_")),
    as.list(pad_numeric(threshold_at_look, max_looks, "enrichment_threshold_")),
    as.list(pad_numeric(ppos_selected_n, max_looks, "ppos_selected_n_")),
    as.list(pad_numeric(ppos_selected, max_looks, "ppos_selected_")),
    as.list(pad_numeric(ppos_largest_n, max_looks, "ppos_largest_n_")),
    as.list(pad_numeric(ppos_largest, max_looks, "ppos_largest_")),
    as.list(pad_numeric(ppos_acceptance_rate, max_looks, "ppos_acceptance_rate_")),
    as.list(pad_numeric(ppos_projection_pool, max_looks, "ppos_projection_pool_")),
    as.list(pad_numeric(ppos_projection_accepted, max_looks, "ppos_projection_accepted_")),
    as.list(pad_character(ssr_decision, max_looks, "ssr_decision_")),
    as.list(pad_numeric(ssr_candidate_count, max_looks, "ssr_candidate_count_")),
    as.list(pad_numeric(ssr_feasible_candidate_count, max_looks, "ssr_feasible_candidate_count_")),
    flatten_count_list(enrollment_count_list, "enrollments", max_looks,
                       biomarker_type = biomarker_type, n_bins = n_bins),
    flatten_count_list(screening_count_list, "screenings", max_looks,
                       biomarker_type = biomarker_type, n_bins = n_bins)
  )

  out
}

# ---- Design-specific wrappers ------------------------------------------------

gsd <- function(trial_data, biomarker_type, max_screen, max_sample_size,
                info_rates, alpha, beta, type_alpha, type_beta, ...) {
  run_one_biomarker_design(
    trial_data = trial_data,
    biomarker_type = biomarker_type,
    design_type = "GSD",
    max_screen = max_screen,
    max_sample_size = max_sample_size,
    info_rates = info_rates,
    alpha = alpha,
    beta = beta,
    type_alpha = type_alpha,
    type_beta = type_beta,
    ...
  )
}

gse <- function(trial_data, biomarker_type, max_screen, max_sample_size,
                info_rates, alpha, beta, type_alpha, type_beta,
                epsilon, d, g, ...) {
  run_one_biomarker_design(
    trial_data = trial_data,
    biomarker_type = biomarker_type,
    design_type = "GSE",
    max_screen = max_screen,
    max_sample_size = max_sample_size,
    info_rates = info_rates,
    alpha = alpha,
    beta = beta,
    type_alpha = type_alpha,
    type_beta = type_beta,
    epsilon = epsilon,
    d = d,
    g = g,
    ...
  )
}

gsd_ssr <- function(trial_data, biomarker_type, max_screen, max_sample_size,
                    info_rates, alpha, beta, type_alpha, type_beta,
                    epsilon, d, g, ...) {
  run_one_biomarker_design(
    trial_data = trial_data,
    biomarker_type = biomarker_type,
    design_type = "GSD-SSR",
    max_screen = max_screen,
    max_sample_size = max_sample_size,
    info_rates = info_rates,
    alpha = alpha,
    beta = beta,
    type_alpha = type_alpha,
    type_beta = type_beta,
    epsilon = epsilon,
    d = d,
    g = g,
    ...
  )
}

gse_ssr <- function(trial_data, biomarker_type, max_screen, max_sample_size,
                    info_rates, alpha, beta, type_alpha, type_beta,
                    epsilon, d, g, ...) {
  run_one_biomarker_design(
    trial_data = trial_data,
    biomarker_type = biomarker_type,
    design_type = "GSE-SSR",
    max_screen = max_screen,
    max_sample_size = max_sample_size,
    info_rates = info_rates,
    alpha = alpha,
    beta = beta,
    type_alpha = type_alpha,
    type_beta = type_beta,
    epsilon = epsilon,
    d = d,
    g = g,
    ...
  )
}

# ---- Data-generating scenarios ----------------------------------------------

true_effect_binary <- function(x, scenario_index, effect_size, multiplier = 1.5) {
  if (scenario_index == 1) return(rep(0, length(x)))
  if (scenario_index == 2) return(rep(effect_size, length(x)))
  if (scenario_index == 3) return(effect_size * as.numeric(x == 1))
  if (scenario_index == 4) return(effect_size * ifelse(x == 1, 1, 0.5))
  if (scenario_index == 5) return(effect_size * ifelse(x == 0, 1, 0))
  stop("scenario_index must be 1, 2, 3, 4, or 5.")
}

true_effect_continuous <- function(x, scenario_index, effect_size, multiplier = 1.2) {
  u_shape <- 4 * (x - 0.5)^2
  inv_u <- 4 * x * (1 - x)

  if (scenario_index == 1) return(rep(0, length(x)))
  if (scenario_index == 2) return(rep(effect_size, length(x)))
  if (scenario_index == 3) return(effect_size * x * multiplier)
  if (scenario_index == 4) return(effect_size * u_shape * multiplier)
  if (scenario_index == 5) return(effect_size * inv_u * multiplier)
  stop("scenario_index must be 1, 2, 3, 4, or 5.")
}

simulate_trial_data <- function(sim_index,
                                biomarker_type,
                                max_screen,
                                prev_rate = 0.5,
                                effect_size,
                                scenario_index,
                                binary_multiplier = 1.5,
                                continuous_multiplier = 1.2,
                                sigma_binary = 1,
                                sigma_continuous = 0.1) {
  biomarker_type <- assert_biomarker_type(biomarker_type)
  set.seed(sim_index)

  group_num <- stats::rbinom(max_screen, size = 1, prob = 0.5)

  if (biomarker_type == "binary") {
    x <- stats::rbinom(max_screen, size = 1, prob = prev_rate)
    true_effect <- true_effect_binary(
      x = x,
      scenario_index = scenario_index,
      effect_size = effect_size,
      multiplier = binary_multiplier
    )
    noise <- stats::rnorm(max_screen, mean = 0, sd = sigma_binary)
    y_control <- 5 + 0.1 * x + noise
  } else {
    x <- stats::runif(max_screen, min = 0, max = 1)
    true_effect <- true_effect_continuous(
      x = x,
      scenario_index = scenario_index,
      effect_size = effect_size,
      multiplier = continuous_multiplier
    )
    noise <- stats::rnorm(max_screen, mean = 0, sd = sigma_continuous)

    if (scenario_index %in% c(1, 2, 3)) {
      y_control <- 5 + 0.5 * x + noise
    } else if (scenario_index == 4) {
      y_control <- 3 + 2 * (4 * (x - 0.5)^2) + noise
    } else {
      y_control <- 3 + 2 * (4 * x * (1 - x)) + noise
    }
  }

  y_treat <- y_control + true_effect
  y <- ifelse(group_num == 1, y_treat, y_control)

  data.frame(
    subject_id = seq_len(max_screen),
    x = x,
    group = factor(group_num, levels = c(0, 1),
                   labels = c("Control", "Treatment")),
    y = as.numeric(y),
    true_effect = as.numeric(true_effect)
  )
}

trial_simulation <- function(sim_index,
                             biomarker_type,
                             design_type,
                             prev_rate = 0.5,
                             max_screen,
                             max_sample_size,
                             info_rates,
                             effect_size,
                             alpha,
                             beta,
                             type_alpha,
                             type_beta,
                             epsilon,
                             d,
                             g,
                             scenario_index,
                             analysis_method = "inverse_normal",
                             ssr_min_sample_size = ceiling(max_sample_size * 0.5),
                             ssr_max_sample_size = ceiling(max_sample_size * 1.5),
                             ssr_step_size = 50,
                             allow_sample_size_decrease = TRUE,
                             ssr_min_info_rate = 0,
                             ppos_target = 0.8,
                             ppos_futility = 0.1,
                             projection_pool_size = 1000,
                             model_prior_scale = 10,
                             model_chains = 4,
                             model_iter = 2000,
                             model_cores = 1,
                             n_bins = 20,
                             binary_multiplier = 1.5,
                             continuous_multiplier = 1.2,
                             sigma_binary = 1,
                             sigma_continuous = 0.1) {
  biomarker_type <- assert_biomarker_type(biomarker_type)
  flags <- design_flags(design_type)

  trial_data <- simulate_trial_data(
    sim_index = sim_index,
    biomarker_type = biomarker_type,
    max_screen = max_screen,
    prev_rate = prev_rate,
    effect_size = effect_size,
    scenario_index = scenario_index,
    binary_multiplier = binary_multiplier,
    continuous_multiplier = continuous_multiplier,
    sigma_binary = sigma_binary,
    sigma_continuous = sigma_continuous
  )

  result <- run_one_biomarker_design(
    trial_data = trial_data,
    biomarker_type = biomarker_type,
    design_type = flags$design_type,
    max_screen = max_screen,
    max_sample_size = max_sample_size,
    info_rates = info_rates,
    alpha = alpha,
    beta = beta,
    type_alpha = type_alpha,
    type_beta = type_beta,
    epsilon = epsilon,
    d = d,
    g = g,
    analysis_method = analysis_method,
    ssr_min_sample_size = ssr_min_sample_size,
    ssr_max_sample_size = ssr_max_sample_size,
    ssr_step_size = ssr_step_size,
    allow_sample_size_decrease = allow_sample_size_decrease,
    ssr_min_info_rate = ssr_min_info_rate,
    ppos_target = ppos_target,
    ppos_futility = ppos_futility,
    projection_pool_size = projection_pool_size,
    model_prior_scale = model_prior_scale,
    model_chains = model_chains,
    model_iter = model_iter,
    model_cores = model_cores,
    n_bins = n_bins
  )

  c(
    list(
      sim_index = sim_index,
      biomarker_type = biomarker_type,
      design_type = flags$design_type,
      scenario = scenario_index,
      effect_size = effect_size,
      prev_rate = ifelse(biomarker_type == "binary", prev_rate, NA_real_),
      true_mean_effect = mean(trial_data$true_effect),
      true_sensitive_rate = mean(trial_data$true_effect > epsilon)
    ),
    result
  )
}
