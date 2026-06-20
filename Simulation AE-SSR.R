# =============================================================================
# Simulation driver for the AE + SSR project
#
# Compares:
#   1) GSD
#   2) GSE
#   3) GSD-SSR
#   4) GSE-SSR
#
# Biomarker settings:
#   1) one binary biomarker, prevalence in c(0.3, 0.5, 0.7)
#   2) one continuous biomarker, x ~ Uniform(0, 1)
#
# Output:
#   one CSV per biomarker/design/scenario/effect-size/prevalence setting.
# =============================================================================

library(foreach)
library(doParallel)
library(doRNG)
library(data.table)

# ---- paths ------------------------------------------------------------------

design_fn_path <- "/Users/emily/Documents/Adaptive Enrichment Trial Design/AE and SSR/Code/Design functions AE-SSR.R"
root_out_dir <- "/Users/emily/Documents/Adaptive Enrichment Trial Design/AE and SSR/Simulation results"

env_logical <- function(name, default) {
  raw <- Sys.getenv(name, unset = NA_character_)
  if (is.na(raw) || !nzchar(trimws(raw))) {
    return(default)
  }

  value <- tolower(trimws(raw))
  if (value %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }
  if (value %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }

  stop("Environment variable ", name, " must be true/false, 1/0, or yes/no.")
}

env_csv <- function(name, default) {
  raw <- Sys.getenv(name, unset = NA_character_)
  if (is.na(raw) || !nzchar(trimws(raw))) {
    return(default)
  }

  values <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  values[nzchar(values)]
}

env_numeric <- function(name, default) {
  raw <- Sys.getenv(name, unset = NA_character_)
  if (is.na(raw) || !nzchar(trimws(raw))) {
    return(default)
  }

  value <- suppressWarnings(as.numeric(trimws(raw)))
  if (is.na(value)) {
    stop("Environment variable ", name, " must be numeric.")
  }
  value
}

# ---- global simulation settings ---------------------------------------------

num_sim <- 1000
chunk_size <- 5

biomarker_types <- c("binary", "continuous")
all_design_types <- c("GSD", "GSE", "GSD-SSR", "GSE-SSR")
design_types <- env_csv("AE_SSR_DESIGN_TYPES", all_design_types)
unknown_design_types <- setdiff(design_types, all_design_types)
if (length(unknown_design_types) > 0) {
  stop("Unknown AE_SSR_DESIGN_TYPES value(s): ", paste(unknown_design_types, collapse = ", "))
}

prev_grid <- c(0.3, 0.5, 0.7)
effect_size_grid <- c(0.1, 0.2, 0.3)
scenario_indices <- 1:5

max_screen <- 5000
max_sample_size <- 600

info_rates_list <- list(
  c(0.5, 1),
  c(1/3, 2/3, 1),
  c(0.25, 0.5, 0.75, 1)
)

alpha <- 0.025
beta <- 0.2
type_alpha <- "asOF"
type_beta <- "bsOF"

epsilon <- 0.1
d <- 0.5
g <- 0.1

# Analysis method:
#   observed_stratified = original sample-size-weighted statistic
#   inverse_normal     = pre-specified inverse-normal combination statistic
analysis_method <- "inverse_normal"

# SSR tuning.
# This lower bound allows SSR to reduce N when the interim evidence suggests
# the future enriched population has a larger treatment effect than expected.
# The design function also enforces that the selected N is above the number
# already enrolled and that the next planned information-rate milestone is ahead.
ssr_min_sample_size <- ceiling(max_sample_size * 0.5)
ssr_max_sample_size <- 900
ssr_step_size <- 50
allow_sample_size_decrease <- env_logical("AE_SSR_ALLOW_SAMPLE_SIZE_DECREASE", TRUE)
ssr_min_info_rate <- env_numeric("AE_SSR_MIN_INFO_RATE", 0)
ppos_target <- 0.8
ppos_futility <- 0.1
projection_pool_size <- 1000

# rstanarm model settings used for GSE and GSE-SSR.
# Increase chains/iter for final production runs if desired.
model_prior_scale <- 10
model_chains <- 4
model_iter <- 2000
model_cores <- 1

# Outcome scenario multipliers match Scenario plots.R.
binary_multiplier <- 1.5
continuous_multiplier <- 1.2
sigma_binary <- 1
sigma_continuous <- 0.1

base_seed <- 303
skip_existing <- env_logical("AE_SSR_SKIP_EXISTING", FALSE)

default_output_folder <- if (allow_sample_size_decrease) {
  "SSR sample size may decrease"
} else {
  "SSR no sample size decrease"
}
output_folder <- trimws(Sys.getenv("AE_SSR_OUTPUT_FOLDER", unset = default_output_folder))
run_root_out_dir <- if (nzchar(output_folder)) {
  file.path(root_out_dir, output_folder)
} else {
  root_out_dir
}

message("allow_sample_size_decrease: ", allow_sample_size_decrease)
message("ssr_min_info_rate: ", ssr_min_info_rate)
message("design_types: ", paste(design_types, collapse = ", "))
message("skip_existing: ", skip_existing)
message("Simulation output root: ", run_root_out_dir)

# ---- output helpers ----------------------------------------------------------

format_num <- function(x) {
  format(x, trim = TRUE, scientific = FALSE)
}

ssr_timing_label <- function() {
  if (ssr_min_info_rate <= 0) {
    return("")
  }
  paste0(", SSR min info=", format_num(ssr_min_info_rate))
}

make_out_dir <- function(biomarker_type, design_type, info_rates) {
  out_dir <- file.path(
    run_root_out_dir,
    sprintf(
      "Prior N(0, %s), d=%s, g=%s, method=%s, SSR min=%s, SSR max=%s, SSR target=%s, SSR futility=%s%s",
      format_num(model_prior_scale),
      format_num(d),
      format_num(g),
      analysis_method,
      format_num(ssr_min_sample_size),
      format_num(ssr_max_sample_size),
      format_num(ppos_target),
      format_num(ppos_futility),
      ssr_timing_label()
    ),
    sprintf("%d interim analysis", length(info_rates) - 1L),
    biomarker_type,
    design_type
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir
}

make_file_name <- function(out_dir, biomarker_type, prev_rate,
                           effect_size, scenario_index) {
  if (biomarker_type == "binary") {
    file.path(
      out_dir,
      sprintf(
        "prev_%s_efsize_%s_scenario_%s.csv",
        format_num(prev_rate),
        format_num(effect_size),
        scenario_index
      )
    )
  } else {
    file.path(
      out_dir,
      sprintf(
        "efsize_%s_scenario_%s.csv",
        format_num(effect_size),
        scenario_index
      )
    )
  }
}

# ---- parallel setup ----------------------------------------------------------

detected_cores <- parallel::detectCores()
if (is.na(detected_cores)) {
  detected_cores <- 2L
}
num_cores <- max(1L, detected_cores - 1L)
cl <- makeCluster(num_cores)
registerDoParallel(cl)
on.exit(stopCluster(cl), add = TRUE)

clusterExport(cl, varlist = c("design_fn_path"), envir = environment())

clusterEvalQ(cl, {
  library(rstanarm)
  library(rpact)
  library(data.table)
  Sys.setenv(OMP_NUM_THREADS = "1")
  options(mc.cores = 1)
  invisible(source(design_fn_path))
  TRUE
})

worker_packages <- c("rstanarm", "rpact", "data.table")

# ---- run grid ----------------------------------------------------------------

for (info_rates in info_rates_list) {
  n_interim <- length(info_rates) - 1L

  for (biomarker_type in biomarker_types) {
    prev_values <- if (biomarker_type == "binary") prev_grid else NA_real_

    for (design_type in design_types) {
      out_dir <- make_out_dir(biomarker_type, design_type, info_rates)

      for (prev_rate in prev_values) {
        for (scenario_index in scenario_indices) {
          effect_values <- if (scenario_index == 1) 0 else effect_size_grid

          for (effect_size in effect_values) {
            file_name <- make_file_name(
              out_dir = out_dir,
              biomarker_type = biomarker_type,
              prev_rate = prev_rate,
              effect_size = effect_size,
              scenario_index = scenario_index
            )

            if (skip_existing && file.exists(file_name)) {
              message("Skipping existing: ", file_name)
              next
            }

            message(
              sprintf(
                "Running %s | %s | %d interim | prevalence=%s | scenario=%s | effect=%s",
                biomarker_type,
                design_type,
                n_interim,
                ifelse(is.na(prev_rate), "NA", format_num(prev_rate)),
                scenario_index,
                format_num(effect_size)
              )
            )

            n_chunks <- ceiling(num_sim / chunk_size)

            grid_seed <- base_seed +
              1000000L * as.integer(n_interim) +
              100000L * match(biomarker_type, biomarker_types) +
              10000L * match(design_type, all_design_types) +
              100L * as.integer(scenario_index) +
              as.integer(round(effect_size * 100)) +
              ifelse(is.na(prev_rate), 0L, as.integer(prev_rate * 10L))

            res_list <- foreach(
              chunk_id = seq_len(n_chunks),
              .packages = worker_packages,
              .multicombine = TRUE,
              .maxcombine = 50,
              .options.RNG = grid_seed
            ) %dorng% {
              if (!exists("trial_simulation")) {
                source(design_fn_path)
              }

              n_rows <- min(chunk_size, num_sim - (chunk_id - 1L) * chunk_size)
              rows <- vector("list", length = n_rows)
              i0 <- (chunk_id - 1L) * chunk_size

              for (j in seq_along(rows)) {
                sim_index <- i0 + j

                trial_result <- trial_simulation(
                  sim_index = sim_index,
                  biomarker_type = biomarker_type,
                  design_type = design_type,
                  prev_rate = ifelse(is.na(prev_rate), 0.5, prev_rate),
                  max_screen = max_screen,
                  max_sample_size = max_sample_size,
                  info_rates = info_rates,
                  effect_size = effect_size,
                  alpha = alpha,
                  beta = beta,
                  type_alpha = type_alpha,
                  type_beta = type_beta,
                  epsilon = epsilon,
                  d = d,
                  g = g,
                  scenario_index = scenario_index,
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
                  binary_multiplier = binary_multiplier,
                  continuous_multiplier = continuous_multiplier,
                  sigma_binary = sigma_binary,
                  sigma_continuous = sigma_continuous
                )

                dt <- as.data.table(as.list(trial_result))
                dt[, `:=`(
                  n_interim = as.integer(n_interim),
                  info_rates = paste(info_rates, collapse = ","),
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
                  model_iter = model_iter
                )]
                rows[[j]] <- dt
              }

              rbindlist(rows, use.names = TRUE, fill = TRUE)
            }

            sim_results <- rbindlist(res_list, use.names = TRUE, fill = TRUE)

            fwrite(sim_results, file_name)
            message("Wrote: ", file_name)
          }
        }
      }
    }
  }
}
