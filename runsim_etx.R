# ==============================================================================
# PROJECT: CAMO — Missing Data in Multivariate Hypothesis Testing
# SCRIPT:  simulation_main_extensions.R
# PURPOSE: Extension simulation study comparing GPC NP, GPC P, and NPMV under
#          Benchmark, Complete Case (CC), Available Case (AC), and MI
#          across a grid of delta x rho x correlation structure x dropout,
#          for t-distributed at m = 20 and m = 40.
# ==============================================================================

source("simulation_functions_ext.R")  

# ==============================================================================
# MODERN R PACKAGE ECOSYSTEM
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(furrr)        # Parallel map (drop-in replacement for purrr)
library(future)       # Backend for furrr
library(tibble)
library(stringr)
library(forcats)
library(glue)

library(ggplot2)
library(gridExtra)
library(GGally)
library(RColorBrewer)
library(scales)

library(mvtnorm)
library(MASS)
library(mice)
library(naniar)
library(visdat)

library(knitr)
library(kableExtra)
library(DiagrammeR)

# ------------------------------------------------------------------------------
# Parallelisation
# Use all cores minus one — keeps the OS responsive during the run.
# If you hit memory or socket issues, switch to plan(sequential).
# ------------------------------------------------------------------------------
n_workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = n_workers)
cat(sprintf("Running with %d parallel workers\n", n_workers))

start_time <- Sys.time()
cat(sprintf("Simulation started at: %s\n", format(start_time, "%Y-%m-%d %H:%M:%S")))

# ------------------------------------------------------------------------------
# Simulation settings
# ------------------------------------------------------------------------------
results_path <- "results/results_summary_ext_small.csv"
n_sims       <- 1000

# Single master seed — drawn once here, controls all per-replication seeds
set.seed(2026)
#rep_seeds <- sample.int(1e6, n_sims)

# ------------------------------------------------------------------------------
# Scenario grid
# Same delta values across distributions — allows direct comparison
# CS allows rho = 0; AR1 requires rho > 0
# ------------------------------------------------------------------------------
# Build one full main-paper grid for a given (dist, df) combination
build_grid <- function(dist_val, df_val) {
  bind_rows(
    expand.grid(
      n     = 15,
      str   = "cs",
      rho   = c(0, 0.25, 0.50, 0.75),
      psi0  = c(-2, -1, 0),
      psi1  = c(0, 0.5, 1.0),
      delta = c(0, 0.45, 0.65, 0.85)
    ),
     expand.grid(
       n     = 15,
      str   = "ar1",
      rho   = c(0.25, 0.50, 0.75),
      psi0  = c(-2, -1, 0),
      psi1  = c(0, 0.5, 1.0),
      delta = c(0, 0.45, 0.65, 0.85)
     ),
    # N = 50, moderate sample size grid
    # expand.grid(
    # n     = c(50),
    # str   = "cs",
    # rho   = c(0, 0.25, 0.50, 0.75),
    # psi0  = c(0),
    # psi1  = c(0, 0.5, 1.0),
    # delta = c(0, 0.20, 0.40, 0.60)
    # ),
    # expand.grid(
    #  n     = c(50),
    # str   = "ar1",
    #  rho   = c(0.25, 0.50, 0.75),
    #  psi0  = c(0),
    #  psi1  = c(0, 0.5, 1.0),
    #  delta = c(0, 0.20, 0.40, 0.60)
    # ),
    # N = 100, 150 — large sample grid
    # expand.grid(
    #  n     = c(100, 150),
    #  str   = "cs",
    #  rho   = c(0, 0.25, 0.50, 0.75),
    #  psi0  = c(0),
    #  psi1  = c(0, 0.5, 1.0),
    #  delta = c(0, 0.10, 0.20, 0.30)
    #),
    #expand.grid(
    #  n     = c(100, 150),
    #  str   = "ar1",
    #  rho   = c(0.25, 0.50, 0.75),
    #  psi0  = c(0),
    #  psi1  = c(0, 0.5, 1.0),
    #  delta = c(0, 0.10, 0.20, 0.30)
    # )
  ) %>%
    mutate(dist = dist_val, df = df_val) %>%
    as_tibble()
}

# Combine across all (dist, df) combinations and both m values
sim_grid <- bind_rows(
  build_grid("t",      3),
  build_grid("t",      5)
  #build_grid("binary", NA_real_)
) %>%
  crossing(m = c(20, 40))

cat(sprintf("Total scenarios : %d  (CS: %d  |  AR1: %d)\n",
            nrow(sim_grid),
            sum(sim_grid$str == "cs"),
            sum(sim_grid$str == "ar1")))
cat(sprintf("  t df=3 : %d  |  t df=5 : %d\n",
            sum(sim_grid$dist == "t" & sim_grid$df == 3),
            sum(sim_grid$dist == "t" & sim_grid$df == 5)))
            #sum(sim_grid$dist == "binary")

# ------------------------------------------------------------------------------
# Single-scenario runner
# Outer loop: sequential over scenarios
# Inner loop: parallel over replications via future_map_dfr
# ------------------------------------------------------------------------------
run_sim_scenario <- function(grid_row, idx, total) {
  
  idx <- as.integer(idx)
  scenario_start <- Sys.time()
  
  n     <- grid_row$n
  str   <- grid_row$str
  rho   <- grid_row$rho
  psi0  <- grid_row$psi0
  psi1  <- grid_row$psi1
  delta <- grid_row$delta
  dist  <- grid_row$dist
  df_t  <- grid_row$df
  m_val <- grid_row$m
  
  cat(sprintf(
    "Scenario %d/%d: n=%d, str=%s, rho=%.2f, psi0=%.1f, psi1=%.1f, delta=%.2f, dist=%s, df=%s, m=%d\n",
    idx, total, n, str, rho, psi0, psi1, delta, dist,
    ifelse(is.na(df_t), "-", as.character(df_t)), m_val
  ))
  flush.console()
  
  reps <- future_map_dfr(
    seq_len(n_sims),
    function(sim) {
      #set.seed(rep_seeds[sim])   # fixed seed per replication — reproducible
      # regardless of worker count
      #cat(sprintf("  Rep %d/%d | PID %d\n", sim, n_sims, Sys.getpid()))
      
      dat_raw <- if (dist == "t") {
        generate_sim_data(n, rho, str, delta, dist = dist, df = df_t)
      } else {
        generate_sim_data(n, rho, str, delta, dist = dist)
      }
      dat_miss <- apply_dropout(dat_raw, psi0, psi1)
      an       <- run_analysis(dat_miss, dat_complete = dat_raw,
                               m = m_val, dist = dist)
      as_tibble(an)
    },
    .options = furrr_options(seed = TRUE)
  )
  
  scenario_time <- as.numeric(difftime(Sys.time(), scenario_start, units = "secs"))
  elapsed_total <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  eta_seconds   <- (elapsed_total / idx) * (total - idx)
  cat(sprintf(
    "✔ Scenario %d done in %.1fs | ETA: %s (%.1f min remaining)\n\n",
    idx, scenario_time,
    format(Sys.time() + eta_seconds, "%H:%M:%S"),
    eta_seconds / 60
  ))
  flush.console()
  
  reps %>%
    summarise(
      across(ends_with("_est"),  list(mean        = ~mean(., na.rm = TRUE))),
      across(ends_with("_pval"), list(reject_rate = ~mean(. < 0.05, na.rm = TRUE)))
    ) %>%
    mutate(n = n, rho = rho, str = str, psi0 = psi0, psi1 = psi1,
           delta = delta, dist = dist, df = df_t, m = m_val)
}

# ------------------------------------------------------------------------------
# Run all scenarios
# ------------------------------------------------------------------------------
n_total <- nrow(sim_grid)

results_summary <- imap_dfr(
  split(sim_grid, seq_len(n_total)),
  function(grid_row, idx) run_sim_scenario(grid_row[1, ], as.integer(idx), n_total)
)

# ------------------------------------------------------------------------------
# Save and clean up
# ------------------------------------------------------------------------------
if (!dir.exists("results")) dir.create("results")
write.csv(results_summary, results_path, row.names = FALSE)
cat(sprintf("Results saved to: %s\n", results_path))
plan(sequential)   # release workers cleanly when done

results_summary

