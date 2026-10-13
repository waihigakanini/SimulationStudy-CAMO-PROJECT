# ==============================================================================
# PROJECT: CAMO — Missing Data in Multivariate Hypothesis Testing
# SCRIPT:  simulation_main.R
# PURPOSE: Full simulation study comparing GPC NP, GPC P, and NPMV under
#          Benchmark, Complete Case (CC), Available Case (AC), and MI
#          across a grid of delta x rho x correlation structure x dropout
# ==============================================================================

source("simulation_functions_4_2.R")

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
# Plot theme
# ------------------------------------------------------------------------------
theme_report <- function() {
  theme_minimal(base_family = "sans") +
    theme(
      plot.title       = element_text(face = "bold", size = 16, color = "#2d3748"),
      plot.subtitle    = element_text(size = 12, color = "#718096", margin = margin(b = 20)),
      axis.title       = element_text(face = "bold", size = 11, color = "#4a5568"),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#edf2f7"),
      strip.background = element_rect(fill = "#f7fafc", color = "#edf2f7"),
      strip.text       = element_text(face = "bold", color = "#4a5568")
    )
}

# ------------------------------------------------------------------------------
# Colour and shape palettes
# ------------------------------------------------------------------------------
method_colors <- c(
  "GPC_NP_Benchmark" = "#1A202C", "GPC_NP_CC" = "#718096",
  "GPC_NP_AC"        = "#E53E3E", "GPC_NP_MI" = "#38A169",
  "GPC_P_Benchmark"  = "#2D3748", "GPC_P_CC"  = "#A0AEC0",
  "GPC_P_AC"         = "#FC8181", "GPC_P_MI"  = "#68D391",
  "NPMV_Benchmark"   = "#3182CE", "NPMV_CC"   = "#63B3ED",
  "NPMV_AC"          = "#90CDF4", "NPMV_MI"   = "#BEE3F8"
)

method_shapes <- c(
  "Benchmark" = 17, "MI" = 15, "AC" = 1, "CC" = 0
)

# ------------------------------------------------------------------------------
# Simulation settings
# ------------------------------------------------------------------------------
results_path <- "results/simulation_results_40imp.csv"
n_sims       <- 1000   # increase to 500-1000 for the main simulation
M            <- 40     # number of imputations for this run

# Single master seed — drawn once here, controls all per-replication seeds
set.seed(2026)
#rep_seeds <- sample.int(1e6, n_sims)

# ------------------------------------------------------------------------------
# Scenario grid
# CS allows rho = 0; AR1 requires rho > 0
# ------------------------------------------------------------------------------
sim_grid <- bind_rows(
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
  expand.grid(
    n     = c(50),
    str   = "cs",
    rho   = c(0, 0.25, 0.50, 0.75),
    psi0  = c(-2, -1, 0),
    psi1  = c(0, 0.5, 1.0),
    delta = c(0, 0.20, 0.40, 0.60)
  ),
  expand.grid(
    n     = c(50),
    str   = "ar1",
    rho   = c(0.25, 0.50, 0.75),
    psi0  = c(-2, -1, 0),
    psi1  = c(0, 0.5, 1.0),
    delta = c(0, 0.20, 0.40, 0.60)
  ),
  # N = 100, 150 — large sample grid
  expand.grid(
    n     = c(100, 150),
    str   = "cs",
    rho   = c(0, 0.25, 0.50, 0.75),
    psi0  = c(-2, -1, 0),
    psi1  = c(0, 0.5, 1.0),
    delta = c(0, 0.10, 0.20, 0.30)
  ),
  expand.grid(
    n     = c(100, 150),
    str   = "ar1",
    rho   = c(0.25, 0.50, 0.75),
    psi0  = c(-2, -1, 0),
    psi1  = c(0, 0.5, 1.0),
    delta = c(0, 0.10, 0.20, 0.30)
  )
  
) %>% as_tibble()

cat(sprintf("Total scenarios : %d  (CS: %d  |  AR1: %d)\n",
            nrow(sim_grid),
            sum(sim_grid$str == "cs"),
            sum(sim_grid$str == "ar1")))

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
  
  cat(sprintf(
    "Scenario %d/%d: n=%d, str=%s, rho=%.2f, psi0=%.1f, psi1=%.1f, delta=%.1f\n",
    idx, total, n, str, rho, psi0, psi1, delta
  ))
  flush.console()
  
  reps <- future_map_dfr(
    seq_len(n_sims),
    function(sim) {
      #set.seed(rep_seeds[sim])   # fixed seed per replication — reproducible
      # regardless of worker count
      cat(sprintf("  Rep %d/%d | PID %d\n", sim, n_sims, Sys.getpid()))
      
      dat_raw  <- generate_sim_data(n, rho, str, delta)
      dat_miss <- apply_dropout(dat_raw, psi0, psi1)
      an       <- run_analysis(dat_miss, dat_complete = dat_raw, m = M)
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
           delta = delta, m = M)
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
