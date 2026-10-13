# SimulationStudy-CAMO-PROJECT

Code accompanying:

> **Handling missing multivariate outcomes in nonparametric treatment
> comparisons.**
> Kanini L., Geroldinger M., Lauth W., Thiébaut R., Zimmermann G. (2026).
> *Manuscript under review.*

The study compares non-prioritized GPC (GPC-NP), prioritized GPC
(GPC-P) and the nonparametric multivariate ANOVA-type statistic (NPMV)
under four missing-data strategies: complete-data benchmark, complete
case (CC), available case (AC) and multiple imputation (MI). MI
inference uses Rubin's rules for GPC and D₂ pooling (Li *et al.* 1991 /
Grund *et al.* 2016) for the NPMV F-statistics.

---

## Files

| File | What it does |
|------|--------------|
| `simulation_functions.R`     | Core engine for the main (Gaussian) simulation: data generation, GPC-NP, GPC-P, NPMV, and the MI pipeline with Rubin's rules and D2 pooling. Sourced by the runners below. |
| `simulation_functions_ext.R` | Extended engine supporting multivariate t and binary outcomes. Sourced by the heavy-tailed runner. |
| `runsim.R`                   | Main simulation runner (Gaussian outcomes, `m = 20` imputations). Writes `results/results_summary.csv`. |
| `runsim_40imp.R`             | Same simulation at `m = 40` imputations. Writes `results/simulation_results_40imp.csv`. |
| `runsim_etx.R`               | Heavy-tailed robustness extension, small-sample grid (`n = 15`, full `psi0` range, both `m = 20` and `m = 40`). Writes `results/results_summary_ext_small.csv`. |
| `runsim_etx_large_m20.R`     | Heavy-tailed extension, large-sample grid (`n ∈ {50, 100, 150}`, full `psi0` range, `m = 20`). Writes `results/results_summary_ext_large_m20_FINAL.csv`. |
| `runsim_etx_large_m40.R`     | Same large-sample grid at `m = 40`. Writes `results/results_summary_ext_large_m40_FINAL.csv`. |
| `runsim_etx_severe_largen.R` | Heavy-tailed extension at large `n` with severe dropout only (`psi0 = 0`), both `m` values. Writes `results/results_summary_ext_severe_largen.csv`. |

Each runner writes an output CSV that includes an **`m` column**, so
the imputation count can always be read directly from the data without
relying on the filename. The four heavy-tail CSVs are produced as
separate runs and are combined downstream.
| `Peblescode.R`               | Real-data illustration on the PEBLES epidermolysis bullosa registry. |

Each runner script writes a results CSV to a `results/` subfolder. The
analysis / figure code that reads those CSVs is kept outside this
repository.

---

## Running the code

Set the working directory to the folder that holds these files, and
make sure `results/` exists (the runners create it if not). Then:

```r
# main Gaussian simulation
source("runsim.R")                     # m = 20
source("runsim_40imp.R")               # m = 40

# heavy-tailed extension (four runs)
source("runsim_etx.R")                 # n = 15, both m
source("runsim_etx_large_m20.R")       # n in {50, 100, 150}, m = 20
source("runsim_etx_large_m40.R")       # n in {50, 100, 150}, m = 40
source("runsim_etx_severe_largen.R")   # n in {50, 100, 150}, psi0 = 0 only

# real-data illustration
source("Peblescode.R")       # requires the PEBLES CSV (see below)
```

The runners parallelise over scenarios with `furrr`
(`plan(multisession, workers = detectCores() - 1)`); on a modern 8-core
laptop, the full main grid at `n_sims = 1000` takes several hours.
Switch to `plan(sequential)` to debug.

---

## Software requirements

- **R** ≥ 4.3
- CRAN packages:
  `BuyseTest`, `npmv`, `mice`, `mvtnorm`, `MASS`,
  `dplyr`, `tidyr`, `purrr`, `furrr`, `future`, `tibble`,
  `stringr`, `forcats`, `glue`, `lubridate`,
  `ggplot2`, `patchwork`, `naniar`.

Each script calls `sessionInfo()` at the end so the exact versions used
can be checked.

---

## Data availability

The PEBLES epidermolysis bullosa registry extract used in
`Peblescode.R` is **not publicly shareable**. Access requests go
through the PEBLES registry; please contact the corresponding author
of the manuscript.

To reproduce the illustration locally, place `PEBLES_Extract_Loise.csv`
in the path expected by the script and run from the appropriate working
directory.

---

## Funding

This research was funded in part by the Research Foundation – Flanders
(**FWO**) grant **G0A9225N** and the Austrian Science Fund (**FWF**)
[grant number to be added]. The project is part of the FWO–FWF Weave
programme *Covariate Adjustment for Multivariate Outcomes* (CAMO).

---

## Contact

Loise Kanini — Paris Lodron University of Salzburg — loycekarun@gmail.com
