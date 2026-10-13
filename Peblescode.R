# ==============================================================================
# PEBLES REAL-DATA ILLUSTRATION
#
# Procedures:
#   1. Non-prioritized GPC (GPC-NP)
#   2. Prioritized GPC     (GPC-P)
#   3. NPMV
#
# Missing-data strategies:
#   - Complete Case (CC)
#   - Available Case (AC)
#   - Multiple Imputation (MI)
#
# Preliminary outcome vector:
#   Y1 = BEBS_Total
#   Y2 = LIS_sev
#   Y3 = QOLEB_total
#
# Group comparison:
#   Subtype 2 (Intermediate) versus Subtype 1 (Severe)
#
# Analysis occasion:
#   Earliest available review per patient
# ==============================================================================

source("simulation_functions.R")
# ==============================================================================
# 0. PACKAGES
# ==============================================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(mice)
library(BuyseTest)
library(npmv)

set.seed(2026)

# ==============================================================================
# 1. READ DATA
# ==============================================================================
pebles <- read.csv2("../results/PEBLES_Extract_Loise.csv",na.strings = c("", "NA"))
print(dim(pebles))
# View(pebles)   # commented out: requires RStudio; uncomment to inspect
                 # the extract interactively.
# ==============================================================================
# 2. CLEAN VARIABLES
# ==============================================================================

candidate_vars <- c(
  "Background_Pain",
  "BEBS_Total",
  "IScoreEB_Clinician",
  "IScoreEB_Patient",
  "IScoreEB_Total",
  "LIS_freq",
  "LIS_dur",
  "LIS_sev",
  "LIS_cons",
  "LIS_dist",
  "LIS_surf",
  "QOLEB_functioning",
  "QOLEB_emotions",
  "QOLEB_total"
)

pebles <- pebles %>%
  mutate(Review_Date = dmy(Review_Date),
    
    # Treat UNK as missing for exploratory analysis
     Background_Pain = na_if(Background_Pain, "UNK"),
    
    # Convert clinical scores to numeric
    across(all_of(candidate_vars),~ as.numeric(.)),
    
    Review_Age = as.numeric(Review_Age),
    
    Gender = factor(Gender,levels = c(1, 2),labels = c("Male", "Female")),
    
    Subtype = factor(Subtype)
  )

str(pebles[
  c(
    "ID",
    "Review_No",
    "Review_Date",
    "Gender",
    "Subtype",
    "BEBS_Total",
    "LIS_sev",
    "QOLEB_total"
  )
])
# ==============================================================================
# 3. SELECT ONE OBSERVATION PER PATIENT
#    EARLIEST AVAILABLE REVIEW
# ==============================================================================

baseline <- pebles %>%
  arrange(ID, Review_Date) %>%
  group_by(ID) %>%
  slice_head(n = 1) %>%
  ungroup()

print(n_distinct(baseline$ID))
print(table(baseline$Review_No, useNA = "ifany"))


# ==============================================================================
# 4. RESTRICT TO SUBTYPE 1 VS SUBTYPE 2
# ==============================================================================

outcomes <- c("BEBS_Total", "LIS_sev", "QOLEB_total")

dat <- baseline %>%
  filter(Subtype %in% c("1", "2")) %>%
  mutate(Group = factor(Subtype,levels = c("1", "2"),labels = c("Severe", "Intermediate"))) %>%
  select(ID, Group,Gender, Review_Age, all_of(outcomes))

print(table(dat$Group))
print(nrow(dat))


# ==============================================================================
# 5. DESCRIPTIVE STATISTICS BY GROUP
# ==============================================================================

descriptive_summary <- dat %>%
  group_by(Group) %>%
  summarise(
    across(
      all_of(outcomes),
      list(
        n = ~ sum(!is.na(.)),
        mean = ~ mean(., na.rm = TRUE),
        sd = ~ sd(., na.rm = TRUE),
        median = ~ median(., na.rm = TRUE),
        q1 = ~ quantile(., 0.25, na.rm = TRUE),
        q3 = ~ quantile(., 0.75, na.rm = TRUE)
      )
    ),
    .groups = "drop"
  )

print(descriptive_summary)

# ==============================================================================
# BOXPLOTS: distribution of the three outcomes by subtype
# ==============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)   # for arranging the three panels + shared title

plot_long <- dat %>%
  select(Group, BEBS_Total, LIS_sev, QOLEB_total) %>%
  pivot_longer(cols = c(BEBS_Total, LIS_sev, QOLEB_total),
               names_to = "outcome", values_to = "score") %>%
  filter(!is.na(score))

n_labels <- plot_long %>%
  count(outcome, Group, name = "n") %>%
  mutate(label = paste0(Group, "\n(n=", n, ")"))

make_panel <- function(outcome_name, title_text) {
  
  d <- plot_long %>% filter(outcome == outcome_name)
  
  # Get the group labels with n counts for THIS outcome
  lbls <- n_labels %>% filter(outcome == outcome_name)
  label_map <- setNames(lbls$label, as.character(lbls$Group))
  
  ggplot(d, aes(x = Group, y = score, fill = Group)) +
    geom_boxplot(width = 0.55, outlier.shape = 1, outlier.size = 2,
                 colour = "grey25") +
    scale_fill_manual(values = c("Severe"       = "#C56A65",   # salmon-red
                                 "Intermediate" = "#5A8BB0")) + # muted blue
    scale_x_discrete(labels = label_map) +
    labs(title = title_text, x = NULL, y = "Score") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "plain", size = 14),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position  = "none",
      axis.text.x      = element_text(size = 11),
    )
}

p_bebs  <- make_panel("BEBS_Total",  "BEBS (severity)")
p_lis   <- make_panel("LIS_sev",     "LIS_sev (itch severity)")
p_qoleb <- make_panel("QOLEB_total", "QOLEB (quality of life)")

boxplot_fig <- (p_bebs | p_lis | p_qoleb) +
  plot_annotation(
    title = "Distribution of outcomes by subtype (earliest available review)",
    theme = theme(plot.title = element_text(hjust = 0.5,
                                            face = "bold",
                                            size = 16,
                                            margin = margin(b = 10)))
  )

print(boxplot_fig)
ggsave("boxplot_outcomes_by_subtype.png", boxplot_fig,
       width = 12, height = 4.5, dpi = 200, bg = "white")


# ==============================================================================
# 6. MARGINAL MISSINGNESS
# ==============================================================================

missing_summary <- dat %>%
  summarise(across(all_of(outcomes),
      list(
        observed = ~ sum(!is.na(.)),
        missing = ~ sum(is.na(.)),
        pctmissing = ~ 100 * mean(is.na(.))
      )
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to = c("Outcome", ".value"),
    names_pattern = "(.*)_(observed|missing|pctmissing)"
  )
print(missing_summary)

# ==============================================================================
# JOINT MISSINGNESS PATTERNS
#
# O = observed
# M = missing
# ==============================================================================

joint_missingness <- dat %>%
  mutate(
    BEBS_status =
      ifelse(is.na(BEBS_Total), "M", "O"),
    
    LIS_status =
      ifelse(is.na(LIS_sev), "M", "O"),
    
    QOLEB_status =
      ifelse(is.na(QOLEB_total), "M", "O"),
    
    pattern =
      paste0(
        BEBS_status,
        LIS_status,
        QOLEB_status
      )
  ) %>%
  count(pattern) %>%
  mutate(
    percent = 100 * n / sum(n)
  ) %>%
  arrange(desc(n))

print(joint_missingness)


# ==============================================================================
# MISSINGNESS BY GROUP
# ==============================================================================

missing_by_group <- dat %>%
  group_by(Group) %>%
  summarise(across(all_of(outcomes),~ 100 * mean(is.na(.))),
    .groups = "drop")
print(missing_by_group)

# ==============================================================================
# COMPLETE-CASE DATA
# ==============================================================================

cc <- dat %>%filter(complete.cases(across(all_of(outcomes))))
print(table(cc$Group))

# ==============================================================================
# GPC-NP: COMPLETE CASE
#
# Lower scores are assumed favourable.
# ==============================================================================

gpc_np_cc <- BuyseTest(Group ~cont(BEBS_Total,threshold = 0,operator = "<0") +
    cont(LIS_sev,threshold = 0,operator = "<0") +
    cont(QOLEB_total,threshold = 0,operator = "<0"),
    data = cc,
    hierarchical = FALSE,
    method.inference = "u-statistic",
    trace = FALSE
)
print(summary(gpc_np_cc))
print(model.tables(gpc_np_cc))
print(confint(gpc_np_cc))


# ==============================================================================
# GPC-NP COMPLETE CASE USING THE SAME INFERENCE AS THE SIMULATION
# ==============================================================================

# 1. Keep complete cases for the three outcomes
dat_cc <- dat %>%
  filter(
    !is.na(BEBS_Total),
    !is.na(LIS_sev),
    !is.na(QOLEB_total)
  )

table(dat_cc$Group)
nrow(dat_cc)

# Expected:
# Severe        = 15
# Intermediate  = 17
# Total         = 32


# ------------------------------------------------------------------------------
# 2. Put the data into the long format required by gpc_nonprioritized()
# ------------------------------------------------------------------------------

dat_cc_sim <- dat_cc %>%
  mutate(
    # Simulation function treats V as treatment and P as control
    # We want Intermediate = treatment, Severe = control
    Group_sim = ifelse(Group == "Intermediate", "V", "P")
  ) %>%
  select(
    ID,
    Group_sim,
    BEBS_Total,
    LIS_sev,
    QOLEB_total
  ) %>%
  rename(
    Y1 = BEBS_Total,
    Y2 = LIS_sev,
    Y3 = QOLEB_total
  ) %>%
  pivot_longer(
    cols = c(Y1, Y2, Y3),
    names_to = "Time",
    values_to = "Outcome"
  ) %>%
  rename(Group = Group_sim)


# ------------------------------------------------------------------------------
# 3. Run exactly the GPC-NP function used in the simulation
# ------------------------------------------------------------------------------

gpc_np_cc_sim <- gpc_nonprioritized(
  data = dat_cc_sim,
  repeated = c("Y1", "Y2", "Y3"),
  best = "lower",   # lower BEBS, LIS and QOLEB scores are better
  side = 2          # two-sided test
)

gpc_np_cc_sim

# ------------------------------------------------------------------------------
# 4. Extract inference
# ------------------------------------------------------------------------------

estimate <- gpc_np_cc_sim$NB
variance <- gpc_np_cc_sim$Var
se       <- sqrt(variance)
z_value  <- estimate / se
p_value  <- gpc_np_cc_sim$p_value

ci_lower <- estimate - qnorm(0.975) * se
ci_upper <- estimate + qnorm(0.975) * se

results_sim_inference <- data.frame(
  Estimate = estimate,
  Variance = variance,
  SE = se,
  Z = z_value,
  CI_lower = ci_lower,
  CI_upper = ci_upper,
  P_value = p_value
)

print(results_sim_inference)


# ==============================================================================
# 12. GPC-NP: AVAILABLE CASE
#
# Preliminary BuyseTest implementation.
#
# Gehan scoring:
#   comparisons unavailable because of missing outcomes
#   are treated as uninformative.
#
# No redistribution of uninformative comparisons.
# ==============================================================================

gpc_np_ac <- BuyseTest(Group ~cont(BEBS_Total, threshold = 0,operator = "<0") +
    cont(LIS_sev,threshold = 0,operator = "<0") +
      cont(QOLEB_total,threshold = 0,operator = "<0"),
  data = dat,
  hierarchical = FALSE,
  scoring.rule = "Gehan",
  correction.uninf = FALSE,
  method.inference = "u-statistic",
  trace = FALSE
)
print(summary(gpc_np_ac))
print(model.tables(gpc_np_ac))
print(confint(gpc_np_ac))

# ==============================================================================
# GPC-NP AVAILABLE CASE USING THE SAME INFERENCE AS THE SIMULATION
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Use ALL patients -- do NOT remove incomplete outcome vectors
# ------------------------------------------------------------------------------

dat_ac <- dat

table(dat_ac$Group)
nrow(dat_ac)

# Expected:
# Severe        = 21
# Intermediate  = 25
# Total         = 46


# ------------------------------------------------------------------------------
# 2. Convert to the format required by the simulation GPC function
# ------------------------------------------------------------------------------

dat_ac_sim <- dat_ac %>%
  mutate(
    # In the simulation function:
    # V = treatment
    # P = control
    #
    # Here:
    # Intermediate = treatment
    # Severe       = control
    Group_sim = ifelse(Group == "Intermediate", "V", "P")
  ) %>%
  select(
    ID,
    Group_sim,
    BEBS_Total,
    LIS_sev,
    QOLEB_total
  ) %>%
  rename(
    Y1 = BEBS_Total,
    Y2 = LIS_sev,
    Y3 = QOLEB_total
  ) %>%
  pivot_longer(
    cols = c(Y1, Y2, Y3),
    names_to = "Time",
    values_to = "Outcome"
  ) %>%
  rename(Group = Group_sim)


# Check that the missing values have NOT disappeared
dat_ac_sim %>%
  group_by(Time) %>%
  summarise(
    n_total = n(),
    n_observed = sum(!is.na(Outcome)),
    n_missing = sum(is.na(Outcome))
  )

# ------------------------------------------------------------------------------
# 3. Run GPC-NP with the simulation inference
# ------------------------------------------------------------------------------

gpc_np_ac_sim <- gpc_nonprioritized(
  data = dat_ac_sim,
  repeated = c("Y1", "Y2", "Y3"),
  best = "lower",
  side = 2
)

print(gpc_np_ac_sim)

# ------------------------------------------------------------------------------
# 4. Extract estimate, SE, Z statistic, CI and p-value
# ------------------------------------------------------------------------------

estimate <- gpc_np_ac_sim$NB
variance <- gpc_np_ac_sim$Var
se       <- sqrt(variance)

z_value <- estimate / se

p_value <- gpc_np_ac_sim$p_value

ci_lower <- estimate - qnorm(0.975) * se
ci_upper <- estimate + qnorm(0.975) * se

results_ac_sim <- data.frame(
  Method   = "GPC-NP",
  Strategy = "AC",
  N        = nrow(dat_ac),
  Estimate = estimate,
  Variance = variance,
  SE       = se,
  Z        = z_value,
  CI_lower = ci_lower,
  CI_upper = ci_upper,
  P_value  = p_value
)

print(results_ac_sim)
# ==============================================================================
# 13. GPC-P: COMPLETE CASE
#
# PRELIMINARY PRIORITY:
#
#   BEBS > LIS severity > QOLEB
#
# IMPORTANT:
# Thresholds are currently zero.
# This is exploratory only.
#
# Clinical relevance thresholds should be discussed
# before treating this as the final prioritized analysis.
# ==============================================================================

gpc_p_cc <- BuyseTest(Group ~cont(BEBS_Total,threshold = 0,operator = "<0") +
    cont(LIS_sev,threshold = 0,operator = "<0") +
    cont(QOLEB_total,threshold = 0,operator = "<0"),
  
  data = cc,
  hierarchical = TRUE,
  method.inference = "u-statistic",
  trace = FALSE
)

print(summary(gpc_p_cc))
print(model.tables(gpc_p_cc))
print(confint(gpc_p_cc))

# ==============================================================================
# PRIORITISED GPC - COMPLETE CASE
# SAME INFERENCE AS THE SIMULATION
# ==============================================================================

gpc_p_cc_sim <- gpc_prioritized(
  data = dat_cc_sim,
  repeated = c("Y1", "Y2", "Y3"),
  best = "lower",
  side = 2
)

print(gpc_p_cc_sim)

estimate_p_cc <- gpc_p_cc_sim$NB
variance_p_cc <- gpc_p_cc_sim$Var
se_p_cc       <- sqrt(variance_p_cc)
z_p_cc        <- estimate_p_cc / se_p_cc
p_p_cc        <- gpc_p_cc_sim$p_value

ci_lower_p_cc <- estimate_p_cc - qnorm(0.975) * se_p_cc
ci_upper_p_cc <- estimate_p_cc + qnorm(0.975) * se_p_cc

results_p_cc <- data.frame(
  Method   = "GPC-P",
  Strategy = "CC",
  N        = nrow(dat_cc),
  Estimate = estimate_p_cc,
  Variance = variance_p_cc,
  SE       = se_p_cc,
  Z        = z_p_cc,
  CI_lower = ci_lower_p_cc,
  CI_upper = ci_upper_p_cc,
  P_value  = p_p_cc
)

print(results_p_cc)
# ==============================================================================
# 14. GPC-P: AVAILABLE CASE
# ==============================================================================

gpc_p_ac <- BuyseTest(Group ~cont(BEBS_Total, threshold = 0, operator = "<0") +
    cont(LIS_sev,threshold = 0,operator = "<0") +
    cont(QOLEB_total,threshold = 0,operator = "<0"),
  
  data = dat,
  hierarchical = TRUE,
  scoring.rule = "Gehan",
  correction.uninf = FALSE,
  method.inference = "u-statistic",
  trace = FALSE
)

print(summary(gpc_p_ac))
print(model.tables(gpc_p_ac))
print(confint(gpc_p_ac))

# ==============================================================================
# PRIORITISED GPC - AVAILABLE CASE
# SAME INFERENCE AS THE SIMULATION
# ==============================================================================

gpc_p_ac_sim <- gpc_prioritized(
  data = dat_ac_sim,
  repeated = c("Y1", "Y2", "Y3"),
  best = "lower",
  side = 2
)

print(gpc_p_ac_sim)

estimate_p_ac <- gpc_p_ac_sim$NB
variance_p_ac <- gpc_p_ac_sim$Var
se_p_ac       <- sqrt(variance_p_ac)
z_p_ac        <- estimate_p_ac / se_p_ac
p_p_ac        <- gpc_p_ac_sim$p_value

ci_lower_p_ac <- estimate_p_ac - qnorm(0.975) * se_p_ac
ci_upper_p_ac <- estimate_p_ac + qnorm(0.975) * se_p_ac

results_p_ac <- data.frame(
  Method   = "GPC-P",
  Strategy = "AC",
  N        = nrow(dat_ac),
  Estimate = estimate_p_ac,
  Variance = variance_p_ac,
  SE       = se_p_ac,
  Z        = z_p_ac,
  CI_lower = ci_lower_p_ac,
  CI_upper = ci_upper_p_ac,
  P_value  = p_p_ac
)

print(results_p_ac)

# ==============================================================================
# 15. NPMV FUNCTION
#
# Adapted directly from the simulation study.
#
# Returns:
#   p-value
#   F-statistic
#   df1
#   df2
# ==============================================================================

npmv_analysis <- function(
    Y_c,
    Y_t,
    permreps = 10000
) {
  
  tryCatch({
    
    n_c <- nrow(Y_c)
    n_t <- nrow(Y_t)
    
    N <- n_c + n_t
    
    npmv_df <- data.frame(
      Group = factor(
        c(
          rep("P", n_c),
          rep("V", n_t)
        ),
        levels = c("V", "P")
      ),
      
      Y1 = c(
        Y_c[, 1],
        Y_t[, 1]
      ),
      
      Y2 = c(
        Y_c[, 2],
        Y_t[, 2]
      ),
      
      Y3 = c(
        Y_c[, 3],
        Y_t[, 3]
      )
    )
    
    # Same rule as simulation:
    #
    # N < 30:
    # permutation/randomization
    #
    # N >= 30:
    # F approximation
    
    use_perm <- N < 30
    
    npmv_out <- npmv::nonpartest(
      Y1 | Y2 | Y3 ~ Group,
      
      data = npmv_df,
      
      test = c(
        1,
        0,
        0,
        0
      ),
      
      permtest = use_perm,
      
      permreps = permreps,
      
      plots = FALSE
    )
    
    res <- npmv_out$results[1, ]
    
    if (use_perm) {
      
      perm_col <- grep(
        "permutation|perm|P\\.T\\.",
        colnames(res),
        ignore.case = TRUE,
        value = TRUE
      )
      
      p_val <-
        if (length(perm_col) > 0) {
          
          res[[perm_col[1]]]
          
        } else {
          
          res[[ncol(res)]]
          
        }
      
    } else {
      
      p_val <- res[["P-value"]]
      
    }
    
    list(
      p_value = p_val,
      
      f_stat =
        res[["Test Statistic"]],
      
      df1 =
        res[["df1"]],
      
      df2 =
        res[["df2"]],
      
      N = N,
      
      perm_used =
        use_perm
    )
    
  },
  
  error = function(e) {
    
    message(
      "NPMV error: ",
      e$message
    )
    
    list(
      p_value = NA_real_,
      f_stat = NA_real_,
      df1 = NA_real_,
      df2 = NA_real_,
      N = NA_real_,
      perm_used = NA
    )
    
  })
}


# ==============================================================================
# 16. NPMV: COMPLETE CASE
#
# Severe = control
# Intermediate = treatment
# ==============================================================================

Y_c_cc <- cc %>%
  filter(Group == "Severe") %>%
  select(all_of(outcomes)) %>%
  as.matrix()

Y_t_cc <- cc %>%
  filter(Group == "Intermediate") %>%
  select(all_of(outcomes)) %>%
  as.matrix()


npmv_cc <- npmv_analysis(
  Y_c = Y_c_cc,
  Y_t = Y_t_cc
)

print(npmv_cc)

npmv_ac_direct <- npmv::nonpartest(
  BEBS_Total | LIS_sev | QOLEB_total ~ Group,
  data = dat,
  test = c(1, 0, 0, 0),
  permtest = FALSE,
  plots = FALSE
)

print(npmv_ac_direct$results[1, ])

use_perm <- N < 30

npmv_out <- npmv::nonpartest(
  BEBS_Total | LIS_sev | QOLEB_total ~ Group,
  data     = dat,
  test     = c(1, 0, 0, 0),
  permtest = FALSE,
  permreps = permreps,   # only used when use_perm = TRUE
  plots    = FALSE
)

res <- npmv_out$results[1, ]

# ==============================================================================
# 17. NPMV: AVAILABLE CASE
#
# IMPORTANT:
#
# This deliberately follows the simulation implementation:
# incomplete outcome matrices are passed directly to nonpartest().
# ==============================================================================

Y_c_ac <- dat %>%
  filter(Group == "Severe") %>%
  select(all_of(outcomes)) %>%
  as.matrix()

Y_t_ac <- dat %>%
  filter(Group == "Intermediate") %>%
  select(all_of(outcomes)) %>%
  as.matrix()


npmv_ac <- npmv_analysis(
  Y_c = Y_c_ac,
  Y_t = Y_t_ac
)

print(npmv_ac)


# ==============================================================================
# 18. MULTIPLE IMPUTATION DATA
#
# Analysis variable:
#   Group (Subtype: Severe vs Intermediate)
#
# Outcomes:
#   BEBS_Total
#   LIS_sev
#   QOLEB_total
#
# Auxiliary variables:
#   Review_Age
#   Gender
#
# Age and gender are included in the imputation model as auxiliary
# predictors but are not part of the substantive GPC/NPMV analysis.
# ==============================================================================

mi_dat <- dat %>%
  select(
    Group,
    Gender,
    Review_Age,
    BEBS_Total,
    LIS_sev,
    QOLEB_total
  )


# ==============================================================================
# 19. INITIALIZE MICE
# ==============================================================================

ini <- mice(
  mi_dat,
  maxit = 0,
  printFlag = FALSE
)

meth <- ini$method
pred <- ini$predictorMatrix


# ==============================================================================
# 20. IMPUTATION METHODS
# ==============================================================================

# Fully observed variables are not imputed
meth["Group"] <- ""
meth["Gender"] <- ""
meth["Review_Age"] <- ""

# Missing clinical outcomes are imputed using PMM
meth["BEBS_Total"] <- "pmm"
meth["LIS_sev"] <- "pmm"
meth["QOLEB_total"] <- "pmm"


# ==============================================================================
# 21. PREDICTOR MATRIX
#
# Each incomplete outcome is imputed using:
#
#   - Group
#   - Gender
#   - Age
#   - the other clinical outcomes
#
# Group, gender and age themselves are not imputed.
# ==============================================================================

pred[,] <- 0

pred[
  "BEBS_Total",
  c(
    "Group",
    "Gender",
    "Review_Age",
    "LIS_sev",
    "QOLEB_total"
  )
] <- 1

pred[
  "LIS_sev",
  c(
    "Group",
    "Gender",
    "Review_Age",
    "BEBS_Total",
    "QOLEB_total"
  )
] <- 1

pred[
  "QOLEB_total",
  c(
    "Group",
    "Gender",
    "Review_Age",
    "BEBS_Total",
    "LIS_sev"
  )
] <- 1


cat("\n====================================\n")
cat("MI METHOD VECTOR\n")
cat("====================================\n")
print(meth)

cat("\n====================================\n")
cat("MI PREDICTOR MATRIX\n")
cat("====================================\n")
print(pred)


# ==============================================================================
# 22. RUN MULTIPLE IMPUTATION
# ==============================================================================

M <- 40

imp <- mice(
  mi_dat,
  m = M,
  method = meth,
  predictorMatrix = pred,
  maxit = 20,
  seed = 2026,
  printFlag = FALSE
)

print(imp)
# ==============================================================================
# 23. BASIC MI DIAGNOSTICS
# ==============================================================================

print(imp)

plot(imp)

densityplot(
  imp,
  ~ BEBS_Total +
    LIS_sev +
    QOLEB_total
)


# ==============================================================================
# 24. CREATE COMPLETED DATASETS
# ==============================================================================

mi_list <- lapply(
  1:M,
  
  function(i) {
    
    complete(
      imp,
      i
    )
    
  }
)


# ==============================================================================
# 25. RUN GPC-NP ON EACH IMPUTATION
# ==============================================================================

gpc_np_mi_models <- lapply(
  mi_list,
  
  function(d) {
    
    BuyseTest(
      Group ~
        cont(
          BEBS_Total,
          threshold = 0,
          operator = "<0"
        ) +
        cont(
          LIS_sev,
          threshold = 0,
          operator = "<0"
        ) +
        cont(
          QOLEB_total,
          threshold = 0,
          operator = "<0"
        ),
      
      data = d,
      
      hierarchical = FALSE,
      
      method.inference =
        "u-statistic",
      
      trace = FALSE
    )
    
  }
)


# ==============================================================================
# 26. RUN GPC-P ON EACH IMPUTATION
# ==============================================================================

gpc_p_mi_models <- lapply(
  mi_list,
  
  function(d) {
    
    BuyseTest(
      Group ~
        cont(
          BEBS_Total,
          threshold = 0,
          operator = "<0"
        ) +
        cont(
          LIS_sev,
          threshold = 0,
          operator = "<0"
        ) +
        cont(
          QOLEB_total,
          threshold = 0,
          operator = "<0"
        ),
      
      data = d,
      
      hierarchical = TRUE,
      
      method.inference =
        "u-statistic",
      
      trace = FALSE
    )
    
  }
)


# ==============================================================================
# 27. EXTRACT FINAL GLOBAL GPC RESULT
#
# BuyseTest's final cumulative row contains the global result.
# ==============================================================================

extract_gpc <- function(model) {
  
  x <- data.frame(
    confint(model)
  )
  
  x <- x[
    nrow(x),
    ,
    drop = FALSE
  ]
  
  data.frame(
    estimate = x$estimate,
    se = x$se
  )
}


gpc_np_mi_values <- bind_rows(
  lapply(
    gpc_np_mi_models,
    extract_gpc
  )
)


gpc_p_mi_values <- bind_rows(
  lapply(
    gpc_p_mi_models,
    extract_gpc
  )
)


# ==============================================================================
# 28. RUBIN POOLING FOR GPC
#
# This follows the same basic pooling principle as the simulation:
#
#   Qbar = average estimate
#   W    = average within-imputation variance
#   B    = between-imputation variance
#   T    = W + (1 + 1/M)B
#
# We report pooled estimate, SE and p-value.
# ==============================================================================

rubin_gpc <- function(
    estimates,
    ses
) {
  
  valid <-
    is.finite(estimates) &
    is.finite(ses) &
    ses > 0
  
  q <- estimates[valid]
  
  u <- ses[valid]^2
  
  m <- length(q)
  
  Qbar <- mean(q)
  
  W <- mean(u)
  
  B <- var(q)
  
  Tvar <-
    W +
    (1 + 1 / m) * B
  
  SE <- sqrt(Tvar)
  
  if (B <= 1e-12) {
    
    df <- Inf
    
  } else {
    
    r <-
      ((1 + 1 / m) * B) /
      W
    
    df <-
      (m - 1) *
      (1 + 1 / r)^2
    
  }
  
  statistic <-
    Qbar / SE
  
  if (is.infinite(df)) {
    
    p <-
      2 *
      pnorm(
        -abs(statistic)
      )
    
  } else {
    
    p <-
      2 *
      pt(
        -abs(statistic),
        df = df
      )
    
  }
  
  list(
    estimate = Qbar,
    se = SE,
    df = df,
    statistic = statistic,
    p_value = p,
    W = W,
    B = B,
    total_variance = Tvar
  )
}


gpc_np_mi_pooled <- rubin_gpc(
  estimates =
    gpc_np_mi_values$estimate,
  
  ses =
    gpc_np_mi_values$se
)


gpc_p_mi_pooled <- rubin_gpc(
  estimates =
    gpc_p_mi_values$estimate,
  
  ses =
    gpc_p_mi_values$se
)


cat("\n====================================\n")
cat("GPC-NP: POOLED MI\n")
cat("====================================\n")

print(gpc_np_mi_pooled)


cat("\n====================================\n")
cat("GPC-P: POOLED MI\n")
cat("====================================\n")

print(gpc_p_mi_pooled)


# ==============================================================================
# 29. NPMV ON EACH IMPUTED DATASET
# ==============================================================================

npmv_mi_results <- lapply(
  mi_list,
  
  function(d) {
    
    Y_c_imp <- d %>%
      filter(
        Group == "Severe"
      ) %>%
      select(
        all_of(outcomes)
      ) %>%
      as.matrix()
    
    
    Y_t_imp <- d %>%
      filter(
        Group == "Intermediate"
      ) %>%
      select(
        all_of(outcomes)
      ) %>%
      as.matrix()
    
    
    npmv_analysis(
      Y_c = Y_c_imp,
      Y_t = Y_t_imp
    )
    
  }
)


# ==============================================================================
# 30. EXTRACT NPMV MI STATISTICS
# ==============================================================================

MI_NPMV_fstat <- sapply(
  npmv_mi_results,
  function(x) x$f_stat
)

MI_NPMV_df1 <- sapply(
  npmv_mi_results,
  function(x) x$df1
)

MI_NPMV_pvalue <- sapply(
  npmv_mi_results,
  function(x) x$p_value
)


npmv_mi_diagnostics <- data.frame(
  Imputation = 1:M,
  F_statistic = MI_NPMV_fstat,
  df1 = MI_NPMV_df1,
  p_value = MI_NPMV_pvalue
)

cat("\n====================================\n")
cat("NPMV RESULTS IN EACH IMPUTATION\n")
cat("====================================\n")

print(npmv_mi_diagnostics)


# ==============================================================================
# 31. D2 POOLING FOR NPMV MI
#
# Same implementation used in simulation_functions.R
#
# IMPORTANT:
# We pool the test statistics, NOT the p-values.
# ==============================================================================

pool_npmv_D2 <- function(
    f_stats,
    df1_values
) {
  
  valid <-
    is.finite(f_stats) &
    is.finite(df1_values) &
    f_stats > 0 &
    df1_values > 0
  
  if (sum(valid) < 2) {
    
    return(
      list(
        D2 = NA_real_,
        K = NA_real_,
        df2 = NA_real_,
        ARIV2 = NA_real_,
        p_value = NA_real_
      )
    )
    
  }
  
  
  K <-
    mean(
      df1_values[valid]
    )
  
  
  W <-
    K *
    f_stats[valid]
  
  
  m_v <-
    sum(valid)
  
  
  W_bar <-
    mean(W)
  
  
  sqrt_W_bar <-
    mean(
      sqrt(W)
    )
  
  
  ARIV2 <-
    ((m_v + 1) /
       (m_v * (m_v - 1))) *
    sum(
      (
        sqrt(W) -
          sqrt_W_bar
      )^2
    )
  
  
  # If there is essentially no between-imputation
  # variability, the pooled statistic tends toward
  # the ordinary F-statistic.
  
  if (
    !is.finite(ARIV2) ||
    ARIV2 <= 1e-12
  ) {
    
    D2 <-
      W_bar / K
    
    return(
      list(
        D2 = D2,
        K = K,
        df2 = Inf,
        ARIV2 = ARIV2,
        p_value =
          pf(
            D2,
            df1 = K,
            df2 = 1e10,
            lower.tail = FALSE
          )
      )
    )
    
  }
  
  
  D2 <-
    (
      W_bar / K -
        (m_v + 1) /
        (m_v - 1) *
        ARIV2
    ) /
    (
      1 + ARIV2
    )
  
  
  if (
    !is.finite(D2) ||
    D2 < 0
  ) {
    
    return(
      list(
        D2 = NA_real_,
        K = K,
        df2 = NA_real_,
        ARIV2 = ARIV2,
        p_value = NA_real_
      )
    )
    
  }
  
  
  v2 <-
    K^(-3 / m_v) *
    (m_v - 1) *
    (
      1 +
        1 / ARIV2
    )^2
  
  
  p <-
    pf(
      D2,
      df1 = K,
      df2 = v2,
      lower.tail = FALSE
    )
  
  
  list(
    D2 = D2,
    K = K,
    df2 = v2,
    ARIV2 = ARIV2,
    p_value = p
  )
}


npmv_mi_D2 <- pool_npmv_D2(
  f_stats =
    MI_NPMV_fstat,
  
  df1_values =
    MI_NPMV_df1
)


cat("\n====================================\n")
cat("NPMV: MI D2 POOLED RESULT\n")
cat("====================================\n")

print(npmv_mi_D2)


# ==============================================================================
# 32. HELPER TO EXTRACT GPC CC/AC RESULTS
# ==============================================================================

get_final_gpc <- function(model) {
  
  x <- data.frame(
    confint(model)
  )
  
  x <- x[
    nrow(x),
    ,
    drop = FALSE
  ]
  
  list(
    estimate =
      x$estimate,
    
    se =
      x$se,
    
    p_value =
      x$p.value
  )
}


np_cc <- get_final_gpc(
  gpc_np_cc
)

np_ac <- get_final_gpc(
  gpc_np_ac
)

p_cc <- get_final_gpc(
  gpc_p_cc
)

p_ac <- get_final_gpc(
  gpc_p_ac
)


# ==============================================================================
# 33. FINAL P-VALUE COMPARISON TABLE
# ==============================================================================

final_pvalues <- bind_rows(
  
  data.frame(
    Procedure = "GPC-NP",
    Method = "Complete case",
    N = nrow(cc),
    Estimate = np_cc$estimate,
    Statistic = NA_real_,
    P_value = np_cc$p_value
  ),
  
  data.frame(
    Procedure = "GPC-NP",
    Method = "Available case",
    N = nrow(dat),
    Estimate = np_ac$estimate,
    Statistic = NA_real_,
    P_value = np_ac$p_value
  ),
  
  data.frame(
    Procedure = "GPC-NP",
    Method = "Multiple imputation",
    N = nrow(dat),
    Estimate =
      gpc_np_mi_pooled$estimate,
    Statistic =
      gpc_np_mi_pooled$statistic,
    P_value =
      gpc_np_mi_pooled$p_value
  ),
  
  data.frame(
    Procedure = "GPC-P",
    Method = "Complete case",
    N = nrow(cc),
    Estimate = p_cc$estimate,
    Statistic = NA_real_,
    P_value = p_cc$p_value
  ),
  
  data.frame(
    Procedure = "GPC-P",
    Method = "Available case",
    N = nrow(dat),
    Estimate = p_ac$estimate,
    Statistic = NA_real_,
    P_value = p_ac$p_value
  ),
  
  data.frame(
    Procedure = "GPC-P",
    Method = "Multiple imputation",
    N = nrow(dat),
    Estimate =
      gpc_p_mi_pooled$estimate,
    Statistic =
      gpc_p_mi_pooled$statistic,
    P_value =
      gpc_p_mi_pooled$p_value
  ),
  
  data.frame(
    Procedure = "NPMV",
    Method = "Complete case",
    N = npmv_cc$N,
    Estimate = NA_real_,
    Statistic =
      npmv_cc$f_stat,
    P_value =
      npmv_cc$p_value
  ),
  
  data.frame(
    Procedure = "NPMV",
    Method = "Available case",
    N = npmv_ac$N,
    Estimate = NA_real_,
    Statistic =
      npmv_ac$f_stat,
    P_value =
      npmv_ac$p_value
  ),
  
  data.frame(
    Procedure = "NPMV",
    Method = "Multiple imputation (D2)",
    N = nrow(dat),
    Estimate = NA_real_,
    Statistic =
      npmv_mi_D2$D2,
    P_value =
      npmv_mi_D2$p_value
  )
)


cat("\n\n")
cat("============================================================\n")
cat("FINAL PEBLES RESULTS: P-VALUE COMPARISON\n")
cat("============================================================\n")

print(
  final_pvalues,
  digits = 5
)


# ==============================================================================
# 34. ADD SIGNIFICANCE INDICATOR
# ==============================================================================

final_pvalues <- final_pvalues %>%
  mutate(
    Reject_0.05 =
      ifelse(
        is.na(P_value),
        NA,
        P_value < 0.05
      )
  )


cat("\n====================================\n")
cat("INFERENTIAL CONCLUSIONS\n")
cat("====================================\n")

print(
  final_pvalues,
  digits = 5
)


# ==============================================================================
# 35. SAVE TABLES
# ==============================================================================

write.csv(
  final_pvalues,
  "PEBLES_final_pvalue_comparison.csv",
  row.names = FALSE
)

write.csv(
  missing_summary,
  "PEBLES_missingness_summary.csv",
  row.names = FALSE
)

write.csv(
  joint_missingness,
  "PEBLES_joint_missingness.csv",
  row.names = FALSE
)

write.csv(
  missing_by_group,
  "PEBLES_missingness_by_group.csv",
  row.names = FALSE
)

write.csv(
  npmv_mi_diagnostics,
  "PEBLES_NPMV_MI_diagnostics.csv",
  row.names = FALSE
)


# ==============================================================================
# 36. OUTCOME DISTRIBUTION PLOT
# ==============================================================================

plot_data <- dat %>%
  select(
    Group,
    all_of(outcomes)
  ) %>%
  pivot_longer(
    cols =
      all_of(outcomes),
    
    names_to =
      "Outcome",
    
    values_to =
      "Value"
  )


ggplot(
  plot_data,
  aes(
    x = Group,
    y = Value
  )
) +
  geom_boxplot() +
  facet_wrap(
    ~ Outcome,
    scales = "free_y"
  ) +
  labs(
    title =
      "PEBLES outcomes by EB subtype",
    subtitle =
      "Earliest available review",
    x = NULL,
    y = "Score"
  ) +
  theme_minimal()


# ==============================================================================
# 37. P-VALUE PLOT
#
# Plot -log10(p), since the p-values may be very small.
#
# Larger values = stronger evidence against H0.
# ==============================================================================

p_plot <- final_pvalues %>%
  filter(
    !is.na(P_value),
    P_value > 0
  ) %>%
  mutate(
    minus_log10_p =
      -log10(P_value)
  )


ggplot(
  p_plot,
  aes(
    x = Method,
    y = minus_log10_p
  )
) +
  geom_col() +
  geom_hline(
    yintercept =
      -log10(0.05),
    linetype = 2
  ) +
  facet_wrap(
    ~ Procedure,
    scales = "free_x"
  ) +
  labs(
    title =
      "Comparison of inferential results",
    subtitle =
      "Dashed line corresponds to p = 0.05",
    x = NULL,
    y =
      expression(-log[10](p))
  ) +
  theme_minimal() +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )



# ==============================================================================
# MULTIPLE IMPUTATION FOR GPC-NP AND GPC-P
# USING THE SAME GPC INFERENCE AS THE SIMULATION
# ==============================================================================


# ==============================================================================
# 1. PREPARE DATA FOR MULTIPLE IMPUTATION
# ==============================================================================

mi_dat <- dat %>%
  select(
    ID,
    Group,
    Gender,
    Review_Age,
    BEBS_Total,
    LIS_sev,
    QOLEB_total
  )

# Check missingness before imputation
cat("\nMissing values before imputation:\n")
print(colSums(is.na(mi_dat)))


# ==============================================================================
# 2. SPECIFY THE MICE IMPUTATION MODEL
# ==============================================================================

meth <- make.method(mi_dat)
pred <- make.predictorMatrix(mi_dat)

# Start with no predictors
pred[,] <- 0


# ------------------------------------------------------------------------------
# Variables NOT to be imputed
# ------------------------------------------------------------------------------

meth["ID"]         <- ""
meth["Group"]      <- ""
meth["Gender"]     <- ""
meth["Review_Age"] <- ""


# ID should not be used as a predictor
pred[, "ID"] <- 0
pred["ID", ] <- 0


# ------------------------------------------------------------------------------
# Variables to be imputed using predictive mean matching
# ------------------------------------------------------------------------------

meth["BEBS_Total"]  <- "pmm"
meth["LIS_sev"]     <- "pmm"
meth["QOLEB_total"] <- "pmm"


# ------------------------------------------------------------------------------
# Predictor structure
#
# Each outcome is predicted by:
#   - treatment/subtype group
#   - gender
#   - age
#   - the other two outcomes
# ------------------------------------------------------------------------------

pred[
  "BEBS_Total",
  c("Group", "Gender", "Review_Age", "LIS_sev", "QOLEB_total")
] <- 1

pred[
  "LIS_sev",
  c("Group", "Gender", "Review_Age", "BEBS_Total", "QOLEB_total")
] <- 1

pred[
  "QOLEB_total",
  c("Group", "Gender", "Review_Age", "BEBS_Total", "LIS_sev")
] <- 1


# Display imputation specification
cat("\nImputation methods:\n")
print(meth)

cat("\nPredictor matrix:\n")
print(pred)


# ==============================================================================
# 3. RUN MULTIPLE IMPUTATION
# ==============================================================================

imp <- mice(
  data = mi_dat,
  m = 40,
  maxit = 20,
  method = meth,
  predictorMatrix = pred,
  seed = 2026,
  printFlag = FALSE
)

cat("\nMultiple imputation completed:\n")
print(imp)


# Optional diagnostic plot
plot(
  imp,
  c("BEBS_Total", "LIS_sev", "QOLEB_total")
)


# ==============================================================================
# 4. ANALYSE EACH IMPUTED DATASET
# ==============================================================================

m <- imp$m

# Storage for GPC-NP
MI_NP_est <- numeric(m)
MI_NP_var <- numeric(m)

# Storage for GPC-P
MI_P_est <- numeric(m)
MI_P_var <- numeric(m)


for (i in seq_len(m)) {
  
  # Obtain completed dataset i
  comp <- complete(imp, i)
  
  # ------------------------------------------------------------
  # Convert to the structure used in the simulation functions
  #
  # V = treatment = Intermediate
  # P = control   = Severe
  #
  # Y1 = BEBS
  # Y2 = LIS severity
  # Y3 = QOLEB
  # ------------------------------------------------------------
  
  comp_sim <- comp %>%
    mutate(
      Group_sim = ifelse(
        Group == "Intermediate",
        "V",
        "P"
      )
    ) %>%
    select(
      ID,
      Group_sim,
      BEBS_Total,
      LIS_sev,
      QOLEB_total
    ) %>%
    rename(
      Y1 = BEBS_Total,
      Y2 = LIS_sev,
      Y3 = QOLEB_total
    ) %>%
    pivot_longer(
      cols = c(Y1, Y2, Y3),
      names_to = "Time",
      values_to = "Outcome"
    ) %>%
    rename(Group = Group_sim)
  
  
  # ------------------------------------------------------------
  # GPC NON-PRIORITIZED
  # ------------------------------------------------------------
  
  out_np <- gpc_nonprioritized(
    data = comp_sim,
    repeated = c("Y1", "Y2", "Y3"),
    best = "lower",
    side = 2
  )
  
  MI_NP_est[i] <- out_np$NB
  MI_NP_var[i] <- out_np$Var
  
  
  # ------------------------------------------------------------
  # GPC PRIORITIZED
  #
  # Priority:
  # BEBS > LIS > QOLEB
  # ------------------------------------------------------------
  
  out_p <- gpc_prioritized(
    data = comp_sim,
    repeated = c("Y1", "Y2", "Y3"),
    best = "lower",
    side = 2
  )
  
  MI_P_est[i] <- out_p$NB
  MI_P_var[i] <- out_p$Var
}


# ==============================================================================
# 5. FUNCTION TO POOL GPC RESULTS USING RUBIN'S RULES
#    SAME STRUCTURE AS THE SIMULATION CODE
# ==============================================================================

pool_gpc_sim <- function(estimates, variances) {
  
  valid <- is.finite(estimates) & is.finite(variances)
  
  Q <- estimates[valid]
  U <- variances[valid]
  
  m_valid <- length(Q)
  
  if (m_valid < 2) {
    stop("Fewer than two valid imputations.")
  }
  
  # Pooled estimate
  Q_bar <- mean(Q)
  
  # Within-imputation variance
  W <- mean(U)
  
  # Between-imputation variance
  B <- var(Q)
  
  # Total variance
  T_var <- W + (1 + 1 / m_valid) * B
  
  # Standard error
  SE <- sqrt(T_var)
  
  # Relative increase in variance
  r_m <- ((1 + 1 / m_valid) * B) / W
  
  
  # ------------------------------------------------------------
  # Degrees of freedom
  # Same calculation used in simulation code
  # ------------------------------------------------------------
  
  if (!is.finite(r_m) || r_m <= .Machine$double.eps) {
    
    nu <- Inf
    
  } else if (m_valid - 1 > 4) {
    
    nu <- 4 +
      ((m_valid - 1 - 4) *
         (1 - 2 / (m_valid - 1))) /
      r_m^2
    
  } else {
    
    nu <- (m_valid - 1) *
      (1 + 1 / r_m)^2
  }
  
  
  # Test statistic
  test_stat <- Q_bar / SE
  
  
  # ------------------------------------------------------------
  # P-value and confidence interval
  # ------------------------------------------------------------
  
  if (is.infinite(nu)) {
    
    p_value <- 2 * pnorm(-abs(test_stat))
    critical <- qnorm(0.975)
    
  } else {
    
    p_value <- 2 * pt(
      -abs(test_stat),
      df = nu
    )
    
    critical <- qt(
      0.975,
      df = nu
    )
  }
  
  
  CI_lower <- Q_bar - critical * SE
  CI_upper <- Q_bar + critical * SE
  
  
  return(
    data.frame(
      Estimate = Q_bar,
      Within_variance = W,
      Between_variance = B,
      Total_variance = T_var,
      SE = SE,
      df = nu,
      Statistic = test_stat,
      CI_lower = CI_lower,
      CI_upper = CI_upper,
      P_value = p_value,
      m_used = m_valid
    )
  )
}


# ==============================================================================
# 6. POOL GPC-NP RESULTS
# ==============================================================================

gpc_np_mi_sim <- pool_gpc_sim(
  estimates = MI_NP_est,
  variances = MI_NP_var
)


# ==============================================================================
# 7. POOL GPC-P RESULTS
# ==============================================================================

gpc_p_mi_sim <- pool_gpc_sim(
  estimates = MI_P_est,
  variances = MI_P_var
)


# ==============================================================================
# 8. DISPLAY RESULTS
# ==============================================================================

cat("\n============================================================\n")
cat("GPC-NP MULTIPLE IMPUTATION RESULTS\n")
cat("============================================================\n")
print(gpc_np_mi_sim)


cat("\n============================================================\n")
cat("GPC-P MULTIPLE IMPUTATION RESULTS\n")
cat("============================================================\n")
print(gpc_p_mi_sim)


# ==============================================================================
# 9. SIMPLE FINAL TABLE
# ==============================================================================

MI_results <- bind_rows(
  
  data.frame(
    Method = "GPC-NP",
    Strategy = "MI",
    N = nrow(dat),
    Estimate = gpc_np_mi_sim$Estimate,
    SE = gpc_np_mi_sim$SE,
    CI_lower = gpc_np_mi_sim$CI_lower,
    CI_upper = gpc_np_mi_sim$CI_upper,
    P_value = gpc_np_mi_sim$P_value
  ),
  
  data.frame(
    Method = "GPC-P",
    Strategy = "MI",
    N = nrow(dat),
    Estimate = gpc_p_mi_sim$Estimate,
    SE = gpc_p_mi_sim$SE,
    CI_lower = gpc_p_mi_sim$CI_lower,
    CI_upper = gpc_p_mi_sim$CI_upper,
    P_value = gpc_p_mi_sim$P_value
  )
)

cat("\n============================================================\n")
cat("FINAL MI RESULTS TABLE\n")
cat("============================================================\n")
print(MI_results)
# ==============================================================================
# 38. SESSION INFORMATION
# ==============================================================================

cat("\n====================================\n")
cat("SESSION INFORMATION\n")
cat("====================================\n")

sessionInfo()