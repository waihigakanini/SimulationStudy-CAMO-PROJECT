#' Generate Trivariate Data
#'
#' @param n Sample size per group
#' @param rho Correlation coefficient parameter
#' @param structure "cs" (Compound Symmetry) or "ar1" (Auto-Regressive)
#' @param mu_diff Effect size (Delta). Applied to Y1, Y2, Y3.
#' @param dist "normal" (default) or "t" for multivariate t
#' @param df Degrees of freedom for t distribution (ignored if dist = "normal")
#' @return List with Y_c (control matrix) and Y_t (treatment matrix)

generate_sim_data <- function(n, rho, structure = "cs", mu_diff = 0,
                              dist = "normal", df = 5) {
  
  mu_control   <- c(0, 0, 0)
  mu_treatment <- c(mu_diff, mu_diff, mu_diff)
  
  if (structure == "cs") {
    Sigma <- matrix(rho, 3, 3)
    diag(Sigma) <- 1
  } else if (structure == "ar1") {
    Sigma <- matrix(c(1,     rho,   rho^2,
                      rho,   1,     rho,
                      rho^2, rho,   1), 3, 3, byrow = TRUE)
  } else {
    stop("structure must be 'cs' or 'ar1'")
  }
  
  if (dist == "normal") {
    Y_c <- rmvnorm(n, mu_control,   Sigma)
    Y_t <- rmvnorm(n, mu_treatment, Sigma)
    
  } else if (dist == "t") {
    if (df <= 2) stop("df must be > 2")
    Sigma_scaled <- Sigma * (df - 2) / df
    Y_c <- rmvt(n, sigma = Sigma_scaled, df = df, delta = mu_control)
    Y_t <- rmvt(n, sigma = Sigma_scaled, df = df, delta = mu_treatment)
    
  } else if (dist == "binary") {
    # Generate latent normal, then dichotomize at 0
    Y_c_latent <- rmvnorm(n, mu_control,   Sigma)
    Y_t_latent <- rmvnorm(n, mu_treatment, Sigma)
    Y_c <- (Y_c_latent > 0) * 1L
    Y_t <- (Y_t_latent > 0) * 1L
    
  } else {
    stop("dist must be 'normal', 't', or 'binary'")
  }
  
  return(list(Y_c = Y_c, Y_t = Y_t))
}

# ------------------------------------------------------------------------------

#' Apply Monotone Dropout
#'
#' @param Y_list List from generate_sim_data containing Y_c and Y_t
#' @param psi0 Baseline log-odds of dropout
#' @param psi1 Log-odds ratio for association with previous outcome (0 = MCAR)
#' @return List with Y_c, Y_t (with NAs), D_c, D_t (dropout indicators)
apply_dropout <- function(Y_list, psi0, psi1) {
  
  drop_mat <- function(Y) {
    n  <- nrow(Y)
    
    p2 <- 1 / (1 + exp(-(psi0 + psi1 * Y[, 1])))
    d2 <- runif(n) < p2
    
    p3 <- 1 / (1 + exp(-(psi0 + psi1 * Y[, 2])))
    d3 <- runif(n) < p3
    
    Y_miss <- Y
    Y_miss[d2, 2:3]     <- NA
    Y_miss[!d2 & d3, 3] <- NA
    
    D <- rep(4, n)
    D[d2]       <- 2
    D[!d2 & d3] <- 3
    
    list(Y = Y_miss, D = D)
  }
  
  ctrl <- drop_mat(Y_list$Y_c)
  trt  <- drop_mat(Y_list$Y_t)
  
  list(Y_c = ctrl$Y, D_c = ctrl$D,
       Y_t = trt$Y,  D_t = trt$D)
}

wide_to_long <- function(Y_c, Y_t) {
  df <- data.frame(
    Id    = seq_len(nrow(Y_c) + nrow(Y_t)),
    Group = factor(
      c(rep("P", nrow(Y_c)), rep("V", nrow(Y_t))),
      levels = c("V", "P")
    ),
    Y1 = c(Y_c[, 1], Y_t[, 1]),
    Y2 = c(Y_c[, 2], Y_t[, 2]),
    Y3 = c(Y_c[, 3], Y_t[, 3])
  )
  
  df_long <- df %>%
    pivot_longer(cols      = c(Y1, Y2, Y3),
                 names_to  = "Time",
                 values_to = "Outcome") %>%
    mutate(Time = factor(Time, levels = c("Y1", "Y2", "Y3")))
  
  return(df_long)
}


# ==============================================================================
# 2. GPC FUNCTIONS
# ==============================================================================

gpc_nonprioritized <- function(data, repeated = c("Y1", "Y2", "Y3"),
                               best = "higher", side = 2) {
  
  ScoreV <- function(Outcome, Trt) {
    diff_mat <- outer(Outcome, Outcome, function(a, b) a - b)
    S        <- sign(diff_mat)
    na_vec   <- is.na(Outcome)
    na_mask  <- outer(na_vec, na_vec, "|")
    S[na_mask] <- 0
    if (best == "higher") S else -S
  }
  
  data      <- data %>% dplyr::arrange(factor(Time, levels = repeated))
  data$Time <- factor(data$Time, levels = repeated)
  Outcome   <- split(data$Outcome, data$Time)
  db_trt    <- dplyr::filter(data, Time == repeated[1])
  Trt       <- ifelse(db_trt$Group == "V", 1, 0)
  nTest     <- sum(Trt == 1)
  nControl  <- sum(Trt == 0)
  nPatients <- length(Trt)
  
  if (nTest == 0 || nControl == 0)
    return(list(NB = NA, Var = NA, p_value = NA))
  
  list_npD           <- numeric()
  list_npV           <- numeric()
  listnpD_cumulative <- numeric()
  listnpV_cumulative <- numeric()
  Score_npprev       <- matrix(0, nPatients, nPatients)
  
  for (i in seq_along(Outcome)) {
    Score_V      <- ScoreV(unlist(Outcome[[i]]), Trt)
    Score_npD    <- Score_V[which(Trt == 1), which(Trt == 0), drop = FALSE]
    Score_npprev <- Score_npprev + Score_V
    
    list_npD[i]           <- mean(Score_npD)
    listnpD_cumulative[i] <- sum(list_npD[1:i])
    list_npV[i]           <- sum(rowSums(Score_V)^2) /
      (nTest * nControl * nPatients * (nPatients - 1))
    listnpV_cumulative[i] <- sum(rowSums(Score_npprev)^2) /
      (nTest * nControl * nPatients * (nPatients - 1))
  }
  
  npNB     <- listnpD_cumulative[length(Outcome)] / length(Outcome)
  npNB_var <- listnpV_cumulative[length(Outcome)] / length(Outcome)^2
  
  if (is.na(npNB_var) || npNB_var <= 0)
    return(list(NB = npNB, Var = NA, p_value = NA))
  
  p_value <- if (side == 1) {
    pnorm(-npNB / sqrt(npNB_var))
  } else {
    2 * pnorm(-abs(npNB / sqrt(npNB_var)))
  }
  
  return(list(NB = npNB, Var = npNB_var, p_value = p_value))
}

# ------------------------------------------------------------------------------

gpc_prioritized <- function(data, repeated = c("Y1", "Y2", "Y3"),
                            best = "higher", side = 2) {
  
  ScoreV <- function(Outcome, Trt) {
    diff_mat <- outer(Outcome, Outcome, function(a, b) a - b)
    S        <- sign(diff_mat)
    na_vec   <- is.na(Outcome)
    na_mask  <- outer(na_vec, na_vec, "|")
    S[na_mask] <- 0
    if (best == "higher") S else -S
  }
  
  data      <- data %>% dplyr::arrange(factor(Time, levels = repeated))
  data$Time <- factor(data$Time, levels = repeated)
  Outcome   <- split(data$Outcome, data$Time)
  db_trt    <- dplyr::filter(data, Time == repeated[1])
  Trt       <- ifelse(db_trt$Group == "V", 1, 0)
  nTest     <- sum(Trt == 1)
  nControl  <- sum(Trt == 0)
  nPatients <- length(Trt)
  
  list_T           <- numeric()
  list_C           <- numeric()
  list_D           <- numeric()
  listD_cumulative <- numeric()
  list_V           <- numeric()
  listV_cumulative <- numeric()
  Score_prev       <- 0
  
  for (i in seq_along(Outcome)) {
    Score_V    <- ScoreV(unlist(Outcome[[i]]), Trt)
    Score_pV   <- Score_V * (1 - abs(Score_prev))
    Score_D    <- Score_pV[which(Trt == 1), which(Trt == 0)]
    Score_prev <- Score_prev + Score_pV
    
    list_T[i]           <- sum(Score_D[Score_D > 0])
    list_C[i]           <- sum(Score_D[Score_D < 0])
    list_D[i]           <- mean(Score_D)
    listD_cumulative[i] <- sum(list_D[1:i])
    
    list_V[i]           <- sum(rowSums(Score_pV)^2) /
      (nTest * nControl * nPatients * (nPatients - 1))
    listV_cumulative[i] <- sum(rowSums(Score_prev)^2) /
      (nTest * nControl * nPatients * (nPatients - 1))
  }
  
  pNB     <- listD_cumulative[length(Outcome)]
  pNB_var <- listV_cumulative[length(Outcome)]
  
  p_value <- if (side == 1) {
    pnorm(-pNB / sqrt(pNB_var))
  } else {
    2 * pnorm(-abs(pNB / sqrt(pNB_var)))
  }
  
  return(list(NB = pNB, Var = pNB_var, p_value = p_value))
}

# ==============================================================================
# 3. NPMV  (returns F-statistic and df alongside p-value)
# ==============================================================================

#' Run the NPMV ANOVA-type statistic (ATS)
#'
#' Returns the ATS F-statistic and estimated degrees of freedom in addition to
#' the p-value so the MI pooling step has what it needs for D2 pooling.
#'
#' @return List with p_value, f_stat, df1, df2
npmv <- function(Y_c, Y_t, permreps = 10000) {
  
  # Wrap the entire function body in tryCatch so that ANY failure inside
  # nonpartest() — including cases where it returns an atomic vector instead
  # of a list, causing $ operator errors — returns NA gracefully rather than
  # crashing the simulation.
  tryCatch({
    
    n_c <- nrow(Y_c)
    n_t <- nrow(Y_t)
    N   <- n_c + n_t
    
    df <- data.frame(
      Group = factor(c(rep("P", n_c), rep("V", n_t)), levels = c("V", "P")),
      Y1    = c(Y_c[, 1], Y_t[, 1]),
      Y2    = c(Y_c[, 2], Y_t[, 2]),
      Y3    = c(Y_c[, 3], Y_t[, 3])
    )
    
    # ---------------------------------------------------------------
    # Small-sample switch (npmv documentation, Section 1.4):
    #   N <  10 : full permutation test (all N! permutations)
    #   10 <= N < 30 : randomization test (permreps random permutations)
    #   N >= 30 : F-approximation only
    # When use_perm = TRUE the package chooses between full enumeration
    # and random sampling automatically based on N; we just need to
    # pass permtest = TRUE and a sensible permreps value.
    # ---------------------------------------------------------------
    use_perm <- N < 30
    
    npmv_out <- npmv::nonpartest(
      Y1 | Y2 | Y3 ~ Group,
      data     = df,
      test     = c(1, 0, 0, 0),
      permtest = use_perm,
      permreps = permreps,   # only used when use_perm = TRUE
      plots    = FALSE
    )
    
    res <- npmv_out$results[1, ]
    
    # ---------------------------------------------------------------
    # P-value selection:
    #   Small N  -> permutation / randomization p-value (more reliable)
    #   Large N  -> F-approximation p-value
    #
    # The column names returned by nonpartest() are:
    #   "P-value"                  — F-approximation
    #   "Permutation Test p-value" — permutation / randomization
    #
    # We try the named column first; if the name has changed across
    # package versions we fall back to the last column of the result.
    # ---------------------------------------------------------------
    if (use_perm) {
      perm_col <- grep("permutation|perm|P\\.T\\.",
                       colnames(res), ignore.case = TRUE, value = TRUE)
      p_val <- if (length(perm_col) > 0) {
        res[[perm_col[1]]]
      } else {
        res[[ncol(res)]]   # safe fallback: last column
      }
    } else {
      p_val <- res[["P-value"]]
    }
    
    # f_stat / df1 / df2 always come from the F-approximation row and
    # are returned unchanged — they are still needed for D2 pooling in MI.
    list(
      p_value  = p_val,
      f_stat   = res[["Test Statistic"]],
      df1      = res[["df1"]],
      df2      = res[["df2"]],
      N        = N,           # expose N for diagnostics / callers
      perm_used = use_perm    # flag so callers know which branch fired
    )
    
  }, error = function(e) {
    list(p_value   = NA_real_, f_stat    = NA_real_,
         df1       = NA_real_, df2       = NA_real_,
         N         = NA_real_, perm_used = NA)
  })
}

# ==============================================================================
# 4. ANALYSIS PIPELINE
# ==============================================================================

#' Run the full analysis pipeline for one simulated dataset
#'
#' @param dat          List from apply_dropout() with Y_c, Y_t (with NAs), D_c, D_t
#' @param dat_complete List from generate_sim_data() before dropout (for Benchmark)
#' @param m            Number of multiple imputations
#' @param dist         Data type: "normal" (default), "t", or "binary".
#'                     Controls the imputation method:
#'                       - "normal" and "t" use `norm` (Bayesian linear regression)
#'                       - "binary" uses `logreg` (Bayesian logistic regression)
#'                     Under "t", `norm` is used deliberately to test robustness
#'                     of the standard pipeline to heavy tails.
#' @return List of estimates and p-values for all method x strategy combinations
run_analysis <- function(dat, dat_complete = NULL, m = 40, dist = "normal") {
  
  results <- list()
  
  # ----------------------------------------------------------
  # 0. BENCHMARK (complete data, no missingness)
  # ----------------------------------------------------------
  
  if (!is.null(dat_complete)) {
    df_bench <- wide_to_long(dat_complete$Y_c, dat_complete$Y_t)
    
    res <- gpc_nonprioritized(df_bench)
    results$GPC_NP_Benchmark_est  <- res$NB
    results$GPC_NP_Benchmark_pval <- res$p_value
    
    res <- gpc_prioritized(df_bench)
    results$GPC_P_Benchmark_est  <- res$NB
    results$GPC_P_Benchmark_pval <- res$p_value
    
    res <- npmv(Y_c = dat_complete$Y_c, Y_t = dat_complete$Y_t)
    results$NPMV_Benchmark_pval <- res$p_value
    
  } else {
    results$GPC_NP_Benchmark_est  <- NA
    results$GPC_NP_Benchmark_pval <- NA
    results$GPC_P_Benchmark_est   <- NA
    results$GPC_P_Benchmark_pval  <- NA
    results$NPMV_Benchmark_pval   <- NA
  }
  
  # ----------------------------------------------------------
  # 1. COMPLETE CASE (CC)
  # ----------------------------------------------------------
  
  cc_idx_t <- which(dat$D_t == 4)
  cc_idx_c <- which(dat$D_c == 4)
  
  if (length(cc_idx_t) >= 3 && length(cc_idx_c) >= 3) {
    Y_c_cc <- dat$Y_c[cc_idx_c, , drop = FALSE]
    Y_t_cc <- dat$Y_t[cc_idx_t, , drop = FALSE]
    df_cc  <- wide_to_long(Y_c = Y_c_cc, Y_t = Y_t_cc)
    
    res <- gpc_nonprioritized(df_cc)
    results$GPC_NP_CC_est  <- res$NB
    results$GPC_NP_CC_pval <- res$p_value
    
    res <- gpc_prioritized(df_cc)
    results$GPC_P_CC_est  <- res$NB
    results$GPC_P_CC_pval <- res$p_value
    
    res <- npmv(Y_c = Y_c_cc, Y_t = Y_t_cc)
    results$NPMV_CC_pval <- res$p_value
    
  } else {
    results$GPC_NP_CC_est  <- NA
    results$GPC_NP_CC_pval <- NA
    results$GPC_P_CC_est   <- NA
    results$GPC_P_CC_pval  <- NA
    results$NPMV_CC_pval   <- NA
  }
  
  # ----------------------------------------------------------
  # 2. AVAILABLE CASE (AC)
  # ----------------------------------------------------------
  
  df_ac <- wide_to_long(dat$Y_c, dat$Y_t)
  
  res <- gpc_nonprioritized(df_ac)
  results$GPC_NP_AC_est  <- res$NB
  results$GPC_NP_AC_pval <- res$p_value
  
  res <- gpc_prioritized(df_ac)
  results$GPC_P_AC_est  <- res$NB
  results$GPC_P_AC_pval <- res$p_value
  
  res <- tryCatch(
    npmv(dat$Y_c, dat$Y_t),
    error = function(e) list(p_value = NA_real_)
  )
  results$NPMV_AC_pval <- res$p_value
  
  # ----------------------------------------------------------
  # 3. MULTIPLE IMPUTATION (MI)
  # ----------------------------------------------------------
  
  df_miss <- rbind(
    data.frame(grp = 0, y1 = dat$Y_c[, 1], y2 = dat$Y_c[, 2], y3 = dat$Y_c[, 3]),
    data.frame(grp = 1, y1 = dat$Y_t[, 1], y2 = dat$Y_t[, 2], y3 = dat$Y_t[, 3])
  )
  
  # ---- Imputation method depends on data type ----
  #   normal / t  -> `norm` (Bayesian linear regression). Under `t` this is
  #                  a deliberate misspecification: we test whether the
  #                  standard pipeline is robust to heavy tails.
  #   binary      -> `logreg` (Bayesian logistic regression). Requires y1,y2,y3
  #                  to be factors, so we coerce them before passing to mice()
  #                  and convert back to integers after complete().
  if (dist == "binary") {
    imp_method <- "logreg"
    df_miss$y1 <- factor(df_miss$y1, levels = c(0, 1))
    df_miss$y2 <- factor(df_miss$y2, levels = c(0, 1))
    df_miss$y3 <- factor(df_miss$y3, levels = c(0, 1))
  } else {
    imp_method <- "norm"
  }
  
  pred <- make.predictorMatrix(df_miss)
  pred[,]                    <- 0
  pred["y2", "y1"]           <- 1
  pred["y3", c("y1", "y2")]  <- 1
  pred["y2", "grp"]          <- 1
  pred["y3", "grp"]          <- 1
  
  meth        <- make.method(df_miss)
  meth["grp"] <- ""
  meth["y1"]  <- ""
  meth["y2"]  <- imp_method
  meth["y3"]  <- imp_method
  
  imp <- tryCatch(
    mice(df_miss, m = m, method = meth, predictorMatrix = pred,
         maxit = 5, printFlag = FALSE),
    error = function(e) NULL
  )
  
  if (!is.null(imp)) {
    
    # Storage vectors — one entry per imputation
    MI_NP_est     <- numeric(m)
    MI_NP_var     <- numeric(m)
    MI_P_est      <- numeric(m)
    MI_P_var      <- numeric(m)
    MI_NPMV_fstat <- numeric(m)
    MI_NPMV_df1   <- numeric(m)
    
    for (i in 1:m) {
      comp <- complete(imp, i)
      
      # Binary: mice returns factors, convert back to 0/1 integers
      if (dist == "binary") {
        comp$y1 <- as.integer(as.character(comp$y1))
        comp$y2 <- as.integer(as.character(comp$y2))
        comp$y3 <- as.integer(as.character(comp$y3))
      }
      
      Y_c_imp <- as.matrix(comp[comp$grp == 0, c("y1", "y2", "y3")])
      Y_t_imp <- as.matrix(comp[comp$grp == 1, c("y1", "y2", "y3")])
      colnames(Y_c_imp) <- colnames(Y_t_imp) <- c("Y1", "Y2", "Y3")
      
      df_long <- wide_to_long(Y_c = Y_c_imp, Y_t = Y_t_imp)
      
      # GPC Non-Prioritized
      out <- gpc_nonprioritized(df_long)
      MI_NP_est[i] <- out$NB
      MI_NP_var[i] <- out$Var
      
      # GPC Prioritized
      out <- gpc_prioritized(df_long)
      MI_P_est[i] <- out$NB
      MI_P_var[i] <- out$Var
      
      # NPMV — collect F-statistic and df1 for D2 pooling
      out <- npmv(Y_c = Y_c_imp, Y_t = Y_t_imp)
      MI_NPMV_fstat[i] <- out$f_stat
      MI_NPMV_df1[i]   <- out$df1
    }
    
    # ---- Rubin's Rules: GPC Non-Prioritized ----
    Q_bar  <- mean(MI_NP_est, na.rm = TRUE)
    W      <- mean(MI_NP_var, na.rm = TRUE)
    B      <- var(MI_NP_est,  na.rm = TRUE)
    if (is.na(W) || W <= 0) W <- B
    Var_mi <- W + (1 + 1/m) * B
    
    if (!is.na(Var_mi) && Var_mi > 0) {
      t_stat <- Q_bar / sqrt(Var_mi)
      r_m    <- ((1 + 1/m) * B) / W
      if (m - 1 > 4) {
        nu <- 4 + ((m - 1 - 4) * (1 - 2/(m - 1))) / r_m^2
      } else {
        nu <- (m - 1) * (1 + 1/r_m)^2
      }
      results$GPC_NP_MI_est  <- Q_bar
      results$GPC_NP_MI_pval <- 2 * pt(-abs(t_stat), df = nu)
    } else {
      results$GPC_NP_MI_est  <- NA
      results$GPC_NP_MI_pval <- NA
    }
    
    # ---- Rubin's Rules: GPC Prioritized ----
    Q_bar  <- mean(MI_P_est, na.rm = TRUE)
    W      <- mean(MI_P_var, na.rm = TRUE)
    B      <- var(MI_P_est,  na.rm = TRUE)
    if (is.na(W) || W <= 0) W <- B
    Var_mi <- W + (1 + 1/m) * B
    
    if (!is.na(Var_mi) && Var_mi > 0) {
      t_stat <- Q_bar / sqrt(Var_mi)
      r_m    <- ((1 + 1/m) * B) / W
      if (m - 1 > 4) {
        nu <- 4 + ((m - 1 - 4) * (1 - 2/(m - 1))) / r_m^2
      } else {
        nu <- (m - 1) * (1 + 1/r_m)^2
      }
      results$GPC_P_MI_est  <- Q_bar
      results$GPC_P_MI_pval <- 2 * pt(-abs(t_stat), df = nu)
    } else {
      results$GPC_P_MI_est  <- NA
      results$GPC_P_MI_pval <- NA
    }
    
    # ----------------------------------------------------------
    # NPMV MI POOLING — D2 (Enders 2012 / Li et al. 1991)
    # Operates on F-statistics; no p-value pooling needed.
    # ----------------------------------------------------------
    #
    #   K   = effective numerator df of the multivariate ANOVA-type statistic.
    #         In the univariate case (Grund et al.), K = I - 1 is a fixed
    #         design constant. Here the NPMV ANOVA-type statistic is
    #         multivariate with p = 3 outcomes, so its numerator df
    #         f1 = tr(G)^2 / tr(G^2) is data-dependent and reflects both
    #         the group structure and the covariance structure across all
    #         three outcomes. We therefore use mean(f1) across imputations
    #         as the K-equivalent, which is the closest analogue to the
    #         univariate K for the purpose of the chi-squared scaling
    #         W_m = K * F_m.
    #
    #   ARIV2 = ((M+1) / (M*(M-1))) * sum( (sqrt(W_m) - mean(sqrt(W_m)))^2 )
    #         Centering must use mean(sqrt(W_m)), i.e. the average of the
    #         square roots, NOT sqrt(mean(W_m)) (Equation 5, Grund et al. 2016).
    #
    #   v2  = K^(-3/M) * (M-1) * (1 + ARIV2^-1)^2   (Equation 6)
    
    valid <- is.finite(MI_NPMV_fstat) & is.finite(MI_NPMV_df1) &
      MI_NPMV_fstat > 0 & MI_NPMV_df1 > 0
    
    if (sum(valid) >= 2) {
      K     <- mean(MI_NPMV_df1[valid])   # effective numerator df
      W     <- K * MI_NPMV_fstat[valid]   # chi-squared-scale statistics
      m_v   <- sum(valid)
      W_bar <- mean(W)
      
      # ARIV2: centre on mean(sqrt(W_m)), not sqrt(mean(W_m))  [Eq. 5]
      sqrt_W_bar <- mean(sqrt(W))
      ARIV2 <- ((m_v + 1) / (m_v * (m_v - 1))) *
        sum((sqrt(W) - sqrt_W_bar)^2)
      
      # D2 statistic  [Eq. 4]
      D2 <- (W_bar / K - (m_v + 1) / (m_v - 1) * ARIV2) / (1 + ARIV2)
      
      if (!is.finite(D2) || D2 < 0 || !is.finite(ARIV2) || ARIV2 <= 0) {
        results$NPMV_MI_D2_pval <- NA
      } else {
        # v2 denominator df  [Eq. 6]: K^(-3/M) * (M-1) * (1 + ARIV2^-1)^2
        v2 <- K^(-3 / m_v) * (m_v - 1) * (1 + 1 / ARIV2)^2
        results$NPMV_MI_D2_pval <- pf(D2, df1 = K, df2 = v2,
                                      lower.tail = FALSE)
      }
    } else {
      results$NPMV_MI_D2_pval <- NA
    }
    
  } else {
    # Imputation failed entirely
    results$GPC_NP_MI_est   <- NA
    results$GPC_NP_MI_pval  <- NA
    results$GPC_P_MI_est    <- NA
    results$GPC_P_MI_pval   <- NA
    results$NPMV_MI_D2_pval <- NA
  }
  
  return(results)
}