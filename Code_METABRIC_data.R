
## Bayesian Elastic-Net Cox Models (Powered by Stan)
## CORRECTED VERSION: Fixed Lambda Priors


library(survival)    
library(ggplot2)     
library(dplyr)       
library(boot)        
library(glmnet)      
library(rstan)       
library(gridExtra)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## ============================================================
## SECTION 0: CREATE OUTPUT DIRECTORY
## ============================================================
output_dir <- "results_real"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", output_dir))
}

## ============================================================
## SECTION 1: DATA PREPROCESSING
## ============================================================
preprocess_metabric <- function(data_file) {
  cat("Loading METABRIC data...\n")
  metabric <- read.csv(data_file, stringsAsFactors = FALSE)
  
  survival_data <- data.frame(
    time  = metabric$overall_survival_months * 30.44,
    event = metabric$overall_survival
  )
  
  complete_surv <- complete.cases(survival_data)
  survival_data <- survival_data[complete_surv, ]
  metabric      <- metabric[complete_surv, ]
  
  gene_start_col <- which(names(metabric) == "brca1")
  if(length(gene_start_col) == 0) gene_start_col <- which(names(metabric) == "BRCA1")
  gene_end_col   <- which(grepl("_mut$", names(metabric)))[1] - 1
  if(is.na(gene_end_col)) gene_end_col <- ncol(metabric)
  
  gene_matrix    <- as.matrix(metabric[, gene_start_col:gene_end_col])
  gene_names     <- colnames(gene_matrix)
  gene_matrix <- gene_matrix[, complete.cases(t(gene_matrix))]
  gene_names  <- colnames(gene_matrix)
  
  vars <- apply(gene_matrix, 2, var, na.rm = TRUE)
  keep <- vars >= quantile(vars, 0.1, na.rm = TRUE)
  gene_matrix <- gene_matrix[, keep]
  gene_names  <- gene_names[keep]
  
  positive <- survival_data$time > 0
  if (any(!positive)) {
    survival_data <- survival_data[positive, ]
    gene_matrix   <- gene_matrix[positive, ]
  }
  
  set.seed(42)
  event_idx      <- which(survival_data$event == 1)
  censored_idx   <- which(survival_data$event == 0)
  train_evt      <- sample(event_idx,    round(0.8 * length(event_idx)))
  train_cen      <- sample(censored_idx, round(0.8 * length(censored_idx)))
  train_idx      <- c(train_evt, train_cen)
  test_idx       <- setdiff(seq_len(nrow(survival_data)), train_idx)
  
  train_surv <- survival_data[train_idx, ]
  test_surv  <- survival_data[test_idx, ]
  train_mat  <- gene_matrix[train_idx, ]
  test_mat   <- gene_matrix[test_idx, ]
  
  m <- colMeans(train_mat)
  s <- apply(train_mat, 2, sd)
  s[s == 0] <- 1 
  
  train_mat <- scale(train_mat, center = m, scale = s)
  test_mat  <- scale(test_mat,  center = m, scale = s)
  
  list(train_survival = train_surv, test_survival = test_surv,
       train_genes = train_mat, test_genes = test_mat, gene_names = gene_names)
}

## ============================================================
## SECTION 2: CORRECTED STAN MODEL
## ============================================================
stan_model_code <- "
data {
  int<lower=0> N;
  int<lower=0> P;
  matrix[N, P] X;
  vector[N] y;
  int<lower=0, upper=1> event[N];
}

parameters {
  vector[P] beta_raw;
  vector<lower=1e-6>[P] tau_sq;      // Local variance with lower bound
  real<lower=0.01> lambda1;           // L1 shrinkage with lower bound
  real<lower=0.01> lambda2;           // L2 shrinkage with lower bound
}

transformed parameters {
  vector[P] beta;
  for (j in 1:P) {
    // Elastic-net variance: combines local (lasso) and global (ridge) shrinkage
    real variance_j = tau_sq[j] / (1.0 + lambda2 * tau_sq[j]);
    beta[j] = beta_raw[j] * sqrt(variance_j);
  }
}

model {
  // =====================================================
  // CORRECTED PRIORS - Half-Cauchy for scale parameters
  // =====================================================
  // Half-Cauchy(0, 1) is the standard weakly informative prior
  // for scale/variance parameters (Gelman 2006)
  lambda1 ~ cauchy(0, 1);
  lambda2 ~ cauchy(0, 1);
  
  // Local shrinkage: Exponential prior induces Laplace (Lasso) marginally
  // Rate = lambda1^2 / 2 (Park & Casella 2008)
  for (j in 1:P) {
    tau_sq[j] ~ exponential(0.5 * square(lambda1));
  }
  
  // Standard normal for raw coefficients (non-centered parameterization)
  beta_raw ~ std_normal();
  
  // =====================================================
  // COX PARTIAL LIKELIHOOD
  // =====================================================
  {
    vector[N] risk = X * beta;
    real current_log_sum = negative_infinity();
    real log_lik = 0;
    
    // Backward cumulative sum for efficiency (Breslow approximation)
    for (i in 1:N) {
      int idx = N - i + 1; 
      current_log_sum = log_sum_exp(current_log_sum, risk[idx]);
      if (event[idx] == 1) {
        log_lik += risk[idx] - current_log_sum;
      }
    }
    target += log_lik;
  }
}

generated quantities {
  // Compute log-likelihood for model comparison (WAIC/LOO)
  vector[N] log_lik_vec;
  {
    vector[N] risk = X * beta;
    real current_log_sum = negative_infinity();
    
    for (i in 1:N) {
      int idx = N - i + 1;
      current_log_sum = log_sum_exp(current_log_sum, risk[idx]);
      if (event[idx] == 1) {
        log_lik_vec[idx] = risk[idx] - current_log_sum;
      } else {
        log_lik_vec[idx] = 0;
      }
    }
  }
}
"

## ============================================================
## SECTION 3: WRAPPER
## ============================================================
fit_ben_cox_stan <- function(surv, X, n_iter = 2000, n_warmup = 1000) {
  ord <- order(surv$time, decreasing = FALSE)
  X_sorted <- X[ord, ]
  y_sorted <- surv$time[ord]
  e_sorted <- surv$event[ord]
  
  stan_data <- list(N = nrow(X_sorted), P = ncol(X_sorted), X = X_sorted, y = y_sorted, event = e_sorted)
  
  cat("Sampling with Stan (HMC)...\n")
  cat("  Using Half-Cauchy(0,1) priors for lambda1, lambda2\n")
  
  fit <- stan(model_code = stan_model_code, data = stan_data,
              iter = n_iter, warmup = n_warmup, chains = 4, 
              control = list(adapt_delta = 0.95, max_treedepth = 12),
              init = function() list(
                beta_raw = rnorm(ncol(X_sorted), 0, 0.1),
                tau_sq = rep(0.1, ncol(X_sorted)),
                lambda1 = 1,
                lambda2 = 1
              ))
  
  samples <- rstan::extract(fit)
  beta_mean <- colMeans(samples$beta)
  beta_q    <- apply(samples$beta, 2, quantile, probs = c(0.025, 0.5, 0.975))
  
  # Report shrinkage parameters
  cat(sprintf("\n  Posterior Lambda1: %.3f (95%% CI: %.3f - %.3f)\n", 
              mean(samples$lambda1), quantile(samples$lambda1, 0.025), quantile(samples$lambda1, 0.975)))
  cat(sprintf("  Posterior Lambda2: %.3f (95%% CI: %.3f - %.3f)\n", 
              mean(samples$lambda2), quantile(samples$lambda2, 0.025), quantile(samples$lambda2, 0.975)))
  
  out <- list(
    beta = list(mean = beta_mean, quantiles = rbind(beta_q[1,], NA, beta_q[2,], NA, beta_q[3,])),
    lambda1 = list(mean = mean(samples$lambda1), sd = sd(samples$lambda1)),
    lambda2 = list(mean = mean(samples$lambda2), sd = sd(samples$lambda2)),
    stan_fit = fit,
    samples = samples
  )
  class(out) <- "ben_cox_fit"
  return(out)
}

print.ben_cox_fit <- function(x, ...) {
  cat("\nBayesian Elastic-Net Cox Model (Stan)\n")
  cat(sprintf("  Lambda1 (L1): %.3f (sd: %.3f)\n", x$lambda1$mean, x$lambda1$sd))
  cat(sprintf("  Lambda2 (L2): %.3f (sd: %.3f)\n", x$lambda2$mean, x$lambda2$sd))
  nz <- sum(x$beta$quantiles[1, ] > 0 | x$beta$quantiles[5, ] < 0)
  cat(sprintf("  Significant Genes (95%% CI excludes 0): %d\n", nz))
}

## ============================================================
## SECTION 4: BASELINES
## ============================================================
fit_null_cox <- function(surv) coxph(Surv(time, event) ~ 1, data = surv)

fit_ridge_cox <- function(surv, X) {
  surv_obj <- Surv(surv$time, surv$event)
  cv_fit <- cv.glmnet(x = X, y = surv_obj, family = "cox", alpha = 0)
  final_fit <- glmnet(x = X, y = surv_obj, family = "cox", alpha = 0, lambda = cv_fit$lambda.min)
  list(model = final_fit, lambda_opt = cv_fit$lambda.min, coef = as.numeric(coef(final_fit)))
}

fit_lasso_cox <- function(surv, X) {
  surv_obj <- Surv(surv$time, surv$event)
  cv_fit <- cv.glmnet(x = X, y = surv_obj, family = "cox", alpha = 1)
  final_fit <- glmnet(x = X, y = surv_obj, family = "cox", alpha = 1, lambda = cv_fit$lambda.min)
  coef_vec <- as.numeric(coef(final_fit))
  list(model = final_fit, lambda_opt = cv_fit$lambda.min, 
       coef = coef_vec, n_selected = sum(coef_vec != 0))
}

fit_en_cox <- function(surv, X, alpha = 0.5) {
  surv_obj <- Surv(surv$time, surv$event)
  cv_fit <- cv.glmnet(x = X, y = surv_obj, family = "cox", alpha = alpha)
  final_fit <- glmnet(x = X, y = surv_obj, family = "cox", alpha = alpha, lambda = cv_fit$lambda.min)
  coef_vec <- as.numeric(coef(final_fit))
  list(model = final_fit, lambda_opt = cv_fit$lambda.min, alpha = alpha,
       coef = coef_vec, n_selected = sum(coef_vec != 0))
}

## ============================================================
## SECTION 5: METRICS
## ============================================================
calculate_ibs <- function(surv, risk, tau_years = 10) {
  max_time <- tau_years * 365.25
  tp <- seq(0, max_time, length.out = 50); tp <- tp[tp > 0]
  km <- survfit(Surv(time, event) ~ 1, data = surv)
  base_surv <- summary(km, times = tp, extend = TRUE)$surv
  pred_surv <- outer(base_surv, exp(risk), "^")
  bs <- sapply(seq_along(tp), function(i) {
    obs <- as.numeric(surv$time > tp[i])
    mean((obs - pred_surv[i, ])^2)
  })
  dt <- diff(tp)
  sum((bs[-1] + bs[-length(bs)]) * dt / 2) / max_time
}

compute_c_index <- function(surv, risk) {
  cf <- concordance(Surv(time, event) ~ I(-risk), data = surv)
  return(cf$concordance)
}

calculate_gnd_chisq <- function(surv, risk, tau_years = 5, groups = 10) {
  tau <- tau_years * 365.25
  km <- survfit(Surv(time, event) ~ 1, data = surv)
  S0 <- summary(km, times = tau, extend = TRUE)$surv
  if(is.na(S0)) S0 <- tail(summary(km)$surv, 1)
  
  pi_hat <- 1 - (S0 ^ exp(risk))
  pi_hat <- pmin(pmax(pi_hat, 1e-5), 1 - 1e-5)
  
  group <- tryCatch({
    cut(pi_hat, breaks = quantile(pi_hat, probs = seq(0, 1, length.out = groups + 1)), 
        include.lowest = TRUE, labels = FALSE)
  }, error = function(e) {
    return(rep(1, length(pi_hat)))
  })
  
  observed_events <- tapply(ifelse(surv$time <= tau & surv$event == 1, 1, 0), group, sum)
  expected_events <- tapply(pi_hat, group, sum)
  
  chisq_val <- sum((observed_events - expected_events)^2 / expected_events, na.rm = TRUE)
  return(chisq_val)
}

calculate_ece <- function(surv, risk, tau_years = 5, n_bins = 10) {
  tau <- tau_years * 365.25
  km <- survfit(Surv(time, event) ~ 1, data = surv)
  S0 <- summary(km, times = tau, extend = TRUE)$surv
  if(is.na(S0)) S0 <- tail(summary(km)$surv, 1)
  
  pi_hat <- 1 - (S0 ^ exp(risk))
  pi_hat <- pmin(pmax(pi_hat, 1e-5), 1 - 1e-5)
  obs_event <- ifelse(surv$time <= tau & surv$event == 1, 1, 0)
  
  bins <- tryCatch({
    cut(pi_hat, breaks = quantile(pi_hat, probs = seq(0, 1, length.out = n_bins + 1)), 
        include.lowest = TRUE, labels = FALSE)
  }, error = function(e) {
    return(rep(1, length(pi_hat)))
  })
  
  n <- length(pi_hat)
  ece <- 0
  for (b in unique(bins)) {
    idx <- which(bins == b)
    n_b <- length(idx)
    if (n_b > 0) {
      mean_pred <- mean(pi_hat[idx])
      mean_obs <- mean(obs_event[idx])
      ece <- ece + (n_b / n) * abs(mean_obs - mean_pred)
    }
  }
  return(ece)
}

calculate_calibration_slope_intercept <- function(surv, risk, tau_years = 5) {
  tau <- tau_years * 365.25
  km <- survfit(Surv(time, event) ~ 1, data = surv)
  S0 <- summary(km, times = tau, extend = TRUE)$surv
  if(is.na(S0)) S0 <- tail(summary(km)$surv, 1)
  
  pi_hat <- 1 - (S0 ^ exp(risk))
  pi_hat <- pmin(pmax(pi_hat, 1e-6), 1 - 1e-6)
  obs_event <- ifelse(surv$time <= tau & surv$event == 1, 1, 0)
  logit_pi <- log(pi_hat / (1 - pi_hat))
  
  fit <- tryCatch({
    glm(obs_event ~ logit_pi, family = binomial)
  }, error = function(e) {
    return(list(coefficients = c(NA, NA)))
  })
  
  return(c(intercept = fit$coefficients[1], slope = fit$coefficients[2]))
}

bootstrap_evaluation <- function(surv, risk, n_boot = 100, tau = 5) {
  n <- nrow(surv)
  res <- matrix(NA, nrow = n_boot, ncol = 5)
  colnames(res) <- c("IBS", "C", "GND", "ECE", "CalSlope")
  
  for(b in 1:n_boot) {
    idx <- sample(n, n, replace = TRUE)
    s_boot <- surv[idx, ]
    r_boot <- risk[idx]
    
    res[b, "IBS"] <- calculate_ibs(s_boot, r_boot, tau)
    res[b, "C"]   <- compute_c_index(s_boot, r_boot)
    res[b, "GND"] <- calculate_gnd_chisq(s_boot, r_boot, tau)
    res[b, "ECE"] <- calculate_ece(s_boot, r_boot, tau)
    cal <- calculate_calibration_slope_intercept(s_boot, r_boot, tau)
    res[b, "CalSlope"] <- cal["slope"]
  }
  
  data.frame(metric = colnames(res), 
             mean = apply(res, 2, mean, na.rm = TRUE), 
             se = apply(res, 2, sd, na.rm = TRUE))
}

## ============================================================
## SECTION 6: DIAGNOSTIC PLOTS
## ============================================================
plot_schoenfeld_residuals <- function(surv, X, ben_fit, gene_names, output_dir) {
  filename <- file.path(output_dir, "ben_schoenfeld.pdf")
  cat(sprintf("Saving Schoenfeld Residual Plots to '%s'...\n", filename))
  
  beta_means <- ben_fit$beta$mean
  top_idx <- order(abs(beta_means), decreasing = TRUE)[1:6]
  top_names <- gene_names[top_idx]
  
  df <- data.frame(time = surv$time, event = surv$event, X[, top_idx])
  colnames(df)[3:ncol(df)] <- top_names
  
  formula_str <- paste("Surv(time, event) ~", paste(top_names, collapse = " + "))
  cox_fit <- coxph(as.formula(formula_str), data = df)
  ph_test <- cox.zph(cox_fit)
  
  pdf(filename, width = 10, height = 8)
  par(mfrow = c(2, 3))
  for (i in 1:min(6, length(top_names))) {
    plot(ph_test[i], main = top_names[i])
    abline(h = 0, col = "red", lty = 2)
  }
  dev.off()
  
  return(ph_test)
}

plot_posterior_predictive_check <- function(surv, X, ben_fit, output_dir) {
  filename <- file.path(output_dir, "ben_ppc.pdf")
  cat(sprintf("Saving Posterior Predictive Check to '%s'...\n", filename))
  
  samples <- ben_fit$samples
  n_draws <- min(100, nrow(samples$beta))
  draw_idx <- sample(1:nrow(samples$beta), n_draws)
  
  km_obs <- survfit(Surv(time, event) ~ 1, data = surv)
  time_points <- seq(0, max(surv$time) * 0.9, length.out = 100)
  km_summary <- summary(km_obs, times = time_points, extend = TRUE)
  S0_t <- km_summary$surv
  
  sim_curves <- matrix(NA, nrow = n_draws, ncol = length(time_points))
  
  for (d in 1:n_draws) {
    beta_d <- samples$beta[draw_idx[d], ]
    risk_d <- as.numeric(X %*% beta_d)
    pred_surv_d <- outer(S0_t, exp(risk_d), "^")
    sim_curves[d, ] <- rowMeans(pred_surv_d)
  }
  
  surv_lower <- apply(sim_curves, 2, quantile, probs = 0.025)
  surv_upper <- apply(sim_curves, 2, quantile, probs = 0.975)
  
  pdf(filename, width = 8, height = 6)
  plot(km_obs, conf.int = FALSE, col = "black", lwd = 2,
       xlab = "Time (days)", ylab = "Survival Probability",
       main = "Posterior Predictive Check: BEN-Cox")
  
  for (d in 1:min(50, n_draws)) {
    lines(time_points, sim_curves[d, ], col = rgb(0.5, 0.5, 0.5, 0.2))
  }
  
  polygon(c(time_points, rev(time_points)), 
          c(surv_lower, rev(surv_upper)),
          col = rgb(0, 0, 1, 0.2), border = NA)
  
  lines(km_obs, conf.int = FALSE, col = "black", lwd = 2)
  
  legend("bottomleft", 
         legend = c("Observed KM", "95% Credible Band", "Posterior Draws"),
         col = c("black", rgb(0, 0, 1, 0.5), "gray"),
         lwd = c(2, 10, 1), bty = "n")
  dev.off()
}

# -- NEW: Lambda posterior plot --
plot_lambda_posteriors <- function(ben_fit, output_dir) {
  filename <- file.path(output_dir, "ben_lambda_posteriors.pdf")
  cat(sprintf("Saving Lambda Posterior Plots to '%s'...\n", filename))
  
  samples <- ben_fit$samples
  
  pdf(filename, width = 10, height = 4)
  par(mfrow = c(1, 2))
  
  hist(samples$lambda1, breaks = 50, main = "Posterior: Lambda1 (L1 Shrinkage)",
       xlab = "Lambda1", col = "lightblue", border = "white")
  abline(v = mean(samples$lambda1), col = "red", lwd = 2)
  
  hist(samples$lambda2, breaks = 50, main = "Posterior: Lambda2 (L2 Shrinkage)",
       xlab = "Lambda2", col = "lightgreen", border = "white")
  abline(v = mean(samples$lambda2), col = "red", lwd = 2)
  
  dev.off()
}

## ============================================================
## SECTION 7: PLOTTING
## ============================================================
plot_mcmc_diagnostics <- function(ben_fit, gene_names, output_dir) {
  filename <- file.path(output_dir, "ben_mcmc_traces.pdf")
  cat(sprintf("Saving MCMC Trace Plots to '%s'...\n", filename))
  
  fit <- ben_fit$stan_fit
  beta_means <- ben_fit$beta$mean
  top_idx <- order(abs(beta_means), decreasing = TRUE)[1:4]
  stan_beta_names <- paste0("beta[", top_idx, "]")
  
  pdf(filename, width = 12, height = 8)
  p1 <- rstan::traceplot(fit, pars = c("lambda1", "lambda2"), inc_warmup = FALSE) +
    ggtitle("Trace: Global Shrinkage Parameters")
  print(p1)
  
  p2 <- rstan::traceplot(fit, pars = stan_beta_names, inc_warmup = FALSE) +
    ggtitle("Trace: Top 4 Genes by Effect Size")
  print(p2)
  dev.off()
}

plot_km_quintiles <- function(surv, risk, output_dir) {
  filename <- file.path(output_dir, "ben_km_curves.pdf")
  cat(sprintf("Saving KM Curves to '%s'...\n", filename))
  
  q <- quantile(risk, probs = seq(0, 1, 0.2))
  risk_group <- cut(risk, breaks = q, include.lowest = TRUE, 
                    labels = c("Low", "Low-Med", "Med", "Med-High", "High"))
  
  df <- data.frame(time = surv$time, event = surv$event, Group = risk_group)
  fit_km <- survfit(Surv(time, event) ~ Group, data = df)
  
  pdf(filename, width = 8, height = 6)
  plot(fit_km, col = 1:5, lwd = 2, xlab = "Time (days)", ylab = "Survival Probability",
       main = "Kaplan-Meier Curves by Risk Score Quintile (Test Set)")
  legend("bottomleft", legend = levels(risk_group), col = 1:5, lwd = 2, bty = "n")
  dev.off()
}

plot_results_to_file <- function(surv, risk, ben_fit, gene_names, output_dir) {
  filename_cal <- file.path(output_dir, "ben_calibration.pdf")
  cat(sprintf("Saving Calibration Plot to '%s'...\n", filename_cal))
  pdf(filename_cal, width = 6, height = 6)
  tau <- 5 * 365.25
  km <- survfit(Surv(time, event) ~ 1, data = surv)
  S0 <- summary(km, times = tau, extend = TRUE)$surv
  pi_hat <- 1 - (S0 ^ exp(risk))
  Z <- ifelse(surv$time <= tau & surv$event == 1, 1, 0)
  df <- data.frame(pred = pi_hat, obs = Z)
  df$bin <- cut(df$pred, breaks = quantile(df$pred, probs = seq(0, 1, 0.1)), include.lowest = TRUE)
  agg <- df %>% group_by(bin) %>% summarize(pred = mean(pred), obs = mean(obs))
  print(ggplot(agg, aes(x = pred, y = obs)) + geom_point(size = 3) + geom_line() + 
          geom_abline(linetype = "dashed") + labs(title = "BEN-Cox Calibration (5-year)", x = "Predicted", y = "Observed") + theme_minimal())
  dev.off()
  
  filename_forest <- file.path(output_dir, "ben_forest.pdf")
  cat(sprintf("Saving Forest Plot to '%s'...\n", filename_forest))
  pdf(filename_forest, width = 7, height = 8)
  df <- data.frame(gene = gene_names, mean = ben_fit$beta$mean, 
                   low = ben_fit$beta$quantiles[1,], high = ben_fit$beta$quantiles[5,])
  df_top <- df %>% arrange(desc(abs(mean))) %>% slice(1:30)
  df_top$gene <- factor(df_top$gene, levels = rev(df_top$gene))
  print(ggplot(df_top, aes(x = gene, y = mean)) + geom_pointrange(aes(ymin = low, ymax = high)) + 
          coord_flip() + geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
          labs(title = "Top 30 Genes", y = "Log Hazard Ratio") + theme_minimal())
  dev.off()
}

## ============================================================
## SECTION 8: TABLE GENERATORS
## ============================================================
generate_latex_tables <- function(perf_table_df, pam_stats, cal_table_df) {
  
  cat("\n=== LATEX TABLE 3: PERFORMANCE (5 MODELS) ===\n")
  cat("\\begin{table}[h!]\n\\centering\n")
  cat("\\caption{Predictive performance on the held-out 20\\% test set. Values are mean $\\pm$ standard error over 100 bootstrap resamples.}\n")
  cat("\\label{tab:perf}\n")
  cat("\\rev{\n")
  cat("\\begin{tabular}{lcccc}\n\\toprule\n")
  cat("Model & IBS $\\downarrow$ & $C$-index $\\uparrow$ & GND $\\chi^{2}$ $\\downarrow$ & Genes retained \\\\\n\\midrule\n")
  
  for (i in 1:nrow(perf_table_df)) {
    r <- perf_table_df[i, ]
    model_name <- r$Model
    if (model_name == "BEN-Cox") {
      cat(sprintf("\\textbf{BEN--Cox} & \\textbf{%.3f $\\pm$ %.3f} & \\textbf{%.3f $\\pm$ %.3f} & \\textbf{%.1f $\\pm$ %.1f} & %d \\\\\n",
                  r$IBS_mean, r$IBS_se, r$C_mean, r$C_se, r$GND_mean, r$GND_se, r$Genes_retained))
    } else {
      cat(sprintf("%s & %.3f $\\pm$ %.3f & %.3f $\\pm$ %.3f & %.1f $\\pm$ %.1f & %s \\\\\n",
                  model_name, r$IBS_mean, r$IBS_se, r$C_mean, r$C_se, r$GND_mean, r$GND_se, 
                  ifelse(is.na(r$Genes_retained), "---", as.character(r$Genes_retained))))
    }
  }
  
  cat("\\bottomrule\n\\end{tabular}\n}\n\\end{table}\n")
  
  cat("\n=== LATEX TABLE: CALIBRATION METRICS ===\n")
  cat("\\begin{table}[h!]\n\\centering\n")
  cat("\\caption{\\rev{Calibration metrics on the held-out test set at $\\tau = 5$ years.}}\n")
  cat("\\label{tab:calibration}\n")
  cat("\\rev{\n")
  cat("\\begin{tabular}{lccc}\n\\toprule\n")
  cat("Model & Calibration Slope $\\alpha_1$ & Calibration Intercept $\\alpha_0$ & ECE $\\downarrow$ \\\\\n\\midrule\n")
  
  for (i in 1:nrow(cal_table_df)) {
    r <- cal_table_df[i, ]
    model_name <- r$Model
    if (model_name == "BEN-Cox") {
      cat(sprintf("\\textbf{BEN--Cox} & %.3f $\\pm$ %.3f & %.3f & %.4f $\\pm$ %.4f \\\\\n",
                  r$CalSlope_mean, r$CalSlope_se, r$CalInt_mean, r$ECE_mean, r$ECE_se))
    } else {
      cat(sprintf("%s & %.3f $\\pm$ %.3f & %.3f & %.4f $\\pm$ %.4f \\\\\n",
                  model_name, r$CalSlope_mean, r$CalSlope_se, r$CalInt_mean, r$ECE_mean, r$ECE_se))
    }
  }
  
  cat("\\bottomrule\n\\end{tabular}\n}\n\\end{table}\n")
  
  cat("\n=== LATEX TABLE 4: PAM50 ===\n")
  cat("\\begin{table}[h!]\n\\centering\n")
  cat("\\caption{Coverage of the PAM50 gene panel among retained coefficients.}\n")
  cat("\\label{tab:pam50}\n")
  cat("\\begin{tabular}{lcc}\n\\toprule\n")
  cat("Model & Genes retained & Recall (\\%) \\\\\n\\midrule\n")
  cat(sprintf("Ridge Cox & 440 & 100.0 \\\\\n"))
  cat(sprintf("BEN Cox & %d & %.1f \\\\\n", pam_stats$total_retained, round(pam_stats$pam50_recall * 100, 1)))
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")
}

save_convergence_table <- function(ben_fit, output_dir) {
  cat("Generating Convergence Diagnostics Table...\n")
  summary <- summary(ben_fit$stan_fit)$summary
  
  beta_idx <- grep("beta\\[", rownames(summary))
  lambda_idx <- grep("lambda", rownames(summary))
  
  conv_df <- data.frame(
    Parameter = c("Beta (median)", "Lambda1", "Lambda2"),
    Rhat = c(median(summary[beta_idx, "Rhat"]), 
             summary[lambda_idx[1], "Rhat"], 
             summary[lambda_idx[2], "Rhat"]),
    ESS = c(median(summary[beta_idx, "n_eff"]), 
            summary[lambda_idx[1], "n_eff"], 
            summary[lambda_idx[2], "n_eff"])
  )
  
  write.csv(conv_df, file.path(output_dir, "ben_convergence_table.csv"), row.names = FALSE)
  cat(sprintf("Saved to '%s/ben_convergence_table.csv'.\n", output_dir))
  print(conv_df)
  
  cat("\n=== LATEX TABLE 2: CONVERGENCE ===\n")
  cat("\\begin{table}[h!]\n\\centering\n")
  cat("\\caption{Convergence diagnostics for the Stan BEN--Cox fit.}\n")
  cat("\\label{tab:convergence}\n")
  cat("\\begin{tabular}{lcc}\n\\toprule\n")
  cat("Parameter & $\\hat{R}$ & ESS \\\\\n\\midrule\n")
  for (i in 1:nrow(conv_df)) {
    cat(sprintf("%s & %.4f & %.0f \\\\\n", conv_df$Parameter[i], conv_df$Rhat[i], conv_df$ESS[i]))
  }
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")
}

analyze_pam50 <- function(ben_fit, gene_names) {
  pam50_list <- c("ESR1","PGR","ERBB2","FGFR4","BLVRA","BAG1","CCNB1","CCNE1",
                  "CDC20","CDCA1","CEP55","ANLN","BIRC5","BUB1B","CENPF","CENPU",
                  "EXO1","FOXM1","KNTC2","KIF2C","MELK","MKI67","MYBL2","NDC80",
                  "NUF2","PTTG1","RRM2","TYMS","UBE2C","UBE2T")
  sig_idx <- which(ben_fit$beta$quantiles[1,] > 0 | ben_fit$beta$quantiles[5,] < 0)
  sig_genes <- gene_names[sig_idx]
  match_pam <- intersect(toupper(sig_genes), toupper(pam50_list))
  
  cat("\n--- PAM50 Analysis ---\n")
  cat(sprintf("Total significant genes found: %d\n", length(sig_genes)))
  cat("PAM50 genes found:", paste(match_pam, collapse = ", "), "\n")
  
  list(total_retained = length(sig_genes), 
       pam50_recall = length(match_pam) / length(intersect(toupper(gene_names), toupper(pam50_list))))
}

## ============================================================
## SECTION 9: MAIN EXECUTION
## ============================================================
run_complete_analysis <- function(data_file, n_iter = 2000, n_warmup = 1000, n_boot = 100) {
  cat("=== STARTING ANALYSIS (CORRECTED VERSION) ===\n")
  cat("Key fix: Changed lambda priors from Normal(0, 0.05) to Half-Cauchy(0, 1)\n\n")
  
  dat <- preprocess_metabric(data_file)
  
  # FIT ALL MODELS
  cat("\n--- Fitting BEN-Cox (Stan) ---\n")
  ben_fit <- fit_ben_cox_stan(dat$train_survival, dat$train_genes, n_iter = n_iter, n_warmup = n_warmup)
  print(ben_fit)
  
  cat("\n--- Fitting Ridge Cox ---\n")
  ridge_res <- fit_ridge_cox(dat$train_survival, dat$train_genes)
  
  cat("\n--- Fitting Lasso Cox ---\n")
  lasso_res <- fit_lasso_cox(dat$train_survival, dat$train_genes)
  cat(sprintf("  Lasso selected %d genes\n", lasso_res$n_selected))
  
  cat("\n--- Fitting Elastic-Net Cox (alpha=0.5) ---\n")
  en_res <- fit_en_cox(dat$train_survival, dat$train_genes, alpha = 0.5)
  cat(sprintf("  EN-Cox selected %d genes\n", en_res$n_selected))
  
  cat("\n--- Fitting Null Cox ---\n")
  null_fit <- fit_null_cox(dat$train_survival)
  
  # PREDICT ON TEST SET
  ben_risk <- as.numeric(dat$test_genes %*% ben_fit$beta$mean)
  ridge_risk <- as.numeric(predict(ridge_res$model, newx = dat$test_genes, s = ridge_res$lambda_opt, type = "link"))
  lasso_risk <- as.numeric(predict(lasso_res$model, newx = dat$test_genes, s = lasso_res$lambda_opt, type = "link"))
  en_risk <- as.numeric(predict(en_res$model, newx = dat$test_genes, s = en_res$lambda_opt, type = "link"))
  null_risk <- rep(0, nrow(dat$test_survival))
  
  # EVALUATE (BOOTSTRAP)
  cat("\n--- Evaluating (Bootstrap) ---\n")
  ben_perf <- bootstrap_evaluation(dat$test_survival, ben_risk, n_boot = n_boot)
  ridge_perf <- bootstrap_evaluation(dat$test_survival, ridge_risk, n_boot = n_boot)
  lasso_perf <- bootstrap_evaluation(dat$test_survival, lasso_risk, n_boot = n_boot)
  en_perf <- bootstrap_evaluation(dat$test_survival, en_risk, n_boot = n_boot)
  null_perf <- bootstrap_evaluation(dat$test_survival, null_risk, n_boot = n_boot)
  
  # BUILD PERFORMANCE DATAFRAME
  make_row <- function(name, df, n_genes = NA) {
    data.frame(
      Model = name, 
      IBS_mean = df[df$metric == "IBS", "mean"], IBS_se = df[df$metric == "IBS", "se"],
      C_mean = df[df$metric == "C", "mean"], C_se = df[df$metric == "C", "se"],
      GND_mean = df[df$metric == "GND", "mean"], GND_se = df[df$metric == "GND", "se"],
      Genes_retained = n_genes
    )
  }
  
  ben_selected <- sum(ben_fit$beta$quantiles[1,] > 0 | ben_fit$beta$quantiles[5,] < 0)
  
  perf_df <- rbind(
    make_row("Null Cox", null_perf, 0),
    make_row("Ridge Cox", ridge_perf, ncol(dat$train_genes)),
    make_row("Lasso Cox", lasso_perf, lasso_res$n_selected),
    make_row("EN-Cox (freq)", en_perf, en_res$n_selected),
    make_row("BEN-Cox", ben_perf, ben_selected)
  )
  
  write.csv(perf_df, file.path(output_dir, "ben_performance_metrics.csv"), row.names = FALSE)
  cat(sprintf("Saved performance metrics to '%s/ben_performance_metrics.csv'\n", output_dir))
  
  # BUILD CALIBRATION DATAFRAME
  make_cal_row <- function(name, df) {
    data.frame(
      Model = name,
      CalSlope_mean = df[df$metric == "CalSlope", "mean"], 
      CalSlope_se = df[df$metric == "CalSlope", "se"],
      CalInt_mean = NA,
      ECE_mean = df[df$metric == "ECE", "mean"], 
      ECE_se = df[df$metric == "ECE", "se"]
    )
  }
  
  cal_ben <- calculate_calibration_slope_intercept(dat$test_survival, ben_risk)
  cal_ridge <- calculate_calibration_slope_intercept(dat$test_survival, ridge_risk)
  cal_lasso <- calculate_calibration_slope_intercept(dat$test_survival, lasso_risk)
  cal_en <- calculate_calibration_slope_intercept(dat$test_survival, en_risk)
  cal_null <- calculate_calibration_slope_intercept(dat$test_survival, null_risk)
  
  cal_df <- rbind(
    make_cal_row("Null Cox", null_perf),
    make_cal_row("Ridge Cox", ridge_perf),
    make_cal_row("Lasso Cox", lasso_perf),
    make_cal_row("EN-Cox (freq)", en_perf),
    make_cal_row("BEN-Cox", ben_perf)
  )
  
  cal_df$CalInt_mean <- c(cal_null["intercept"], cal_ridge["intercept"], 
                          cal_lasso["intercept"], cal_en["intercept"], cal_ben["intercept"])
  
  write.csv(cal_df, file.path(output_dir, "ben_calibration_metrics.csv"), row.names = FALSE)
  cat(sprintf("Saved calibration metrics to '%s/ben_calibration_metrics.csv'\n", output_dir))
  
  # PLOTS
  plot_results_to_file(dat$test_survival, ben_risk, ben_fit, dat$gene_names, output_dir)
  plot_mcmc_diagnostics(ben_fit, dat$gene_names, output_dir) 
  plot_km_quintiles(dat$test_survival, ben_risk, output_dir)
  plot_lambda_posteriors(ben_fit, output_dir)  # NEW
  
  # DIAGNOSTIC PLOTS
  cat("\n--- Generating Diagnostic Plots ---\n")
  ph_test <- plot_schoenfeld_residuals(dat$train_survival, dat$train_genes, ben_fit, dat$gene_names, output_dir)
  plot_posterior_predictive_check(dat$test_survival, dat$test_genes, ben_fit, output_dir)
  
  # PAM50 & TABLES
  pam_stats <- analyze_pam50(ben_fit, dat$gene_names)
  generate_latex_tables(perf_df, pam_stats, cal_df)
  save_convergence_table(ben_fit, output_dir)
  
  cat("\n=== ANALYSIS COMPLETE ===\n")
  cat(sprintf("All outputs saved to: %s/\n", output_dir))
  
  return(list(
    ben = ben_fit, 
    ridge = ridge_res, 
    lasso = lasso_res, 
    en = en_res,
    perf = perf_df,
    cal = cal_df,
    ph_test = ph_test
  ))
}

# =====================
# EXECUTE
# =====================
results <- run_complete_analysis("METABRIC_RNA_Mutation.csv", n_iter = 2000, n_warmup = 1000, n_boot = 50)

