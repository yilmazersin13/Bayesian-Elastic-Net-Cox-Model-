
## Simulation Study for BEN-Cox Paper (CORRECTED)
## Fixed: Lambda priors changed to Half-Cauchy(0,1)


library(survival)
library(glmnet)
library(rstan)
library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## ============================================================
## SECTION 1: CORRECTED STAN MODEL
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
  vector<lower=1e-6>[P] tau_sq;
  real<lower=0.01> lambda1;
  real<lower=0.01> lambda2;
}

transformed parameters {
  vector[P] beta;
  for (j in 1:P) {
    real variance_j = tau_sq[j] / (1.0 + lambda2 * tau_sq[j]);
    beta[j] = beta_raw[j] * sqrt(variance_j);
  }
}

model {
  // CORRECTED PRIORS - Half-Cauchy for scale parameters
  lambda1 ~ cauchy(0, 1);
  lambda2 ~ cauchy(0, 1);
  
  // Local shrinkage
  for (j in 1:P) {
    tau_sq[j] ~ exponential(0.5 * square(lambda1));
  }
  
  beta_raw ~ std_normal();
  
  // Cox partial likelihood
  {
    vector[N] risk = X * beta;
    real current_log_sum = negative_infinity();
    real log_lik = 0;
    
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
"

cat("Compiling Stan model (with corrected priors)...\n")
stan_model_compiled <- stan_model(model_code = stan_model_code)

## ============================================================
## SECTION 2: DATA GENERATION
## ============================================================
generate_covariance_matrix <- function(p, type = "independent", rho = 0.5, block_size = 10) {
  if (type == "independent") {
    return(diag(p))
  } else if (type == "toeplitz") {
    Sigma <- matrix(0, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        Sigma[i, j] <- rho^abs(i - j)
      }
    }
    return(Sigma)
  } else if (type == "block") {
    Sigma <- diag(p)
    n_blocks <- ceiling(p / block_size)
    for (b in 1:n_blocks) {
      start_idx <- (b - 1) * block_size + 1
      end_idx <- min(b * block_size, p)
      for (i in start_idx:end_idx) {
        for (j in start_idx:end_idx) {
          if (i != j) {
            Sigma[i, j] <- rho
          }
        }
      }
    }
    return(Sigma)
  }
}

generate_survival_data <- function(n, p, s0 = 10, beta_val = 0.5, cov_type = "independent", rho = 0.5) {
  Sigma <- generate_covariance_matrix(p, type = cov_type, rho = rho)
  X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  X <- scale(X)
  
  beta_true <- rep(0, p)
  beta_true[1:s0] <- beta_val
  
  linear_pred <- as.numeric(X %*% beta_true)
  event_times <- rexp(n, rate = exp(linear_pred))
  event_indicator <- rep(1, n)
  
  list(
    X = X,
    time = event_times,
    event = event_indicator,
    beta_true = beta_true,
    Sigma = Sigma
  )
}

## ============================================================
## SECTION 3: MODEL FITTING FUNCTIONS
## ============================================================
fit_null_cox_sim <- function(surv_df) {
  list(risk = rep(0, nrow(surv_df)), coef = NULL, n_selected = 0)
}

fit_ridge_cox_sim <- function(surv_df, X) {
  surv_obj <- Surv(surv_df$time, surv_df$event)
  cv_fit <- cv.glmnet(x = X, y = surv_obj, family = "cox", alpha = 0)
  final_fit <- glmnet(x = X, y = surv_obj, family = "cox", alpha = 0, lambda = cv_fit$lambda.min)
  coef_vec <- as.numeric(coef(final_fit))
  list(model = final_fit, coef = coef_vec, lambda = cv_fit$lambda.min, n_selected = sum(coef_vec != 0))
}

fit_lasso_cox_sim <- function(surv_df, X) {
  surv_obj <- Surv(surv_df$time, surv_df$event)
  cv_fit <- cv.glmnet(x = X, y = surv_obj, family = "cox", alpha = 1)
  final_fit <- glmnet(x = X, y = surv_obj, family = "cox", alpha = 1, lambda = cv_fit$lambda.min)
  coef_vec <- as.numeric(coef(final_fit))
  list(model = final_fit, coef = coef_vec, lambda = cv_fit$lambda.min, n_selected = sum(coef_vec != 0))
}

fit_en_cox_sim <- function(surv_df, X, alpha = 0.5) {
  surv_obj <- Surv(surv_df$time, surv_df$event)
  cv_fit <- cv.glmnet(x = X, y = surv_obj, family = "cox", alpha = alpha)
  final_fit <- glmnet(x = X, y = surv_obj, family = "cox", alpha = alpha, lambda = cv_fit$lambda.min)
  coef_vec <- as.numeric(coef(final_fit))
  list(model = final_fit, coef = coef_vec, lambda = cv_fit$lambda.min, n_selected = sum(coef_vec != 0))
}

fit_ben_cox_sim <- function(surv_df, X, n_iter = 1000, n_warmup = 500) {
  ord <- order(surv_df$time, decreasing = FALSE)
  X_sorted <- X[ord, ]
  y_sorted <- surv_df$time[ord]
  e_sorted <- surv_df$event[ord]
  
  stan_data <- list(N = nrow(X_sorted), P = ncol(X_sorted), X = X_sorted, y = y_sorted, event = e_sorted)
  
  fit <- sampling(stan_model_compiled, data = stan_data,
                  iter = n_iter, warmup = n_warmup, chains = 2, 
                  control = list(adapt_delta = 0.9, max_treedepth = 12),
                  init = function() list(
                    beta_raw = rnorm(ncol(X_sorted), 0, 0.1),
                    tau_sq = rep(0.1, ncol(X_sorted)),
                    lambda1 = 1,
                    lambda2 = 1
                  ),
                  refresh = 0, show_messages = FALSE)
  
  samples <- rstan::extract(fit)
  beta_mean <- colMeans(samples$beta)
  beta_q <- apply(samples$beta, 2, quantile, probs = c(0.025, 0.975))
  
  selected <- (beta_q[1, ] > 0) | (beta_q[2, ] < 0)
  
  list(coef = beta_mean, lower = beta_q[1, ], upper = beta_q[2, ], n_selected = sum(selected), selected = selected,
       lambda1 = mean(samples$lambda1), lambda2 = mean(samples$lambda2))
}

## ============================================================
## SECTION 4: EVALUATION METRICS
## ============================================================
compute_estimation_error <- function(beta_hat, beta_true) {
  sqrt(sum((beta_hat - beta_true)^2))
}

compute_sensitivity <- function(selected, beta_true, threshold = 1e-6) {
  true_nonzero <- abs(beta_true) > threshold
  if (sum(true_nonzero) == 0) return(NA)
  sum(selected & true_nonzero) / sum(true_nonzero)
}

compute_specificity <- function(selected, beta_true, threshold = 1e-6) {
  true_zero <- abs(beta_true) <= threshold
  if (sum(true_zero) == 0) return(NA)
  sum(!selected & true_zero) / sum(true_zero)
}

compute_c_index_sim <- function(surv_df, risk) {
  cf <- concordance(Surv(time, event) ~ I(-risk), data = surv_df)
  return(cf$concordance)
}

compute_ibs_sim <- function(surv_df, risk, tau_quantile = 0.9) {
  max_time <- quantile(surv_df$time, tau_quantile)
  tp <- seq(0, max_time, length.out = 30)
  tp <- tp[tp > 0]
  
  km <- survfit(Surv(time, event) ~ 1, data = surv_df)
  base_surv <- summary(km, times = tp, extend = TRUE)$surv
  pred_surv <- outer(base_surv, exp(risk), "^")
  
  bs <- sapply(seq_along(tp), function(i) {
    obs <- as.numeric(surv_df$time > tp[i])
    mean((obs - pred_surv[i, ])^2)
  })
  
  dt <- diff(tp)
  sum((bs[-1] + bs[-length(bs)]) * dt / 2) / max_time
}

## ============================================================
## SECTION 5: SINGLE REPLICATE EVALUATION
## ============================================================
run_single_replicate <- function(n, p, s0 = 10, beta_val = 0.5, cov_type = "independent", 
                                 rho = 0.5, n_test = 300, ben_iter = 1000, ben_warmup = 500) {
  
  train_data <- generate_survival_data(n, p, s0, beta_val, cov_type, rho)
  train_df <- data.frame(time = train_data$time, event = train_data$event)
  
  test_data <- generate_survival_data(n_test, p, s0, beta_val, cov_type, rho)
  test_df <- data.frame(time = test_data$time, event = test_data$event)
  
  beta_true <- train_data$beta_true
  
  # Fit all models
  null_fit <- fit_null_cox_sim(train_df)
  ridge_fit <- fit_ridge_cox_sim(train_df, train_data$X)
  lasso_fit <- fit_lasso_cox_sim(train_df, train_data$X)
  en_fit <- fit_en_cox_sim(train_df, train_data$X, alpha = 0.5)
  ben_fit <- fit_ben_cox_sim(train_df, train_data$X, n_iter = ben_iter, n_warmup = ben_warmup)
  
  # Compute test risks
  null_risk <- rep(0, n_test)
  ridge_risk <- as.numeric(predict(ridge_fit$model, newx = test_data$X, s = ridge_fit$lambda, type = "link"))
  lasso_risk <- as.numeric(predict(lasso_fit$model, newx = test_data$X, s = lasso_fit$lambda, type = "link"))
  en_risk <- as.numeric(predict(en_fit$model, newx = test_data$X, s = en_fit$lambda, type = "link"))
  ben_risk <- as.numeric(test_data$X %*% ben_fit$coef)
  
  # Variable selection
  ridge_selected <- rep(TRUE, p)
  lasso_selected <- lasso_fit$coef != 0
  en_selected <- en_fit$coef != 0
  ben_selected <- ben_fit$selected
  
  # Compute metrics
  results <- data.frame(
    Model = c("Null Cox", "Ridge Cox", "Lasso Cox", "EN-Cox", "BEN-Cox"),
    EstError = c(
      compute_estimation_error(rep(0, p), beta_true),
      compute_estimation_error(ridge_fit$coef, beta_true),
      compute_estimation_error(lasso_fit$coef, beta_true),
      compute_estimation_error(en_fit$coef, beta_true),
      compute_estimation_error(ben_fit$coef, beta_true)
    ),
    Sensitivity = c(
      NA, NA,
      compute_sensitivity(lasso_selected, beta_true),
      compute_sensitivity(en_selected, beta_true),
      compute_sensitivity(ben_selected, beta_true)
    ),
    Specificity = c(
      NA, NA,
      compute_specificity(lasso_selected, beta_true),
      compute_specificity(en_selected, beta_true),
      compute_specificity(ben_selected, beta_true)
    ),
    C_index = c(
      compute_c_index_sim(test_df, null_risk),
      compute_c_index_sim(test_df, ridge_risk),
      compute_c_index_sim(test_df, lasso_risk),
      compute_c_index_sim(test_df, en_risk),
      compute_c_index_sim(test_df, ben_risk)
    ),
    IBS = c(
      compute_ibs_sim(test_df, null_risk),
      compute_ibs_sim(test_df, ridge_risk),
      compute_ibs_sim(test_df, lasso_risk),
      compute_ibs_sim(test_df, en_risk),
      compute_ibs_sim(test_df, ben_risk)
    ),
    N_selected = c(0, p, lasso_fit$n_selected, en_fit$n_selected, ben_fit$n_selected)
  )
  
  return(results)
}

## ============================================================
## SECTION 6: RUN FULL SIMULATION
## ============================================================
run_simulation_study <- function(n_reps = 100, 
                                 scenarios = NULL,
                                 ben_iter = 1000, 
                                 ben_warmup = 500,
                                 n_test = 300,
                                 s0 = 10,
                                 beta_val = 0.5) {
  
  if (is.null(scenarios)) {
    scenarios <- expand.grid(
      n = c(150, 300),
      p = c(50, 150),
      cov_type = c("independent", "toeplitz", "block"),
      rho = 0.5,
      stringsAsFactors = FALSE
    )
  }
  
  all_results <- list()
  
  for (sc in 1:nrow(scenarios)) {
    n <- scenarios$n[sc]
    p <- scenarios$p[sc]
    cov_type <- scenarios$cov_type[sc]
    rho <- scenarios$rho[sc]
    
    cat(sprintf("\n=== Scenario %d/%d: n=%d, p=%d, cov=%s ===\n", 
                sc, nrow(scenarios), n, p, cov_type))
    
    scenario_results <- list()
    
    for (rep in 1:n_reps) {
      if (rep %% 10 == 0) cat(sprintf("  Replicate %d/%d\n", rep, n_reps))
      
      tryCatch({
        res <- run_single_replicate(n, p, s0, beta_val, cov_type, rho, n_test, ben_iter, ben_warmup)
        res$replicate <- rep
        res$n <- n
        res$p <- p
        res$cov_type <- cov_type
        res$rho <- rho
        scenario_results[[rep]] <- res
      }, error = function(e) {
        cat(sprintf("    Error in rep %d: %s\n", rep, e$message))
      })
    }
    
    all_results[[sc]] <- do.call(rbind, scenario_results)
  }
  
  final_results <- do.call(rbind, all_results)
  return(final_results)
}

## ============================================================
## SECTION 7: SUMMARIZE AND GENERATE TABLES
## ============================================================
summarize_simulation_results <- function(results) {
  summary_df <- results %>%
    group_by(n, p, cov_type, rho, Model) %>%
    summarize(
      EstError_mean = mean(EstError, na.rm = TRUE),
      EstError_se = sd(EstError, na.rm = TRUE) / sqrt(n()),
      Sensitivity_mean = mean(Sensitivity, na.rm = TRUE),
      Sensitivity_se = sd(Sensitivity, na.rm = TRUE) / sqrt(sum(!is.na(Sensitivity))),
      Specificity_mean = mean(Specificity, na.rm = TRUE),
      Specificity_se = sd(Specificity, na.rm = TRUE) / sqrt(sum(!is.na(Specificity))),
      C_index_mean = mean(C_index, na.rm = TRUE),
      C_index_se = sd(C_index, na.rm = TRUE) / sqrt(n()),
      IBS_mean = mean(IBS, na.rm = TRUE),
      IBS_se = sd(IBS, na.rm = TRUE) / sqrt(n()),
      N_selected_mean = mean(N_selected, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(summary_df)
}

generate_simulation_latex_table <- function(summary_df, n_val, p_val, cov_val) {
  df <- summary_df %>% filter(n == n_val, p == p_val, cov_type == cov_val)
  
  if (nrow(df) == 0) {
    cat(sprintf("No data for n=%d, p=%d, cov=%s\n", n_val, p_val, cov_val))
    return()
  }
  
  cat(sprintf("\n=== LATEX TABLE: n=%d, p=%d, %s correlation ===\n", n_val, p_val, cov_val))
  cat("\\begin{table}[h!]\n\\centering\n")
  cat(sprintf("\\caption{\\rev{Simulation results for $(n, p) = (%d, %d)$ with %s correlation. Values are mean $\\pm$ SE over 100 replicates.}}\n", n_val, p_val, cov_val))
  cat("\\label{tab:sim_main}\n")
  cat("\\rev{\n")
  cat("\\begin{tabular}{lccccc}\n\\toprule\n")
  cat("Model & $\\|\\widehat{\\bm{\\beta}} - \\bm{\\beta}^\\star\\|_2$ $\\downarrow$ & Sens. $\\uparrow$ & Spec. $\\uparrow$ & C-index $\\uparrow$ & IBS $\\downarrow$ \\\\\n\\midrule\n")
  
  for (i in 1:nrow(df)) {
    r <- df[i, ]
    sens_str <- ifelse(is.na(r$Sensitivity_mean), "---", sprintf("%.2f $\\pm$ %.2f", r$Sensitivity_mean, r$Sensitivity_se))
    spec_str <- ifelse(is.na(r$Specificity_mean), "---", sprintf("%.2f $\\pm$ %.2f", r$Specificity_mean, r$Specificity_se))
    
    if (r$Model == "BEN-Cox") {
      cat(sprintf("\\textbf{BEN--Cox} & %.2f $\\pm$ %.2f & %s & %s & %.3f $\\pm$ %.3f & %.3f $\\pm$ %.3f \\\\\n",
                  r$EstError_mean, r$EstError_se, sens_str, spec_str, r$C_index_mean, r$C_index_se, r$IBS_mean, r$IBS_se))
    } else {
      cat(sprintf("%s & %.2f $\\pm$ %.2f & %s & %s & %.3f $\\pm$ %.3f & %.3f $\\pm$ %.3f \\\\\n",
                  r$Model, r$EstError_mean, r$EstError_se, sens_str, spec_str, r$C_index_mean, r$C_index_se, r$IBS_mean, r$IBS_se))
    }
  }
  
  cat("\\bottomrule\n\\end{tabular}\n}\n\\end{table}\n")
}

## ============================================================
## SECTION 8: PLOTTING
## ============================================================
plot_simulation_boxplots <- function(results, output_dir) {
  filename <- file.path(output_dir, "sim_boxplots.pdf")
  cat(sprintf("Saving simulation boxplots to '%s'...\n", filename))
  
  plot_data <- results %>% filter(n == 300, p == 150)
  
  if (nrow(plot_data) == 0) {
    plot_data <- results %>% filter(n == max(results$n), p == max(results$p))
  }
  
  p1 <- ggplot(plot_data, aes(x = Model, y = EstError, fill = cov_type)) +
    geom_boxplot() +
    labs(title = sprintf("Estimation Error (n=%d, p=%d)", plot_data$n[1], plot_data$p[1]),
         y = expression(paste("||", hat(beta), " - ", beta^"*", "||")[2]),
         x = "", fill = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(plot_data, aes(x = Model, y = C_index, fill = cov_type)) +
    geom_boxplot() +
    labs(title = sprintf("C-index (n=%d, p=%d)", plot_data$n[1], plot_data$p[1]),
         y = "C-index", x = "", fill = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  pdf(filename, width = 12, height = 5)
  grid.arrange(p1, p2, ncol = 2)
  dev.off()
}

plot_simulation_by_pn_ratio <- function(results, output_dir) {
  filename <- file.path(output_dir, "sim_pn_ratio.pdf")
  cat(sprintf("Saving p/n ratio plots to '%s'...\n", filename))
  
  plot_data <- results %>%
    mutate(pn_ratio = p / n) %>%
    filter(cov_type == "toeplitz")
  
  if (nrow(plot_data) == 0) {
    plot_data <- results %>% mutate(pn_ratio = p / n)
  }
  
  summary_data <- plot_data %>%
    group_by(pn_ratio, Model) %>%
    summarize(
      EstError_mean = mean(EstError, na.rm = TRUE),
      EstError_se = sd(EstError, na.rm = TRUE) / sqrt(n()),
      C_index_mean = mean(C_index, na.rm = TRUE),
      C_index_se = sd(C_index, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  p1 <- ggplot(summary_data, aes(x = pn_ratio, y = EstError_mean, color = Model, group = Model)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = EstError_mean - EstError_se, ymax = EstError_mean + EstError_se), width = 0.02) +
    labs(title = "Estimation Error vs p/n Ratio",
         x = "p/n Ratio", y = "Estimation Error") +
    theme_minimal()
  
  p2 <- ggplot(summary_data, aes(x = pn_ratio, y = C_index_mean, color = Model, group = Model)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = C_index_mean - C_index_se, ymax = C_index_mean + C_index_se), width = 0.02) +
    labs(title = "C-index vs p/n Ratio",
         x = "p/n Ratio", y = "C-index") +
    theme_minimal()
  
  pdf(filename, width = 12, height = 5)
  grid.arrange(p1, p2, ncol = 2)
  dev.off()
}

plot_sensitivity_specificity <- function(results, output_dir) {
  filename <- file.path(output_dir, "sim_sens_spec.pdf")
  cat(sprintf("Saving sensitivity/specificity plot to '%s'...\n", filename))
  
  plot_data <- results %>%
    filter(Model %in% c("Lasso Cox", "EN-Cox", "BEN-Cox"))
  
  # Use largest n and p available
  max_n <- max(plot_data$n)
  max_p <- max(plot_data$p)
  plot_data <- plot_data %>% filter(n == max_n, p == max_p)
  
  p1 <- ggplot(plot_data, aes(x = Model, y = Sensitivity, fill = cov_type)) +
    geom_boxplot() +
    labs(title = sprintf("Sensitivity (n=%d, p=%d)", max_n, max_p),
         y = "Sensitivity", x = "", fill = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, 1)
  
  p2 <- ggplot(plot_data, aes(x = Model, y = Specificity, fill = cov_type)) +
    geom_boxplot() +
    labs(title = sprintf("Specificity (n=%d, p=%d)", max_n, max_p),
         y = "Specificity", x = "", fill = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, 1)
  
  pdf(filename, width = 12, height = 5)
  grid.arrange(p1, p2, ncol = 2)
  dev.off()
}

## ============================================================
## SECTION 9: MAIN EXECUTION (WITH output_dir INSIDE FUNCTION)
## ============================================================
run_full_simulation_study <- function(n_reps = 100, ben_iter = 1000, ben_warmup = 500) {
  
  # === CREATE OUTPUT DIRECTORY INSIDE FUNCTION ===
  output_dir <- "results_simulation"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(sprintf("Created output directory: %s\n", output_dir))
  }
  # ===============================================
  
  cat("=== STARTING SIMULATION STUDY (CORRECTED VERSION) ===\n")
  cat("Key fix: Changed lambda priors from Normal(0, 0.05) to Half-Cauchy(0, 1)\n\n")
  cat(sprintf("Replicates: %d\n", n_reps))
  cat(sprintf("BEN-Cox iterations: %d (warmup: %d)\n", ben_iter, ben_warmup))
  cat(sprintf("Output directory: %s\n", output_dir))
  
  # REDUCED scenarios: 2 n x 2 p x 3 correlations = 12 scenarios
  scenarios <- expand.grid(
    n = c(150, 300),
    p = c(50, 150),
    cov_type = c("independent", "toeplitz", "block"),
    rho = 0.5,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("\nTotal scenarios: %d\n", nrow(scenarios)))
  print(scenarios)
  
  # Run simulation
  results <- run_simulation_study(
    n_reps = n_reps,
    scenarios = scenarios,
    ben_iter = ben_iter,
    ben_warmup = ben_warmup
  )
  
  # Save raw results
  write.csv(results, file.path(output_dir, "simulation_raw_results.csv"), row.names = FALSE)
  cat(sprintf("\nSaved raw results to '%s/simulation_raw_results.csv'\n", output_dir))
  
  # Summarize
  summary_df <- summarize_simulation_results(results)
  write.csv(summary_df, file.path(output_dir, "simulation_summary.csv"), row.names = FALSE)
  cat(sprintf("Saved summary to '%s/simulation_summary.csv'\n", output_dir))
  
  # Generate LaTeX tables (main scenarios)
  generate_simulation_latex_table(summary_df, 300, 150, "toeplitz")
  generate_simulation_latex_table(summary_df, 300, 150, "block")
  generate_simulation_latex_table(summary_df, 150, 150, "toeplitz")
  
  # Generate plots
  plot_simulation_boxplots(results, output_dir)
  plot_simulation_by_pn_ratio(results, output_dir)
  plot_sensitivity_specificity(results, output_dir)
  
  cat("\n=== SIMULATION STUDY COMPLETE ===\n")
  cat(sprintf("All outputs saved to: %s/\n", output_dir))
  cat("Output files:\n")
  cat("  - simulation_raw_results.csv\n")
  cat("  - simulation_summary.csv\n")
  cat("  - sim_boxplots.pdf\n")
  cat("  - sim_pn_ratio.pdf\n")
  cat("  - sim_sens_spec.pdf\n")
  
  return(list(results = results, summary = summary_df))
}

# =====================
# EXECUTE
# =====================
# For testing: n_reps = 5, ben_iter = 100, ben_warmup = 50
# For final paper: n_reps = 100, ben_iter = 1000, ben_warmup = 500

sim_results <- run_full_simulation_study(n_reps = 50, ben_iter = 1000, ben_warmup = 500)
