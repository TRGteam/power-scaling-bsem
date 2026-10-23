# ================================================================
# SIMULASI UTAMA - STRICT POWER-SCALING (OPSI B)
# 12 kondisi × 50 replikasi = 600 fit
# Sample = 2000, Burnin = 1000
# ================================================================

rm(list = ls())
set.seed(2026)

# ================================================================
# KONFIGURASI PARALEL — MAKSIMAL
# ================================================================
NCORES <- 12
R_REP  <- 50

Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)

suppressPackageStartupMessages({
  library(blavaan); library(future); library(future.apply)
  library(dplyr); library(tidyr); library(ggplot2)
  library(MASS); library(rstan)
})

options(parallelly.maxWorkers.localhost = Inf)
options(mc.cores = 1)
rstan_options(auto_write = TRUE)
options(future.globals.maxSize = 12 * 1024^3)
plan(multisession, workers = NCORES)

cat("[OK] Parallel backend:", NCORES, "workers\n")

setwd("~/SMIBSEMREVISI6")
dir.create("Simulasi/HasilSimulasi_V6c", showWarnings = FALSE, recursive = TRUE)
setwd("Simulasi")

# ================================================================
# TRUE VALUES (MODERATE)
# ================================================================
true_params <- list(
  beta_z1 = -0.25, gamma11 = 0.15, gamma12 = 0.40,
  loading = 0.70, error_var = 1 - 0.70^2
)

true_vals <- c(
  beta_z1   = true_params$beta_z1,
  gamma11   = true_params$gamma11,
  gamma12   = true_params$gamma12,
  indirect1 = true_params$gamma11 * true_params$beta_z1,
  indirect2 = true_params$gamma12 * true_params$beta_z1
)

cat("\n=== TRUE VALUES (MODERATE) ===\n"); print(round(true_vals, 3))

# ================================================================
# GENERATE DATA
# ================================================================
generate_data <- function(n, rho, dist = "normal", true_params, seed) {
  set.seed(seed)
  Phi <- matrix(c(1, rho, rho, 1), 2, 2)
  ksi <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Phi)
  ksi1 <- ksi[, 1]; ksi2 <- ksi[, 2]
  
  resid_var_eta1 <- 1 - (true_params$gamma11^2 + true_params$gamma12^2 +
                           2 * true_params$gamma11 * true_params$gamma12 * rho)
  if (resid_var_eta1 <= 0.05) return(NULL)
  resid_var_Z <- 1 - true_params$beta_z1^2
  
  eta1 <- true_params$gamma11 * ksi1 +
    true_params$gamma12 * ksi2 + rnorm(n, 0, sqrt(resid_var_eta1))
  Z <- true_params$beta_z1 * eta1 + rnorm(n, 0, sqrt(resid_var_Z))
  
  lam <- true_params$loading; err <- sqrt(true_params$error_var)
  if (dist == "normal") e <- function() rnorm(n, 0, err)
  else if (dist == "t") {
    df_t <- 5; scale_t <- err * sqrt((df_t - 2) / df_t)
    e <- function() scale_t * rt(n, df_t)
  }
  
  x11 <- lam * ksi1 + e(); x21 <- lam * ksi1 + e(); x31 <- lam * ksi1 + e()
  x12 <- lam * ksi2 + e(); x22 <- lam * ksi2 + e(); x32 <- lam * ksi2 + e()
  y11 <- lam * eta1 + e(); y21 <- lam * eta1 + e(); y31 <- lam * eta1 + e()
  
  data.frame(x11, x21, x31, x12, x22, x32, y11, y21, y31, Z)
}

# ================================================================
# ANALISIS — STRICT
# ================================================================
analyse_bsem <- function(data, tau_theta = 1, tau_beta = 1,
                         tau_lambda = 1, tau_psi = 1, seed) {
  set.seed(seed)
  data_std <- scale(data) |> as.data.frame()
  
  prior_theta  <- paste0("gamma(", tau_theta * 0 + 1, ", ", 0.5 * tau_theta, ")")
  prior_psi    <- paste0("gamma(", 3 * tau_psi + 1, ", ", 1 * tau_psi, ")")
  prior_beta   <- paste0("normal(0, ", 0.5 / sqrt(tau_beta), ")")
  prior_lambda <- paste0("normal(1, ", 0.5 / sqrt(tau_lambda), ")")
  
  dp <- dpriors(lambda = prior_lambda, beta = prior_beta,
                psi = prior_psi, theta = prior_theta)
  
  model <- '
    Z ~ beta_z1*eta1
    eta1 ~ gamma11*ksi1 + gamma12*ksi2
    eta1 =~ y11 + y21 + y31
    ksi1 =~ x11 + x21 + x31
    ksi2 =~ x12 + x22 + x32
    ksi1 ~~ phi21*ksi2
  '
  
  fit <- tryCatch(
    bsem(model, data = data_std, dp = dp,
         burnin = 1000, sample = 2000, save.lvs = FALSE,
         target = "stanclassic",
         bcontrol = list(adapt_delta = 0.99, max_treedepth = 12)),
    error = function(e) NULL
  )
  
  na_result <- c(beta_z1 = NA, gamma11 = NA, gamma12 = NA,
                 indirect1 = NA, indirect2 = NA,
                 ci_lo_beta = NA, ci_hi_beta = NA,
                 ci_lo_g11 = NA, ci_hi_g11 = NA,
                 ci_lo_g12 = NA, ci_hi_g12 = NA,
                 ci_lo_ind1 = NA, ci_hi_ind1 = NA,
                 ci_lo_ind2 = NA, ci_hi_ind2 = NA,
                 rhat_max = NA, ess_min = NA, converged = NA)
  
  if (is.null(fit)) return(na_result)
  post <- tryCatch(as.matrix(blavInspect(fit, "mcmc")), error = function(e) NULL)
  if (is.null(post)) return(na_result)
  if (!all(c("beta_z1","gamma11","gamma12") %in% colnames(post))) return(na_result)
  
  b <- post[, "beta_z1"]; g1 <- post[, "gamma11"]; g2 <- post[, "gamma12"]
  ind1 <- g1 * b; ind2 <- g2 * b
  q <- function(x, p) quantile(x, p, names = FALSE)
  rhat_vals <- tryCatch(blavInspect(fit, "psrf"), error = function(e) NA)
  ess_vals  <- tryCatch(blavInspect(fit, "neff"), error = function(e) NA)
  
  c(beta_z1 = mean(b), gamma11 = mean(g1), gamma12 = mean(g2),
    indirect1 = mean(ind1), indirect2 = mean(ind2),
    ci_lo_beta = q(b, 0.025), ci_hi_beta = q(b, 0.975),
    ci_lo_g11 = q(g1, 0.025), ci_hi_g11 = q(g1, 0.975),
    ci_lo_g12 = q(g2, 0.025), ci_hi_g12 = q(g2, 0.975),
    ci_lo_ind1 = q(ind1, 0.025), ci_hi_ind1 = q(ind1, 0.975),
    ci_lo_ind2 = q(ind2, 0.025), ci_hi_ind2 = q(ind2, 0.975),
    rhat_max = ifelse(all(is.na(rhat_vals)), NA, max(rhat_vals, na.rm = TRUE)),
    ess_min  = ifelse(all(is.na(ess_vals)), NA, min(ess_vals, na.rm = TRUE)),
    converged = ifelse(is.na(max(rhat_vals, na.rm = TRUE)), NA,
                       max(rhat_vals, na.rm = TRUE) < 1.05))
}

# ================================================================
# DESAIN
# ================================================================
Design <- expand.grid(
  n        = c(75, 119, 250),
  rho      = c(0.30, 0.80),
  dist     = c("normal", "t"),
  tau_beta = c(1),
  stringsAsFactors = FALSE
)

cat("\n=== DESAIN SIMULASI ===\n")
cat("Kondisi:", nrow(Design), "| Reps:", R_REP,
    "| Total:", nrow(Design) * R_REP, "\n")
cat("Sample/chain: 2000 | Burnin: 1000\n")
cat("Estimasi: ~18 jam\n\n")

# ================================================================
# JALANKAN
# ================================================================
run_condition <- function(row_idx, Design, R_REP, true_params, true_vals) {
  cond <- Design[row_idx, ]
  cat(sprintf("\n>>> Kondisi %d/%d: n=%d, rho=%.2f, dist=%s\n",
              row_idx, nrow(Design), cond$n, cond$rho, cond$dist))
  t_start <- Sys.time()
  
  results <- future_lapply(1:R_REP, function(r) {
    seed_val <- as.integer(row_idx * 10000 + r) %% 2147483647L
    dat <- generate_data(n = cond$n, rho = cond$rho, dist = cond$dist,
                         true_params = true_params, seed = seed_val)
    if (is.null(dat)) return(rep(NA, 17))
    analyse_bsem(data = dat, tau_beta = cond$tau_beta, seed = seed_val)
  }, future.seed = TRUE)
  
  results_df <- do.call(rbind, lapply(results, function(x) as.data.frame(t(x))))
  
  est_names <- c("beta_z1", "gamma11", "gamma12", "indirect1", "indirect2")
  out <- data.frame(n = cond$n, rho = cond$rho, dist = cond$dist,
                    tau_beta = cond$tau_beta)
  
  for (p in est_names) {
    est <- results_df[[p]]; tv <- true_vals[[p]]; ok <- !is.na(est)
    out[[paste0("bias_", p)]] <- mean(est[ok] - tv, na.rm = TRUE)
    out[[paste0("rmse_", p)]] <- sqrt(mean((est[ok] - tv)^2, na.rm = TRUE))
    out[[paste0("sd_", p)]]   <- sd(est[ok], na.rm = TRUE)
  }
  
  cover <- function(lo, hi, tv) mean(lo <= tv & tv <= hi, na.rm = TRUE)
  out$cover_beta_z1   <- cover(results_df$ci_lo_beta, results_df$ci_hi_beta, true_vals[["beta_z1"]])
  out$cover_gamma11   <- cover(results_df$ci_lo_g11,  results_df$ci_hi_g11,  true_vals[["gamma11"]])
  out$cover_gamma12   <- cover(results_df$ci_lo_g12,  results_df$ci_hi_g12,  true_vals[["gamma12"]])
  out$cover_indirect1 <- cover(results_df$ci_lo_ind1, results_df$ci_hi_ind1, true_vals[["indirect1"]])
  out$cover_indirect2 <- cover(results_df$ci_lo_ind2, results_df$ci_hi_ind2, true_vals[["indirect2"]])
  
  out$mean_rhat <- mean(results_df$rhat_max, na.rm = TRUE)
  out$prop_rhat_ok <- mean(results_df$rhat_max < 1.05, na.rm = TRUE)
  out$mean_ess <- mean(results_df$ess_min, na.rm = TRUE)
  out$prop_conv <- mean(results_df$converged, na.rm = TRUE)
  out$n_fail <- sum(is.na(results_df$beta_z1))
  
  t_el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  cat(sprintf("  Selesai %.1f menit. Rhat=%.3f | ESS=%.0f | fail=%d\n",
              t_el, out$mean_rhat, out$mean_ess, out$n_fail))
  
  list(summary = out, results = results_df)
}

t0 <- Sys.time()
all_summaries <- list(); all_results <- list()

for (i in 1:nrow(Design)) {
  res_i <- run_condition(i, Design, R_REP, true_params, true_vals)
  all_summaries[[i]] <- res_i$summary
  all_results[[i]]   <- res_i$results
  saveRDS(all_summaries, "HasilSimulasi_V6c/progress_summaries.rds")
  saveRDS(all_results,   "HasilSimulasi_V6c/progress_results.rds")
}

t_total <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
cat("\n>>> SELESAI dalam:", round(t_total, 2), "jam\n")

res_df <- do.call(rbind, all_summaries)
write.csv(res_df, "HasilSimulasi_V6c/sim_summary_v6c.csv", row.names = FALSE)

# ================================================================
# VERIFIKASI ESS
# ================================================================
cat("\n=== VERIFIKASI ===\n")
cat("ESS min:", round(min(res_df$mean_ess), 1), "\n")
cat("ESS < 400 di", sum(res_df$mean_ess < 400), "kondisi\n")
cat("Rhat max:", round(max(res_df$mean_rhat), 4), "\n")
cat("Total fail:", sum(res_df$n_fail), "\n")

# ================================================================
# TABEL 2 UNTUK NASKAH
# ================================================================
table_2 <- res_df %>%
  dplyr::select(n, rho, dist,
                bias_beta_z1, rmse_beta_z1, cover_beta_z1,
                bias_gamma11, rmse_gamma11, cover_gamma11,
                bias_gamma12, rmse_gamma12, cover_gamma12,
                mean_rhat, mean_ess, prop_conv, n_fail) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

write.csv(table_2, "HasilSimulasi_V6c/table_2_main_simulation.csv",
          row.names = FALSE)
cat("[OK] table_2_main_simulation.csv\n")

# ================================================================
# FIGURE 2: BIAS & RMSE
# ================================================================
res_long <- res_df %>%
  pivot_longer(cols = starts_with(c("bias_", "rmse_", "cover_")),
               names_to = c("metric", "param"),
               names_pattern = "(bias|rmse|cover)_(.+)",
               values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = c("bias","rmse","cover")),
    param  = factor(param,
                    levels = c("beta_z1","gamma11","gamma12","indirect1","indirect2"),
                    labels = c("beta[z1]~(eta1%->%Z)",
                               "gamma[11]~(ksi1%->%eta1)",
                               "gamma[12]~(ksi2%->%eta1)",
                               "Indirect[1]",
                               "Indirect[2]"))
  )

p1_data <- res_long %>%
  filter(metric %in% c("bias","rmse"),
         param %in% c("beta[z1]~(eta1%->%Z)",
                      "gamma[11]~(ksi1%->%eta1)",
                      "gamma[12]~(ksi2%->%eta1)"))

p_fig2 <- ggplot(p1_data, aes(x = factor(n), y = value, fill = dist)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_grid(metric ~ param, scales = "free_y", labeller = label_parsed) +
  scale_fill_manual(values = c("normal"="#4C72B0","t"="#DD8452"),
                    name = "Error dist.") +
  labs(x = "Sample size (n)", y = "Value",
       title = "Bias and RMSE across sample sizes") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 8))

ggsave("HasilSimulasi_V6c/Fig2_Bias_RMSE.png", p_fig2,
       width = 10, height = 6, dpi = 300)
cat("[OK] Fig2_Bias_RMSE.png\n")

# ================================================================
# FIGURE 3: COVERAGE
# ================================================================
p2_data <- res_long %>%
  filter(metric == "cover",
         param %in% c("beta[z1]~(eta1%->%Z)",
                      "gamma[11]~(ksi1%->%eta1)",
                      "gamma[12]~(ksi2%->%eta1)",
                      "Indirect[2]"))

p_fig3 <- ggplot(p2_data, aes(x = factor(n), y = value,
                              color = dist, group = dist)) +
  geom_hline(yintercept = 0.95, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  geom_point(size = 2.5, position = position_dodge(0.4)) +
  geom_line(position = position_dodge(0.4), linewidth = 0.8) +
  facet_wrap(~ param, labeller = label_parsed, nrow = 1) +
  scale_color_manual(values = c("normal"="#4C72B0","t"="#DD8452"),
                     name = "Error dist.") +
  labs(x = "Sample size (n)", y = "Coverage of 95% CI",
       title = "Coverage probability") +
  ylim(0.7, 1.0) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("HasilSimulasi_V6c/Fig3_Coverage.png", p_fig3,
       width = 10, height = 4.5, dpi = 300)
cat("[OK] Fig3_Coverage.png\n")

# ================================================================
# FIGURE 4: MCMC CONVERGENCE DIAGNOSTICS
# ================================================================
conv_data <- res_df %>%
  dplyr::select(n, rho, dist, mean_rhat, mean_ess, prop_rhat_ok, prop_conv) %>%
  pivot_longer(
    cols = c(mean_rhat, mean_ess, prop_rhat_ok, prop_conv),
    names_to = "diagnostic", values_to = "value"
  ) %>%
  mutate(diagnostic = factor(diagnostic,
                             levels = c("mean_rhat", "mean_ess",
                                        "prop_rhat_ok", "prop_conv"),
                             labels = c("Mean Rhat", "Mean ESS",
                                        "Prop. Rhat < 1.05",
                                        "Prop. converged")))

p_fig4 <- ggplot(conv_data, aes(x = factor(n), y = value, fill = dist)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  facet_wrap(~ diagnostic, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("normal"="#4C72B0","t"="#DD8452"),
                    name = "Error dist.") +
  labs(x = "Sample size (n)", y = "Value",
       title = "MCMC convergence diagnostics") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("HasilSimulasi_V6c/Fig4_Convergence.png", p_fig4,
       width = 12, height = 4, dpi = 300)
cat("[OK] Fig4_Convergence.png\n")

# ================================================================
# FIGURE S: HEATMAP RMSE
# ================================================================
heat_data <- res_df %>%
  dplyr::select(n, rho, dist, rmse_beta_z1, rmse_gamma12) %>%
  pivot_longer(cols = c(rmse_beta_z1, rmse_gamma12),
               names_to = "param", values_to = "rmse") %>%
  mutate(
    param = factor(param,
                   levels = c("rmse_beta_z1", "rmse_gamma12"),
                   labels = c("beta[z1]", "gamma[12]")),
    n   = factor(n),
    rho = factor(rho)
  )

p_heat <- ggplot(heat_data, aes(x = n, y = rho, fill = rmse)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", rmse)),
            color = "white", size = 3) +
  facet_grid(param ~ dist, labeller = label_parsed) +
  scale_fill_viridis_c(option = "C", direction = -1, name = "RMSE") +
  labs(x = "Sample size (n)",
       y = expression(rho~"(exogenous correlation)"),
       title = "RMSE heatmap across design conditions") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank())

ggsave("HasilSimulasi_V6c/FigS_Heatmap.png", p_heat,
       width = 10, height = 6, dpi = 300)
cat("[OK] FigS_Heatmap.png\n")

cat("\n========================================\n")
cat("SIMULASI UTAMA SELESAI\n")
cat("========================================\n")