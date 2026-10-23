# ================================================================
# SIMULASI SUPLEMEN - STRICT POWER-SCALING (DENGAN n=35)
# True values: empiris n=35
# Desain: 6 kondisi × 50 replikasi = 300 fit
#   - n = 35  → sample = 3000, adapt_delta = 0.999
#   - n = 119 → sample = 2000, adapt_delta = 0.99
#   - n = 250 → sample = 2000, adapt_delta = 0.99
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
dir.create("Simulasi/HasilSimulasi_V7c", showWarnings = FALSE, recursive = TRUE)
setwd("Simulasi")

# ================================================================
# TRUE VALUES: EMPIRIS n=35
# ================================================================
true_params <- list(
  beta_z1   = -0.587,
  gamma11   = -0.648,
  gamma12   =  0.421,
  loading   =  0.70,
  error_var = 1 - 0.70^2
)

true_vals <- c(
  beta_z1   = true_params$beta_z1,
  gamma11   = true_params$gamma11,
  gamma12   = true_params$gamma12,
  indirect1 = true_params$gamma11 * true_params$beta_z1,
  indirect2 = true_params$gamma12 * true_params$beta_z1
)

cat("\n=== TRUE VALUES (EMPIRIS n=35) ===\n")
print(round(true_vals, 3))

# ================================================================
# VERIFIKASI RESIDUAL VARIANCE
# ================================================================
cat("\n=== VERIFIKASI RESIDUAL VARIANCE ===\n")
for (rho_val in c(0.30, 0.80)) {
  rv <- 1 - (true_params$gamma11^2 + true_params$gamma12^2 +
               2 * true_params$gamma11 * true_params$gamma12 * rho_val)
  status <- if (rv > 0.05) "OK" else "SKIP"
  cat(sprintf("  rho = %.2f → resid_var = %.4f %s\n", rho_val, rv, status))
}

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
# ANALISIS — STRICT POWER-SCALING
# KONFIGURASI ADAPTIF BERDASARKAN n
# ================================================================
analyse_bsem <- function(data, tau_theta = 1, tau_beta = 1,
                         tau_lambda = 1, tau_psi = 1, seed) {
  set.seed(seed)
  data_std <- scale(data) |> as.data.frame()
  
  # --- Konfigurasi adaptif berdasarkan ukuran sampel ---
  n_obs <- nrow(data_std)
  if (n_obs <= 50) {
    # n=35: butuh sample lebih besar dan adapt_delta lebih tinggi
    sample_size  <- 3000
    adapt_delta_use <- 0.999
    max_td_use   <- 15
  } else {
    # n=119, 250: konfigurasi standar Opsi B
    sample_size  <- 2000
    adapt_delta_use <- 0.99
    max_td_use   <- 12
  }
  
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
         burnin = 1000,
         sample = sample_size,
         save.lvs = FALSE,
         target = "stanclassic",
         bcontrol = list(adapt_delta = adapt_delta_use,
                         max_treedepth = max_td_use)),
    error = function(e) NULL
  )
  
  na <- c(beta_z1 = NA, gamma11 = NA, gamma12 = NA,
          indirect1 = NA, indirect2 = NA,
          ci_lo_beta = NA, ci_hi_beta = NA,
          ci_lo_g11 = NA, ci_hi_g11 = NA,
          ci_lo_g12 = NA, ci_hi_g12 = NA,
          ci_lo_ind1 = NA, ci_hi_ind1 = NA,
          ci_lo_ind2 = NA, ci_hi_ind2 = NA,
          rhat_max = NA, ess_min = NA, converged = NA)
  
  if (is.null(fit)) return(na)
  post <- tryCatch(as.matrix(blavInspect(fit, "mcmc")), error = function(e) NULL)
  if (is.null(post)) return(na)
  if (!all(c("beta_z1","gamma11","gamma12") %in% colnames(post))) return(na)
  
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
# DESAIN SUPLEMEN — 6 KONDISI × 50 REP = 300 FIT
# ================================================================
Design <- expand.grid(
  n        = c(35, 119, 250),    # ← n=35 DITAMBAHKAN
  rho      = c(0.30, 0.80),
  dist     = c("normal"),
  tau_beta = c(1),
  stringsAsFactors = FALSE
)

cat("\n=== DESAIN SUPLEMEN (DENGAN n=35) ===\n")
cat("Kondisi   :", nrow(Design), "\n")
cat("Replikasi :", R_REP, "\n")
cat("Total fit :", nrow(Design) * R_REP, "\n")
cat("Sample    : n=35 → 3000 | n=119,250 → 2000\n")
cat("Estimasi  : ~5-6 jam\n\n")
print(Design)

# ================================================================
# JALANKAN
# ================================================================
t0 <- Sys.time()
all_summaries <- list()
all_results   <- list()

for (i in 1:nrow(Design)) {
  cond <- Design[i, ]
  cat(sprintf("\n>>> Kondisi %d/%d: n=%d, rho=%.2f, dist=%s\n",
              i, nrow(Design), cond$n, cond$rho, cond$dist))
  
  t_start <- Sys.time()
  
  results <- future_lapply(1:R_REP, function(r) {
    seed_val <- as.integer(i * 10000 + r)
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
  out$mean_ess  <- mean(results_df$ess_min, na.rm = TRUE)
  out$prop_conv <- mean(results_df$converged, na.rm = TRUE)
  out$n_fail    <- sum(is.na(results_df$beta_z1))
  
  t_el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  cat(sprintf("  Selesai %.1f menit. Rhat=%.3f | ESS=%.0f | fail=%d\n",
              t_el, out$mean_rhat, out$mean_ess, out$n_fail))
  
  all_summaries[[i]] <- out
  all_results[[i]]   <- results_df
  saveRDS(all_summaries, "HasilSimulasi_V7c/progress_v7c.rds")
  saveRDS(all_results,   "HasilSimulasi_V7c/progress_results_v7c.rds")
  cat(sprintf("  [CHECKPOINT] Kondisi %d tersimpan\n", i))
}

t_total <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
cat(sprintf("\n>>> SELESAI dalam %.2f jam\n", t_total))

res_df <- do.call(rbind, all_summaries)
write.csv(res_df, "HasilSimulasi_V7c/sim_summary_v7c.csv", row.names = FALSE)

# ================================================================
# VERIFIKASI
# ================================================================
cat("\n========================================\n")
cat("VERIFIKASI SIMULASI SUPLEMEN\n")
cat("========================================\n")

cat("Jumlah kondisi :", nrow(res_df), "\n")
cat("Total fail     :", sum(res_df$n_fail), "\n")
cat("Rhat max       :", round(max(res_df$mean_rhat), 4), "\n")
cat("ESS range      :", round(min(res_df$mean_ess), 1), "-",
    round(max(res_df$mean_ess), 1), "\n")
cat("ESS < 400 di   :", sum(res_df$mean_ess < 400), "kondisi\n")

cat("\n--- Ringkasan per kondisi ---\n")
res_df %>%
  mutate(kondisi = paste0("n=", n, ", rho=", rho)) %>%
  dplyr::select(kondisi, mean_rhat, mean_ess, prop_conv, n_fail,
                bias_beta_z1, cover_beta_z1,
                bias_gamma11, cover_gamma11,
                bias_gamma12, cover_gamma12) %>%
  print(digits = 3)

# ================================================================
# TABEL 3 UNTUK NASKAH
# ================================================================
table_3 <- res_df %>%
  dplyr::select(n, rho,
                bias_beta_z1, rmse_beta_z1, cover_beta_z1,
                bias_gamma11, rmse_gamma11, cover_gamma11,
                bias_gamma12, rmse_gamma12, cover_gamma12,
                mean_rhat, mean_ess, prop_conv, n_fail) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

write.csv(table_3, "HasilSimulasi_V7c/table_3_supplementary.csv",
          row.names = FALSE)
cat("[OK] table_3_supplementary.csv\n")

# ================================================================
# FIGURE S1: BIAS & RMSE (SUPPLEMENTARY)
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

# --- FigS1: Bias & RMSE ---
p_S1 <- ggplot(res_long %>%
                 filter(metric %in% c("bias","rmse"),
                        param %in% c("beta[z1]~(eta1%->%Z)",
                                     "gamma[11]~(ksi1%->%eta1)",
                                     "gamma[12]~(ksi2%->%eta1)")),
               aes(x = factor(n), y = value, fill = factor(rho))) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_grid(metric ~ param, scales = "free_y",
             labeller = label_parsed) +
  scale_fill_manual(values = c("0.3"="#4C72B0", "0.8"="#DD8452"),
                    name = expression(rho)) +
  labs(x = "Sample size (n)", y = "Value",
       title = "Supplementary Simulation: Bias and RMSE (with n=35)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("HasilSimulasi_V7c/FigS1_Bias_RMSE_Supp.png", p_S1,
       width = 10, height = 6, dpi = 300)
cat("[OK] FigS1_Bias_RMSE_Supp.png\n")

# --- FigS2: Coverage ---
p_S2 <- ggplot(res_long %>%
                 filter(metric == "cover",
                        param %in% c("beta[z1]~(eta1%->%Z)",
                                     "gamma[11]~(ksi1%->%eta1)",
                                     "gamma[12]~(ksi2%->%eta1)")),
               aes(x = factor(n), y = value,
                   color = factor(rho), group = factor(rho))) +
  geom_hline(yintercept = 0.95, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  geom_point(size = 3, position = position_dodge(0.4)) +
  geom_line(position = position_dodge(0.4), linewidth = 0.8) +
  facet_wrap(~ param, labeller = label_parsed, nrow = 1) +
  scale_color_manual(values = c("0.3"="#4C72B0", "0.8"="#DD8452"),
                     name = expression(rho)) +
  labs(x = "Sample size (n)", y = "Coverage of 95% CI",
       title = "Supplementary Simulation: Coverage (with n=35)") +
  ylim(0, 1.0) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("HasilSimulasi_V7c/FigS2_Coverage_Supp.png", p_S2,
       width = 10, height = 4.5, dpi = 300)
cat("[OK] FigS2_Coverage_Supp.png\n")

# ================================================================
# FIGURE S3: KONTINUM n (35 → 119 → 250)
# ================================================================
# Visualisasi khusus untuk menunjukkan efek ukuran sampel
continuum_data <- res_df %>%
  pivot_longer(cols = c(cover_beta_z1, cover_gamma11, cover_gamma12),
               names_to = "param", values_to = "coverage") %>%
  mutate(
    param = factor(param,
                   levels = c("cover_beta_z1", "cover_gamma11", "cover_gamma12"),
                   labels = c("beta[z1]", "gamma[11]", "gamma[12]")),
    n_factor = factor(n, levels = c(35, 119, 250))
  )

p_S3 <- ggplot(continuum_data,
               aes(x = n_factor, y = coverage,
                   color = factor(rho), group = factor(rho))) +
  geom_hline(yintercept = 0.95, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  geom_point(size = 3.5, position = position_dodge(0.4)) +
  geom_line(position = position_dodge(0.4), linewidth = 1) +
  facet_wrap(~ param, labeller = label_parsed, nrow = 1) +
  scale_color_manual(values = c("0.3"="#4C72B0", "0.8"="#DD8452"),
                     name = expression(rho)) +
  labs(x = "Sample size (n)", y = "Coverage of 95% CI",
       title = "Sample Size Continuum: Coverage from n=35 to n=250") +
  ylim(0, 1.0) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("HasilSimulasi_V7c/FigS3_SampleSize_Continuum.png", p_S3,
       width = 10, height = 4.5, dpi = 300)
cat("[OK] FigS3_SampleSize_Continuum.png\n")

# ================================================================
# ANALISIS KHUSUS n=35
# ================================================================
cat("\n========================================\n")
cat("ANALISIS KHUSUS n=35\n")
cat("========================================\n")

n35_results <- res_df %>% filter(n == 35)

cat("\n=== HASIL n=35 ===\n")
print(n35_results, digits = 3)

cat("\n=== PERBANDINGAN COVERAGE γ₁₁ ===\n")
cat("n=35  :", round(mean(n35_results$cover_gamma11), 3), "\n")
cat("n=119 :", round(mean(res_df$cover_gamma11[res_df$n == 119]), 3), "\n")
cat("n=250 :", round(mean(res_df$cover_gamma11[res_df$n == 250]), 3), "\n")

cat("\n=== INTERPRETASI ===\n")
cat("Jika coverage γ₁₁ di n=35 rendah, ini mengkonfirmasi\n")
cat("bahwa small sample size menyebabkan sensitivity γ₁₁,\n")
cat("bahkan ketika ρ moderat (0.30).\n")

# Simpan hasil n=35 khusus
write.csv(n35_results, "HasilSimulasi_V7c/results_n35_only.csv",
          row.names = FALSE)
cat("[OK] results_n35_only.csv\n")

# ================================================================
# RINGKASAN AKHIR
# ================================================================
cat("\n========================================\n")
cat("SIMULASI SUPLEMEN SELESAI\n")
cat("Output: HasilSimulasi_V7c/\n")
cat("Total waktu:", round(t_total, 2), "jam\n")
cat("========================================\n")

cat("\n=== FILE YANG DIHASILKAN ===\n")
files <- list.files("HasilSimulasi_V7c", pattern = "\\.(csv|png|rds)$")
for (f in files) cat("  -", f, "\n")