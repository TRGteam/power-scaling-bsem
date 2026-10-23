# ================================================================
# BAYESIAN SEM - ANALISIS EMPIRIS (MODEL BARU) - FINAL
# Konstruk: Socioeconomic Vulnerability, WSH, Reproductive Behavior
# Data: Jawa Tengah (n = 35)
# Strategi Dua-Track:
#   - STRICT power-scaling (utama, 11 model)
#   - PARAMETER-SCALED (tambahan untuk θ, ψ, 4 model)
# Total: 15 model
# ================================================================

rm(list = ls())
set.seed(2026)

library(blavaan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(coda)
library(rstan)

# ================================================================
# KONFIGURASI PARALEL — MAKSIMAL UNTUK 12 PHYSICAL CORES
# ================================================================
n_physical <- parallel::detectCores(logical = FALSE)
cat("[OK] Physical cores:", n_physical, "\n")

# 6 chains per model (setengah dari physical cores)
options(mc.cores = 6)
rstan_options(auto_write = TRUE)
options(future.globals.maxSize = 8 * 1024^3)

cat("  Chains per model:", getOption("mc.cores"), "\n")
cat("  Strategi: sequential per model\n\n")

# ================================================================
# OUTPUT
# ================================================================
setwd("~/SMIBSEMREVISI6")
dir.create("HasilFinal4", showWarnings = FALSE)
sink("HasilFinal4/full_output.txt", split = TRUE)
on.exit(sink(), add = TRUE)

cat("\n========================================================\n")
cat("ANALISIS BAYESIAN SEM - DUA-TRACK POWER-SCALING\n")
cat("Waktu:", date(), "\n")
cat("========================================================\n\n")

# ================================================================
# BAGIAN 1: LOAD DATA
# ================================================================
stunting <- read_excel("Jateng.xlsx", sheet = "Sheet3")

kolom_wajib <- c("x11", "x21", "x31",
                 "x12", "x22", "x32",
                 "y11", "y21", "y31", "Z")

stunting.num <- stunting[, kolom_wajib]
stunting.std <- scale(stunting.num) |> as.data.frame()

cat("\nObservasi:", nrow(stunting.std), "| Indikator:", ncol(stunting.std), "\n")

# ================================================================
# BAGIAN 2: SPESIFIKASI MODEL
# ================================================================
model <- '
  Z ~ beta_z1*eta1
  eta1 ~ gamma11*ksi1 + gamma12*ksi2
  eta1 =~ y11 + y21 + y31
  ksi1 =~ x11 + x21 + x31
  ksi2 =~ x12 + x22 + x32
  ksi1 ~~ phi21*ksi2
'

# ================================================================
# BAGIAN 3: FUNGSI — STRICT POWER-SCALING
# ================================================================
run_model_strict <- function(tau_theta = 1, tau_beta = 1,
                             tau_lambda = 1, tau_psi = 1,
                             sample = 30000, seed = 2026) {
  set.seed(seed)
  
  prior_theta  <- paste0("gamma(", tau_theta * 0 + 1, ", ",
                         0.5 * tau_theta, ")")
  prior_psi    <- paste0("gamma(", 3 * tau_psi + 1, ", ",
                         1 * tau_psi, ")")
  prior_beta   <- paste0("normal(0, ", 0.5 / sqrt(tau_beta), ")")
  prior_lambda <- paste0("normal(1, ", 0.5 / sqrt(tau_lambda), ")")
  
  dp <- dpriors(lambda = prior_lambda, beta = prior_beta,
                psi = prior_psi, theta = prior_theta)
  
  cat("\n[STRICT] τ_θ =", tau_theta, "| τ_β =", tau_beta,
      "| τ_λ =", tau_lambda, "| τ_ψ =", tau_psi, "\n")
  
  fit <- bsem(model, data = stunting.std, dp = dp,
              burnin = 30000, sample = sample, save.lvs = FALSE)
  return(fit)
}

# ================================================================
# BAGIAN 4: FUNGSI — PARAMETER-SCALED
# ================================================================
run_model_ps <- function(tau_theta = 1, tau_beta = 1,
                         tau_lambda = 1, tau_psi = 1,
                         sample = 30000, seed = 2026) {
  set.seed(seed)
  
  prior_theta  <- paste0("gamma(", tau_theta, ", ",
                         0.5 * tau_theta, ")")
  prior_psi    <- paste0("gamma(", 4 * tau_psi, ", ",
                         1 * tau_psi, ")")
  prior_beta   <- paste0("normal(0, ", 0.5 / sqrt(tau_beta), ")")
  prior_lambda <- paste0("normal(1, ", 0.5 / sqrt(tau_lambda), ")")
  
  dp <- dpriors(lambda = prior_lambda, beta = prior_beta,
                psi = prior_psi, theta = prior_theta)
  
  cat("\n[PS] τ_θ =", tau_theta, "| τ_β =", tau_beta,
      "| τ_λ =", tau_lambda, "| τ_ψ =", tau_psi, "\n")
  
  fit <- bsem(model, data = stunting.std, dp = dp,
              burnin = 30000, sample = sample, save.lvs = FALSE)
  return(fit)
}

# ================================================================
# BAGIAN 5: JALANKAN MODEL
# ================================================================
all_fits <- list()

cat("\n--- STRICT: BASELINE ---\n")
all_fits[["strict_baseline"]] <- run_model_strict(
  tau_theta = 1, tau_beta = 1, tau_lambda = 1, tau_psi = 1,
  sample = 100000)

cat("\n--- STRICT: SKENARIO ---\n")
for (t in c(0.5, 0.25, 2.0)) {
  all_fits[[paste0("strict_theta_", t)]] <-
    run_model_strict(tau_theta = t, sample = 30000)
}
for (t in c(0.5, 2.0)) {
  all_fits[[paste0("strict_beta_", t)]] <-
    run_model_strict(tau_beta = t, sample = 30000)
}
for (t in c(0.5, 2.0)) {
  all_fits[[paste0("strict_lambda_", t)]] <-
    run_model_strict(tau_lambda = t, sample = 30000)
}
for (t in c(0.5, 2.0)) {
  all_fits[[paste0("strict_psi_", t)]] <-
    run_model_strict(tau_psi = t, sample = 30000)
}
all_fits[["strict_extreme"]] <- run_model_strict(
  tau_theta = 0.25, tau_beta = 0.5,
  tau_lambda = 0.5, tau_psi = 0.5, sample = 30000)

cat("\n--- PARAMETER-SCALED: TAMBAHAN ---\n")
for (t in c(0.5, 2.0)) {
  all_fits[[paste0("ps_theta_", t)]] <-
    run_model_ps(tau_theta = t, sample = 30000)
}
for (t in c(0.5, 2.0)) {
  all_fits[[paste0("ps_psi_", t)]] <-
    run_model_ps(tau_psi = t, sample = 30000)
}

cat("\n SEMUA MODEL SELESAI!\n")
save(all_fits, stunting.std, model, file = "HasilFinal4/all_fits.RData")

# ================================================================
# BAGIAN 6: VERIFIKASI KONVERGENSI
# ================================================================
cat("\n========================================\n")
cat("VERIFIKASI KONVERGENSI\n")
cat("========================================\n")

conv_summary <- data.frame()
for (name in names(all_fits)) {
  fit <- all_fits[[name]]
  rhat <- tryCatch(blavInspect(fit, "psrf"), error = function(e) NA)
  ess  <- tryCatch(blavInspect(fit, "neff"), error = function(e) NA)
  
  max_rhat <- ifelse(all(is.na(rhat)), NA, max(rhat, na.rm = TRUE))
  min_ess  <- ifelse(all(is.na(ess)), NA, min(ess, na.rm = TRUE))
  
  cat(sprintf("%-22s | Rhat = %.4f | ESS = %.0f\n",
              name, max_rhat, min_ess))
  
  conv_summary <- rbind(conv_summary,
                        data.frame(Model = name,
                                   max_Rhat = max_rhat,
                                   min_ESS = min_ess))
}
write.csv(conv_summary, "HasilFinal4/convergence_summary.csv", row.names = FALSE)

# ================================================================
# BAGIAN 6b: ESS PER PARAMETER
# ================================================================
cat("\n=== ESS PER PARAMETER (BASELINE) ===\n")

ess_detail <- tidy(all_fits[["strict_baseline"]], conf.int = TRUE) %>%
  dplyr::select(term, op, ess, rhat) %>%
  filter(!is.na(ess))

cat("ESS minimum:", round(min(ess_detail$ess, na.rm = TRUE), 1), "\n")
cat("Rhat maximum:", round(max(ess_detail$rhat, na.rm = TRUE), 4), "\n")
cat("Param dengan ESS < 400:", sum(ess_detail$ess < 400, na.rm = TRUE), "\n")

write.csv(ess_detail, "HasilFinal4/ess_per_parameter.csv", row.names = FALSE)

# ================================================================
# BAGIAN 7: EKSTRAKSI KOEFISIEN
# ================================================================
get_coefs <- function(fit) {
  pe <- tidy(fit, conf.int = TRUE)
  pe_sub <- pe[pe$op == "~", c("term","estimate","conf.low","conf.high")]
  target <- c("Z~eta1", "eta1~ksi1", "eta1~ksi2")
  out <- pe_sub[pe_sub$term %in% target, ]
  if (nrow(out) == 0) warning("Tidak ada parameter target.")
  names(out) <- c("Parameter","Estimate","Lower","Upper")
  rownames(out) <- NULL
  return(out)
}

all_results <- data.frame()
for (name in names(all_fits)) {
  df <- get_coefs(all_fits[[name]])
  if (nrow(df) > 0) {
    df$Fit <- name
    all_results <- rbind(all_results, df)
  }
}
write.csv(all_results, "HasilFinal4/all_results.csv", row.names = FALSE)

# ================================================================
# BAGIAN 7b: STANDARDIZED ESTIMATES
# ================================================================
cat("\n=== STANDARDIZED ESTIMATES (BASELINE) ===\n")

pe_std <- tryCatch(
  tidy(all_fits[["strict_baseline"]], conf.int = TRUE, standardized = TRUE),
  error = function(e) tidy(all_fits[["strict_baseline"]], conf.int = TRUE)
)

pe_std_struct <- pe_std %>%
  filter(op == "~", term %in% c("Z~eta1", "eta1~ksi1", "eta1~ksi2")) %>%
  dplyr::select(term, estimate, std.all, conf.low, conf.high)

print(pe_std_struct)
write.csv(pe_std_struct, "HasilFinal4/standardized_estimates.csv",
          row.names = FALSE)

# ================================================================
# BAGIAN 8: EFEK MEDIASI (BASELINE STRICT)
# ================================================================
mcmc_samples <- as.matrix(blavInspect(all_fits[["strict_baseline"]], "mcmc"))

ind1 <- mcmc_samples[, "gamma11"] * mcmc_samples[, "beta_z1"]
ind2 <- mcmc_samples[, "gamma12"] * mcmc_samples[, "beta_z1"]

cat("\n=== EFEK MEDIASI (STRICT BASELINE) ===\n")
cat("\nξ₁ → η₁ → Z:  mean =", round(mean(ind1), 3),
    "| CI = [", round(quantile(ind1, 0.025), 3), ",",
    round(quantile(ind1, 0.975), 3), "]",
    "| P(<0) =", round(mean(ind1 < 0), 3), "\n")
cat("ξ₂ → η₁ → Z:  mean =", round(mean(ind2), 3),
    "| CI = [", round(quantile(ind2, 0.025), 3), ",",
    round(quantile(ind2, 0.975), 3), "]",
    "| P(<0) =", round(mean(ind2 < 0), 3), "\n")

write.csv(data.frame(indirect = ind1),
          "HasilFinal4/indirect_samples_ksi1.csv", row.names = FALSE)
write.csv(data.frame(indirect = ind2),
          "HasilFinal4/indirect_samples_ksi2.csv", row.names = FALSE)

# ================================================================
# BAGIAN 9: R-SQUARED DAN PARAMETER TAMBAHAN
# ================================================================
rsq <- blavInspect(all_fits[["strict_baseline"]], "rsquare")
cat("\nR-Squared:\n"); print(rsq)

pe_full <- tidy(all_fits[["strict_baseline"]], conf.int = TRUE)
phi21 <- pe_full[pe_full$op == "~~" & grepl("ksi", pe_full$term), ]
cat("\nPhi21:\n"); print(phi21)

write.csv(pe_full, "HasilFinal4/parameter_full.csv", row.names = FALSE)

# ================================================================
# BAGIAN 9b: MODEL FIT (PPP, BRMSEA)
# ================================================================
cat("\n=== MODEL FIT (BASELINE STRICT) ===\n")

fit_measures <- tryCatch({
  s <- summary(all_fits[["strict_baseline"]], fit.measures = TRUE)
  s$FIT
}, error = function(e) {
  cat("[WARN] Tidak dapat ekstrak fit measures:", e$message, "\n")
  NULL
})

if (!is.null(fit_measures)) {
  print(fit_measures)
  write.csv(as.data.frame(fit_measures),
            "HasilFinal4/model_fit_strict.csv", row.names = FALSE)
}

# ================================================================
# BAGIAN 10: KLASIFIKASI ROBUSTNESS (STRICT)
# ================================================================
strict_results <- all_results %>% filter(grepl("^strict_", Fit))

classif_strict <- strict_results %>%
  group_by(Parameter) %>%
  summarise(
    n          = n(),
    mean_est   = round(mean(Estimate), 3),
    min_est    = round(min(Estimate), 3),
    max_est    = round(max(Estimate), 3),
    range_est  = round(max(Estimate) - min(Estimate), 3),
    all_ci_neg = all(Upper < 0),
    all_ci_pos = all(Lower > 0),
    status = case_when(
      all_ci_neg ~ "ROBUST (Negatif)",
      all_ci_pos ~ "ROBUST (Positif)",
      TRUE       ~ "Sensitive"
    ),
    .groups = "drop"
  )

print(classif_strict)
write.csv(classif_strict, "HasilFinal4/classification_strict.csv",
          row.names = FALSE)

# ================================================================
# BAGIAN 11: PERBANDINGAN STRICT vs PS
# ================================================================
comparison <- data.frame()
for (param in c("theta", "psi")) {
  for (t in c(0.5, 2.0)) {
    strict_name <- paste0("strict_", param, "_", t)
    ps_name     <- paste0("ps_", param, "_", t)
    
    if (strict_name %in% names(all_fits) && ps_name %in% names(all_fits)) {
      strict_df <- get_coefs(all_fits[[strict_name]])
      ps_df     <- get_coefs(all_fits[[ps_name]])
      
      for (p in unique(strict_df$Parameter)) {
        comparison <- rbind(comparison, data.frame(
          Kelas    = param,
          Tau      = t,
          Parameter = p,
          Strict   = round(strict_df$Estimate[strict_df$Parameter == p], 4),
          PS       = round(ps_df$Estimate[ps_df$Parameter == p], 4),
          Diff     = round(ps_df$Estimate[ps_df$Parameter == p] -
                             strict_df$Estimate[strict_df$Parameter == p], 4)
        ))
      }
    }
  }
}

print(comparison)
write.csv(comparison, "HasilFinal4/comparison_strict_ps.csv", row.names = FALSE)

# ================================================================
# BAGIAN 11b: MODEL FIT COMPARISON (STRICT vs PS)
# ================================================================
cat("\n=== MODEL FIT COMPARISON ===\n")

get_fit <- function(fit) {
  tryCatch({
    s <- summary(fit, fit.measures = TRUE)
    # s$FIT adalah vektor, bukan data.frame
    ppp_val <- ifelse("PPP" %in% names(s$FIT), 
                      as.numeric(s$FIT["PPP"]), NA)
    # BRMSEA, BCFI, BTLI mungkin tidak tersedia di semua versi blavaan
    brmsea_val <- ifelse("BRMSEA" %in% names(s$FIT), 
                         as.numeric(s$FIT["BRMSEA"]), NA)
    
    data.frame(
      PPP    = ppp_val,
      BRMSEA = brmsea_val
    )
  }, error = function(e) {
    data.frame(PPP = NA, BRMSEA = NA)
  })
}

fit_comp <- data.frame()
for (name in c("strict_baseline", "ps_theta_0.5", "ps_theta_2",
               "ps_psi_0.5", "ps_psi_2")) {
  if (name %in% names(all_fits)) {
    fit_row <- get_fit(all_fits[[name]])
    fit_row$Model <- name
    fit_comp <- rbind(fit_comp, fit_row)
  }
}

print(fit_comp)
write.csv(fit_comp, "HasilFinal4/model_fit_comparison.csv", row.names = FALSE)

# ================================================================
# BAGIAN 12: PLOTTING — SENSITIVITY (STRICT)
# ================================================================
param_labels <- c(
  "Z~eta1"    = "beta[z1]~(eta1%->%Z)",
  "eta1~ksi1" = "gamma[11]~(ksi1%->%eta1)",
  "eta1~ksi2" = "gamma[12]~(ksi2%->%eta1)"
)

strict_results <- strict_results %>%
  mutate(
    TauClass = case_when(
      grepl("theta_",  Fit) ~ "theta",
      grepl("beta_",   Fit) ~ "beta",
      grepl("lambda_", Fit) ~ "lambda",
      grepl("psi_",    Fit) ~ "psi",
      Fit == "strict_baseline" ~ "baseline",
      Fit == "strict_extreme"  ~ "extreme",
      TRUE ~ NA_character_
    ),
    TauValue = suppressWarnings(as.numeric(gsub("[^0-9.]", "", Fit)))
  )
strict_results$TauValue[strict_results$Fit == "strict_baseline"] <- 1

baseline_row <- strict_results %>% filter(Fit == "strict_baseline")
baseline_exp <- do.call(rbind, lapply(c("theta","beta","lambda","psi"),
                                      function(tc) {
                                        baseline_row %>% mutate(TauClass = tc, TauValue = 1)
                                      }))

plot_data <- strict_results %>%
  filter(TauClass %in% c("theta","beta","lambda","psi")) %>%
  bind_rows(baseline_exp) %>%
  mutate(TauClass = factor(TauClass,
                           levels = c("theta","beta","lambda","psi")))

# --- Figure 7: Sensitivity Plot ---
p_sens <- ggplot(plot_data, aes(x = TauValue, y = Estimate,
                                color = Parameter, fill = Parameter)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_grid(Parameter ~ TauClass, scales = "free",
             labeller = labeller(Parameter = as_labeller(param_labels,
                                                         label_parsed))) +
  scale_x_reverse() +
  labs(x = expression(tau~"(power-scaling factor)"),
       y = "Posterior estimate",
       title = "Prior sensitivity (strict power-scaling)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 9))

ggsave("HasilFinal4/Figure7_Sensitivity_Strict.png", p_sens,
       width = 12, height = 10, dpi = 300)
cat("[OK] Figure7_Sensitivity_Strict.png\n")

# --- Figure 10: Comparison Strict vs PS ---
comparison_long <- comparison %>%
  pivot_longer(cols = c("Strict", "PS"),
               names_to = "Formulation", values_to = "Estimate")

p_comp <- ggplot(comparison_long,
                 aes(x = factor(Tau), y = Estimate, fill = Formulation)) +
  geom_col(position = position_dodge(0.8), alpha = 0.85) +
  facet_grid(Parameter ~ Kelas, scales = "free_y",
             labeller = labeller(Parameter = as_labeller(param_labels,
                                                         label_parsed))) +
  scale_fill_manual(values = c("Strict" = "#4C72B0", "PS" = "#DD8452")) +
  labs(x = expression(tau), y = "Posterior estimate",
       title = "Strict vs Parameter-Scaled formulation") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

ggsave("HasilFinal4/Figure10_Comparison_Strict_PS.png", p_comp,
       width = 12, height = 8, dpi = 300)
cat("[OK] Figure10_Comparison_Strict_PS.png\n")

# ================================================================
# BAGIAN 13: FIGURE 6 — POSTERIOR INDIRECT EFFECT
# ================================================================
p_ind2 <- ggplot(data.frame(indirect = ind2),
                 aes(x = indirect)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "steelblue", color = "white", alpha = 0.7) +
  geom_density(color = "darkblue", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "red", linewidth = 0.8) +
  geom_vline(xintercept = mean(ind2), linetype = "solid",
             color = "darkgreen", linewidth = 1) +
  annotate("text", x = mean(ind2), y = Inf,
           label = paste0("Mean = ", round(mean(ind2), 3)),
           vjust = 2, hjust = -0.1, color = "darkgreen", size = 4) +
  annotate("text", x = 0, y = Inf,
           label = "0", vjust = 2, hjust = -0.1, color = "red", size = 4) +
  labs(x = expression("Indirect path coefficient ("*
                        xi[2] %->% eta[1] %->% Z*")"),
       y = "Density",
       title = "Posterior distribution of indirect path") +
  theme_bw(base_size = 11)

ggsave("HasilFinal4/Figure6_Posterior_Indirect.png", p_ind2,
       width = 8, height = 5, dpi = 300)
cat("[OK] Figure6_Posterior_Indirect.png\n")

# ================================================================
# BAGIAN 14: DELTA PLOT
# ================================================================
cat("\n========================================\n")
cat("DELTA PLOT vs BASELINE\n")
cat("========================================\n")

baseline_df <- strict_results %>%
  filter(Fit == "strict_baseline") %>%
  dplyr::select(Parameter,
                Baseline_Est   = Estimate,
                Baseline_Lower = Lower,
                Baseline_Upper = Upper)

delta_df <- strict_results %>%
  filter(Fit != "strict_baseline", Fit != "strict_extreme") %>%
  left_join(baseline_df, by = "Parameter") %>%
  mutate(
    Delta_Est   = Estimate - Baseline_Est,
    Delta_Lower = Lower    - Baseline_Lower,
    Delta_Upper = Upper    - Baseline_Upper,
    TauClass = case_when(
      grepl("strict_theta_",  Fit) ~ "theta",
      grepl("strict_beta_",   Fit) ~ "beta",
      grepl("strict_lambda_", Fit) ~ "lambda",
      grepl("strict_psi_",    Fit) ~ "psi",
      TRUE ~ NA_character_
    ),
    TauValue = suppressWarnings(as.numeric(gsub("[^0-9.]", "", Fit)))
  ) %>%
  filter(!is.na(TauClass))

delta_summary <- delta_df %>%
  group_by(Parameter, TauClass) %>%
  summarise(
    Delta_min    = round(min(Delta_Est), 4),
    Delta_max    = round(max(Delta_Est), 4),
    Delta_range  = round(max(Delta_Est) - min(Delta_Est), 4),
    CI_shift_max = round(max(abs(Delta_Upper)), 4),
    .groups = "drop"
  )

print(delta_summary)
write.csv(delta_summary, "HasilFinal4/delta_summary.csv", row.names = FALSE)

# --- Figure 8: Delta Absolute ---
plot_delta <- delta_df %>%
  mutate(TauClass = factor(TauClass,
                           levels = c("theta", "beta", "lambda", "psi")))

p_delta <- ggplot(plot_delta,
                  aes(x = factor(TauValue), y = Delta_Est,
                      fill = Parameter)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  geom_col(position = position_dodge(0.8), alpha = 0.85) +
  geom_errorbar(aes(ymin = Delta_Lower, ymax = Delta_Upper),
                position = position_dodge(0.8),
                width = 0.25, linewidth = 0.7) +
  facet_grid(Parameter ~ TauClass, scales = "free_y",
             labeller = labeller(Parameter = as_labeller(param_labels,
                                                         label_parsed))) +
  scale_fill_manual(values = c("Z~eta1"    = "#E41A1C",
                               "eta1~ksi1" = "#377EB8",
                               "eta1~ksi2" = "#4DAF4A"),
                    guide = "none") +
  labs(x = expression(tau~"(power-scaling factor)"),
       y = expression(Delta~"estimate vs baseline"),
       title = "Shift in posterior mean and 95% CI vs baseline") +
  theme_bw(base_size = 10) +
  theme(strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 9),
        panel.grid.minor = element_blank())

ggsave("HasilFinal4/Figure8_Delta_Absolute.png", p_delta,
       width = 12, height = 9, dpi = 300)
cat("[OK] Figure8_Delta_Absolute.png\n")

# --- Figure 9: Delta Range ---
p_range <- ggplot(delta_summary,
                  aes(x = TauClass, y = Delta_range, fill = TauClass)) +
  geom_col(alpha = 0.85) +
  facet_wrap(~ Parameter,
             labeller = as_labeller(param_labels, label_parsed),
             nrow = 1) +
  scale_fill_manual(values = c("theta"  = "#4C72B0",
                               "beta"   = "#DD8452",
                               "lambda" = "#55A868",
                               "psi"    = "#C44E52"),
                    guide = "none") +
  labs(x = "Power-scaling class",
       y = expression("Range of "*Delta~"estimate"),
       title = "Magnitude of prior sensitivity by parameter class") +
  theme_bw(base_size = 11)

ggsave("HasilFinal4/Figure9_Delta_Range.png", p_range,
       width = 10, height = 4, dpi = 300)
cat("[OK] Figure9_Delta_Range.png\n")

# ================================================================
# BAGIAN 14.7: TABEL 5b — DELTA RANGE SUMMARY
# ================================================================
cat("\n=== TABEL 5b: DELTA RANGE ===\n")

table_5b <- delta_summary %>%
  dplyr::select(Parameter, TauClass, Delta_range) %>%
  pivot_wider(names_from = TauClass,
              values_from = Delta_range,
              names_prefix = "range_") %>%
  mutate(across(starts_with("range_"), ~ round(.x, 4)))

print(table_5b)
write.csv(table_5b, "HasilFinal4/table_5b_delta_summary.csv",
          row.names = FALSE)

# ================================================================
# RINGKASAN AKHIR
# ================================================================
cat("\n========================================\n")
cat("ANALISIS SELESAI!\n")
cat("Output: HasilFinal4/\n")
cat("Catatan: Figure 1 & 5 dibuat manual.\n")
cat("========================================\n")

cat("\n=== FILE OUTPUT ===\n")
files <- list.files("HasilFinal4", pattern = "\\.(csv|png)$")
for (f in files) cat("  -", f, "\n")