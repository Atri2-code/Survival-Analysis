# =============================================================================
# Advanced Survival Analysis Pipeline
# =============================================================================
# Comprehensive time-to-event analysis pipeline for patient-level RCT data
# implementing methods used in NICE health technology assessment submissions.
#
# Methods implemented:
#   1. Kaplan-Meier analysis with log-rank test
#   2. Cox proportional hazards (univariable + multivariable)
#   3. Proportional hazards assumption testing (Schoenfeld residuals)
#   4. Parametric survival models (exponential, Weibull, log-normal,
#      log-logistic, Gompertz, generalised gamma)
#   5. Model selection via AIC/BIC
#   6. Restricted Mean Survival Time (RMST)
#   7. Landmark analysis
#   8. Subgroup forest plot
#   9. Extrapolation plots for HTA submissions
#  10. Statistical Analysis Plan (SAP) document
#
# Author: Atrija Haldar
# =============================================================================

# ── 0. Setup ──────────────────────────────────────────────────────────────────

required_packages <- c(
  "survival", "survminer", "ggplot2", "dplyr", "tidyr",
  "flexsurv", "broom", "gridExtra", "knitr", "scales",
  "survRM2", "ggfortify"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
}

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(flexsurv)
library(broom)
library(gridExtra)
library(scales)
library(survRM2)

set.seed(42)
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)

cat("Advanced Survival Analysis Pipeline\n")
cat(strrep("=", 50), "\n\n")


# ── 1. Simulate patient-level RCT data ────────────────────────────────────────

simulate_trial <- function(
    n           = 500,
    hr_os       = 0.65,
    hr_pfs      = 0.55,
    median_os_c = 14,
    median_pfs_c = 6,
    dropout_rate = 0.08
) {
  cat("Simulating trial data...\n")

  n_arm <- n / 2
  lambda_os_c  <- log(2) / median_os_c
  lambda_os_t  <- lambda_os_c  * hr_os
  lambda_pfs_c <- log(2) / median_pfs_c
  lambda_pfs_t <- lambda_pfs_c * hr_pfs

  # Overall survival times
  os_c <- rexp(n_arm, lambda_os_c)
  os_t <- rexp(n_arm, lambda_os_t)

  # PFS times (always <= OS)
  pfs_c <- pmin(rexp(n_arm, lambda_pfs_c), os_c)
  pfs_t <- pmin(rexp(n_arm, lambda_pfs_t), os_t)

  # Administrative censoring at 36 months
  admin_censor <- 36

  # Random dropout
  dropout_c <- rexp(n_arm, rate = dropout_rate)
  dropout_t <- rexp(n_arm, rate = dropout_rate)

  censor_c <- pmin(dropout_c, admin_censor)
  censor_t <- pmin(dropout_t, admin_censor)

  data.frame(
    patient_id   = 1:n,
    arm          = c(rep("Control", n_arm), rep("Treatment", n_arm)),
    os_time      = c(pmin(os_c, censor_c), pmin(os_t, censor_t)),
    os_event     = c(as.integer(os_c <= censor_c),
                     as.integer(os_t <= censor_t)),
    pfs_time     = c(pmin(pfs_c, censor_c), pmin(pfs_t, censor_t)),
    pfs_event    = c(as.integer(pfs_c <= censor_c),
                     as.integer(pfs_t <= censor_t)),
    age          = round(rnorm(n, mean = 63, sd = 9)),
    ecog_ps      = sample(0:2, n, replace = TRUE,
                          prob = c(0.40, 0.42, 0.18)),
    prior_lines  = sample(1:3, n, replace = TRUE,
                          prob = c(0.48, 0.36, 0.16)),
    pdl1_pos     = sample(c(TRUE, FALSE), n, replace = TRUE,
                          prob = c(0.45, 0.55)),
    biomarker    = rnorm(n, mean = 0, sd = 1),
    weight_kg    = round(rnorm(n, mean = 72, sd = 14))
  )
}

trial <- simulate_trial(n = 500)
write.csv(trial, file.path(output_dir, "trial_data.csv"), row.names = FALSE)

cat(sprintf("  N = %d patients (%d events OS, %d events PFS)\n\n",
            nrow(trial), sum(trial$os_event), sum(trial$pfs_event)))


# ── 2. Baseline characteristics table ─────────────────────────────────────────

cat("Generating baseline characteristics table...\n")

baseline <- trial %>%
  group_by(arm) %>%
  summarise(
    n               = n(),
    median_age      = round(median(age), 1),
    sd_age          = round(sd(age), 1),
    ecog_0_pct      = round(mean(ecog_ps == 0) * 100, 1),
    ecog_1_pct      = round(mean(ecog_ps == 1) * 100, 1),
    ecog_2_pct      = round(mean(ecog_ps == 2) * 100, 1),
    pdl1_pos_pct    = round(mean(pdl1_pos) * 100, 1),
    prior_1_pct     = round(mean(prior_lines == 1) * 100, 1),
    os_events_pct   = round(mean(os_event) * 100, 1),
    pfs_events_pct  = round(mean(pfs_event) * 100, 1),
    .groups = "drop"
  )

write.csv(baseline, file.path(output_dir, "baseline_table.csv"),
          row.names = FALSE)
cat("  Saved: output/baseline_table.csv\n\n")


# ── 3. Kaplan-Meier analysis ──────────────────────────────────────────────────

cat("Fitting Kaplan-Meier curves...\n")

# OS
km_os  <- survfit(Surv(os_time, os_event) ~ arm, data = trial)
# PFS
km_pfs <- survfit(Surv(pfs_time, pfs_event) ~ arm, data = trial)

# Log-rank tests
lr_os  <- survdiff(Surv(os_time, os_event) ~ arm, data = trial)
lr_pfs <- survdiff(Surv(pfs_time, pfs_event) ~ arm, data = trial)

p_os  <- round(1 - pchisq(lr_os$chisq, df = 1), 4)
p_pfs <- round(1 - pchisq(lr_pfs$chisq, df = 1), 4)

cat(sprintf("  OS  log-rank p = %.4f\n", p_os))
cat(sprintf("  PFS log-rank p = %.4f\n\n", p_pfs))

# KM plots
km_os_plot <- ggsurvplot(
  km_os, data = trial,
  pval = TRUE, pval.method = TRUE,
  conf.int = TRUE, risk.table = TRUE,
  surv.median.line = "hv",
  palette = c("#1F4E79", "#D85A30"),
  title = "Overall Survival by Treatment Arm",
  xlab = "Time (months)", ylab = "Survival probability",
  legend.title = "Arm",
  legend.labs = c("Control", "Treatment"),
  xlim = c(0, 36), break.time.by = 6,
  ggtheme = theme_bw(base_size = 12)
)

km_pfs_plot <- ggsurvplot(
  km_pfs, data = trial,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE,
  surv.median.line = "hv",
  palette = c("#1F4E79", "#D85A30"),
  title = "Progression-Free Survival by Treatment Arm",
  xlab = "Time (months)", ylab = "Survival probability",
  legend.title = "Arm",
  legend.labs = c("Control", "Treatment"),
  xlim = c(0, 36), break.time.by = 6,
  ggtheme = theme_bw(base_size = 12)
)

ggsave(file.path(output_dir, "km_os.png"),
       print(km_os_plot), width = 10, height = 8, dpi = 300)
ggsave(file.path(output_dir, "km_pfs.png"),
       print(km_pfs_plot), width = 10, height = 8, dpi = 300)
cat("  Saved: output/km_os.png, output/km_pfs.png\n\n")


# ── 4. Cox proportional hazards ───────────────────────────────────────────────

cat("Fitting Cox proportional hazards models...\n")

# Univariable
cox_uni <- coxph(Surv(os_time, os_event) ~ arm, data = trial)

# Multivariable
cox_multi <- coxph(
  Surv(os_time, os_event) ~ arm + age + factor(ecog_ps) +
    prior_lines + pdl1_pos + biomarker,
  data = trial
)

uni_tidy   <- tidy(cox_uni,   exponentiate = TRUE, conf.int = TRUE)
multi_tidy <- tidy(cox_multi, exponentiate = TRUE, conf.int = TRUE)

hr_uni   <- uni_tidy[uni_tidy$term == "armTreatment", ]
hr_multi <- multi_tidy[multi_tidy$term == "armTreatment", ]

cat(sprintf("  Unadjusted HR:  %.3f (95%% CI: %.3f-%.3f), p=%.4f\n",
            hr_uni$estimate, hr_uni$conf.low,
            hr_uni$conf.high, hr_uni$p.value))
cat(sprintf("  Adjusted HR:    %.3f (95%% CI: %.3f-%.3f), p=%.4f\n\n",
            hr_multi$estimate, hr_multi$conf.low,
            hr_multi$conf.high, hr_multi$p.value))

write.csv(multi_tidy,
          file.path(output_dir, "cox_results.csv"),
          row.names = FALSE)


# ── 5. PH assumption testing ──────────────────────────────────────────────────

cat("Testing proportional hazards assumption...\n")

ph_test <- cox.zph(cox_multi)
cat("\n  Schoenfeld residuals:\n")
print(ph_test$table)

ph_violated <- any(ph_test$table[, "p"] < 0.05, na.rm = TRUE)
cat(sprintf("\n  PH assumption: %s\n\n",
            ifelse(!ph_violated,
                   "Not violated (p > 0.05 for all covariates)",
                   "Potentially violated for some covariates")))

png(file.path(output_dir, "ph_test.png"),
    width = 1400, height = 900, res = 150)
par(mfrow = c(2, 3))
plot(ph_test)
dev.off()
cat("  Saved: output/ph_test.png\n\n")


# ── 6. Parametric survival models ─────────────────────────────────────────────

cat("Fitting parametric survival models (for HTA extrapolation)...\n")

distributions <- c("exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma")
dist_names    <- c("Exponential", "Weibull", "Log-normal",
                   "Log-logistic", "Gompertz", "Gen. Gamma")

param_fits  <- list()
param_stats <- data.frame()

for (i in seq_along(distributions)) {
  d <- distributions[i]
  tryCatch({
    fit <- flexsurvreg(
      Surv(os_time, os_event) ~ arm,
      data = trial, dist = d
    )
    param_fits[[d]] <- fit
    param_stats <- rbind(param_stats, data.frame(
      Distribution = dist_names[i],
      AIC          = round(AIC(fit), 2),
      BIC          = round(BIC(fit), 2),
      LogLik       = round(fit$loglik, 2)
    ))
    cat(sprintf("  %-15s AIC=%.1f  BIC=%.1f\n",
                dist_names[i], AIC(fit), BIC(fit)))
  }, error = function(e) {
    cat(sprintf("  %-15s FAILED: %s\n", dist_names[i], e$message))
  })
}

param_stats <- param_stats %>% arrange(AIC)
best_dist   <- param_stats$Distribution[1]
cat(sprintf("\n  Best fit (lowest AIC): %s\n\n", best_dist))

write.csv(param_stats,
          file.path(output_dir, "parametric_model_fit.csv"),
          row.names = FALSE)

# Extrapolation plot — key HTA output
if (length(param_fits) > 0) {
  png(file.path(output_dir, "survival_extrapolation.png"),
      width = 1400, height = 800, res = 150)

  plot(km_os$strata[[1]], conf.int = FALSE,
       col = "#1F4E79", lwd = 2,
       xlim = c(0, 60), ylim = c(0, 1),
       xlab = "Time (months)",
       ylab = "Overall Survival",
       main = "Parametric Model Extrapolation (Control Arm)",
       las = 1)

  colours <- c("#E41A1C", "#377EB8", "#4DAF4A",
               "#984EA3", "#FF7F00", "#A65628")

  for (i in seq_along(param_fits)) {
    d   <- names(param_fits)[i]
    fit <- param_fits[[d]]
    tryCatch({
      t_seq <- seq(0, 60, by = 0.5)
      s_hat <- summary(fit, t = t_seq, type = "survival",
                       newdata = data.frame(arm = "Control"))
      lines(s_hat[[1]]$time, s_hat[[1]]$est,
            col = colours[i], lwd = 1.5, lty = 2)
    }, error = function(e) NULL)
  }

  legend("topright",
         legend = c("KM (observed)", dist_names[seq_along(param_fits)]),
         col    = c("#1F4E79", colours[seq_along(param_fits)]),
         lwd    = c(2, rep(1.5, length(param_fits))),
         lty    = c(1, rep(2, length(param_fits))),
         cex    = 0.8)
  dev.off()
  cat("  Saved: output/survival_extrapolation.png\n\n")
}


# ── 7. Restricted Mean Survival Time (RMST) ───────────────────────────────────

cat("Computing Restricted Mean Survival Time (RMST)...\n")

rmst_result <- tryCatch({
  rmst2(
    time   = trial$os_time,
    status = trial$os_event,
    arm    = as.integer(trial$arm == "Treatment"),
    tau    = 24  # Restrict to 24 months
  )
}, error = function(e) NULL)

if (!is.null(rmst_result)) {
  cat(sprintf("  RMST at 24 months:\n"))
  cat(sprintf("    Control:   %.2f months\n",
              rmst_result$RMST.arm0$rmst[1]))
  cat(sprintf("    Treatment: %.2f months\n",
              rmst_result$RMST.arm1$rmst[1]))
  cat(sprintf("    Difference: %.2f months (p=%.4f)\n\n",
              rmst_result$unadjusted.result[1, 1],
              rmst_result$unadjusted.result[1, 4]))
}


# ── 8. Landmark analysis ──────────────────────────────────────────────────────

cat("Running landmark analysis (6-month landmark)...\n")

landmark_time <- 6
landmark_data <- trial %>%
  filter(os_time >= landmark_time) %>%
  mutate(
    os_time_lm  = os_time - landmark_time,
    os_event_lm = os_event
  )

km_landmark <- survfit(
  Surv(os_time_lm, os_event_lm) ~ arm,
  data = landmark_data
)

lr_landmark <- survdiff(
  Surv(os_time_lm, os_event_lm) ~ arm,
  data = landmark_data
)
p_landmark <- round(1 - pchisq(lr_landmark$chisq, df = 1), 4)

cat(sprintf("  Landmark analysis (t=%d): n=%d, log-rank p=%.4f\n\n",
            landmark_time, nrow(landmark_data), p_landmark))

lm_plot <- ggsurvplot(
  km_landmark, data = landmark_data,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE,
  palette = c("#1F4E79", "#D85A30"),
  title = sprintf("Landmark Analysis: OS from Month %d", landmark_time),
  xlab = sprintf("Time from %d-month landmark (months)", landmark_time),
  ylab = "Survival probability",
  ggtheme = theme_bw(base_size = 12)
)

ggsave(file.path(output_dir, "landmark_analysis.png"),
       print(lm_plot), width = 10, height = 8, dpi = 300)
cat("  Saved: output/landmark_analysis.png\n\n")


# ── 9. Subgroup forest plot ───────────────────────────────────────────────────

cat("Running subgroup analysis...\n")

subgroups <- list(
  "All patients"       = trial,
  "Age < 65"           = trial[trial$age < 65, ],
  "Age \u2265 65"      = trial[trial$age >= 65, ],
  "ECOG PS 0"          = trial[trial$ecog_ps == 0, ],
  "ECOG PS 1-2"        = trial[trial$ecog_ps >= 1, ],
  "PD-L1 positive"     = trial[trial$pdl1_pos == TRUE, ],
  "PD-L1 negative"     = trial[trial$pdl1_pos == FALSE, ],
  "1 prior line"       = trial[trial$prior_lines == 1, ],
  "2+ prior lines"     = trial[trial$prior_lines >= 2, ]
)

forest_data <- do.call(rbind, lapply(names(subgroups), function(sg) {
  d <- subgroups[[sg]]
  if (nrow(d) < 15 || sum(d$os_event) < 5) return(NULL)
  fit <- tryCatch(
    coxph(Surv(os_time, os_event) ~ arm, data = d),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  r <- tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  r <- r[r$term == "armTreatment", ]
  data.frame(
    subgroup = sg,
    n        = nrow(d),
    events   = sum(d$os_event),
    hr       = round(r$estimate, 2),
    lower    = round(r$conf.low, 2),
    upper    = round(r$conf.high, 2),
    p_value  = round(r$p.value, 3)
  )
}))

forest_data$subgroup <- factor(forest_data$subgroup,
                                levels = rev(forest_data$subgroup))
forest_data$label    <- sprintf("%.2f (%.2f\u2013%.2f)",
                                  forest_data$hr,
                                  forest_data$lower,
                                  forest_data$upper)

fp <- ggplot(forest_data,
             aes(x = hr, y = subgroup, xmin = lower, xmax = upper)) +
  geom_point(size = 3, colour = "#1F4E79") +
  geom_errorbarh(height = 0.3, colour = "#1F4E79") +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey50", linewidth = 0.8) +
  geom_text(aes(label = label, x = 2.8),
            size = 3.2, hjust = 0) +
  geom_text(aes(label = paste0("n=", n), x = 0.18),
            size = 3, hjust = 0, colour = "grey40") +
  scale_x_continuous(limits = c(0.15, 3.5),
                     trans  = "log",
                     breaks = c(0.25, 0.5, 1, 2)) +
  labs(
    title   = "Subgroup Analysis: Hazard Ratio (Treatment vs Control)",
    x       = "Hazard Ratio (log scale)",
    y       = NULL,
    caption = "HR < 1 favours treatment. Error bars = 95% CI."
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(output_dir, "forest_plot.png"),
       fp, width = 12, height = 7, dpi = 300)

write.csv(forest_data,
          file.path(output_dir, "subgroup_analysis.csv"),
          row.names = FALSE)
cat("  Saved: output/forest_plot.png\n\n")


# ── 10. Statistical Analysis Plan ─────────────────────────────────────────────

source("generate_sap.R")


# ── Summary ───────────────────────────────────────────────────────────────────

cat(strrep("=", 55), "\n", sep = "")
cat("  ANALYSIS SUMMARY\n")
cat(strrep("=", 55), "\n", sep = "")
cat(sprintf("  Population:         %d patients (%d OS events)\n",
            nrow(trial), sum(trial$os_event)))
cat(sprintf("  OS log-rank p:      %.4f\n", p_os))
cat(sprintf("  Unadjusted HR:      %.3f (%.3f\u2013%.3f)\n",
            hr_uni$estimate, hr_uni$conf.low, hr_uni$conf.high))
cat(sprintf("  Adjusted HR:        %.3f (%.3f\u2013%.3f)\n",
            hr_multi$estimate, hr_multi$conf.low, hr_multi$conf.high))
cat(sprintf("  Best parametric fit: %s\n", best_dist))
cat(sprintf("  PH assumption:      %s\n",
            ifelse(!ph_violated, "Not violated", "Review required")))
cat(strrep("=", 55), "\n\n", sep = "")
cat("All outputs saved to /output\n")
cat("Analysis complete.\n")
