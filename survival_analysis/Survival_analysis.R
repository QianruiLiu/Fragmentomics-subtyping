#!/usr/bin/env Rscript

# ==== Dependencies ====
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(survival)
  library(survminer)
  library(readr)
  library(janitor)
})

# ==== CLI options ====
option_list <- list(
  make_option("--dev_clin",  type = "character", default = "lpWGS_clinical_ubc.txt",
              help = "Development cohort: clinical file (must include OS time 'os', event 'dead', and a ctF column). Default: %default"),
  make_option("--dev_pred",  type = "character", default = "combined_LLR_softTF_predictions.tsv",
              help = "Development cohort: prediction file (columns: sample, .pred_NEPC). Default: %default"),
  make_option("--clinic_pred", type = "character", default = "prediction_clinicaldata.txt",
              help = "Deployment cohort: prediction file (must include sample_id and p_nepc/.pred_nepc). Default: %default"),
  make_option("--clinic_tf", type = "character", default = "clinic_TF.txt",
              help = "Deployment cohort: ctF file (must include sample_id and a ct_fraction synonym). Default: %default"),
  make_option("--times", type = "character", default = "12,24,36",
              help = "Time points for absolute survival projection (same unit as 'os'), comma-separated. Default: %default"),
  make_option("--out_prefix", type = "character", default = "",
              help = "Optional output filename prefix. Default: empty"),
  make_option("--seed", type = "integer", default = 42,
              help = "Random seed. Default: %default")
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

# Helper for prefixed outputs
out_path <- function(name){
  if (is.null(opt$out_prefix) || opt$out_prefix == "") return(name)
  dir.create(dirname(file.path(opt$out_prefix)), showWarnings = FALSE, recursive = TRUE)
  paste0(opt$out_prefix, "_", name)
}

# ==== A) Development cohort: ctF-aware probability calibration (unlabeled residualization) + Cox ====

# 1) Load development cohort clinical and predictions
# Clinical must contain: sample_id (or sample), OS time 'os', event 'dead', and a ctF column (see synonyms below).
clin <- read.delim(opt$dev_clin, sep = "", quote = "\"", check.names = FALSE)
if ("sample" %in% names(clin) && !("sample_id" %in% names(clin))) {
  clin <- clin %>% rename(sample_id = sample)
}

pred <- readr::read_tsv(opt$dev_pred,
                        col_types = readr::cols(
                          sample = readr::col_character(),
                          .pred_ARPC = readr::col_double(),
                          .pred_NEPC = readr::col_double()
                        )) %>%
  transmute(sample_id = sample, p_NEPC = `.pred_NEPC`)

# Auto-detect a ctF column (first match wins)
cand <- c("est_ctfrac_targeted","ct_fraction","ichor_tfx","tfx","TF","tf","tumour_fraction","tumor_fraction")
ct_col <- intersect(cand, names(clin))
stopifnot("No ctF column found in development clinical file" = length(ct_col) >= 1)
clin <- clin %>% mutate(ct_fraction = .data[[ct_col[1]]])

# Merge and build logits
dat <- clin %>%
  inner_join(pred, by = "sample_id") %>%
  select(sample_id, os, dead, p_NEPC, ct_fraction) %>%
  drop_na(os, dead, p_NEPC, ct_fraction) %>%
  mutate(
    p_NEPC = pmin(pmax(p_NEPC, 1e-6), 1-1e-6),
    logit_p = qlogis(p_NEPC)
  )

# 2) Unlabeled residualization: logit(p_NEPC) ~ ct_fraction
m_resid <- lm(logit_p ~ ct_fraction, data = dat)

# Shift residuals to preserve cohort-wide mean, then inverse-logit to calibrated probability
dat <- dat %>%
  mutate(
    logit_resid = resid(m_resid),
    intercept_shift = mean(logit_p) - mean(logit_resid),
    logit_adj = logit_resid + intercept_shift,
    p_NEPC_cal = plogis(logit_adj),
    p_NEPC_cal_10 = p_NEPC_cal * 10
  )

# 3) Cox PH (univariable): main predictor is p_NEPC_cal × 10
fit_cox <- coxph(Surv(os, dead) ~ p_NEPC_cal_10, data = dat)
print(summary(fit_cox))
print(cox.zph(fit_cox))  # PH assumption test

# Extract baseline cumulative hazard and beta for deployment projections
bh <- basehaz(fit_cox, centered = FALSE)
beta <- coef(fit_cox)[["p_NEPC_cal_10"]]

# 4) Visualization only: KM by top 10% of p_NEPC_cal (High10 vs Other)
thr_cal <- quantile(dat$p_NEPC_cal, 0.90, na.rm = TRUE)
dat <- dat %>%
  mutate(lpWGS_top10 = factor(if_else(p_NEPC_cal >= thr_cal, "High10", "Other"),
                              levels = c("Other","High10")))

fit_km_top10 <- survfit(Surv(os, dead) ~ lpWGS_top10, data = dat)

p_top10 <- ggsurvplot(
  fit_km_top10,
  data = dat,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "Group",
  legend.labs = c("Other (≤90th pct)", "Top10% (≥90th pct)"),
  title = sprintf("OS by p_NEPC_cal Top 10%% (threshold = %.3f)", thr_cal),
  xlab = "Time",
  ylab = "Survival probability"
)

ggsave(out_path("OS_by_pNEPCcal_Top10_lpWGS.png"), p_top10$plot, width = 7, height = 5, dpi = 300)
ggsave(out_path("OS_by_pNEPCcal_Top10_lpWGS_risktable.png"), p_top10$table, width = 7, height = 2.8, dpi = 300)

# ==== B) Deployment cohort: individualized absolute survival probabilities ====

# 5) Load deployment cohort (predictions + ctF), minimally harmonize column names
clin_pred <- readr::read_tsv(opt$clinic_pred, show_col_types = FALSE) %>% clean_names()
if (!"sample_id" %in% names(clin_pred)) {
  id_hits <- intersect(c("sample_id","sample","id"), names(clin_pred))
  stopifnot("Deployment prediction file lacks sample_id/sample/id" = length(id_hits) >= 1)
  clin_pred <- clin_pred %>% rename(sample_id = !!dplyr::all_of(id_hits[1]))
}
if (!"p_nepc" %in% names(clin_pred)) {
  p_hits <- intersect(c("p_nepc",".pred_nepc","pred_nepc","pnepc"), names(clin_pred))
  stopifnot("Deployment prediction file lacks p_nepc/.pred_nepc" = length(p_hits) >= 1)
  clin_pred <- clin_pred %>% rename(p_nepc = !!dplyr::all_of(p_hits[1]))
}

clin_tf <- read.delim(opt$clinic_tf, sep = "", quote = "\"", check.names = FALSE) %>% clean_names()
if (!"sample_id" %in% names(clin_tf)) {
  id_hits_tf <- intersect(c("sample_id","sample","id"), names(clin_tf))
  stopifnot("Deployment ctF file lacks sample_id/sample/id" = length(id_hits_tf) >= 1)
  clin_tf <- clin_tf %>% rename(sample_id = !!dplyr::all_of(id_hits_tf[1]))
}
tf_hits <- intersect(c("ct_fraction","tumor_fraction","tumour_fraction","ichor_tfx","tfx","tf"), names(clin_tf))
stopifnot("Deployment ctF file lacks a ct_fraction synonym" = length(tf_hits) >= 1)
clin_tf <- clin_tf %>% rename(ct_fraction = !!dplyr::all_of(tf_hits[1]))

# Merge and apply the same ctF-aware calibration (using development residual model and intercept shift)
new_dat <- clin_pred %>%
  inner_join(clin_tf %>% select(sample_id, ct_fraction), by = "sample_id") %>%
  mutate(
    p_nepc = pmin(pmax(p_nepc, 1e-6), 1-1e-6),
    logit_p = qlogis(p_nepc)
  ) %>%
  mutate(
    lp_hat = predict(m_resid, newdata = tibble(ct_fraction = ct_fraction)),
    logit_resid = logit_p - lp_hat,
    intercept_shift = mean(dat$logit_p) - mean(residuals(m_resid)),  # keep the same mean scale as development set
    p_NEPC_cal = plogis(logit_resid + intercept_shift),
    p_NEPC_cal_10 = p_NEPC_cal * 10
  )

# 6) Absolute survival probabilities:
# S(t|x) = exp( - H0(t) * exp(β * p_NEPC_cal_10) ), at prespecified times
times_num <- as.numeric(strsplit(opt$times, ",")[[1]])
H0_vec <- approx(bh$time, bh$hazard, xout = times_num, rule = 2)$y

new_dat <- new_dat %>%
  mutate(lp = as.numeric(beta * p_NEPC_cal_10))

for (i in seq_along(times_num)) {
  nm <- paste0("S_", times_num[i], "m")
  new_dat[[nm]] <- exp(- H0_vec[i] * exp(new_dat$lp))
}

# 7) Export deployment results
out <- new_dat %>%
  select(sample_id, ct_fraction, p_NEPC = p_nepc, p_NEPC_cal, p_NEPC_cal_10, risk_lp = lp,
         starts_with("S_"))
readr::write_tsv(out, out_path("clinic_survival_predictions.tsv"))
message(sprintf("Wrote: %s", out_path("clinic_survival_predictions.tsv")))
