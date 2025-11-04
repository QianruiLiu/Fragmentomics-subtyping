# =======================
# Environment & CLI
# =======================
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(janitor)
  library(tidymodels)
  library(dials)
  library(glmnet)
  library(recipes)
  library(tidyselect)
  library(broom)
  library(ggplot2)
  library(yardstick)
  library(stringr)
  library(scales)
  library(purrr)
  library(tidyr)
  library(readr)
  library(dplyr)
})

# ---- CLI options ----
option_list <- list(
  make_option("--train", type = "character", help = "Training feature matrix TSV (LuCaP).", metavar = "FILE"),
  make_option("--test", type = "character", default = "", help = "External feature matrix TSV (lpWGS). Optional.", metavar = "FILE"),
  make_option("--labels", type = "character", help = "Subtype label table (must contain pc_phenotype or PC.phenotype).", metavar = "FILE"),
  make_option("--tf_train", type = "character", help = "Training ctF table (columns: sample, tumor_fraction).", metavar = "FILE"),
  make_option("--tf_test", type = "character", default = "", help = "External ctF table (columns: sample, tumor_fraction). Optional.", metavar = "FILE"),
  make_option("--healthy", type = "character", help = "Healthy feature matrix TSV for baseline (TF ~ 0).", metavar = "FILE"),
  make_option("--gamma", type = "double", default = 0.5, help = "Soft TF-correction gamma (default: 0.5)."),
  make_option("--arpc_low_tf_q", type = "double", default = 0.30, help = "Quantile for low-TF ARPC to approximate Healthy when healthy_df is NULL (default: 0.30)."),
  make_option("--folds", type = "integer", default = NA, help = "CV folds (v). If NA, use v = min(5, #NEPC), lower-bounded at 2."),
  make_option("--repeats", type = "integer", default = 3, help = "CV repeats (default: 3)."),
  make_option("--grid_size", type = "integer", default = 25, help = "Latin hypercube grid size (default: 25)."),
  make_option("--seed", type = "integer", default = 200, help = "Seed for reproducibility (default: 200)."),
  make_option("--out_prefix", type = "character", default = "", help = "If set, prefix all outputs as <out_prefix>_<name>.*")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---- Basic checks ----
stopifnot(!is.null(opt$train), opt$train != "")
stopifnot(!is.null(opt$labels), opt$labels != "")
stopifnot(!is.null(opt$tf_train), opt$tf_train != "")
stopifnot(!is.null(opt$healthy), opt$healthy != "")

# optional files can be empty string
has_test <- !is.null(opt$test) && nchar(opt$test) > 0
has_tf_te  <- !is.null(opt$tf_test) && nchar(opt$tf_test) > 0

set.seed(opt$seed)
tidymodels_prefer()

# Helper: output path with optional prefix
out_path <- function(name) {
  if (is.null(opt$out_prefix) || opt$out_prefix == "") return(name)
  dir.create(dirname(file.path(opt$out_prefix)), showWarnings = FALSE, recursive = TRUE)
  paste0(opt$out_prefix, "_", name)
}

# =======================
# 0A) Custom step: LLR features (TF-aware interpolation)
# =======================

step_llr_tf <- function(recipe, features, subtype_col = "subtype", tf_col = "tumor_fraction", healthy_df = NULL,
                        arpc_low_tf_q = 0.30, sigma_floor = 1e-6, role = "predictor", skip = FALSE,
                        id = rand_id("llr_tf")) {
  add_step(recipe,
    step_llr_tf_new(
      terms = enquo(features),
      subtype_col = subtype_col,
      tf_col = tf_col,
      healthy_df = healthy_df,
      arpc_low_tf_q = arpc_low_tf_q,
      sigma_floor = sigma_floor,
      role = role,
      skip = skip,
      trained = FALSE,
      feats = NULL,
      mu_ref = NULL,
      sd_ref = NULL,
      id = id
    )
  )
}

step_llr_tf_new <- function(terms, subtype_col, tf_col, healthy_df, arpc_low_tf_q, sigma_floor,
                            role, skip, trained, feats, mu_ref, sd_ref, id) {
  step(
    subclass = "llr_tf",
    terms = terms,
    role = role,
    skip = skip,
    trained = trained,
    subtype_col = subtype_col,
    tf_col = tf_col,
    healthy_df = healthy_df,
    arpc_low_tf_q = arpc_low_tf_q,
    sigma_floor = sigma_floor,
    feats = feats,
    mu_ref = mu_ref,
    sd_ref = sd_ref,
    id = id
  )
}

prep.step_llr_tf <- function(x, training, info = NULL, ...) {
  # Select feature columns from recipe terms
  feats <- names(tidyselect::eval_select(x$terms, data = training))
  if (!x$tf_col %in% names(training)) stop("`tf_col` not found in training data.")
  if (!x$subtype_col %in% names(training)) stop("`subtype_col` not found in training data.")
  
  tr <- training
  ar <- tr %>% filter(.data[[x$subtype_col]] == "ARPC")
  ne <- tr %>% filter(.data[[x$subtype_col]] == "NEPC")
  if (nrow(ar) == 0 || nrow(ne) == 0) stop("Both ARPC and NEPC are required in training.")
  
  # Healthy reference: use provided healthy_df; otherwise low-ctF ARPC as proxy
  if (!is.null(x$healthy_df)) {
    he_raw <- x$healthy_df
    feats_he <- intersect(feats, names(he_raw))
    if (length(feats_he) == 0) stop("`healthy_df` has no overlap with training features.")
    he <- he_raw[, feats_he, drop = FALSE]
    if (length(setdiff(feats, feats_he)) > 0) {
      he <- tibble::add_column(
        he,
        !!!setNames(rep(list(NA_real_), length(setdiff(feats, feats_he))), setdiff(feats, feats_he))
      )
      he <- he[, feats, drop = FALSE]
    }
  } else {
    q <- quantile(ar[[x$tf_col]], probs = x$arpc_low_tf_q, na.rm = TRUE)
    he <- ar %>% filter(.data[[x$tf_col]] <= q) %>% select(all_of(feats))
    if (nrow(he) == 0) stop("No low-TF ARPC available to approximate Healthy.")
  }
  
  # Class-conditional means/sds per feature
  mu_ref <- tibble(
    feature = feats,
    Healthy = sapply(feats, function(f) mean(he[[f]], na.rm = TRUE)),
    ARPC = sapply(feats, function(f) mean(ar[[f]], na.rm = TRUE)),
    NEPC = sapply(feats, function(f) mean(ne[[f]], na.rm = TRUE))
  )
  sd_ref <- tibble(
    feature = feats,
    Healthy = pmax(sapply(feats, function(f) stats::sd(he[[f]], na.rm = TRUE)), x$sigma_floor),
    ARPC = pmax(sapply(feats, function(f) stats::sd(ar[[f]], na.rm = TRUE)), x$sigma_floor),
    NEPC = pmax(sapply(feats, function(f) stats::sd(ne[[f]], na.rm = TRUE)), x$sigma_floor)
  )
  
  mu_ref <- column_to_rownames(mu_ref, "feature")
  sd_ref <- column_to_rownames(sd_ref, "feature")
  
  step_llr_tf_new(
    terms = x$terms,
    subtype_col = x$subtype_col,
    tf_col = x$tf_col,
    healthy_df = NULL,
    arpc_low_tf_q = x$arpc_low_tf_q,
    sigma_floor = x$sigma_floor,
    role = x$role,
    skip = x$skip,
    trained = TRUE,
    feats = feats,
    mu_ref = mu_ref,
    sd_ref = sd_ref,
    id = x$id
  )
}

bake.step_llr_tf <- function(object, new_data, ...) {
  # Compute per-row log-likelihoods under ARPC/NEPC mixtures given ctF (alpha)
  feats <- intersect(object$feats, names(new_data))
  if (!object$tf_col %in% names(new_data)) stop("`tf_col` not found in new_data during bake().")
  
  alpha <- pmin(pmax(new_data[[object$tf_col]], 0), 1)
  
  muH <- object$mu_ref[feats, "Healthy"]; vH <- object$sd_ref[feats, "Healthy"]^2
  muA <- object$mu_ref[feats, "ARPC"];    vA <- object$sd_ref[feats, "ARPC"]^2
  muN <- object$mu_ref[feats, "NEPC"];    vN <- object$sd_ref[feats, "NEPC"]^2
  
  llA <- llN <- numeric(nrow(new_data))
  for (i in seq_len(nrow(new_data))) {
    a <- ifelse(is.na(alpha[i]), 0, alpha[i])
    x  <- as.numeric(new_data[i, feats, drop = TRUE])
    # Linear interpolation of means/variances: alpha·Tumor + (1-alpha)·Healthy
    muA_p <- a * muA + (1 - a) * muH
    muN_p <- a * muN + (1 - a) * muH
    vA_p  <- a * vA  + (1 - a) * vH
    vN_p  <- a * vN  + (1 - a) * vH
    
    logpdf <- function(x, mu, v) { -0.5 * ( ((x - mu)^2 / v) + log(2*pi*v) ) }
    
    llA[i] <- sum(logpdf(x, muA_p, vA_p), na.rm = TRUE)
    llN[i] <- sum(logpdf(x, muN_p, vN_p), na.rm = TRUE)
  }
  new_data$ll_arpc <- llA
  new_data$ll_nepc <- llN
  new_data$llr <- llA - llN
  new_data
}

print.step_llr_tf <- function(x, width = max(20, options()$width - 30), ...) {
  cat("LLR on ", length(x$feats), " features\n", sep = "")
  invisible(x)
}

# =======================
# 0B) Custom step: Soft TF correction
# =======================

step_tf_softcorr <- function(recipe,
                             features,
                             tf_col = "tumor_fraction",
                             healthy_df,
                             gamma = 0.5,
                             role = "predictor",
                             skip = FALSE,
                             id = rand_id("tf_softcorr")) {
  add_step(
    recipe,
    step_tf_softcorr_new(
      terms = enquo(features),
      tf_col = tf_col,
      healthy_df = healthy_df,
      gamma = gamma,
      role = role,
      skip = skip,
      trained = FALSE,
      feats = NULL,
      muH = NULL,
      xbar = NULL,
      id = id
    )
  )
}

step_tf_softcorr_new <- function(terms, tf_col, healthy_df, gamma,
                                 role, skip, trained, feats, muH, xbar, id) {
  step(
    subclass = "tf_softcorr",
    terms = terms,
    role = role,
    skip = skip,
    trained = trained,
    tf_col = tf_col,
    healthy_df = healthy_df,
    gamma = gamma,
    feats = feats,
    muH = muH,
    xbar = xbar,
    id = id
  )
}

prep.step_tf_softcorr <- function(x, training, info = NULL, ...) {
  # Select numeric features to be corrected
  feats <- names(tidyselect::eval_select(x$terms, data = training))
  if (!x$tf_col %in% names(training)) stop("`tf_col` not found in training data.")
  
  feats <- feats[vapply(training[feats], is.numeric, TRUE)]
  if (length(feats) == 0) stop("No numeric features found for TF correction.")
  
  if (is.null(x$healthy_df)) stop("`healthy_df` must be provided.")
  feats_he <- intersect(feats, names(x$healthy_df))
  if (length(feats_he) == 0) stop("`healthy_df` has no overlap with training features.")
  muH <- setNames(rep(NA_real_, length(feats)), feats)
  muH[feats_he] <- vapply(feats_he, function(f) mean(x$healthy_df[[f]], na.rm = TRUE), numeric(1))
  
  xbar <- vapply(feats, function(f) mean(training[[f]], na.rm = TRUE), numeric(1))
  
  step_tf_softcorr_new(
    terms = x$terms,
    tf_col = x$tf_col,
    healthy_df = NULL,
    gamma = x$gamma,
    role = x$role,
    skip = x$skip,
    trained = TRUE,
    feats = feats,
    muH = muH,
    xbar = xbar,
    id = x$id
  )
}

bake.step_tf_softcorr <- function(object, new_data, ...) {
  # Apply soft correction row-wise using the provided ctF column
  if (!object$tf_col %in% names(new_data)) stop("`tf_col` not found in new_data during bake().")
  
  feats <- intersect(object$feats, names(new_data))
  tfv <- new_data[[object$tf_col]]
  tfv <- pmin(pmax(tfv, 0), 1)
  gamma <- object$gamma
  
  for (f in feats) {
    if (is.numeric(new_data[[f]])) {
      new_data[[f]] <- new_data[[f]] - gamma * (1 - tfv) * (object$muH[[f]] - object$xbar[[f]])
    }
  }
  new_data
}

print.step_tf_softcorr <- function(x, width = max(20, options()$width - 30), ...) {
  cat("Soft TF correction on ", length(x$feats), " features; gamma=", x$gamma, "\n", sep = "")
  invisible(x)
}

# =======================
# 1) Load data
# =======================
train <- read_tsv(opt$train, show_col_types = FALSE) |>
  clean_names() |>
  rename(sample = sample_id)

test <- if (has_test) {
  read_tsv(opt$test, show_col_types = FALSE) |>
    clean_names() |>
    rename(sample = sample_id)
} else {
  tibble()
}

# =======================
# 2) Load TF
# =======================
tf_tr <- read_tsv(opt$tf_train, show_col_types = FALSE) |>
  clean_names() |>
  select(sample, tumor_fraction)

tf_te <- if (has_tf_te) {
  read_tsv(opt$tf_test, show_col_types = FALSE) |>
    clean_names() |>
    select(sample, tumor_fraction)
} else {
  tibble()
}

# =======================
# 3) Load labels
# =======================
lab_raw <- read_tsv(opt$labels, show_col_types = FALSE) |> clean_names()
if ("pc_phenotype" %in% names(lab_raw)) {
  label_col <- "pc_phenotype"
} else if ("PC.phenotype" %in% names(lab_raw)) {
  lab_raw <- lab_raw |> rename(pc_phenotype = `PC.phenotype`)
  label_col <- "pc_phenotype"
} else {
  stop("Label file must contain pc_phenotype / PC.phenotype.")
}

# =======================
# 4) Sample name normalization
# =======================
norm_sample_from_label <- function(x){
  x <- as.character(x)
  x <- gsub("-", "_", x)
  need_prefix <- grepl("^(\\d+(_\\d+)?|\\d+CR(_rep\\d+)?)$", x)
  x <- ifelse(need_prefix & !grepl("^LuCaP_", x, ignore.case = TRUE),
              paste0("LuCaP_", x), x)
  x
}
to_base <- function(x){
  x <- as.character(x)
  x <- sub("_rep[0-9]+$", "", x, perl = TRUE)
  x <- sub("_recal$", "", x, perl = TRUE)
  x
}

# =======================
# 5) Build 2-class labels
# =======================
labels <- lab_raw |>
  mutate(
    sample = as.character(sample),
    sample_base = norm_sample_from_label(sample),
    subtype_raw = .data[[label_col]]
  ) |>
  filter(!is.na(subtype_raw), subtype_raw != "") |>
  filter(subtype_raw %in% c("ARPC", "NEPC")) |>
  transmute(
    sample_base,
    subtype = factor(subtype_raw, levels = c("ARPC","NEPC"))
  )
if (nrow(labels) == 0) stop("No ARPC/NEPC labels found.")

# =======================
# 6) Join TF + labels
# =======================
tf_tr <- tf_tr |>
  mutate(sample_base = to_base(sample)) |>
  select(sample_base, tumor_fraction)

train <- train |>
  mutate(sample_base = to_base(sample)) |>
  left_join(tf_tr, by = "sample_base") |>
  left_join(labels, by = "sample_base") |>
  filter(!is.na(subtype))

test <- if (nrow(test) > 0) {
  if (nrow(tf_te) == 0) stop("External set provided but --tf_test is missing.")
  test |>
    mutate(sample_base = to_base(sample)) |>
    left_join(tf_te, by = "sample")
} else test

cat("Final training n:", nrow(train), "\n")
cat("Class distribution:\n"); print(table(train$subtype))

# =======================
# 7) Prepare feature sets
# =======================
nm <- names(train)

# Define feature families for LLR (avoid PCs)
occ_cols <- nm[grepl("coverage|atac", nm, ignore.case = TRUE)]
clea_cols <- nm[grepl("cleavage", nm, ignore.case = TRUE)]
len_cols <- nm[grepl("fraglen", nm, ignore.case = TRUE)]
motif_cols <- nm[grepl("motif", nm, ignore.case = TRUE)]
wps_cols <- nm[grepl("wps", nm, ignore.case = TRUE)]

dedup <- function(x, taken) setdiff(x, taken)
taken <- character(0)
occ_cols <- dedup(occ_cols, taken); taken <- union(taken, occ_cols)
clea_cols <- dedup(clea_cols, taken); taken <- union(taken, clea_cols)
len_cols <- dedup(len_cols, taken); taken <- union(taken, len_cols)
motif_cols <- dedup(motif_cols, taken); taken <- union(taken, motif_cols)
wps_cols <- dedup(wps_cols, taken); taken <- union(taken, wps_cols)

drop_non_predictor <- function(cols, dat){
  cols <- intersect(cols, names(dat))
  cols <- cols[vapply(dat[cols], is.numeric, TRUE)]
  setdiff(cols, c("sample","sample_base","subtype","tumor_fraction"))
}

occ_all <- drop_non_predictor(occ_cols,   train)
clea_all <- drop_non_predictor(clea_cols,  train)
len_all <- drop_non_predictor(len_cols,   train)
motif_all <- drop_non_predictor(motif_cols, train)
wps_all <- drop_non_predictor(wps_cols,   train)

# LLR features
llr_feats <- unique(c(occ_all, clea_all, len_all, wps_all))
if (length(llr_feats) == 0) {
  message("[LLR] Family regex matched 0 columns; falling back to ALL numeric predictors.")
  num_cols <- names(train)[vapply(train, is.numeric, TRUE)]
  llr_feats <- setdiff(num_cols, c("tumor_fraction","sample","sample_base"))
}

# Soft correction features (all numeric predictors)
all_numeric_feats <- names(train)[vapply(train, is.numeric, TRUE)]
all_numeric_feats <- setdiff(all_numeric_feats, c("tumor_fraction","sample","sample_base"))

cat("LLR features:", length(llr_feats), "\n")
cat("Soft correction features:", length(all_numeric_feats), "\n")

# Heavy-tailed candidates for asinh transform
heavy_candidates <- llr_feats[grepl("total_fragments|n_intervals|coverage|counts",
                                    llr_feats, ignore.case = TRUE)]

# =======================
# 7.5) Load external Healthy baseline
# =======================
healthy_df <- read_tsv(opt$healthy, show_col_types = FALSE) |>
  clean_names()

# =======================
# 8) Recipe: impute -> LLR -> Soft TF correction -> asinh -> ZV -> normalize
# =======================
rec <- recipe(subtype ~ ., data = train) %>%
  update_role(all_of(c("sample","sample_base","tumor_fraction")), new_role = "id") %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_llr_tf(features = all_of(llr_feats),
              subtype_col = "subtype",
              tf_col = "tumor_fraction",
              healthy_df = healthy_df,
              arpc_low_tf_q = opt$arpc_low_tf_q,
              sigma_floor = 1e-6) %>%
  step_tf_softcorr(features = all_of(all_numeric_feats),
                   tf_col = "tumor_fraction",
                   healthy_df = healthy_df,
                   gamma = opt$gamma) %>%
  step_mutate(across(all_of(heavy_candidates), ~ asinh(.x))) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_numeric_predictors())

# =======================
# 9) Cross-validation
# =======================
train <- train %>% mutate(subtype = fct_relevel(subtype, "ARPC","NEPC"))
options(yardstick.event_first = FALSE)

n_nepc <- sum(train$subtype == "NEPC")
v_auto <- max(2, min(5, n_nepc))
v <- if (is.na(opt$folds)) v_auto else opt$folds
set.seed(2025)
folds <- vfold_cv(train, v = v, repeats = opt$repeats, strata = subtype)
metrics_used <- metric_set(roc_auc, pr_auc, accuracy)

# =======================
# 10) Model & tuning: elastic net + 1-SE rule
# =======================
logit_en <- logistic_reg(penalty = tune(), mixture = tune()) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

wf <- workflow() %>% add_recipe(rec) %>% add_model(logit_en)

param <- extract_parameter_set_dials(wf) %>%
  update(
    penalty = penalty(range = c(-4, 1)),
    mixture = mixture(range = c(0, 1))
  )

ctrl <- control_grid(
  allow_par = FALSE,          
  save_pred = TRUE,
  save_workflow = TRUE,
  verbose = TRUE
)
set.seed(2025)
grid <- grid_latin_hypercube(param, size = opt$grid_size)

cat("Tuning...\n")
res <- tune_grid(wf, folds, grid = grid, metrics = metrics_used, control = ctrl)

# Choose by ROC-AUC 1-SE rule (simpler model preference)
best <- select_by_one_std_err(res, metric = "roc_auc", desc(penalty))
cat("Chosen by 1-SE:\n"); print(best)

final_wf <- finalize_workflow(wf, best)

# =======================
# 11) Fit on full LuCaP
# =======================
final_fit <- fit(final_wf, train)

# CV summary
perf <- collect_metrics(res)
cat("\n=== CV Performance ===\n")
print(perf)

# =======================
# 12) External predictions (lpWGS)
# =======================
if (nrow(test) > 0) {
  test_pred <- predict(final_fit, test, type = "prob")
  test_class <- predict(final_fit, test, type = "class")
  
  results <- bind_cols(
    test %>% select(sample),
    test_pred,
    test_class
  )
  
  cat("\n=== External predictions (first 20) ===\n")
  print(head(results, 20))
  
  # Keep original filename; also write with out_prefix if provided
  write_tsv(results, out_path("combined_LLR_softTF_predictions.tsv"))
  cat("\nSaved to ", out_path("combined_LLR_softTF_predictions.tsv"), "\n", sep = "")
} else {
  cat("No external test data.\n")
}

# =======================
# 13) Sanity checks
# =======================
cat("\n=== Training TF missingness ===\n")
print(table(is.na(train$tumor_fraction)))

cat("\n=== Top 20 coefficients (by absolute value) ===\n")
coef_tbl <- tidy(extract_fit_parsnip(final_fit)) %>%
  filter(term != "(Intercept)") %>%
  arrange(desc(abs(estimate)))
print(head(coef_tbl, 20))

# Top 10 samples with highest NEPC probability (if external provided)
if (nrow(test) > 0) {
  cat("\n=== Top 10 samples with highest NEPC probability ===\n")
  top10_nepc <- bind_cols(
    test %>% select(sample),
    predict(final_fit, test, type = "prob"),
    predict(final_fit, test, type = "class")
  ) %>%
    arrange(desc(.pred_NEPC)) %>%
    slice_head(n = 10)
  
  print(top10_nepc)
  write_tsv(top10_nepc, out_path("top10_nepc_prob_samples.tsv"))
  cat("\nSaved top 10 to ", out_path("top10_nepc_prob_samples.tsv"), "\n", sep = "")
}

# ----------------------------------------------------------
# CV predictions restricted to chosen config
# ----------------------------------------------------------
pred_cv <- collect_predictions(res) %>%
  filter(.config == best$.config)
options(yardstick.event_first = FALSE)

auc_by_resample_corr <- pred_cv %>%
  group_by(id, id2, .config) %>%
  roc_auc(truth = subtype, .pred_NEPC, event_level = "second") %>%
  ungroup() %>%
  mutate(model = "ctF-aware (LLR+softTF)") %>%
  select(model, id, id2, auc = .estimate)

pr_by_resample_corr <- pred_cv %>%
  group_by(id, id2, .config) %>%
  pr_auc(truth = subtype, .pred_NEPC, event_level = "second") %>%
  ungroup() %>%
  mutate(model = "ctF-aware (LLR+softTF)") %>%
  select(model, id, id2, pr_auc = .estimate)

write_tsv(auc_by_resample_corr, out_path("cv_auc_per_resample_corrected.tsv"))
write_tsv(pr_by_resample_corr,  out_path("cv_pr_auc_per_resample_corrected.tsv"))

# ----------------------------------------------------------
# Figure A: family-level importance (sum of |coef|)
# ----------------------------------------------------------
coef_final <- tidy(extract_fit_parsnip(final_fit)) %>%
  filter(term != "(Intercept)") %>%
  mutate(abs_est = abs(estimate))

family_of <- function(x){
  x <- tolower(x)
  case_when(
    str_detect(x, "wps") ~ "WPS",
    str_detect(x, "cleavage") ~ "Cleavage",
    str_detect(x, "fraglen|fragment") ~ "Length",
    str_detect(x, "motif|entropy") ~ "Motif",
    str_detect(x, "coverage|atac") ~ "Coverage",
    str_detect(x, "^llr$|^ll_") ~ "LLR-derived",
    TRUE ~ "Other"
  )
}

fam_imp <- coef_final %>%
  mutate(family = family_of(term)) %>%
  group_by(family) %>%
  summarise(importance = sum(abs_est), .groups = "drop") %>%
  arrange(desc(importance))

fig5a_fam <- ggplot(fam_imp, aes(x = reorder(family, importance), y = importance)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = "Sum of |coefficients|",
       title = "Family-level importance (glmnet coefficients)") +
  theme_bw()

# ----------------------------------------------------------
# Figure B: glmnet coefficient paths (lambda path)
# ----------------------------------------------------------
fit_eng <- final_fit %>% extract_fit_parsnip() %>% extract_fit_engine()
path <- broom::tidy(fit_eng)  # columns: term, step, estimate, lambda

top_terms <- coef_final %>%
  arrange(desc(abs_est)) %>%
  slice_head(n = 15) %>%
  pull(term)

path_sel <- path %>%
  filter(term %in% top_terms, term != "(Intercept)")

fig5b_path <- ggplot(path_sel, aes(x = lambda, y = estimate, group = term, color = term)) +
  geom_line() +
  scale_x_log10() +
  labs(x = "Lambda (penalty, log scale)",
       y = "Coefficient",
       title = "Coefficient paths of top features") +
  theme_bw() +
  theme(legend.position = "right")

print(fig5a_fam)
print(fig5b_path)

ggsave(out_path("FigA_FAM.png"), fig5a_fam, width = 5, height = 4, dpi = 300)
ggsave(out_path("FigB_PATH.png"), fig5b_path, width = 6, height = 4, dpi = 300)

cat("\nAll done.\n")
