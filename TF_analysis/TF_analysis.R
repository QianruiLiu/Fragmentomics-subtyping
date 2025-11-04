# =======================
# Environment & CLI
# =======================
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(janitor)
  library(yardstick)
  library(rsample)
  library(broom)
  library(rlang)
  library(ggrepel)
  library(patchwork)
})

# ---- CLI options ----
option_list <- list(
  make_option("--subtypes", type = "character", help = "Subtype label file (must contain sample & pc_phenotype).", metavar = "FILE"),
  make_option("--tf_file", type = "character", help = "Tumor fraction file (columns: sample, tumor_fraction).", metavar = "FILE"),
  make_option("--sites_dir", type = "character", default = ".", help = "Root directory to search for *_sites folders (default: current dir)."),
  make_option("--topN", type = "integer", default = 20, help = "Number of top TFs to plot in violin plots (default: 20)."),
  make_option("--topK", type = "integer", default = 50, help = "Number of top TFs for bar/dot discriminability plots (default: 50)."),
  make_option("--seed", type = "integer", default = 42, help = "Random seed for reproducibility (default: 42)."),
  make_option("--out_prefix", type = "character", default = "", help = "Optional output prefix for all files.")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Basic checks
stopifnot(!is.null(opt$subtypes), opt$subtypes != "")
stopifnot(!is.null(opt$tf_file), opt$tf_file != "")

set.seed(opt$seed)
options(warn = 1)

log_info <- function(...){ message(format(Sys.time(), "[%H:%M:%S] "), ...) }

# Helper: output path with optional prefix
out_path <- function(name) {
  if (is.null(opt$out_prefix) || opt$out_prefix == "") return(name)
  dir.create(dirname(file.path(opt$out_prefix)), showWarnings = FALSE, recursive = TRUE)
  paste0(opt$out_prefix, "_", name)
}

# =======================
# Unified sample naming: keep Griffin style (e.g., 70CR_rep2_LuCaP)
# =======================
clean_griffin_id <- function(x){
  y <- x
  y <- gsub("\\s+", "", y)
  y <- gsub("-", "_", y)
  y
}

# =======================
# 1) Load subtypes and TFX
# =======================
log_info("[S1] Loading subtypes and TFX ...")
subtypes <- readr::read_tsv(opt$subtypes, show_col_types = FALSE) %>%
  clean_names()

nm <- names(subtypes)
col_sample <- nm[grepl("^sample$", nm)]
col_type <- nm[grepl("pc_phenotype", nm)]
if (!(length(col_sample)==1 && length(col_type)==1)) {
  abort(paste0(
    "[S1] Cannot recognize columns in subtypes file.\n",
    "  Required: sample, pc_phenotype\n",
    "  Found: ", paste(nm, collapse = ", ")
  ))
}

subtypes <- subtypes %>%
  transmute(
    sample_std = clean_griffin_id(!!sym(col_sample)),
    subtype = !!sym(col_type)
  ) %>%
  filter(subtype %in% c("ARPC","NEPC")) %>%
  mutate(subtype = factor(subtype, levels = c("ARPC","NEPC")))

log_info("[S1] Subtypes rows: ", nrow(subtypes))
log_info("[S1] Subtypes levels: ", paste(levels(subtypes$subtype), collapse = ", "))

# TFX (optional)
tfx <- readr::read_tsv(opt$tf_file, show_col_types = FALSE) %>%
  clean_names()
if (!"sample" %in% names(tfx)) {
  nm_tfx <- names(tfx); col_s <- nm_tfx[grepl("^sample(_std)?$", nm_tfx)][1]
  if (!is.na(col_s)) names(tfx)[names(tfx)==col_s] <- "sample" else {
    log_info("[S1] TF file has no 'sample' column, continuing with empty TFX. Can be ignored.")
    tfx$sample <- character(0)
  }
}
tfx <- tfx %>%
  mutate(sample_std = clean_griffin_id(sample)) %>%
  select(-sample)

log_info("[S1] TFX rows: ", nrow(tfx))
if (nrow(tfx) > 0) log_info("[S1] TFX columns: ", paste(names(tfx), collapse = ", "))

# =======================
# 2) Scan and read coverage.tsv files
# =======================
log_info("[S2] Scanning *_sites directories in: ", opt$sites_dir)
site_dirs <- list.dirs(opt$sites_dir, recursive = FALSE, full.names = TRUE) %>%
  keep(~ grepl("_sites$", .x, ignore.case = TRUE))
log_info("[S2] Found sites directories: ", length(site_dirs))
if (length(site_dirs) > 0) log_info("[S2] Example sites dir: ", site_dirs[[1]])

cov_files <- purrr::map(site_dirs, ~ list.files(.x,
                                                pattern = "GC_corrected\\.coverage\\.tsv$",
                                                recursive = TRUE, full.names = TRUE)) %>% unlist()
log_info("[S2] Matched coverage files: ", length(cov_files))
if (length(cov_files) == 0) {
  log_info("[S2][NOTICE] No coverage files matched. Check structure/suffix/scan depth.")
} else {
  log_info("[S2] Example coverage file: ", cov_files[[1]])
}

read_one <- function(fp){
  tf_folder <- basename(dirname(dirname(fp)))  # e.g., A_sites
  sample_dir <- basename(dirname(fp))           # e.g., 70CR_rep2_LuCaP
  
  df_raw <- readr::read_tsv(fp, show_col_types = FALSE)
  orig_names <- names(df_raw)
  
  # Pure numeric columns (with optional minus sign)
  pos_cols_orig <- orig_names[grepl("^[-]?\\d+$", orig_names)]
  if (length(pos_cols_orig) == 0) {
    log_info("[read_one] Skip (no numeric position columns): ", fp); return(NULL)
  }
  
  df <- df_raw %>% clean_names()
  pos_cols_clean <- janitor::make_clean_names(pos_cols_orig)
  
  has_site_name <- "site_name" %in% names(df)
  base_cols <- intersect(c("mean_coverage","central_coverage","amplitude"), names(df))
  
  # Derive shoulder mean from curve: |pos| ∈ [90,180]
  pos_num <- as.numeric(gsub("^x_?", "", gsub("_", "", pos_cols_clean)))
  shoulder_idx <- which(abs(pos_num) >= 90 & abs(pos_num) <= 180)
  if (length(shoulder_idx) == 0) {
    log_info("[read_one] Skip (no shoulder indices 90-180bp): ", fp); return(NULL)
  }
  
  shoulder_mean <- df %>%
    select(all_of(pos_cols_clean[shoulder_idx])) %>%
    mutate(.row_id = dplyr::row_number()) %>%
    tidyr::pivot_longer(-.row_id, values_to = "val") %>%
    dplyr::group_by(.row_id) %>%
    dplyr::summarise(shoulder_mean = mean(val, na.rm = TRUE), .groups = "drop") %>%
    dplyr::pull(shoulder_mean)
  
  out <- df %>%
    mutate(
      tf_folder = tf_folder,
      sample_std = clean_griffin_id(sample_dir),
      shoulder_mean = shoulder_mean
    )
  
  # If central_coverage not in table, approximate with |pos| <= 30bp
  if (!"central_coverage" %in% names(out)) {
    center_idx <- which(abs(pos_num) <= 30)
    if (length(center_idx) == 0) {
      log_info("[read_one] Skip (no center indices <=30bp): ", fp); return(NULL)
    }
    central_cov <- df %>%
      select(all_of(pos_cols_clean[center_idx])) %>%
      mutate(.row_id = dplyr::row_number()) %>%
      tidyr::pivot_longer(-.row_id, values_to = "val") %>%
      dplyr::group_by(.row_id) %>%
      dplyr::summarise(central_cov = mean(val, na.rm = TRUE), .groups = "drop") %>%
      dplyr::pull(central_cov)
    out$central_coverage <- central_cov
  }
  
  # Metrics: center_over_shoulder & contrast
  out <- out %>%
    mutate(
      center_over_shoulder = central_coverage / pmax(shoulder_mean, 1e-9),
      contrast = (shoulder_mean - central_coverage) / pmax(shoulder_mean, 1e-9)
    )
  
  keep_cols <- c(
    if (has_site_name) "site_name",
    "sample","sample_std","tf_folder",
    base_cols,"shoulder_mean","center_over_shoulder","contrast"
  )
  
  res <- intersect(keep_cols, names(out))
  if (length(res) == 0) { log_info("[read_one] Skip (empty keep cols): ", fp); return(NULL) }
  out[, res]
}

log_info("[S2] Batch reading coverage files ...")
dat_raw <- purrr::map_dfr(
  cov_files,
  ~ tryCatch(
    read_one(.x),
    error = function(e) {
      log_info("[read_one][ERROR] File failed: ", .x, " | ", conditionMessage(e))
      NULL
    }
  )
) %>%
  clean_names()

log_info("[S2] dat_raw rows: ", nrow(dat_raw), " | cols: ", length(names(dat_raw)))
log_info("[S2] dat_raw columns: ", paste(names(dat_raw), collapse = ", "))

# Standardize TF name column
if (!"site_name" %in% names(dat_raw)) {
  abort(paste0("[S2] Missing column 'site_name', cannot standardize TF name.\n",
               "  Current dat_raw columns: ", paste(names(dat_raw), collapse = ", ")))
}
dat_raw <- dat_raw %>% rename(tf_name = site_name)
log_info("[S2] Renamed site_name to tf_name.")

# =======================
# 2.5) Pre-merge mismatch list
# =======================
raw_keys <- dat_raw %>% distinct(sample_std) %>% mutate(src = "dat_raw")
sub_keys <- subtypes %>% distinct(sample_std) %>% mutate(src = "subtypes")

miss_in_subtypes <- anti_join(raw_keys, sub_keys, by = "sample_std")
miss_in_dat_raw <- anti_join(sub_keys, raw_keys, by = "sample_std")

log_info("[CHECK-pre] dat_raw unique samples: ", nrow(raw_keys))
log_info("[CHECK-pre] subtypes unique samples: ", nrow(sub_keys))
log_info("[CHECK-pre] raw not in subtypes: ", nrow(miss_in_subtypes))
log_info("[CHECK-pre] subtypes not in raw: ", nrow(miss_in_dat_raw))

if (nrow(miss_in_subtypes) > 0) {
  log_info("[CHECK-pre] Example (raw→subtypes): ",
           paste(head(miss_in_subtypes$sample_std, 10), collapse = ", "))
  readr::write_tsv(miss_in_subtypes, out_path("unmatched_in_subtypes.tsv"))
}
if (nrow(miss_in_dat_raw) > 0) {
  log_info("[CHECK-pre] Example (subtypes→raw): ",
           paste(head(miss_in_dat_raw$sample_std, 10), collapse = ", "))
  readr::write_tsv(miss_in_dat_raw, out_path("unmatched_in_dat_raw.tsv"))
}

# =======================
# 3) Merge subtypes / TFX
# =======================
log_info("[S3] Merging subtypes/TFX ...")
dat <- dat_raw %>%
  inner_join(subtypes, by = "sample_std") %>%
  left_join(tfx,      by = "sample_std") %>%
  relocate(tf_name, tf_folder, sample_std, subtype)

log_info("[S3] dat after merge rows: ", nrow(dat), " | cols: ", length(names(dat)))
if (nrow(dat) == 0) {
  abort("[S3] Data empty after merge. Check unmatched_in_* lists and sample_std normalization.")
}

# =======================
# Helper: AUC from arbitrary real-valued scores (Mann-Whitney U)
# =======================
calc_auc_from_scores <- function(score, truth, event = "NEPC"){
  ok <- is.finite(score) & !is.na(truth)
  s <- score[ok]
  y <- as.integer(truth[ok] == event)
  
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0 || length(unique(s)) < 2) return(NA_real_)
  
  r <- rank(s, ties.method = "average")
  R1 <- sum(r[y == 1])
  U1 <- R1 - n1*(n1 + 1)/2
  auc <- U1 / (n1 * n0)
  auc
}

# =======================
# 4) Bidirectional AUC + directional score
# =======================
log_info("[S4] Scoring TFs with contrast (bidirectional AUC) ...")

score_tf <- function(df) {
  df2 <- df %>% filter(is.finite(contrast))
  x_ar <- df2$contrast[df2$subtype == "ARPC"]
  x_ne <- df2$contrast[df2$subtype == "NEPC"]
  n1 <- length(x_ar); n2 <- length(x_ne)
  
  AUC_AR <- tryCatch(calc_auc_from_scores(df2$contrast, df2$subtype, event = "ARPC"),
                      error = function(e) NA_real_)
  AUC_NE <- tryCatch(calc_auc_from_scores(df2$contrast, df2$subtype, event = "NEPC"),
                      error = function(e) NA_real_)
  sAUC_AR <- if (!is.na(AUC_AR)) if (AUC_AR >= 0.5) AUC_AR else (AUC_AR - 1) else NA_real_
  sAUC_NE <- if (!is.na(AUC_NE)) if (AUC_NE >= 0.5) AUC_NE else (AUC_NE - 1) else NA_real_
  
  mu_ar <- mean(x_ar, na.rm = TRUE); mu_ne <- mean(x_ne, na.rm = TRUE)
  direction <- if (is.finite(mu_ar) && is.finite(mu_ne)) {
    if (mu_ar > mu_ne) "ARPC_high" else "NEPC_high"
  } else NA_character_
  
  sd1 <- stats::sd(x_ar, na.rm = TRUE); sd2 <- stats::sd(x_ne, na.rm = TRUE)
  denom <- sqrt((sd1^2 + sd2^2) / 2)
  d <- if (n1 > 0 && n2 > 0 && is.finite(denom) && denom > 0) {
    (mu_ar - mu_ne) / denom
  } else NA_real_
  
  p <- if (n1 > 0 && n2 > 0) {
    tryCatch(stats::wilcox.test(x_ar, x_ne)$p.value, error = function(e) NA_real_)
  } else NA_real_
  
  tibble::tibble(
    AUC_AR  = AUC_AR,  sAUC_AR  = sAUC_AR,
    AUC_NE  = AUC_NE,  sAUC_NE  = sAUC_NE,
    mean_AR = mu_ar,   mean_NE  = mu_ne,
    direction = direction, d = d, p = p,
    n_ARPC = n1, n_NEPC = n2
  )
}

ranked <- dat %>%
  dplyr::group_by(tf_name) %>%
  dplyr::group_modify(~ score_tf(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    # Unified directional score for single-column sorting and readability
    sAUC_dir = dplyr::case_when(
      direction == "ARPC_high" ~ sAUC_AR,
      direction == "NEPC_high" ~ sAUC_NE,
      TRUE ~ NA_real_
    ),
    p_adj = p.adjust(p, method = "BH")
  ) %>%
  dplyr::arrange(dplyr::desc(abs(sAUC_dir)), p_adj)

readr::write_tsv(ranked, out_path("TF_ranking_contrast.tsv"))
log_info("[S4] Written: TF_ranking_contrast.tsv (bidirectional AUC)")

# =======================
# 5) Take TopN & plot (by |sAUC_dir|)
# =======================
log_info("[S5] Taking TopN and plotting ...")
topN_names <- ranked %>% 
  filter(is.finite(sAUC_dir)) %>% 
  slice_max(order_by = abs(sAUC_dir), n = opt$topN, with_ties = FALSE) %>% 
  pull(tf_name)

plot_df <- dat %>%
  filter(tf_name %in% topN_names) %>%
  transmute(tf_name, subtype, value = contrast)

p <- ggplot(plot_df, aes(x = subtype, y = value)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = .1, height = 0, alpha = .6) +
  facet_wrap(~ tf_name, scales = "free_y") +
  labs(x = NULL,
       y = "contrast = (shoulder - center) / shoulder",
       title = sprintf("Top %d TFs by |sAUC_dir| (direction-aware)", opt$topN)) +
  theme_bw()

ggsave(out_path("TF_topN_contrast.png"), p, width = 12, height = 8, dpi = 300)
ranked %>% 
  filter(tf_name %in% topN_names) %>%
  arrange(desc(abs(sAUC_dir)), p_adj) %>%
  readr::write_tsv(out_path("TF_topN_contrast.tsv"))
log_info("[S5] Written: TF_topN_contrast.tsv, TF_topN_contrast.png")

# =======================
# 6) Export AR / NE axis (sorted by respective direction)
# =======================
log_info("[S6] Exporting AR axis / NE axis TFs ...")
ar_axis <- c("AR","FOXA1","HOXB13","GATA2")
ne_axis <- c("ASCL1","NEUROD1","BRN2","INSM1","SOX2")

ranked_AR <- ranked %>%
  filter(tf_name %in% ar_axis) %>%
  arrange(desc(sAUC_AR), p_adj)
ranked_NE <- ranked %>%
  filter(tf_name %in% ne_axis) %>%
  arrange(desc(sAUC_NE), p_adj)

readr::write_tsv(ranked_AR, out_path("TF_ranking_contrast_AR_axis.tsv"))
readr::write_tsv(ranked_NE, out_path("TF_ranking_contrast_NE_axis.tsv"))
log_info("[S6] Written: TF_ranking_contrast_AR_axis.tsv, TF_ranking_contrast_NE_axis.tsv")

# =======================
# 7) Top K discriminability plots (bar & dot)
# =======================
log_info("[S7] Generating Top-K discriminability plots ...")

# Calculate direction-agnostic magnitude
ranked_plot <- ranked %>%
  mutate(
    auc_mag = abs(sAUC_dir),
    direction_label = dplyr::case_when(
      direction == "ARPC_high" ~ "ARPC-high",
      direction == "NEPC_high" ~ "NEPC-high",
      TRUE ~ NA_character_
    )
  )

ranked_K <- ranked_plot %>%
  filter(is.finite(auc_mag)) %>%
  arrange(desc(auc_mag), desc(abs(d)), p_adj) %>%
  slice_head(n = opt$topK)

# For within-direction ordering (strong to weak)
ranked_K <- ranked_K %>%
  group_by(direction_label) %>%
  mutate(tf_order = factor(tf_name, levels = rev(tf_name))) %>%
  ungroup()

# Bar plot (faceted by direction)
p_bar <- ggplot(ranked_K,
                aes(x = auc_mag, y = tf_order, fill = direction_label)) +
  geom_col(width = 0.8) +
  facet_wrap(~ direction_label, scales = "free_y", nrow = 1) +
  geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.3) +
  labs(
    title = sprintf("Top %d TFs by discriminability |sAUC_dir|", nrow(ranked_K)),
    x = "Discriminability (|sAUC_dir|, 0–1)",
    y = NULL,
    fill = "Higher in"
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA)
  )

ggsave(out_path("TF_topK_bar_by_direction.png"), p_bar, width = 12, height = 10, dpi = 300)

# Dot plot (more compact)
p_dot <- ggplot(ranked_K,
                aes(x = auc_mag, y = tf_order, color = direction_label)) +
  geom_point(size = 2) +
  facet_wrap(~ direction_label, scales = "free_y", nrow = 1) +
  geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.3) +
  labs(
    title = sprintf("Top %d TFs by discriminability |sAUC_dir| (dot plot)", nrow(ranked_K)),
    x = "Discriminability (|sAUC_dir|, 0–1)",
    y = NULL,
    color = "Higher in"
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA)
  )

ggsave(out_path("TF_topK_dot_by_direction.png"), p_dot, width = 12, height = 10, dpi = 300)

# Export table for supplementary materials
ranked_K %>%
  select(tf_name, direction = direction_label, auc_mag, sAUC_dir, d, p, p_adj,
         AUC_AR, sAUC_AR, AUC_NE, sAUC_NE, mean_AR, mean_NE, n_ARPC, n_NEPC) %>%
  arrange(desc(auc_mag), p_adj) %>%
  readr::write_tsv(out_path("TF_topK_discriminability.tsv"))

log_info("[S7] Written: TF_topK_bar_by_direction.png, TF_topK_dot_by_direction.png")
log_info("[S7] Written: TF_topK_discriminability.tsv")

log_info("[ALL DONE] Bidirectional AUC + directional score analysis complete.")

# =======================
# 8) Heatmap plotting
# =======================

log_info("[S8] Building heatmap (fixed TF order; faceted by subtype) ...")

# ---- (1) choose exactly the TFs and their order you want to display ----
tf_order <- c("NKX3-1","AR","ERF","ASCL1","REST","HNF4G",
              "HOXB13","KLF6","NEUROG2","TCF4","MYOD1","MYOG")

# keep only those TFs present in your data (and warn if some are missing)
missing_tfs <- setdiff(tf_order, unique(dat$tf_name))
if (length(missing_tfs) > 0) {
  log_info("[S8][WARN] Missing in data and will be dropped: ", paste(missing_tfs, collapse = ", "))
}
tf_keep <- intersect(tf_order, unique(dat$tf_name))

# ---- (2) wide matrix TF x sample from contrast, then row-wise z-score ----
mat_df <- dat %>%
  filter(tf_name %in% tf_keep) %>%
  select(tf_name, sample_std, subtype, contrast) %>%
  group_by(tf_name, sample_std, subtype) %>%
  summarise(contrast = mean(contrast, na.rm = TRUE), .groups = "drop")

# wide for z-score calc across ALL samples (both subtypes together)
mat_wide <- mat_df %>%
  select(tf_name, sample_std, contrast) %>%
  tidyr::pivot_wider(names_from = sample_std, values_from = contrast)

if (nrow(mat_wide) == 0L) {
  log_info("[S8] Skip heatmap: no data.")
} else {
  mat <- as.matrix(mat_wide[,-1, drop = FALSE])
  rownames(mat) <- mat_wide$tf_name
  
  # row-wise z-score (mean 0, sd 1); set non-finite to 0
  mat_z <- t(scale(t(mat), center = TRUE, scale = TRUE))
  mat_z[!is.finite(mat_z)] <- 0
  
  # back to long + put subtype back for faceting
  hm_df <- as.data.frame(mat_z) %>%
    tibble::rownames_to_column("tf_name") %>%
    tidyr::pivot_longer(-tf_name, names_to = "sample_std", values_to = "z") %>%
    left_join(distinct(mat_df, sample_std, subtype), by = "sample_std")
  
  # ---- (3) ordering: fixed TF order; samples grouped by subtype ----
  hm_df <- hm_df %>%
    mutate(tf_name = factor(tf_name, levels = rev(tf_order)))  # reverse so top appears first
  
  # order samples within each subtype by their appearance in data (or alphabetically)
  samp_order_ar <- hm_df %>% filter(subtype == "ARPC")  %>% distinct(sample_std) %>% arrange(sample_std) %>% pull(sample_std)
  samp_order_ne <- hm_df %>% filter(subtype == "NEPC")  %>% distinct(sample_std) %>% arrange(sample_std) %>% pull(sample_std)
  # build a combined order so facets each use their own x-order cleanly
  samp_order_all <- c(samp_order_ar, samp_order_ne)
  
  hm_df <- hm_df %>%
    mutate(sample_std = factor(sample_std, levels = samp_order_all))
  
  # ---- (4) plot: facet by subtype----
  p_hm <- ggplot(hm_df, aes(x = sample_std, y = tf_name, fill = z)) +
    geom_tile() +
    facet_grid(. ~ subtype, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(name = "z-score", low = "#2166ac", mid = "white", high = "#b2182b",
                         limits = c(-2.5, 2.5), oob = scales::squish) +
    labs(
      title = "Top-N TFs: sample × TF z-score (contrast)",
      x = "Samples (grouped by subtype)", y = NULL
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank(),
      legend.title = element_text(),
      legend.position = "right"
    )
  
  ggsave("TF_heatmap.png", p_hm, width = 14, height = 7, dpi = 300)
  
  # export the z-score matrix used (optional)
  as_tibble(mat_z, rownames = "tf_name") %>%
    write_tsv("TF_heatmap_zscore.tsv")
  
  log_info("[S8] Written: TF_heatmap_fixed_faceted.png, TF_heatmap_fixed_faceted_zscore.tsv")
}
