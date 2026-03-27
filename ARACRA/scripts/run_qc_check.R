#!/usr/bin/env Rscript
# ============================================================
#  run_qc_check.R — PCA-based QC outlier detection only
#
#  Runs filtering + (optional) batch correction + PCA,
#  writes JSON for the app to display to the user.
#
#  Does NOT run DESeq2 or DRomics — that happens in Phase 2
#  after the user reviews and confirms sample exclusions.
#
#  Batch correction here is used SOLELY for PCA visualisation
#  so that batch effects don't obscure outlier detection.
#  The corrected matrix is never passed downstream.
# ============================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--counts",      type = "character", help = "Count matrix file"),
  make_option("--metadata",    type = "character", help = "Metadata file (CSV/XLSX)"),
  make_option("--treatment",   type = "character", help = "Treatment name"),
  make_option("--control",     type = "character", help = "Control name"),
  make_option("--outdir",      type = "character", default = "."),
  make_option("--read_thresh", type = "integer",   default = 1000000L),
  make_option("--mode",        type = "character", default = "deseq2",
              help = "deseq2 or dromics — affects grouping column")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Source shared utilities
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
utils_path <- file.path(script_dir, "utils.R")
if (!file.exists(utils_path)) utils_path <- file.path(getwd(), "utils.R")
if (!file.exists(utils_path)) stop("Cannot find utils.R — must be in scripts/ directory")
source(utils_path)

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

cat("=== QC Outlier Check ===\n")
cat("Treatment:", opt$treatment, "\n")
cat("Control:  ", opt$control, "\n\n")

# ── Load data ──────────────────────────────────────────────────────────────────
counts_mat    <- read_count_matrix(opt$counts)
metadata_full <- load_metadata(opt$metadata)
rownames(metadata_full) <- metadata_full$Sample_Name

if (opt$mode == "dromics") {
  meta_sub <- metadata_full %>%
    filter(Treatment %in% c(opt$treatment, opt$control)) %>%
    mutate(
      dose_numeric = if_else(Type == "control", 0, as.numeric(as.character(Dose))),
      condition    = if_else(Type == "control", "Control",
                             paste0("Dose_", as.character(Dose)))
    )
} else {
  meta_sub <- metadata_full %>%
    filter(Treatment %in% c(opt$treatment, opt$control)) %>%
    mutate(
      group     = if_else(Type == "control", "Control",
                          paste0("Treated_", gsub("\\.", "p", as.character(Dose)), "uM")),
      condition = if_else(Type == "control", "Control", "Treated")
    )
}

# ── Step 1: align & filter on RAW counts (no exclusions at QC stage) ──────────
aligned  <- align_samples(counts_mat, meta_sub, exclude_samples_str = NULL)
filtered <- filter_samples_and_genes(
  aligned$counts, aligned$meta,
  read_thresh = opt$read_thresh,
  group_col   = if (opt$mode == "dromics") "condition" else "group"
)

# ── Step 2: batch correction FOR PCA VISUALISATION ONLY ───────────────────────
# The corrected counts are used exclusively to produce a cleaner PCA for
# outlier detection. They are discarded after this block.
has_batch <- "Batch" %in% colnames(filtered$meta) &&
              nlevels(as.factor(filtered$meta$Batch)) > 1

if (has_batch) {
  bc <- apply_batch_correction(
    filtered$counts, filtered$meta,
    group_col = if (opt$mode == "dromics") "condition" else "Type"
  )
  # Use batch-corrected counts only for PCA — better outlier visibility
  counts_for_pca  <- bc$counts
  batch_corrected <- bc$corrected
  cat("Batch correction applied for PCA visualisation only.\n")
  cat("Raw counts will be passed to DESeq2/DRomics in Phase 2.\n\n")
} else {
  counts_for_pca  <- filtered$counts
  batch_corrected <- FALSE
  cat("No batch column detected — PCA on raw VST counts.\n\n")
}

# ── Step 3: PCA outlier detection — writes pca_outlier_flag.json ──────────────
pca_result <- detect_pca_outliers(
  counts_for_pca, filtered$meta, opt$outdir,
  group_col = "condition"
)

# ── Summary JSON ───────────────────────────────────────────────────────────────
summary_out <- list(
  stage              = "qc_check",
  mode               = opt$mode,
  total_samples      = ncol(filtered$counts),   # raw count dimensions
  total_genes        = nrow(filtered$counts),
  batch_corrected    = FALSE,                    # downstream scripts get raw counts
  batch_corrected_pca = batch_corrected,         # batch correction used for PCA only
  pca_outliers       = length(pca_result$outliers),
  flagged_samples    = as.list(pca_result$outliers),
  all_samples        = colnames(filtered$counts) # raw sample list for exclusion UI
)
write_json(summary_out, file.path(opt$outdir, "qc_summary.json"),
           pretty = TRUE, auto_unbox = TRUE)

cat("\n=== QC Check Complete ===\n")
cat("Samples:", ncol(filtered$counts), "| Genes:", nrow(filtered$counts), "\n")
cat("Batch correction for PCA:", batch_corrected, "\n")
cat("Outliers flagged:", length(pca_result$outliers), "\n")
if (length(pca_result$outliers) > 0)
  cat("  Flagged:", paste(pca_result$outliers, collapse = ", "), "\n")
cat("Output:", opt$outdir, "\n")
