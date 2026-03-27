#!/usr/bin/env Rscript
# ============================================================
#  run_deseq2.R — DESeq2/R-ODAF analysis (Phase 2)
#  Runs AFTER user reviews QC and confirms exclusions.
#
#  IMPORTANT — batch correction policy:
#    Batch correction is applied ONLY for PCA visualisation.
#    DESeq2 always receives RAW integer counts.
#    If a Batch column is present it is included in the DESeq2
#    design formula as a covariate instead.
# ============================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--counts",          type = "character", help = "Count matrix file"),
  make_option("--metadata",        type = "character", help = "Metadata file (CSV/XLSX)"),
  make_option("--treatment",       type = "character", help = "Treatment name"),
  make_option("--control",         type = "character", help = "Control name"),
  make_option("--outdir",          type = "character", default = "."),
  make_option("--fdr_strict",      type = "double",    default = 0.01),
  make_option("--fdr_relaxed",     type = "double",    default = 0.05),
  make_option("--log2fc",          type = "double",    default = 1.0),
  make_option("--read_thresh",     type = "integer",   default = 1000000L),
  make_option("--exclude_samples", type = "character", default = NULL)
)

opt <- parse_args(OptionParser(option_list = option_list))

# Source shared utilities (try scripts/ dir, then current dir)
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
utils_path <- file.path(script_dir, "utils.R")
if (!file.exists(utils_path)) utils_path <- file.path(getwd(), "utils.R")
if (!file.exists(utils_path)) stop("Cannot find utils.R")
source(utils_path)

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(org.Hs.eg.db)
})

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

cat("=== DESeq2/R-ODAF Analysis ===\n")
cat("Treatment:", opt$treatment, "\n")
cat("Control:  ", opt$control, "\n")
cat("Output:   ", opt$outdir, "\n\n")

# ── Load & prepare ─────────────────────────────────────────────────────────────
counts_mat    <- read_count_matrix(opt$counts)
metadata_full <- load_metadata(opt$metadata)
rownames(metadata_full) <- metadata_full$Sample_Name

meta_sub <- metadata_full %>%
  filter(Treatment %in% c(opt$treatment, opt$control)) %>%
  mutate(
    group     = if_else(Type == "control", "Control",
                        paste0("Treated_", gsub("\\.", "p", as.character(Dose)), "uM")),
    condition = if_else(Type == "control", "Control", "Treated")
  )

# ── Step 1: align & filter on RAW counts ──────────────────────────────────────
aligned  <- align_samples(counts_mat, meta_sub, opt$exclude_samples)
filtered <- filter_samples_and_genes(aligned$counts, aligned$meta,
                                      read_thresh = opt$read_thresh,
                                      group_col   = "group")

# ── Step 2: batch correction FOR PCA VISUALISATION ONLY ───────────────────────
# bc$counts must NEVER be passed to DESeq2.
# If a Batch column exists we include it as a covariate in the design instead.
has_batch     <- "Batch" %in% colnames(filtered$meta) &&
                  nlevels(as.factor(filtered$meta$Batch)) > 1
batch_corrected_pca <- FALSE

if (has_batch) {
  bc <- apply_batch_correction(filtered$counts, filtered$meta, group_col = "Type")
  batch_corrected_pca <- bc$corrected
  # PCA outlier detection on batch-corrected data (visualisation clarity only)
  pca_result <- detect_pca_outliers(bc$counts, filtered$meta, opt$outdir,
                                     group_col = "condition")
  cat("Note: PCA used batch-corrected counts for visualisation only.\n")
  cat("DESeq2 will receive RAW counts with Batch as a design covariate.\n\n")
} else {
  # No batch — PCA on raw counts directly
  pca_result <- detect_pca_outliers(filtered$counts, filtered$meta, opt$outdir,
                                     group_col = "condition")
}

# ── Step 3: RAW counts → DESeq2 ───────────────────────────────────────────────
final_counts <- filtered$counts   # always raw integer counts
final_meta   <- filtered$meta
cat("Final:", ncol(final_counts), "samples |", nrow(final_counts), "genes\n\n")

# ── Final PCA plot ────────────────────────────────────────────────────────────
# Uses batch-corrected counts (if batch exists) for cleaner visualisation.
# DESeq2 still receives RAW counts — this PCA is for display only.
if (has_batch && batch_corrected_pca) {
  cat("Final PCA: using batch-corrected VST counts (visualisation only)\n")
  pca_counts_for_plot <- bc$counts  # batch-corrected from Step 2
  pca_title <- "PCA \u2014 batch-corrected VST (visualisation only)"
} else {
  cat("Final PCA: using raw VST counts (no batch correction)\n")
  pca_counts_for_plot <- final_counts
  pca_title <- "PCA \u2014 VST counts"
}

pca2     <- prcomp(t(assay(vst(
  DESeqDataSetFromMatrix(pca_counts_for_plot, final_meta, ~1), blind = TRUE))))
pca_df           <- as.data.frame(pca2$x)
pca_df$condition <- final_meta$condition
pca_df$batch     <- as.factor(if (has_batch) final_meta$Batch else rep("1", nrow(pca_df)))
pct <- round(100 * pca2$sdev^2 / sum(pca2$sdev^2))

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = batch,
                              label = rownames(pca_df))) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 15) +
  labs(title = pca_title,
       subtitle = if (has_batch && batch_corrected_pca)
         "Note: DESeq2 receives raw counts with Batch as design covariate" else NULL,
       x = paste0("PC1: ", pct[1], "%"),
       y = paste0("PC2: ", pct[2], "%")) +
  theme_bw()
ggsave(file.path(opt$outdir, "Final_PCA_Plot.png"), p_pca,
       width = 12, height = 8, dpi = 300)

# ── DESeq2 design — include Batch as covariate if present ─────────────────────
# IMPORTANT: DESeq2 always receives final_counts (RAW integer counts), never bc$counts
cat("Confirming: DESeq2 input is RAW counts (", nrow(final_counts), "genes x",
    ncol(final_counts), "samples)\n")
final_meta$condition <- relevel(as.factor(final_meta$condition), ref = "Control")

if (has_batch) {
  final_meta$Batch <- as.factor(final_meta$Batch)
  cat("DESeq2 design: ~ Batch + condition\n")
  dds <- DESeqDataSetFromMatrix(countData = final_counts, colData = final_meta,
                                 design = ~ Batch + condition)
} else {
  cat("DESeq2 design: ~ condition\n")
  dds <- DESeqDataSetFromMatrix(countData = final_counts, colData = final_meta,
                                 design = ~ condition)
}

dds <- DESeq(dds)
res <- results(dds, name = "condition_Treated_vs_Control")

# Gene symbol mapping
all_gene_symbols <- resolve_gene_symbols(rownames(final_counts), count_matrix = counts_mat)

res_df             <- as.data.frame(res)
res_df$ensembl_id  <- rownames(res_df)
res_df             <- res_df[!is.na(res_df$padj), ]
res_df$gene_symbol <- all_gene_symbols[res_df$ensembl_id]
res_df             <- res_df[order(res_df$padj), ]

sig_df  <- res_df[res_df$padj < opt$fdr_relaxed & abs(res_df$log2FoldChange) > opt$log2fc, ]
up_df   <- sig_df[sig_df$log2FoldChange >  opt$log2fc, ]
down_df <- sig_df[sig_df$log2FoldChange < -opt$log2fc, ]

cat("DEGs (FDR<", opt$fdr_relaxed, ", |log2FC|>", opt$log2fc, "):", nrow(sig_df), "\n")
cat("Up:", nrow(up_df), "| Down:", nrow(down_df), "\n\n")

write.csv(res_df,  file.path(opt$outdir, "All_Results.csv"),          row.names = FALSE)
write.csv(sig_df,  file.path(opt$outdir, "Custom_Filtered_DEGs.csv"), row.names = FALSE)
write.csv(up_df,   file.path(opt$outdir, "Upregulated_Genes.csv"),    row.names = FALSE)
write.csv(down_df, file.path(opt$outdir, "Downregulated_Genes.csv"),  row.names = FALSE)

# ── Volcano ───────────────────────────────────────────────────────────────────
res_df$significance <- "Not significant"
res_df$significance[res_df$padj < opt$fdr_relaxed & res_df$log2FoldChange >  opt$log2fc] <- "Upregulated"
res_df$significance[res_df$padj < opt$fdr_relaxed & res_df$log2FoldChange < -opt$log2fc] <- "Downregulated"

top_genes <- res_df %>%
  filter(significance != "Not significant") %>% arrange(padj) %>% head(15)

p_vol <- ggplot(res_df, aes(log2FoldChange, -log10(padj), color = significance)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c(Upregulated     = "red",
                                 Downregulated   = "blue",
                                 "Not significant" = "grey70")) +
  geom_vline(xintercept = c(-opt$log2fc, opt$log2fc), linetype = "dashed") +
  geom_hline(yintercept = -log10(opt$fdr_relaxed),     linetype = "dashed") +
  geom_text_repel(data = top_genes,
                  aes(label = ifelse(!is.na(gene_symbol), gene_symbol, ensembl_id)),
                  size = 3, max.overlaps = 10) +
  labs(title = paste("Volcano:", opt$treatment, "vs", opt$control),
       x = "log2FC", y = "-log10(padj)") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(opt$outdir, "Volcano_Plot.png"), p_vol,
       width = 10, height = 8, dpi = 300)

# ── Dose-response DEG ─────────────────────────────────────────────────────────
unique_doses <- sort(unique(final_meta$Dose[final_meta$Type == "treatment"]))
if (length(unique_doses) > 1) {
  cat("Dose-response DEG...\n")
  final_meta$group <- relevel(as.factor(final_meta$group), ref = "Control")

  if (has_batch) {
    dds_dose <- DESeqDataSetFromMatrix(final_counts, final_meta, ~ Batch + group)
  } else {
    dds_dose <- DESeqDataSetFromMatrix(final_counts, final_meta, ~ group)
  }
  dds_dose <- DESeq(dds_dose)

  for (comp in resultsNames(dds_dose)[grepl("^group_Treated", resultsNames(dds_dose))]) {
    r             <- as.data.frame(results(dds_dose, name = comp))
    r$ensembl_id  <- rownames(r)
    r$gene_symbol <- all_gene_symbols[r$ensembl_id]
    r             <- r[!is.na(r$padj), ]
    clean_name    <- gsub("group_|_vs_Control", "", comp)
    write.csv(r, file.path(opt$outdir, paste0("Dose_DEG_", clean_name, ".csv")),
              row.names = FALSE)
  }
}

# ── Summary JSON ──────────────────────────────────────────────────────────────
write_json(list(
  treatment              = opt$treatment,
  control                = opt$control,
  total_genes            = nrow(res_df),
  sig_strict             = sum(res_df$padj < opt$fdr_strict, na.rm = TRUE),
  sig_relaxed            = nrow(sig_df),
  upregulated            = nrow(up_df),
  downregulated          = nrow(down_df),
  samples_analyzed       = ncol(final_counts),
  genes_analyzed         = nrow(final_counts),
  batch_corrected        = FALSE,          # DESeq2 always gets raw counts
  batch_as_covariate     = has_batch,      # TRUE if Batch included in design
  batch_corrected_pca    = batch_corrected_pca,
  pca_outliers_flagged   = length(pca_result$outliers),
  doses_analyzed         = length(unique_doses)
), file.path(opt$outdir, "analysis_summary.json"), pretty = TRUE, auto_unbox = TRUE)

cat("\n=== DESeq2 Complete ===\n")
