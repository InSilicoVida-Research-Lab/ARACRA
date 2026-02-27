#!/usr/bin/env Rscript
# ============================================================
#  run_deseq2.R — DESeq2/R-ODAF analysis
#  Metadata columns: Sample_Name, Batch, Treatment, Dose, Type
# ============================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--counts",      type="character", help="featureCounts output"),
  make_option("--metadata",    type="character", help="Metadata file (CSV/XLSX)"),
  make_option("--treatment",   type="character", help="Treatment name"),
  make_option("--control",     type="character", help="Control name"),
  make_option("--outdir",      type="character", default="."),
  make_option("--fdr_strict",  type="double",    default=0.01),
  make_option("--fdr_relaxed", type="double",    default=0.05),
  make_option("--log2fc",      type="double",    default=1.0),
  make_option("--read_thresh", type="integer",   default=1000000L),
  make_option("--combat_seq",  type="character", default=NULL),
  make_option("--helper_seq",  type="character", default=NULL)
)

opt <- parse_args(OptionParser(option_list=option_list))

# Source ComBat_seq helpers
if (!is.null(opt$combat_seq) && file.exists(opt$combat_seq)) {
  source(opt$combat_seq)
  cat("ComBat_seq loaded\n")
} else {
  cat("ComBat_seq not found — will use sva::ComBat\n")
}
if (!is.null(opt$helper_seq) && file.exists(opt$helper_seq)) {
  source(opt$helper_seq)
}

suppressPackageStartupMessages({
  library(DESeq2); library(edgeR); library(sva)
  library(dplyr);  library(readxl); library(ggplot2)
  library(ggrepel); library(stringr); library(org.Hs.eg.db)
  library(jsonlite)
})

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive=TRUE)

cat("=== DESeq2/R-ODAF Analysis ===\n")
cat("Treatment:", opt$treatment, "\n")
cat("Control:  ", opt$control,   "\n")
cat("Output:   ", opt$outdir,    "\n\n")

# ── Load metadata ──────────────────────────────────────────────────────────────
load_metadata <- function(path) {
  if (grepl("\\.xlsx?$", path, ignore.case=TRUE))
    as.data.frame(readxl::read_excel(path))
  else
    read.csv(path, stringsAsFactors=FALSE)
}

# ── Load featureCounts ─────────────────────────────────────────────────────────
read_featurecounts <- function(path) {
  cat("Reading counts:", path, "\n")
  
  # Detect format: featureCounts .out vs cleaned CSV/TSV
  first_lines <- readLines(path, n=3)
  is_fc_format <- any(grepl("^#", first_lines)) || any(grepl("Geneid.*Chr.*Start.*End.*Strand.*Length", first_lines))
  
  if (is_fc_format) {
    # Standard featureCounts output (tab-separated with annotation columns)
    skip <- sum(grepl("^#", first_lines))
    raw <- tryCatch(
      read.table(path, header=TRUE, sep="\t", comment.char="",
                 stringsAsFactors=FALSE, skip=skip, check.names=FALSE),
      error=function(e)
        read.table(path, header=TRUE, sep="\t", comment.char="",
                   stringsAsFactors=FALSE, skip=skip, fill=TRUE, check.names=FALSE)
    )
    len_idx <- if ("Length" %in% colnames(raw)) which(colnames(raw)=="Length") else 6
    mat     <- raw[, (len_idx+1):ncol(raw), drop=FALSE]
    rownames(mat) <- if ("Geneid" %in% colnames(raw)) raw$Geneid else raw[,1]
  } else {
    # Cleaned CSV/TSV (gene_id + sample columns)
    sep <- if (grepl("\\.csv$", path, ignore.case=TRUE)) "," else "\t"
    raw <- read.table(path, header=TRUE, sep=sep, stringsAsFactors=FALSE, 
                      check.names=FALSE, row.names=1)
    mat <- raw
  }
  
  # Clean column names to SRR IDs
  clean <- str_extract(colnames(mat), "SRR[0-9]+|ERR[0-9]+|DRR[0-9]+")
  still_na <- is.na(clean)
  if (any(still_na))
    clean[still_na] <- gsub("\\.[^.]*$","", basename(colnames(mat)[still_na]))
  colnames(mat) <- clean
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  cat("Genes:", nrow(mat), "| Samples:", ncol(mat), "\n")
  mat
}

# ── Batch correction ───────────────────────────────────────────────────────────
apply_batch_correction <- function(counts, metadata) {
  # Your metadata uses 'Batch' column
  if (!"Batch" %in% colnames(metadata)) {
    cat("No Batch column — skipping batch correction\n")
    return(list(counts=counts, corrected=FALSE))
  }
  batch <- as.factor(metadata$Batch)
  if (nlevels(batch) < 2) {
    cat("Only 1 batch — skipping\n")
    return(list(counts=counts, corrected=FALSE))
  }
  # Use Type column for group
  grp <- as.factor(metadata$Type)
  tryCatch({
    if (exists("ComBat_seq")) {
      corrected <- ComBat_seq(counts=counts, batch=batch, group=grp, shrink=FALSE)
      cat("ComBat_seq correction applied\n")
    } else {
      cpm_log   <- edgeR::cpm(counts, log=TRUE, prior.count=1)
      corr_log  <- sva::ComBat(dat=cpm_log, batch=batch,
                                mod=model.matrix(~grp), par.prior=TRUE)
      corrected <- round(pmax(2^corr_log - 1, 0))
      cat("sva::ComBat correction applied\n")
    }
    list(counts=corrected, corrected=TRUE)
  }, error=function(e) {
    cat("Batch correction failed:", e$message, "— using original\n")
    list(counts=counts, corrected=FALSE)
  })
}

# ── Main ───────────────────────────────────────────────────────────────────────
counts_mat    <- read_featurecounts(opt$counts)
metadata_full <- load_metadata(opt$metadata)

# Your metadata uses Sample_Name as the ID column
rownames(metadata_full) <- metadata_full$Sample_Name

# Subset to treatment vs control using Treatment and Type columns
meta_sub <- metadata_full %>%
  filter(Treatment %in% c(opt$treatment, opt$control)) %>%
  mutate(
    group     = if_else(Type == "control", "Control",
                        paste0("Treated_", gsub("\\.", "p", as.character(Dose)), "uM")),
    condition = if_else(Type == "control", "Control", "Treated")
  )

# Align samples between count matrix and metadata
common       <- intersect(colnames(counts_mat), rownames(meta_sub))
cat("Common samples between counts and metadata:", length(common), "\n")
if (length(common) == 0) stop("No common samples found! Check Sample_Name vs count matrix column names.")

counts_sub   <- counts_mat[, common]
meta_sub     <- meta_sub[common, ]

# Low-read filter
lib_size <- colSums(counts_sub)
keep_samp <- lib_size >= opt$read_thresh
cat("Samples passing read filter (>=", opt$read_thresh, "):", sum(keep_samp), "/", length(keep_samp), "\n")
counts_sub <- counts_sub[, keep_samp]
meta_sub   <- meta_sub[keep_samp, ]

# CPM gene filter
cpm_mat    <- edgeR::cpm(counts_sub)
keep_genes <- rep(FALSE, nrow(cpm_mat))
meta_sub$group <- as.factor(meta_sub$group)
for (grp in levels(meta_sub$group)) {
  samps <- intersect(rownames(meta_sub[meta_sub$group==grp,]), colnames(cpm_mat))
  if (!length(samps)) next
  prop  <- rowSums(cpm_mat[,samps,drop=FALSE] >= 1) / length(samps)
  keep_genes[prop >= 0.75] <- TRUE
}
cat("Genes after CPM filter:", sum(keep_genes), "/", nrow(counts_sub), "\n")
counts_sub <- counts_sub[keep_genes, ]

# Batch correction
bc <- apply_batch_correction(counts_sub, meta_sub)

# PCA + outlier detection
dds_qc <- DESeqDataSetFromMatrix(countData=bc$counts, colData=meta_sub, design=~1)
vst_qc <- vst(dds_qc, blind=TRUE)

# Outlier detection per condition group
pca      <- prcomp(t(assay(vst_qc)))
pca_data <- as.data.frame(pca$x[,1:2])
outliers <- c()
for (grp in unique(meta_sub$condition)) {
  samps <- rownames(meta_sub[meta_sub$condition==grp,])
  if (length(samps) < 3) next
  gdat  <- pca_data[samps,,drop=FALSE]
  dists <- sqrt(rowSums(sweep(gdat, 2, colMeans(gdat), '-')^2))
  thr   <- mean(dists) + 2.0*sd(dists)
  outliers <- c(outliers, names(dists[dists > thr]))
}
cat("Outliers removed:", length(outliers), "\n")

keep_samp2  <- setdiff(rownames(meta_sub), outliers)
final_counts <- bc$counts[, keep_samp2]
final_meta   <- meta_sub[keep_samp2, ]

cat("Final: ", ncol(final_counts), "samples |", nrow(final_counts), "genes\n\n")

# PCA plot
pca2     <- prcomp(t(assay(vst(DESeqDataSetFromMatrix(
              final_counts, final_meta, ~1), blind=TRUE))))
pca_df   <- as.data.frame(pca2$x)
pca_df$condition <- final_meta$condition
pca_df$batch     <- as.factor(final_meta$Batch)
pct      <- round(100*pca2$sdev^2/sum(pca2$sdev^2))
p_pca <- ggplot(pca_df, aes(PC1, PC2, color=condition, shape=batch,
                             label=rownames(pca_df))) +
  geom_point(size=4, alpha=0.8) +
  geom_text_repel(size=3, max.overlaps=15) +
  labs(title="PCA — Final QC",
       x=paste0("PC1: ",pct[1],"% variance"),
       y=paste0("PC2: ",pct[2],"% variance")) +
  theme_bw()
ggsave(file.path(opt$outdir,"Final_PCA_Plot.png"), p_pca, width=12, height=8, dpi=300)
cat("PCA plot saved\n")

# DESeq2
final_meta$condition <- relevel(as.factor(final_meta$condition), ref="Control")
dds <- DESeqDataSetFromMatrix(countData=final_counts,
                               colData=final_meta, design=~condition)
dds <- DESeq(dds)

res    <- results(dds, name="condition_Treated_vs_Control")
res_df <- as.data.frame(res)
res_df$ensembl_id <- rownames(res_df)
res_df <- res_df[!is.na(res_df$padj),]

# Gene symbol mapping
gene_symbols <- tryCatch(
  mapIds(org.Hs.eg.db, keys=res_df$ensembl_id, column="SYMBOL",
         keytype="ENSEMBL", multiVals="first"),
  error=function(e) setNames(rep(NA_character_, nrow(res_df)), res_df$ensembl_id)
)
res_df$gene_symbol <- gene_symbols[res_df$ensembl_id]
res_df <- res_df[order(res_df$padj),]

# Filter DEGs
sig_df  <- res_df[!is.na(res_df$padj) &
                  res_df$padj < opt$fdr_relaxed &
                  abs(res_df$log2FoldChange) > opt$log2fc, ]
up_df   <- sig_df[sig_df$log2FoldChange >  opt$log2fc, ]
down_df <- sig_df[sig_df$log2FoldChange < -opt$log2fc, ]

cat("Total genes tested:", nrow(res_df), "\n")
cat("Significant DEGs (FDR<", opt$fdr_relaxed, "& |log2FC|>", opt$log2fc, "):", nrow(sig_df), "\n")
cat("Upregulated:", nrow(up_df), "\n")
cat("Downregulated:", nrow(down_df), "\n\n")

# Save CSVs
write.csv(res_df,   file.path(opt$outdir,"All_Results.csv"),           row.names=FALSE)
write.csv(sig_df,   file.path(opt$outdir,"Custom_Filtered_DEGs.csv"),  row.names=FALSE)
write.csv(up_df,    file.path(opt$outdir,"Upregulated_Genes.csv"),     row.names=FALSE)
write.csv(down_df,  file.path(opt$outdir,"Downregulated_Genes.csv"),   row.names=FALSE)

# Volcano plot
res_df$significance <- "Not significant"
res_df$significance[!is.na(res_df$padj) &
                    res_df$padj < opt$fdr_relaxed &
                    res_df$log2FoldChange >  opt$log2fc] <- "Upregulated"
res_df$significance[!is.na(res_df$padj) &
                    res_df$padj < opt$fdr_relaxed &
                    res_df$log2FoldChange < -opt$log2fc] <- "Downregulated"

top_genes <- res_df %>%
  filter(significance %in% c("Upregulated","Downregulated")) %>%
  arrange(padj) %>% head(15)

p_vol <- ggplot(res_df, aes(log2FoldChange, -log10(padj), color=significance)) +
  geom_point(alpha=0.5, size=1) +
  scale_color_manual(values=c(Upregulated="red", Downregulated="blue",
                               "Not significant"="grey70")) +
  geom_vline(xintercept=c(-opt$log2fc, opt$log2fc), linetype="dashed", alpha=0.6) +
  geom_hline(yintercept=-log10(opt$fdr_relaxed),    linetype="dashed", alpha=0.6) +
  geom_text_repel(data=top_genes,
                  aes(label=ifelse(!is.na(gene_symbol), gene_symbol, ensembl_id)),
                  size=3, max.overlaps=10) +
  labs(title=paste("Volcano:", opt$treatment, "vs", opt$control),
       x="log2 Fold Change", y="-log10(adj. p-value)", color="") +
  theme_bw() + theme(legend.position="bottom")
ggsave(file.path(opt$outdir,"Volcano_Plot.png"), p_vol, width=10, height=8, dpi=300)
cat("Volcano plot saved\n")

# Dose-response DEG (each dose vs control)
unique_doses <- sort(unique(final_meta$Dose[final_meta$Type=="treatment"]))
if (length(unique_doses) > 1) {
  cat("\nRunning dose-response DEG analysis...\n")
  final_meta$group <- relevel(as.factor(final_meta$group), ref="Control")
  dds_dose <- DESeqDataSetFromMatrix(final_counts, final_meta, ~group)
  dds_dose <- DESeq(dds_dose)
  for (comp in resultsNames(dds_dose)[grepl("^group_Treated", resultsNames(dds_dose))]) {
    r   <- as.data.frame(results(dds_dose, name=comp))
    r$ensembl_id  <- rownames(r)
    r$gene_symbol <- gene_symbols[r$ensembl_id]
    r <- r[!is.na(r$padj),]
    clean_name <- gsub("group_|_vs_Control","",comp)
    write.csv(r, file.path(opt$outdir, paste0("Dose_DEG_",clean_name,".csv")), row.names=FALSE)
    cat("  Saved:", clean_name, "— sig:", sum(r$padj < opt$fdr_relaxed, na.rm=TRUE), "\n")
  }
}

# Save JSON summary
summary <- list(
  treatment=opt$treatment, control=opt$control,
  total_genes=nrow(res_df),
  sig_strict=sum(!is.na(res_df$padj) & res_df$padj < opt$fdr_strict),
  sig_relaxed=nrow(sig_df),
  upregulated=nrow(up_df), downregulated=nrow(down_df),
  samples_analyzed=ncol(final_counts),
  genes_analyzed=nrow(final_counts),
  batch_corrected=bc$corrected,
  outliers_removed=length(outliers),
  doses_analyzed=length(unique_doses)
)
write_json(summary, file.path(opt$outdir,"analysis_summary.json"), pretty=TRUE)

cat("\n=== DESeq2 Complete ===\n")
cat("Results saved to:", opt$outdir, "\n")
