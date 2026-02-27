#!/usr/bin/env Rscript
# ============================================================
#  run_dromics.R — DRomics/BMD dose-response analysis
#  Metadata columns: Sample_Name, Batch, Treatment, Dose, Type
# ============================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--counts",     type="character", help="featureCounts output"),
  make_option("--metadata",   type="character", help="Metadata file"),
  make_option("--treatment",  type="character", help="Treatment name"),
  make_option("--control",    type="character", help="Control name"),
  make_option("--outdir",     type="character", default="."),
  make_option("--fdr",        type="double",    default=0.05),
  make_option("--criterion",  type="character", default="AIC"),
  make_option("--bmd",        type="logical",   default=TRUE),
  make_option("--bootstrap",  type="logical",   default=FALSE),
  make_option("--combat_seq", type="character", default=NULL),
  make_option("--helper_seq", type="character", default=NULL)
)

opt <- parse_args(OptionParser(option_list=option_list))

if (!is.null(opt$combat_seq) && file.exists(opt$combat_seq)) source(opt$combat_seq)
if (!is.null(opt$helper_seq) && file.exists(opt$helper_seq)) source(opt$helper_seq)

suppressPackageStartupMessages({
  library(DESeq2); library(edgeR); library(sva)
  library(dplyr);  library(readxl); library(ggplot2)
  library(ggrepel); library(stringr); library(org.Hs.eg.db)
  library(DRomics); library(jsonlite); library(parallel)
})

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive=TRUE)

cat("=== DRomics/BMD Analysis ===\n")
cat("Treatment:", opt$treatment, "\n")
cat("Control:  ", opt$control, "\n\n")

# ── Helpers ────────────────────────────────────────────────────────────────────
load_metadata <- function(path) {
  if (grepl("\\.xlsx?$", path, ignore.case=TRUE))
    as.data.frame(readxl::read_excel(path))
  else read.csv(path, stringsAsFactors=FALSE)
}

read_featurecounts <- function(path) {
  cat("Reading counts:", path, "\n")
  
  first_lines <- readLines(path, n=3)
  is_fc_format <- any(grepl("^#", first_lines)) || any(grepl("Geneid.*Chr.*Start.*End.*Strand.*Length", first_lines))
  
  if (is_fc_format) {
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
    sep <- if (grepl("\\.csv$", path, ignore.case=TRUE)) "," else "\t"
    raw <- read.table(path, header=TRUE, sep=sep, stringsAsFactors=FALSE,
                      check.names=FALSE, row.names=1)
    mat <- raw
  }
  
  clean <- str_extract(colnames(mat), "SRR[0-9]+|ERR[0-9]+|DRR[0-9]+")
  still_na <- is.na(clean)
  if (any(still_na))
    clean[still_na] <- gsub("\\.[^.]*$","",basename(colnames(mat)[still_na]))
  colnames(mat) <- clean
  mat <- as.matrix(mat); mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  cat("Genes:", nrow(mat), "| Samples:", ncol(mat), "\n")
  mat
}

apply_batch_correction <- function(counts, metadata) {
  if (!"Batch" %in% colnames(metadata))
    return(list(counts=counts, corrected=FALSE))
  batch <- as.factor(metadata$Batch)
  if (nlevels(batch) < 2)
    return(list(counts=counts, corrected=FALSE))
  dose_fac <- as.factor(metadata$Dose)
  tryCatch({
    if (exists("ComBat_seq")) {
      corrected <- ComBat_seq(counts=counts, batch=batch, group=dose_fac, shrink=FALSE)
    } else {
      cpm_log  <- edgeR::cpm(counts, log=TRUE, prior.count=1)
      corr_log <- sva::ComBat(dat=cpm_log, batch=batch,
                               mod=model.matrix(~dose_fac), par.prior=TRUE)
      corrected <- round(pmax(2^corr_log-1, 0))
    }
    list(counts=corrected, corrected=TRUE)
  }, error=function(e) {
    cat("Batch correction failed:", e$message, "\n")
    list(counts=counts, corrected=FALSE)
  })
}

# ── Load data ──────────────────────────────────────────────────────────────────
counts_mat    <- read_featurecounts(opt$counts)
metadata_full <- load_metadata(opt$metadata)
rownames(metadata_full) <- metadata_full$Sample_Name

# Subset to treatment + control using Treatment column
meta_sub <- metadata_full %>%
  filter(Treatment %in% c(opt$treatment, opt$control)) %>%
  mutate(dose_numeric = if_else(Type=="control", 0, as.numeric(as.character(Dose))))

common     <- intersect(colnames(counts_mat), rownames(meta_sub))
cat("Common samples:", length(common), "\n")
if (length(common) == 0) stop("No common samples! Check Sample_Name vs count matrix columns.")

counts_sub <- counts_mat[, common]
meta_sub   <- meta_sub[common, ]

cat("Dose levels:", paste(sort(unique(meta_sub$dose_numeric)), collapse=", "), "\n\n")

# QC: low read filter
keep <- colSums(counts_sub) >= 1000000
counts_sub <- counts_sub[, keep]; meta_sub <- meta_sub[keep, ]

# QC: CPM filter per dose group
cpm_mat    <- edgeR::cpm(counts_sub)
keep_genes <- rep(FALSE, nrow(cpm_mat))
for (d in unique(meta_sub$dose_numeric)) {
  samps <- intersect(rownames(meta_sub[meta_sub$dose_numeric==d,]), colnames(cpm_mat))
  if (!length(samps)) next
  prop  <- rowSums(cpm_mat[,samps,drop=FALSE] >= 1) / length(samps)
  keep_genes[prop >= 0.75] <- TRUE
}
cat("Genes after CPM filter:", sum(keep_genes), "\n")
counts_sub <- counts_sub[keep_genes, ]

# Batch correction
bc           <- apply_batch_correction(counts_sub, meta_sub)
final_counts <- bc$counts
final_meta   <- meta_sub

cat("Final:", ncol(final_counts), "samples |", nrow(final_counts), "genes\n\n")

# Gene symbol mapping
gene_ids     <- rownames(final_counts)
gene_symbols <- tryCatch(
  mapIds(org.Hs.eg.db, keys=gene_ids, column="SYMBOL",
         keytype="ENSEMBL", multiVals="first"),
  error=function(e) setNames(gene_ids, gene_ids)
)
gene_symbols[is.na(gene_symbols)] <- names(gene_symbols)[is.na(gene_symbols)]
gene_map <- data.frame(ensembl_id=names(gene_symbols),
                        gene_symbol=gene_symbols, stringsAsFactors=FALSE)

# ── DRomics ────────────────────────────────────────────────────────────────────
dose_vec   <- final_meta$dose_numeric
samp_names <- rownames(final_meta)

dromics_data <- formatdata4DRomics(signalmatrix=final_counts,
                                    dose=dose_vec, samplenames=samp_names)

transfo <- if (ncol(final_counts) < 30) "rlog" else "vst"
cat("Transformation:", transfo, "\n")

o <- RNAseqdata(file=dromics_data, check=TRUE,
                transfo.method=transfo, transfo.blind=TRUE, round.counts=FALSE)

# QC plots
png(file.path(opt$outdir,"DRomics_QC.png"), width=800, height=600)
plot(o); dev.off()
png(file.path(opt$outdir,"DRomics_PCA.png"), width=800, height=600)
PCAdataplot(o, label=TRUE); dev.off()

# Gene selection
s <- itemselect(omicdata=o, select.method="quadratic", FDR=opt$fdr)
n_selected <- length(s$selectindex)
cat("Dose-responsive genes:", n_selected, "\n")

bmd_performed <- FALSE
bmd_summary   <- list()
n_fitted      <- 0

if (n_selected > 0) {
  # Save selected genes
  sel_genes <- rownames(o$data)[s$selectindex]
  sel_df    <- data.frame(
    ensembl_id  = sel_genes,
    gene_symbol = gene_map$gene_symbol[match(sel_genes, gene_map$ensembl_id)]
  )
  write.csv(sel_df, file.path(opt$outdir,"selected_responding_genes.csv"), row.names=FALSE)

  # Fit models
  f <- drcfit(itemselect=s, information.criterion=opt$criterion,
               progressbar=FALSE, parallel="no")
  n_fitted <- nrow(f$fitres)
  cat("Models fitted:", n_fitted, "\n")

  if (n_fitted > 0) {
    f$fitres$gene_symbol <- gene_map$gene_symbol[match(f$fitres$id, gene_map$ensembl_id)]
    write.csv(f$fitres, file.path(opt$outdir,"dose_response_models.csv"), row.names=FALSE)

    # Plot curves
    png(file.path(opt$outdir,"dose_response_curves.png"), width=1200, height=800)
    plot(f, items=min(12, n_fitted)); dev.off()

    # ── BMD ──────────────────────────────────────────────────────────────────
    if (opt$bmd) {
      cat("\n=== BMD Analysis ===\n")
      r_bmd <- tryCatch(bmdcalc(f=f, z=1.0, x=10), error=function(e) NULL)

      if (!is.null(r_bmd)) {
        if (opt$bootstrap) {
          cat("Bootstrap (500 iterations)...\n")
          n_cores <- min(parallel::detectCores()-1, 4)
          cl <- b_bmd <- NULL
          tryCatch({
            cl    <- parallel::makeCluster(n_cores)
            parallel::clusterEvalQ(cl, library(DRomics))
            parallel::clusterExport(cl, "r_bmd", envir=environment())
            b_bmd <- bmdboot(r=r_bmd, niter=500, progressbar=FALSE, cl=cl)
          }, error=function(e) cat("Parallel bootstrap failed\n"),
          finally = { if (!is.null(cl)) parallel::stopCluster(cl) })
          if (is.null(b_bmd))
            b_bmd <- tryCatch(bmdboot(r=r_bmd, niter=500, parallel="no"),
                               error=function(e) NULL)
          filtered <- if (!is.null(b_bmd))
            bmdfilter(b_bmd$res, BMDfilter="definedCI",  BMDtype="zSD")
          else
            bmdfilter(r_bmd$res, BMDfilter="definedBMD", BMDtype="zSD")
        } else {
          filtered <- bmdfilter(r_bmd$res, BMDfilter="definedBMD", BMDtype="zSD")
        }

        if (nrow(filtered) > 0) {
          filtered$gene_symbol <- gene_map$gene_symbol[match(filtered$id, gene_map$ensembl_id)]
          filtered$gene_symbol[is.na(filtered$gene_symbol)] <- filtered$id[is.na(filtered$gene_symbol)]
          filtered <- filtered[order(filtered$BMD.zSD, na.last=TRUE),]

          valid <- filtered[!is.na(filtered$BMD.zSD) & filtered$BMD.zSD > 0,]
          if (nrow(valid) > 0) {
            bmd_performed <- TRUE
            bmd_summary   <- list(
              total_bmds = nrow(valid),
              min_bmd    = min(valid$BMD.zSD),
              max_bmd    = max(valid$BMD.zSD),
              median_bmd = median(valid$BMD.zSD),
              mean_bmd   = mean(valid$BMD.zSD)
            )
            cat("Genes with valid BMDs:", nrow(valid), "\n")
            cat("Median BMD:", round(bmd_summary$median_bmd,4), "μM\n")

            write.csv(filtered,        file.path(opt$outdir,"bmd_results.csv"),              row.names=FALSE)
            write.csv(head(valid,50),  file.path(opt$outdir,"top50_sensitive_genes.csv"),    row.names=FALSE)

            # BMD distribution plot
            p_bmd <- ggplot(valid, aes(x=BMD.zSD)) +
              geom_histogram(bins=30, fill="#1a6b3a", alpha=0.8, color="white") +
              geom_vline(xintercept=median(valid$BMD.zSD), color="red", linetype="dashed", size=1) +
              labs(title="BMD Distribution", subtitle=paste(nrow(valid),"genes"),
                   x="BMD (μM)", y="Count") +
              theme_minimal()
            ggsave(file.path(opt$outdir,"bmd_distribution.png"), p_bmd, width=10, height=6, dpi=300)

            # Top sensitive genes bar plot
            if (nrow(valid) >= 10) {
              top25 <- head(valid, 25)
              p_sens <- ggplot(top25, aes(x=reorder(gene_symbol,-BMD.zSD), y=BMD.zSD)) +
                geom_col(fill="#2196F3", alpha=0.8) +
                coord_flip() +
                labs(title=paste("Top 25 Most Sensitive Genes —", opt$treatment),
                     x="Gene", y="BMD (μM)") +
                theme_minimal()
              ggsave(file.path(opt$outdir,"top25_sensitive_genes.png"), p_sens,
                     width=10, height=8, dpi=300)
            }
          }
        }
      }
    }
  }
}

# Save JSON summary
summary_stats <- list(
  treatment=opt$treatment, control=opt$control,
  total_genes_final=nrow(final_counts),
  total_samples_final=ncol(final_counts),
  dose_levels=length(unique(dose_vec)),
  genes_selected=n_selected, models_fitted=n_fitted,
  batch_corrected=bc$corrected,
  bmd_performed=bmd_performed,
  bmd_summary=if(bmd_performed) bmd_summary else NULL
)
write_json(summary_stats, file.path(opt$outdir,"dromics_summary.json"), pretty=TRUE)

cat("\n=== DRomics Complete ===\n")
cat("Output:", opt$outdir, "\n")
