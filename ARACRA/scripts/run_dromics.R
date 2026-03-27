#!/usr/bin/env Rscript
# ============================================================
#  run_dromics.R — DRomics/BMD dose-response analysis (Phase 2)
#  Runs AFTER user reviews QC and confirms exclusions.
#
#  IMPORTANT — batch correction policy:
#    Batch correction is applied ONLY for PCA visualisation.
#    DRomics always receives RAW integer counts and performs
#    its own VST/rlog transformation internally via RNAseqdata().
# ============================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--counts",          type = "character", help = "Count matrix file"),
  make_option("--metadata",        type = "character", help = "Metadata file"),
  make_option("--treatment",       type = "character", help = "Treatment name"),
  make_option("--control",         type = "character", help = "Control name"),
  make_option("--outdir",          type = "character", default = "."),
  make_option("--fdr",             type = "double",    default = 0.05),
  make_option("--criterion",       type = "character", default = "AICc"),
  make_option("--bmd",             type = "logical",   default = TRUE),
  make_option("--bootstrap",       type = "logical",   default = FALSE),
  make_option("--select_method",   type = "character", default = "quadratic",
              help = "Gene selection method: quadratic, linear, or ANOVA"),
  make_option("--transfo_method",  type = "character", default = "auto",
              help = "Transformation: auto (vst/rlog by sample count), vst, or rlog"),
  make_option("--bmd_z",           type = "double",    default = 1.0,
              help = "BMD z-value (signal-to-noise ratio threshold)"),
  make_option("--bmd_x",           type = "double",    default = 10.0,
              help = "BMD x-value (percent change threshold)"),
  make_option("--niter",           type = "integer",   default = 20000L,
              help = "Bootstrap iterations for BMD confidence intervals"),
  make_option("--exclude_samples", type = "character", default = NULL),
  make_option("--bmdu_bmdl_ratio", type = "double",    default = 40,
              help = "BMDU/BMDL ratio threshold (NTP=40, EFSA=50)"),
  make_option("--fold_change_min", type = "double",    default = 0,
              help = "Min absolute VST/rlog difference (approx log2FC) for gene pre-filter (0=disabled, NTP recommends ~1.5 FC = 0.585)"),
  make_option("--bmd_max_dose_filter", type = "logical", default = TRUE,
              help = "Remove genes with BMD > highest tested dose (NTP 2018)"),
  make_option("--bmd_extrap_factor",   type = "double",  default = 10,
              help = "Flag genes with BMD < lowest_dose / this factor (NTP 2018)"),
  make_option("--read_thresh",         type = "integer",  default = 0L,
              help = "Min total reads per sample (0=disabled; default 0 for TempO-Seq compatibility)")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Source shared utilities
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
utils_path <- file.path(script_dir, "utils.R")
if (!file.exists(utils_path)) utils_path <- file.path(getwd(), "utils.R")
if (!file.exists(utils_path)) stop("Cannot find utils.R")
source(utils_path)

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(org.Hs.eg.db)
  library(DRomics); library(parallel); library(AnnotationDbi)
  if (requireNamespace("GO.db", quietly = TRUE)) library(GO.db)
})

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

cat("=== DRomics/BMD Analysis ===\n")
cat("Treatment:", opt$treatment, "\n")
cat("Control:  ", opt$control, "\n\n")

# ── Load & prepare ─────────────────────────────────────────────────────────────
counts_mat    <- read_count_matrix(opt$counts)
metadata_full <- load_metadata(opt$metadata)
rownames(metadata_full) <- metadata_full$Sample_Name

meta_sub <- metadata_full %>%
  filter(Treatment %in% c(opt$treatment, opt$control)) %>%
  mutate(
    dose_numeric = if_else(Type == "control", 0, as.numeric(as.character(Dose))),
    dose_group   = paste0("dose_", dose_numeric),
    condition    = if_else(Type == "control", "Control",
                           paste0("Dose_", as.character(Dose)))
  )

# ── Step 1: align & filter on RAW counts ──────────────────────────────────────
aligned  <- align_samples(counts_mat, meta_sub, opt$exclude_samples)
filtered <- filter_samples_and_genes(aligned$counts, aligned$meta,
                                      read_thresh = opt$read_thresh,
                                      group_col = "dose_group")

# ── Step 2: batch correction FOR PCA VISUALISATION ONLY ───────────────────────
# bc$counts must NEVER be passed to DRomics.
# RNAseqdata() expects raw integer counts and applies VST/rlog internally.
has_batch           <- "Batch" %in% colnames(filtered$meta) &&
                        nlevels(as.factor(filtered$meta$Batch)) > 1
batch_corrected_pca <- FALSE

if (has_batch) {
  bc <- apply_batch_correction(filtered$counts, filtered$meta,
                                group_col = "dose_group")
  batch_corrected_pca <- bc$corrected
  # PCA outlier detection on batch-corrected data (visualisation clarity only)
  pca_result <- detect_pca_outliers(bc$counts, filtered$meta, opt$outdir,
                                     group_col = "condition")
  cat("Note: PCA used batch-corrected counts for visualisation only.\n")
  cat("DRomics will receive RAW counts.\n\n")
} else {
  # No batch — PCA on raw counts directly
  pca_result <- detect_pca_outliers(filtered$counts, filtered$meta, opt$outdir,
                                     group_col = "condition")
}

# ── Step 3: RAW counts → DRomics ──────────────────────────────────────────────
final_counts <- filtered$counts   # always raw integer counts
final_meta   <- filtered$meta

cat("Final:", ncol(final_counts), "samples |", nrow(final_counts), "genes\n")
cat("Doses:", paste(sort(unique(final_meta$dose_numeric)), collapse = ", "), "\n")
cat("Confirming: DRomics input is RAW counts (", nrow(final_counts), "genes x",
    ncol(final_counts), "samples)\n")
cat("DRomics will apply its own VST/rlog transformation internally.\n\n")

# Gene symbol mapping
gene_ids     <- rownames(final_counts)
gene_symbols <- resolve_gene_symbols(gene_ids, count_matrix = counts_mat)
gene_map     <- data.frame(ensembl_id  = gene_ids,
                            gene_symbol = gene_symbols,
                            stringsAsFactors = FALSE)

# ── DRomics ────────────────────────────────────────────────────────────────────
# Verify sample ordering is consistent between counts and metadata
stopifnot("Column/row order mismatch before DRomics" =
            identical(colnames(final_counts), rownames(final_meta)))

dose_vec   <- as.numeric(as.character(final_meta$dose_numeric))
samp_names <- rownames(final_meta)

cat("dose_vec length:", length(dose_vec), "| ncol:", ncol(final_counts), "\n")
stopifnot(length(dose_vec) == ncol(final_counts))

dromics_data <- formatdata4DRomics(signalmatrix = final_counts,
                                    dose = dose_vec, samplenames = samp_names)

# Auto-select transformation based on sample count (30 is DESeq2's own threshold)
transfo <- if (opt$transfo_method == "auto") {
  if (ncol(final_counts) < 30) "rlog" else "vst"
} else {
  opt$transfo_method
}
cat("Transformation:", transfo, "\n")
cat("Selection method:", opt$select_method, "\n")

# RNAseqdata() applies VST/rlog to RAW counts internally
o <- RNAseqdata(file = dromics_data, check = TRUE,
                transfo.method = transfo, transfo.blind = TRUE,
                round.counts = FALSE)

png(file.path(opt$outdir, "DRomics_QC.png"),  width = 800, height = 600); plot(o);             dev.off()
png(file.path(opt$outdir, "DRomics_PCA.png"), width = 800, height = 600); PCAdataplot(o, label = TRUE); dev.off()

s <- itemselect(omicdata = o, select.method = opt$select_method, FDR = opt$fdr)
n_selected <- length(s$selectindex)
cat("Dose-responsive genes (trend test):", n_selected, "\n")

# ── Gap 7: Compute fold-change for post-fit filtering (NTP 2018) ──────────
# NTP recommends combining trend test with fold-change filter for reproducibility.
# We compute FC here but apply the filter AFTER drcfit() to avoid modifying
# the itemselect object (which has internal consistency requirements).
fc_filter_genes <- NULL  # will hold gene IDs that PASS the FC filter
n_fc_filtered <- 0
if (opt$fold_change_min > 0 && n_selected > 0) {
  cat("Computing fold-change for post-fit filter: min |log2FC| >=", opt$fold_change_min, "\n")
  transfo_data <- o$data  # transformed expression matrix (genes x samples)
  ctrl_cols <- which(dose_vec == 0)
  if (length(ctrl_cols) > 0) {
    ctrl_mean <- rowMeans(transfo_data[, ctrl_cols, drop = FALSE])
    non_zero_doses <- sort(unique(dose_vec[dose_vec > 0]))
    # Compute max |log2FC| for each selected gene
    sel_indices <- s$selectindex
    max_abs_fc <- sapply(sel_indices, function(idx) {
      gene_vals <- transfo_data[idx, ]
      max(sapply(non_zero_doses, function(d) {
        dose_cols <- which(dose_vec == d)
        abs(mean(gene_vals[dose_cols]) - ctrl_mean[idx])
      }))
    })
    fc_pass <- max_abs_fc >= opt$fold_change_min
    fc_filter_genes <- rownames(transfo_data)[sel_indices[fc_pass]]
    n_fc_filtered <- sum(!fc_pass)
    cat("  Genes passing FC filter:", sum(fc_pass), "of", length(sel_indices),
        "| Will be applied after model fitting\n")
  }
}

bmd_performed <- FALSE
bmd_summary   <- list()
n_fitted      <- 0

if (n_selected > 0) {
  sel_genes <- rownames(o$data)[s$selectindex]
  sel_df    <- data.frame(
    ensembl_id  = sel_genes,
    gene_symbol = gene_map$gene_symbol[match(sel_genes, gene_map$ensembl_id)]
  )
  write.csv(sel_df, file.path(opt$outdir, "selected_responding_genes.csv"), row.names = FALSE)

  n_cores <- max(1, parallel::detectCores() - 1)
  cat("Using", n_cores, "cores for parallel processing\n")

  f <- drcfit(itemselect = s, information.criterion = opt$criterion,
               progressbar = FALSE, parallel = "snow", ncpus = n_cores)
  n_fitted <- nrow(f$fitres)
  cat("Models fitted (before FC filter):", n_fitted, "\n")

  # Apply fold-change filter on fitted results
  if (!is.null(fc_filter_genes) && n_fitted > 0) {
    fc_keep <- f$fitres$id %in% fc_filter_genes
    n_fc_removed_from_fits <- sum(!fc_keep)
    if (n_fc_removed_from_fits > 0) {
      f$fitres <- f$fitres[fc_keep, ]
      if (!is.null(f$fits)) f$fits <- f$fits[fc_keep]
      n_fitted <- nrow(f$fitres)
      cat("  FC filter removed:", n_fc_removed_from_fits, "fitted genes |",
          n_fitted, "remaining\n")
    }
  }

  if (n_fitted > 0) {
    f$fitres$gene_symbol <- gene_map$gene_symbol[match(f$fitres$id, gene_map$ensembl_id)]
    write.csv(f$fitres, file.path(opt$outdir, "dose_response_models.csv"), row.names = FALSE)

    # ── Model quality metrics: pseudo-R² and SNR ───────────────────────────
    cat("Computing model quality metrics...\n")
    transfo_data <- o$data
    quality_metrics <- data.frame(
      id          = f$fitres$id,
      gene_symbol = f$fitres$gene_symbol,
      model       = f$fitres$model,
      trend       = if ("trend" %in% colnames(f$fitres)) f$fitres$trend else NA,
      SDres       = f$fitres$SDres,
      stringsAsFactors = FALSE
    )

    # Compute total variance per gene and pseudo-R²
    gene_var_total <- sapply(f$fitres$id, function(gid) {
      if (gid %in% rownames(transfo_data)) {
        var(as.numeric(transfo_data[gid, ]))
      } else NA
    })
    quality_metrics$var_total  <- round(gene_var_total, 6)
    quality_metrics$pseudo_R2  <- round(1 - (f$fitres$SDres^2 / gene_var_total), 4)
    quality_metrics$pseudo_R2[quality_metrics$pseudo_R2 < 0] <- 0  # floor at 0

    # Signal-to-noise ratio = maxychange / SDres
    if ("maxychange" %in% colnames(f$fitres)) {
      quality_metrics$maxychange <- round(f$fitres$maxychange, 4)
      quality_metrics$SNR        <- round(abs(f$fitres$maxychange) / f$fitres$SDres, 3)
    } else {
      quality_metrics$maxychange <- NA
      quality_metrics$SNR        <- NA
    }

    # yrange as fraction of y0 (effect size relative to control)
    if (all(c("yrange", "y0") %in% colnames(f$fitres))) {
      quality_metrics$effect_pct <- round(100 * f$fitres$yrange / abs(f$fitres$y0), 2)
    }

    write.csv(quality_metrics, file.path(opt$outdir, "model_quality_metrics.csv"),
              row.names = FALSE)

    # Aggregate quality stats for JSON summary
    model_quality_summary <- list(
      n_fitted             = n_fitted,
      median_pseudo_R2     = round(median(quality_metrics$pseudo_R2, na.rm = TRUE), 3),
      mean_pseudo_R2       = round(mean(quality_metrics$pseudo_R2, na.rm = TRUE), 3),
      median_SNR           = round(median(quality_metrics$SNR, na.rm = TRUE), 3),
      mean_SNR             = round(mean(quality_metrics$SNR, na.rm = TRUE), 3),
      pct_R2_above_0.5     = round(100 * mean(quality_metrics$pseudo_R2 >= 0.5, na.rm = TRUE), 1),
      pct_R2_above_0.3     = round(100 * mean(quality_metrics$pseudo_R2 >= 0.3, na.rm = TRUE), 1),
      pct_SNR_above_1      = round(100 * mean(quality_metrics$SNR >= 1, na.rm = TRUE), 1),
      pct_SNR_above_2      = round(100 * mean(quality_metrics$SNR >= 2, na.rm = TRUE), 1),
      median_SDres         = round(median(f$fitres$SDres, na.rm = TRUE), 4),
      model_distribution   = as.list(table(f$fitres$model)),
      trend_distribution   = as.list(table(if ("trend" %in% colnames(f$fitres)) f$fitres$trend else c()))
    )

    cat("  Median pseudo-R²:", model_quality_summary$median_pseudo_R2, "\n")
    cat("  Median SNR:", model_quality_summary$median_SNR, "\n")
    cat("  Genes with R² >= 0.5:", model_quality_summary$pct_R2_above_0.5, "%\n")
    cat("  Genes with SNR >= 2:", model_quality_summary$pct_SNR_above_2, "%\n")

    # ── Relabel ENSG IDs with gene symbols before the DRomics default plot ──
    # This keeps the correct parametric curves while showing gene names
    rownames(o$data) <- ifelse(
      !is.na(gene_map$gene_symbol[match(rownames(o$data), gene_map$ensembl_id)]) &
        nchar(gene_map$gene_symbol[match(rownames(o$data), gene_map$ensembl_id)]) > 0,
      gene_map$gene_symbol[match(rownames(o$data), gene_map$ensembl_id)],
      rownames(o$data)
    )
    # Also update f so plot() uses the new labels
    # IMPORTANT: preserve ENSEMBL→symbol mapping BEFORE overwriting $id,
    # so pathway-level tPOD can map BMD gene IDs back to ENSEMBL.
    symbol_to_ensembl <- setNames(gene_map$ensembl_id, gene_map$gene_symbol)
    fitres_original_ids <- f$fitres$id  # preserve ENSEMBL IDs before overwriting
    f$fitres$id <- gene_map$gene_symbol[match(f$fitres$id, gene_map$ensembl_id)]
    f$fitres$id[is.na(f$fitres$id)] <- fitres_original_ids[is.na(f$fitres$id)]

    # ── Publishable dose-response plot ──────────────────────────────────────
    # Strategy: extract fitted curves from DRomics internal model objects
    # (f$fits), predict on a linear dose grid from 0→max, plot with
    # gene symbol labels on a linear x-axis.
    cat("Generating publishable dose-response plot...\n")

    # Save DRomics default as reference (now with gene symbols)
    # Use linear dose scale so dose=0 (controls) are visible
    png(file.path(opt$outdir, "dose_response_curves_dromics.png"),
        width = 1200, height = 800)
    plot(f, items = min(12, n_fitted), dose_log_transfo = FALSE); dev.off()

    tryCatch({
      transfo_data <- o$data
      dose_raw     <- final_meta$dose_numeric

      n_plot      <- min(12, n_fitted)
      plot_ids    <- f$fitres$id[1:n_plot]
      plot_models <- f$fitres$model[1:n_plot]
      plot_trends <- if ("trend" %in% colnames(f$fitres)) f$fitres$trend[1:n_plot] else rep("", n_plot)

      # ── Observation data on linear dose scale ──
      obs_list <- list()
      for (k in seq_len(n_plot)) {
        gid <- plot_ids[k]
        if (!(gid %in% rownames(transfo_data))) next
        trend_txt <- if (nchar(plot_trends[k]) > 0) paste0(", ", plot_trends[k]) else ""
        obs_list[[k]] <- data.frame(
          gene_label = paste0(gid, " (", plot_models[k], trend_txt, ")"),
          dose       = dose_raw,
          signal     = as.numeric(transfo_data[gid, ]),
          is_control = (dose_raw == 0),
          stringsAsFactors = FALSE
        )
      }
      obs_data <- do.call(rbind, obs_list)
      obs_data$gene_label <- factor(obs_data$gene_label, levels = unique(obs_data$gene_label))

      # ── Extract fitted curves from DRomics model objects ──
      dose_min <- min(dose_raw[dose_raw > 0])
      dose_max <- max(dose_raw)
      # Linear dose grid from near-0 to max (500 points for smooth curves)
      dose_grid_nz <- seq(dose_min / 5, dose_max * 1.02, length.out = 500)
      # DRomics models work on log(dose), so compute log-dose grid
      log_dose_grid <- log(dose_grid_nz)

      curve_list <- list()
      for (k in seq_len(n_plot)) {
        gid   <- plot_ids[k]
        model <- as.character(plot_models[k])
        trend_txt <- if (nchar(plot_trends[k]) > 0) paste0(", ", plot_trends[k]) else ""

        pred_y <- tryCatch({
          # f$fits is a list of model objects (nls or lm)
          # DRomics models use 'x' as the predictor (= log(dose))
          fit_obj <- f$fits[[k]]
          predict(fit_obj, newdata = data.frame(x = log_dose_grid))
        }, error = function(e) NULL)

        if (!is.null(pred_y) && length(pred_y) == length(dose_grid_nz)) {
          # Extrapolate to dose=0 using the value at lowest prediction point
          y_at_zero <- pred_y[1]  # value at dose_min/5

          curve_list[[k]] <- data.frame(
            gene_label = paste0(gid, " (", model, trend_txt, ")"),
            dose       = c(0, dose_grid_nz),
            signal     = c(y_at_zero, pred_y),
            stringsAsFactors = FALSE
          )
        }
      }

      has_curves <- !sapply(curve_list, is.null)
      if (any(has_curves)) {
        curve_data <- do.call(rbind, curve_list[has_curves])
        curve_data$gene_label <- factor(curve_data$gene_label,
                                         levels = levels(obs_data$gene_label))
      } else {
        curve_data <- NULL
      }

      # ── Build plot ──
      dose_breaks <- sort(unique(c(0, dose_raw[dose_raw > 0])))
      jitter_w <- dose_max * 0.012

      p_dr <- ggplot() +
        geom_point(data = obs_data,
                   aes(x = dose, y = signal, shape = is_control),
                   size = 2, alpha = 0.7,
                   position = position_jitter(width = jitter_w, seed = 42)) +
        scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 16),
                           labels = c("TRUE" = "Control", "FALSE" = "Treated"),
                           name = NULL)

      if (!is.null(curve_data) && nrow(curve_data) > 0) {
        p_dr <- p_dr +
          geom_line(data = curve_data, aes(x = dose, y = signal),
                    color = "#cc0000", linewidth = 0.9)
      } else {
        # Fallback: smooth curve through data
        p_dr <- p_dr +
          geom_smooth(data = obs_data, aes(x = dose, y = signal),
                      method = "loess", formula = y ~ x, se = FALSE,
                      color = "#cc0000", linewidth = 0.8, span = 0.75)
      }

      p_dr <- p_dr +
        facet_wrap(~gene_label, scales = "free_y", ncol = 4) +
        scale_x_continuous(breaks = dose_breaks, labels = dose_breaks,
                           expand = expansion(mult = c(0.02, 0.05))) +
        labs(
          x = expression(paste("Dose (", mu, "M)")),
          y = "Normalized Expression (VST/rlog)",
          title = paste("Dose-Response Curves \u2014", opt$treatment),
          subtitle = paste0(n_selected, " responsive genes (FDR<", opt$fdr,
                            ") | top ", n_plot, " shown")
        ) +
        theme_bw(base_size = 11) +
        theme(
          strip.text       = element_text(size = 8.5, face = "bold"),
          strip.background = element_rect(fill = "#f0f4f8", color = "grey80"),
          panel.grid.minor = element_blank(),
          plot.title       = element_text(size = 13, face = "bold"),
          plot.subtitle    = element_text(size = 9, color = "grey40"),
          axis.text.x      = element_text(size = 8),
          legend.position  = "bottom"
        )

      ggsave(file.path(opt$outdir, "dose_response_curves.png"), p_dr,
             width = 14, height = max(6, 3.5 * ceiling(n_plot / 4)), dpi = 300)
      cat("Publishable dose-response plot saved (linear axis from 0)\n")

    }, error = function(e) {
      cat("Custom plot failed:", e$message, "\n")
      cat("Falling back to DRomics default plot\n")
      file.copy(file.path(opt$outdir, "dose_response_curves_dromics.png"),
                file.path(opt$outdir, "dose_response_curves.png"),
                overwrite = TRUE)
    })

    # ── Model fit summary table (with quality metrics) ─────────────────────
    model_summary <- merge(
      f$fitres[, intersect(c("id", "gene_symbol", "model", "SDres", "trend"),
                            colnames(f$fitres))],
      quality_metrics[, intersect(c("id", "pseudo_R2", "SNR", "effect_pct"),
                                   colnames(quality_metrics))],
      by = "id", all.x = TRUE
    )
    write.csv(model_summary, file.path(opt$outdir, "model_fit_summary.csv"), row.names = FALSE)

    cat("\nModel type distribution:\n"); print(table(f$fitres$model))
    cat("\nTrend distribution:\n");      print(table(f$fitres$trend))

    # ── BMD Analysis ────────────────────────────────────────────────────────
    if (opt$bmd) {
      cat("\n=== BMD Analysis ===\n")
      cat("z =", opt$bmd_z, "| x =", opt$bmd_x, "| bootstrap:", opt$bootstrap, "\n")
      r_bmd <- tryCatch(bmdcalc(f = f, z = opt$bmd_z, x = opt$bmd_x), error = function(e) NULL)

      if (!is.null(r_bmd)) {
        if (opt$bootstrap) {
          cat("Bootstrap (", opt$niter, " iterations) on", n_cores, "cores...\n")
          cl <- b_bmd <- NULL
          tryCatch({
            cl <- parallel::makeCluster(n_cores)
            parallel::clusterEvalQ(cl, library(DRomics))
            parallel::clusterExport(cl, "r_bmd", envir = environment())
            b_bmd <- bmdboot(r = r_bmd, niter = opt$niter, progressbar = FALSE, cl = cl)
            cat("Bootstrap completed successfully\n")
          }, error = function(e) cat("Parallel bootstrap failed:", e$message, "\n"),
          finally = { if (!is.null(cl)) parallel::stopCluster(cl) })
          if (is.null(b_bmd)) {
            cat("Retrying bootstrap without parallelization...\n")
            b_bmd <- tryCatch(bmdboot(r = r_bmd, niter = opt$niter, parallel = "no"),
                               error = function(e) NULL)
          }
          filtered_bmd <- if (!is.null(b_bmd))
            bmdfilter(b_bmd$res, BMDfilter = "definedCI",  BMDtype = "zSD")
          else
            bmdfilter(r_bmd$res, BMDfilter = "definedBMD", BMDtype = "zSD")
        } else {
          filtered_bmd <- bmdfilter(r_bmd$res, BMDfilter = "definedBMD", BMDtype = "zSD")
        }

        if (nrow(filtered_bmd) > 0) {
          # filtered_bmd$id may contain gene symbols (from relabelling) or ENSEMBL IDs
          # Try ENSEMBL match first, then check if ID is already a symbol
          filtered_bmd$gene_symbol <- gene_map$gene_symbol[match(filtered_bmd$id, gene_map$ensembl_id)]
          # For IDs that are already symbols (post-relabelling), use them directly
          still_na <- is.na(filtered_bmd$gene_symbol)
          if (any(still_na)) {
            is_known_symbol <- filtered_bmd$id[still_na] %in% gene_map$gene_symbol
            filtered_bmd$gene_symbol[still_na][is_known_symbol] <-
              filtered_bmd$id[still_na][is_known_symbol]
            # Final fallback: use the ID as-is
            filtered_bmd$gene_symbol[is.na(filtered_bmd$gene_symbol)] <-
              filtered_bmd$id[is.na(filtered_bmd$gene_symbol)]
          }
          filtered_bmd <- filtered_bmd[order(filtered_bmd$BMD.zSD, na.last = TRUE), ]
          valid <- filtered_bmd[!is.na(filtered_bmd$BMD.zSD) & filtered_bmd$BMD.zSD > 0, ]

          if (nrow(valid) > 0) {
            bmd_performed <- TRUE

            # ── Dose range for NTP quality filters ───────────────────────────
            max_dose <- max(dose_vec[dose_vec > 0])
            min_dose <- min(dose_vec[dose_vec > 0])

            # ── Gap 2: Remove BMD > highest tested dose (NTP 2018) ───────────
            n_above_max <- 0
            if (opt$bmd_max_dose_filter) {
              above_max <- valid$BMD.zSD > max_dose
              n_above_max <- sum(above_max, na.rm = TRUE)
              if (n_above_max > 0) {
                cat("NTP max-dose filter: removed", n_above_max, "genes with BMD >",
                    round(max_dose, 4), "(highest tested dose)\n")
                valid <- valid[!above_max, ]
              }
            }

            # ── Gap 3: Flag BMD < lowest_dose / extrap_factor (NTP 2018) ─────
            extrap_threshold <- min_dose / opt$bmd_extrap_factor
            valid$bmd_extrapolated <- valid$BMD.zSD < extrap_threshold
            n_extrapolated <- sum(valid$bmd_extrapolated, na.rm = TRUE)
            if (n_extrapolated > 0) {
              cat("NTP extrapolation flag:", n_extrapolated, "genes with BMD <",
                  round(extrap_threshold, 6), "(lowest dose /", opt$bmd_extrap_factor, ")\n")
              cat("  These are flagged but NOT removed (used in gene-level tPOD,\n")
              cat("  but pathway tPODs driven by flagged genes will be noted)\n")
            }

            # ── Gap 1: BMDU/BMDL ratio filter (configurable threshold) ───────
            # NTP 2018: ratio < 40 | EFSA 2022: ratio < 50
            has_ci <- "BMD.zSD.lower" %in% colnames(valid) &
                       "BMD.zSD.upper" %in% colnames(valid)
            ratio_threshold <- opt$bmdu_bmdl_ratio
            n_before_ratio <- nrow(valid)
            n_ratio_filtered <- 0

            if (has_ci) {
              valid$BMDL <- valid$BMD.zSD.lower
              valid$BMDU <- valid$BMD.zSD.upper
              valid$BMDU_BMDL_ratio <- ifelse(
                !is.na(valid$BMDL) & !is.na(valid$BMDU) &
                is.finite(valid$BMDL) & is.finite(valid$BMDU) &
                valid$BMDL > 0,
                valid$BMDU / valid$BMDL, NA)

              ratio_fail <- !is.na(valid$BMDU_BMDL_ratio) &
                             valid$BMDU_BMDL_ratio > ratio_threshold
              n_ratio_filtered <- sum(ratio_fail, na.rm = TRUE)

              cat("BMDU/BMDL ratio filter (>", ratio_threshold, "):",
                  n_ratio_filtered, "of", n_before_ratio, "genes removed\n")

              valid$ratio_pass <- !ratio_fail | is.na(ratio_fail)
              valid_efsa <- valid[valid$ratio_pass, ]
            } else {
              cat("No bootstrap CI available — skipping BMDU/BMDL ratio filter\n")
              valid$BMDL <- NA
              valid$BMDU <- NA
              valid$BMDU_BMDL_ratio <- NA
              valid$ratio_pass <- TRUE
              valid_efsa <- valid
            }

            # ── NTP quality filter summary ───────────────────────────────────
            ntp_summary <- list(
              max_dose_filter    = opt$bmd_max_dose_filter,
              n_removed_maxdose  = n_above_max,
              max_dose           = max_dose,
              extrap_factor      = opt$bmd_extrap_factor,
              extrap_threshold   = extrap_threshold,
              n_extrapolated     = n_extrapolated,
              ratio_threshold    = ratio_threshold,
              n_removed_ratio    = n_ratio_filtered,
              has_bootstrap_ci   = has_ci,
              n_fold_change_removed = n_fc_filtered,
              fold_change_min    = opt$fold_change_min,
              select_method      = opt$select_method
            )

            # ── Distribution-based tPOD methods ──────────────────────────────
            # Use EFSA-filtered genes for tPOD computation
            efsa_bmds <- sort(valid_efsa$BMD.zSD)
            n_efsa <- length(efsa_bmds)

            tpod_methods <- list()

            # Method 1: 25th ranked gene BMD (Reardon et al. 2021)
            if (n_efsa >= 25) {
              tpod_methods$rank25 <- list(
                value = round(efsa_bmds[25], 6),
                label = "25th Ranked Gene BMD",
                ref   = "Reardon et al. 2021",
                n_genes_used = n_efsa
              )
            }

            # Method 2: 5th percentile (Farmahin et al. 2017)
            if (n_efsa >= 25) {
              tpod_methods$perc05 <- list(
                value = round(quantile(efsa_bmds, 0.05), 6),
                label = "5th Percentile BMD",
                ref   = "Farmahin et al. 2017",
                n_genes_used = n_efsa
              )
            }

            # Method 3: 10th percentile
            if (n_efsa >= 10) {
              tpod_methods$perc10 <- list(
                value = round(quantile(efsa_bmds, 0.10), 6),
                label = "10th Percentile BMD",
                ref   = "Johnson et al. 2022",
                n_genes_used = n_efsa
              )
            }

            # Method 4: First mode (ECDF max curvature) — simplified kneedle
            if (n_efsa >= 25) {
              log_bmds <- log10(efsa_bmds)
              ecdf_y   <- seq_along(log_bmds) / n_efsa
              # Normalize to [0,1]
              xn <- (log_bmds - min(log_bmds)) / max(1e-12, max(log_bmds) - min(log_bmds))
              yn <- ecdf_y
              # Kneedle: max distance from diagonal
              dists <- abs(yn - xn)
              knee_idx <- which.max(dists)
              tpod_methods$first_mode <- list(
                value = round(efsa_bmds[knee_idx], 6),
                label = "First Mode (Max Curvature)",
                ref   = "Pagé-Larivière et al. 2019",
                n_genes_used = n_efsa
              )
            }

            # Method 5: Median BMD of all genes
            tpod_methods$median_all <- list(
              value = round(median(efsa_bmds), 6),
              label = "Median BMD (All Genes)",
              ref   = "Thomas et al. 2013",
              n_genes_used = n_efsa
            )

            cat("Distribution-based tPOD methods computed:", length(tpod_methods), "\n")
            for (nm in names(tpod_methods)) {
              cat("  ", tpod_methods[[nm]]$label, ":", tpod_methods[[nm]]$value, "uM\n")
            }

            bmd_summary   <- list(
              total_bmds  = nrow(valid),
              total_bmds_filtered = n_efsa,
              ntp_quality_filters = ntp_summary,
              has_bootstrap_ci = has_ci,
              min_bmd     = min(valid$BMD.zSD),
              max_bmd     = max(valid$BMD.zSD),
              median_bmd  = median(valid$BMD.zSD),
              mean_bmd    = mean(valid$BMD.zSD),
              tpod_methods = tpod_methods
            )
            cat("Valid BMDs:", nrow(valid),
                "| After quality filters:", n_efsa,
                "| Median:", round(bmd_summary$median_bmd, 4), "uM\n")

            write.csv(filtered_bmd, file.path(opt$outdir, "bmd_results.csv"),        row.names = FALSE)
            write.csv(head(valid, 50), file.path(opt$outdir, "top50_sensitive_genes.csv"), row.names = FALSE)

            p_bmd <- ggplot(valid, aes(x = BMD.zSD)) +
              geom_histogram(bins = 30, fill = "#1a6b3a", alpha = 0.8, color = "white") +
              geom_vline(xintercept = median(valid$BMD.zSD), color = "red", linetype = "dashed") +
              labs(title = "BMD Distribution",
                   subtitle = paste0("Median BMD = ", round(bmd_summary$median_bmd, 3), " \u03bcM"),
                   x = expression(paste("BMD (", mu, "M)")), y = "Count") +
              theme_minimal()
            ggsave(file.path(opt$outdir, "bmd_distribution.png"), p_bmd,
                   width = 10, height = 6, dpi = 300)

            if (nrow(valid) >= 10) {
              top25   <- head(valid, 25)
              p_sens  <- ggplot(top25, aes(x = reorder(gene_symbol, -BMD.zSD), y = BMD.zSD)) +
                geom_col(fill = "#2196F3", alpha = 0.8) + coord_flip() +
                labs(title    = paste("Top 25 Most Sensitive Genes \u2014", opt$treatment),
                     subtitle = "Ordered by BMD (lower = more sensitive)",
                     x = "", y = expression(paste("BMD (", mu, "M)"))) +
                theme_minimal()
              ggsave(file.path(opt$outdir, "top25_sensitive_genes.png"), p_sens,
                     width = 10, height = 8, dpi = 300)
            }
          }
        }
      }
    }
  }
}

# ══════════════════════════════════════════════════════════════════════════════
#  PATHWAY-LEVEL tPOD (NTP 2018 / EPA 2024 gene set-based approach)
# ══════════════════════════════════════════════════════════════════════════════

tpod_performed <- FALSE
tpod_summary   <- list()
n_pathways_pass <- 0

if (bmd_performed && nrow(valid_efsa) >= 5) {
  cat("\n=== Pathway-level tPOD Analysis (NTP 2018) ===\n")
  cat("Valid BMD genes (EFSA-filtered):", nrow(valid_efsa), "\n")

  # valid$id may contain gene symbols (relabelled for plotting).
  # Map them back to ENSEMBL IDs for annotation DB lookups.
  if (exists("symbol_to_ensembl")) {
    bmd_gene_ids <- ifelse(
      valid_efsa$id %in% gene_map$ensembl_id,
      valid_efsa$id,  # already ENSEMBL
      ifelse(!is.na(symbol_to_ensembl[valid_efsa$id]),
             symbol_to_ensembl[valid_efsa$id],
             valid_efsa$id)
    )
  } else {
    bmd_gene_ids <- valid_efsa$id
  }
  cat("ENSEMBL IDs resolved for pathway mapping:", sum(grepl("^ENSG", bmd_gene_ids)),
      "/", length(bmd_gene_ids), "\n")

  # ── Step 1: Map ALL human ENSEMBL genes to GO BP (for pathway sizes) ──────
  # Then subset to our BMD genes. This avoids unreliable reverse lookups.
  cat("Step 1: Mapping genes to GO Biological Process...\n")

  all_go <- tryCatch({
    # Get ALL human ENSEMBL → GO mappings (this is the reliable direction)
    all_ensembl <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
    cat("  Total ENSEMBL IDs in org.Hs.eg.db:", length(all_ensembl), "\n")

    raw <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys    = all_ensembl,
                                  columns = c("ENSEMBL", "GO", "ONTOLOGY"),
                                  keytype = "ENSEMBL")
    # Keep only Biological Process with valid GO IDs
    bp <- raw[!is.na(raw$GO) & !is.na(raw$ONTOLOGY) & raw$ONTOLOGY == "BP",
              c("ENSEMBL", "GO")]
    bp <- unique(bp)
    cat("  GO BP mappings:", nrow(bp), "(", length(unique(bp$GO)), "GO terms )\n")
    bp
  }, error = function(e) {
    cat("  GO mapping FAILED:", e$message, "\n")
    NULL
  })

  # ── Step 1b: KEGG mapping ─────────────────────────────────────────────────
  cat("Step 1b: Mapping genes to KEGG pathways...\n")

  all_kegg <- tryCatch({
    all_ensembl <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
    raw_k <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys    = all_ensembl,
                                    columns = c("ENSEMBL", "PATH"),
                                    keytype = "ENSEMBL")
    kp <- raw_k[!is.na(raw_k$PATH), c("ENSEMBL", "PATH")]
    kp$PATH <- paste0("hsa", kp$PATH)
    kp <- unique(kp)
    cat("  KEGG mappings:", nrow(kp), "(", length(unique(kp$PATH)), "pathways )\n")
    kp
  }, error = function(e) {
    cat("  KEGG mapping FAILED:", e$message, "\n")
    NULL
  })

  # ── Step 1c: MSigDB Hallmark (H) + C2 (curated) gene sets ────────────────
  # Uses msigdbr package for programmatic access to MSigDB with ENSEMBL mappings
  all_hallmark <- NULL
  all_c2       <- NULL

  if (requireNamespace("msigdbr", quietly = TRUE)) {
    cat("Step 1c: Loading MSigDB gene sets via msigdbr...\n")

    # Hallmark (H) — 50 curated, non-redundant gene sets
    tryCatch({
      h_sets <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
      h_map  <- unique(h_sets[, c("ensembl_gene", "gs_name")])
      h_map  <- h_map[nchar(h_map$ensembl_gene) > 0, ]
      colnames(h_map) <- c("ENSEMBL", "PATHWAY")
      all_hallmark <- h_map
      cat("  Hallmark mappings:", nrow(h_map), "(",
          length(unique(h_map$PATHWAY)), "gene sets )\n")
    }, error = function(e) {
      cat("  Hallmark loading FAILED:", e$message, "\n")
    })

    # C2 (curated gene sets) — canonical pathways + chemical/genetic perturbations
    tryCatch({
      c2_sets <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2")
      c2_map  <- unique(c2_sets[, c("ensembl_gene", "gs_name", "gs_subcat")])
      c2_map  <- c2_map[nchar(c2_map$ensembl_gene) > 0, ]
      colnames(c2_map) <- c("ENSEMBL", "PATHWAY", "SUBCAT")
      all_c2 <- c2_map
      n_cp  <- length(unique(c2_map$PATHWAY[c2_map$SUBCAT %in%
                  c("CP", "CP:BIOCARTA", "CP:KEGG", "CP:PID",
                    "CP:REACTOME", "CP:WIKIPATHWAYS")]))
      n_cgp <- length(unique(c2_map$PATHWAY[c2_map$SUBCAT == "CGP"]))
      cat("  C2 mappings:", nrow(c2_map), "(",
          length(unique(c2_map$PATHWAY)), "gene sets:",
          n_cp, "CP +", n_cgp, "CGP )\n")
    }, error = function(e) {
      cat("  C2 loading FAILED:", e$message, "\n")
    })
  } else {
    cat("  msigdbr package not available — skipping MSigDB collections\n")
    cat("  Install with: install.packages('msigdbr')\n")
  }

  # ── Step 2: Compute pathway sizes from ALL genes ──────────────────────────
  pathway_sizes <- list()
  if (!is.null(all_go)) {
    go_sizes <- table(all_go$GO)
    for (nm in names(go_sizes)) pathway_sizes[[nm]] <- as.integer(go_sizes[nm])
    cat("  GO pathway sizes computed:", length(go_sizes), "terms\n")
  }
  if (!is.null(all_kegg)) {
    kegg_sizes <- table(all_kegg$PATH)
    for (nm in names(kegg_sizes)) pathway_sizes[[nm]] <- as.integer(kegg_sizes[nm])
    cat("  KEGG pathway sizes computed:", length(kegg_sizes), "pathways\n")
  }
  if (!is.null(all_hallmark)) {
    h_sizes <- table(all_hallmark$PATHWAY)
    for (nm in names(h_sizes)) pathway_sizes[[nm]] <- as.integer(h_sizes[nm])
    cat("  Hallmark pathway sizes computed:", length(h_sizes), "sets\n")
  }
  if (!is.null(all_c2)) {
    c2_sizes <- table(all_c2$PATHWAY)
    for (nm in names(c2_sizes)) pathway_sizes[[nm]] <- as.integer(c2_sizes[nm])
    cat("  C2 pathway sizes computed:", length(c2_sizes), "sets\n")
  }

  # ── Step 3: Build gene-to-pathway mapping for BMD genes only ──────────────
  gene2pathway <- data.frame(gene_id = character(), pathway_id = character(),
                              category = character(), stringsAsFactors = FALSE)
  if (!is.null(all_go)) {
    bmd_go <- all_go[all_go$ENSEMBL %in% bmd_gene_ids, ]
    if (nrow(bmd_go) > 0) {
      bmd_go_df <- data.frame(gene_id = bmd_go$ENSEMBL, pathway_id = bmd_go$GO,
                               category = "GO:BP", stringsAsFactors = FALSE)
      gene2pathway <- rbind(gene2pathway, bmd_go_df)
    }
  }
  if (!is.null(all_kegg)) {
    bmd_kegg <- all_kegg[all_kegg$ENSEMBL %in% bmd_gene_ids, ]
    if (nrow(bmd_kegg) > 0) {
      bmd_kegg_df <- data.frame(gene_id = bmd_kegg$ENSEMBL, pathway_id = bmd_kegg$PATH,
                                 category = "KEGG", stringsAsFactors = FALSE)
      gene2pathway <- rbind(gene2pathway, bmd_kegg_df)
    }
  }
  if (!is.null(all_hallmark)) {
    bmd_h <- all_hallmark[all_hallmark$ENSEMBL %in% bmd_gene_ids, ]
    if (nrow(bmd_h) > 0) {
      bmd_h_df <- data.frame(gene_id = bmd_h$ENSEMBL, pathway_id = bmd_h$PATHWAY,
                              category = "Hallmark", stringsAsFactors = FALSE)
      gene2pathway <- rbind(gene2pathway, bmd_h_df)
    }
  }
  if (!is.null(all_c2)) {
    bmd_c2 <- all_c2[all_c2$ENSEMBL %in% bmd_gene_ids, ]
    if (nrow(bmd_c2) > 0) {
      # Label CP subcollections as "C2:CP" and CGP as "C2:CGP"
      bmd_c2$cat_label <- ifelse(
        bmd_c2$SUBCAT %in% c("CP", "CP:BIOCARTA", "CP:KEGG", "CP:PID",
                              "CP:REACTOME", "CP:WIKIPATHWAYS"),
        "C2:CP", "C2:CGP")
      bmd_c2_df <- data.frame(gene_id = bmd_c2$ENSEMBL, pathway_id = bmd_c2$PATHWAY,
                               category = bmd_c2$cat_label, stringsAsFactors = FALSE)
      gene2pathway <- rbind(gene2pathway, bmd_c2_df)
    }
  }

  gene2pathway <- unique(gene2pathway)
  cat("BMD gene-pathway mappings:", nrow(gene2pathway),
      "(", length(unique(gene2pathway$pathway_id)), "unique pathways )\n")

  if (nrow(gene2pathway) > 0) {

    # ── Step 4: Merge BMD values ────────────────────────────────────────────
    # gene2pathway$gene_id contains ENSEMBL IDs; valid$id may be gene symbols.
    # Build lookups keyed by ENSEMBL ID.
    bmd_lookup <- setNames(valid_efsa$BMD.zSD, bmd_gene_ids)
    sym_lookup <- setNames(valid_efsa$gene_symbol, bmd_gene_ids)
    gene2pathway$bmd <- bmd_lookup[gene2pathway$gene_id]
    gene2pathway$gene_symbol <- sym_lookup[gene2pathway$gene_id]
    gene2pathway <- gene2pathway[!is.na(gene2pathway$bmd), ]

    # ── Step 5: Compute per-pathway statistics ──────────────────────────────
    pathway_ids <- unique(gene2pathway$pathway_id)
    cat("Computing pathway-level BMD summaries for", length(pathway_ids), "pathways...\n")

    pathway_results <- do.call(rbind, lapply(pathway_ids, function(pid) {
      members <- gene2pathway[gene2pathway$pathway_id == pid, ]
      # Deduplicate: one gene might map to same pathway multiple times
      members <- members[!duplicated(members$gene_id), ]
      n_genes <- nrow(members)
      total_genes <- pathway_sizes[[pid]]
      if (is.null(total_genes)) total_genes <- 0L
      coverage <- if (total_genes > 0) n_genes / total_genes else 0

      data.frame(
        pathway_id   = pid,
        category     = members$category[1],
        n_genes_bmd  = n_genes,
        total_genes  = total_genes,
        coverage_pct = round(coverage * 100, 1),
        median_bmd   = round(median(members$bmd), 6),
        mean_bmd     = round(mean(members$bmd), 6),
        q25_bmd      = round(quantile(members$bmd, 0.25), 6),
        min_bmd      = round(min(members$bmd), 6),
        gene_list    = paste(unique(members$gene_symbol), collapse = ";"),
        stringsAsFactors = FALSE
      )
    }))

    cat("Total pathways computed:", nrow(pathway_results), "\n")

    # ── Step 6: Apply NTP 2018 filter criteria ──────────────────────────────
    ntp_pass <- pathway_results[
      pathway_results$n_genes_bmd >= 3 &
      (pathway_results$coverage_pct >= 5 | pathway_results$total_genes == 0), ]

    # Remove very large generic GO terms (>2000 genes) — only for GO:BP
    ntp_pass <- ntp_pass[
      !(ntp_pass$category == "GO:BP" & ntp_pass$total_genes > 2000), ]
    # Remove very small pathways (< 10 total genes) — only for GO:BP and KEGG
    ntp_pass <- ntp_pass[
      !(ntp_pass$category %in% c("GO:BP", "KEGG") &
        ntp_pass$total_genes > 0 & ntp_pass$total_genes < 10), ]

    ntp_pass <- ntp_pass[order(ntp_pass$median_bmd), ]
    n_pathways_pass <- nrow(ntp_pass)

    cat("Pathways passing NTP criteria:", n_pathways_pass, "\n")

    if (n_pathways_pass > 0) {

      # ── Step 7: Get pathway names ────────────────────────────────────────
      go_names <- c()
      tryCatch({
        go_ids <- ntp_pass$pathway_id[ntp_pass$category == "GO:BP"]
        if (length(go_ids) > 0 && requireNamespace("GO.db", quietly = TRUE)) {
          terms <- AnnotationDbi::select(GO.db::GO.db, keys = go_ids,
                                          columns = "TERM", keytype = "GOID")
          go_names <- setNames(terms$TERM, terms$GOID)
          cat("  GO term names resolved:", length(go_names), "\n")
        }
      }, error = function(e) {
        cat("  GO name lookup failed:", e$message, "\n")
      })

      kegg_names <- c()
      tryCatch({
        kegg_ids <- ntp_pass$pathway_id[ntp_pass$category == "KEGG"]
        if (length(kegg_ids) > 0 && requireNamespace("KEGGREST", quietly = TRUE)) {
          kinfo <- KEGGREST::keggList("pathway", "hsa")
          name_map <- setNames(sub(" - Homo sapiens \\(human\\)", "", kinfo),
                               names(kinfo))
          matched <- name_map[paste0("path:", kegg_ids)]
          kegg_names <- setNames(ifelse(is.na(matched), kegg_ids, matched), kegg_ids)
        }
      }, error = function(e) {
        cat("  KEGG name lookup failed:", e$message, "\n")
      })

      all_names <- c(go_names, kegg_names)
      ntp_pass$pathway_name <- ntp_pass$pathway_id  # default: use ID
      for (i in seq_len(nrow(ntp_pass))) {
        nm <- all_names[ntp_pass$pathway_id[i]]
        if (!is.null(nm) && !is.na(nm) && nchar(nm) > 0) {
          ntp_pass$pathway_name[i] <- nm
        } else if (ntp_pass$category[i] %in% c("Hallmark", "C2:CP", "C2:CGP")) {
          # MSigDB names: clean up HALLMARK_P53_PATHWAY → P53 Pathway
          clean <- ntp_pass$pathway_id[i]
          clean <- sub("^HALLMARK_", "", clean)
          clean <- gsub("_", " ", clean)
          clean <- paste0(toupper(substr(clean, 1, 1)),
                          tolower(substr(clean, 2, nchar(clean))))
          ntp_pass$pathway_name[i] <- clean
        }
      }

      # ── Per-collection tPOD summary ──────────────────────────────────────
      collection_tpod <- list()
      for (cat_name in unique(ntp_pass$category)) {
        cat_rows <- ntp_pass[ntp_pass$category == cat_name, ]
        if (nrow(cat_rows) > 0) {
          top_row <- cat_rows[1, ]  # already sorted by median_bmd
          collection_tpod[[cat_name]] <- list(
            category     = cat_name,
            tpod_value   = top_row$median_bmd,
            pathway_name = top_row$pathway_name,
            pathway_id   = top_row$pathway_id,
            n_genes      = top_row$n_genes_bmd,
            coverage_pct = top_row$coverage_pct,
            n_pathways   = nrow(cat_rows)
          )
          cat("  ", cat_name, ": tPOD =", round(top_row$median_bmd, 4), "uM (",
              top_row$pathway_name, ",", nrow(cat_rows), "pathways passed)\n")
        }
      }

      # ── Step 8: tPOD = lowest pathway median BMD ─────────────────────────
      tpod_row    <- ntp_pass[1, ]
      tpod_value  <- tpod_row$median_bmd
      tpod_performed <- TRUE

      # Add pathway method to the distribution-based tpod_methods
      tpod_methods$pathway_ntp <- list(
        value = round(tpod_value, 6),
        label = "Lowest Pathway Median BMD (NTP 2018)",
        ref   = "NTP 2018; Thomas et al. 2013",
        n_genes_used = tpod_row$n_genes_bmd,
        pathway_name = tpod_row$pathway_name,
        pathway_id   = tpod_row$pathway_id,
        coverage_pct = tpod_row$coverage_pct
      )

      tpod_summary <- list(
        tpod_value       = tpod_value,
        tpod_pathway     = tpod_row$pathway_name,
        tpod_pathway_id  = tpod_row$pathway_id,
        tpod_category    = tpod_row$category,
        tpod_n_genes     = tpod_row$n_genes_bmd,
        tpod_coverage    = tpod_row$coverage_pct,
        tpod_method      = "NTP2018_gene_set",
        n_pathways_tested = n_pathways_pass,
        gene_level_median = bmd_summary$median_bmd,
        gene_level_q25    = round(quantile(valid_efsa$BMD.zSD, 0.25), 6),
        tpod_methods     = tpod_methods,
        collection_tpod  = collection_tpod,
        top5_pathways    = lapply(seq_len(min(5, nrow(ntp_pass))), function(i) {
          list(name = ntp_pass$pathway_name[i], id = ntp_pass$pathway_id[i],
               category = ntp_pass$category[i], median_bmd = ntp_pass$median_bmd[i],
               n_genes = ntp_pass$n_genes_bmd[i], coverage = ntp_pass$coverage_pct[i])
        })
      )

      cat("\n--- tPOD Result ---\n")
      cat("  tPOD:", tpod_value, "uM\n")
      cat("  Pathway:", tpod_row$pathway_name, "(", tpod_row$pathway_id, ")\n")
      cat("  Genes:", tpod_row$n_genes_bmd, "| Coverage:", tpod_row$coverage_pct, "%\n")

      # Save pathway results
      out_cols <- intersect(c("pathway_id", "pathway_name", "category", "n_genes_bmd",
                               "total_genes", "coverage_pct", "median_bmd", "mean_bmd",
                               "q25_bmd", "min_bmd", "gene_list"), colnames(ntp_pass))
      write.csv(ntp_pass[, out_cols],
                file.path(opt$outdir, "pathway_bmd_summary.csv"), row.names = FALSE)

      # Pathway sensitivity bar plot
      n_plot_pw <- min(20, nrow(ntp_pass))
      top_pw    <- ntp_pass[seq_len(n_plot_pw), ]
      top_pw$plot_label <- ifelse(nchar(top_pw$pathway_name) > 45,
        paste0(substr(top_pw$pathway_name, 1, 42), "..."), top_pw$pathway_name)
      top_pw$plot_label <- paste0(top_pw$plot_label, " (", top_pw$n_genes_bmd, ")")
      top_pw$is_tpod    <- seq_len(n_plot_pw) == 1

      p_pathway <- ggplot(top_pw,
                           aes(x = reorder(plot_label, -median_bmd),
                               y = median_bmd, fill = is_tpod)) +
        geom_col(alpha = 0.85, show.legend = FALSE) +
        scale_fill_manual(values = c("FALSE" = "#2196F3", "TRUE" = "#FF5722")) +
        coord_flip() +
        geom_hline(yintercept = tpod_value, linetype = "dashed",
                   color = "#FF5722", linewidth = 0.7) +
        labs(title = paste("Pathway Sensitivity \u2014", opt$treatment),
             subtitle = paste0("NTP 2018 | tPOD = ", round(tpod_value, 4),
                               " \u03bcM (", tpod_row$pathway_name, ")"),
             x = "", y = expression(paste("Median BMD (", mu, "M)"))) +
        theme_bw(base_size = 11) +
        theme(plot.title = element_text(size = 13, face = "bold"),
              plot.subtitle = element_text(size = 9, color = "grey40"),
              axis.text.y = element_text(size = 8))
      ggsave(file.path(opt$outdir, "pathway_sensitivity.png"), p_pathway,
             width = 12, height = max(6, 0.35 * n_plot_pw), dpi = 300)

      write_json(tpod_summary, file.path(opt$outdir, "tpod_summary.json"),
                 pretty = TRUE, auto_unbox = TRUE)
      cat("tPOD outputs saved\n")

    } else {
      cat("No pathways passed NTP criteria\n")
    }
  } else {
    cat("No gene-pathway mappings found\n")
  }
} else if (bmd_performed) {
  cat("Too few valid BMDs for pathway analysis\n")
}
write_json(list(
  treatment              = opt$treatment,
  control                = opt$control,
  total_genes_final      = nrow(final_counts),
  total_samples_final    = ncol(final_counts),
  dose_levels            = length(unique(dose_vec)),
  genes_selected         = n_selected,
  models_fitted          = n_fitted,
  batch_corrected        = FALSE,
  batch_corrected_pca    = batch_corrected_pca,
  pca_outliers_flagged   = length(pca_result$outliers),
  bmd_performed          = bmd_performed,
  bmd_summary            = if (bmd_performed) bmd_summary else NULL,
  model_quality          = if (exists("model_quality_summary")) model_quality_summary else NULL,
  tpod_performed         = tpod_performed,
  tpod_summary           = if (tpod_performed) tpod_summary else NULL,
  n_pathways_tested      = n_pathways_pass
), file.path(opt$outdir, "dromics_summary.json"), pretty = TRUE, auto_unbox = TRUE)

cat("\n=== DRomics Complete ===\n")
