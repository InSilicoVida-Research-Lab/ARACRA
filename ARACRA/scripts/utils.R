#!/usr/bin/env Rscript
# ============================================================
#  utils.R — Shared functions for DESeq2 and DRomics scripts
#  Source this from run_deseq2.R and run_dromics.R
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2); library(edgeR); library(sva)
  library(dplyr); library(readxl); library(stringr)
  library(jsonlite); library(org.Hs.eg.db)
})

# ── Load metadata (CSV or XLSX) ───────────────────────────────────────────────
load_metadata <- function(path) {
  cat("Reading metadata:", path, "\n")
  if (grepl("\\.xlsx?$", path, ignore.case = TRUE))
    as.data.frame(readxl::read_excel(path))
  else
    read.csv(path, stringsAsFactors = FALSE)
}

# ── Load count matrix (featureCounts, CSV, TSV) ───────────────────────────────
read_count_matrix <- function(path) {
  cat("Reading counts:", path, "\n")
  first_lines <- readLines(path, n = 3)

  is_fc_format <- any(grepl("^#", first_lines)) ||
    any(grepl("Geneid.*Chr.*Start.*End.*Strand.*Length", first_lines))

  if (is_fc_format) {
    skip <- sum(grepl("^#", first_lines))
    raw <- tryCatch(
      read.table(path, header = TRUE, sep = "\t", comment.char = "",
                 stringsAsFactors = FALSE, skip = skip, check.names = FALSE),
      error = function(e)
        read.table(path, header = TRUE, sep = "\t", comment.char = "",
                   stringsAsFactors = FALSE, skip = skip, fill = TRUE,
                   check.names = FALSE)
    )
    len_idx <- if ("Length" %in% colnames(raw)) which(colnames(raw) == "Length") else 6
    mat <- raw[, (len_idx + 1):ncol(raw), drop = FALSE]
    rownames(mat) <- if ("Geneid" %in% colnames(raw)) raw$Geneid else raw[, 1]
  } else {
    sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
    # Auto-detect separator if ambiguous
    if (!grepl("\\.(csv|tsv|out)$", path, ignore.case = TRUE)) {
      if (grepl("\t", first_lines[1])) sep <- "\t" else sep <- ","
    }
    raw <- read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE,
                      check.names = FALSE, row.names = 1)
    mat <- raw
  }

  # Clean column names: extract SRR/ERR/DRR IDs or use as-is
  clean <- str_extract(colnames(mat), "SRR[0-9]+|ERR[0-9]+|DRR[0-9]+")
  still_na <- is.na(clean)
  if (any(still_na)) {
    clean[still_na] <- ifelse(
      grepl("[/\\\\]", colnames(mat)[still_na]),
      gsub("\\.[^.]*$", "", basename(colnames(mat)[still_na])),
      colnames(mat)[still_na]
    )
  }
  colnames(mat) <- clean

  # Convert to integer matrix (critical for DESeq2 — salmon counts are floats)
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  mat <- round(mat)
  storage.mode(mat) <- "integer"

  # ── Parse gene IDs from rownames ──
  # Handle pipe-delimited salmon/GENCODE IDs:
  #   ENST...|ENSG...|OTTHUMG...|...|SYMBOL-NNN|SYMBOL|LENGTH|biotype|
  # Extract gene_symbol and ENSG from these, then aggregate to gene level
  raw_ids <- rownames(mat)

  if (any(grepl("\\|", raw_ids))) {
    cat("Detected pipe-delimited IDs (salmon/GENCODE format) — parsing...\n")
    # Parse: field 2 = ENSG, field 6 = gene symbol (in GENCODE format)
    parts <- strsplit(raw_ids, "\\|")
    n_fields <- sapply(parts, length)

    # Extract ENSG ID (field 2) and gene symbol (field 6)
    ensg_ids <- sapply(parts, function(p) if (length(p) >= 2) p[2] else p[1])
    gene_syms <- sapply(parts, function(p) if (length(p) >= 6) p[6] else "")

    # Strip version numbers from ENSG IDs
    ensg_ids <- sub("\\.\\d+$", "", ensg_ids)

    # Aggregate transcript counts to gene level using ENSG IDs
    cat("Aggregating", nrow(mat), "transcripts to gene level...\n")
    gene_counts <- list()
    for (i in seq_len(nrow(mat))) {
      gid <- ensg_ids[i]
      if (is.null(gene_counts[[gid]])) {
        gene_counts[[gid]] <- mat[i, , drop = FALSE]
      } else {
        gene_counts[[gid]] <- gene_counts[[gid]] + mat[i, , drop = FALSE]
      }
    }
    agg_mat <- do.call(rbind, gene_counts)
    rownames(agg_mat) <- names(gene_counts)
    storage.mode(agg_mat) <- "integer"

    # Build gene symbol lookup (ENSG -> symbol)
    sym_lookup <- setNames(gene_syms, ensg_ids)
    sym_lookup <- sym_lookup[!duplicated(names(sym_lookup))]
    # Remove empty/NA symbols
    sym_lookup <- sym_lookup[nchar(sym_lookup) > 0 & sym_lookup != ""]

    cat("Aggregated:", nrow(agg_mat), "genes from", nrow(mat), "transcripts\n")
    mat <- agg_mat
    attr(mat, "gene_symbol_lookup") <- sym_lookup
  } else if (any(grepl("^ENST", raw_ids))) {
    # Transcript-level ENST IDs without pipe delimiter — strip versions
    cat("Detected ENST IDs — stripping version numbers\n")
    rownames(mat) <- sub("\\.\\d+$", "", raw_ids)
  } else if (any(grepl("^ENSG", raw_ids))) {
    # Gene-level ENSG IDs — strip versions
    rownames(mat) <- sub("\\.\\d+$", "", raw_ids)
  }

  cat("Final matrix:", nrow(mat), "genes x", ncol(mat), "samples\n")
  mat
}

# ── Resolve gene symbols from IDs ─────────────────────────────────────────────
resolve_gene_symbols <- function(gene_ids, count_matrix = NULL) {
  # First try: use the lookup embedded in the count matrix (from pipe-delimited parsing)
  sym_lookup <- NULL
  if (!is.null(count_matrix)) {
    sym_lookup <- attr(count_matrix, "gene_symbol_lookup")
  }

  symbols <- rep(NA_character_, length(gene_ids))
  names(symbols) <- gene_ids

  # Apply embedded lookup first
  if (!is.null(sym_lookup)) {
    matched <- sym_lookup[gene_ids]
    valid <- !is.na(matched) & nchar(matched) > 0
    symbols[valid] <- matched[valid]
    cat("Gene symbols from embedded lookup:", sum(valid), "/", length(gene_ids), "\n")
  }

  # For remaining NAs, try org.Hs.eg.db
  still_na <- is.na(symbols)
  if (any(still_na)) {
    ensg_ids <- gene_ids[still_na]
    db_symbols <- tryCatch(
      mapIds(org.Hs.eg.db, keys = ensg_ids, column = "SYMBOL",
             keytype = "ENSEMBL", multiVals = "first"),
      error = function(e) {
        cat("org.Hs.eg.db lookup failed:", e$message, "\n")
        setNames(rep(NA_character_, length(ensg_ids)), ensg_ids)
      }
    )
    matched2 <- !is.na(db_symbols)
    symbols[still_na][matched2] <- db_symbols[matched2]
    cat("Gene symbols from org.Hs.eg.db:", sum(matched2), "/", sum(still_na), "\n")
  }

  # Final fallback: use the ID itself
  symbols[is.na(symbols)] <- gene_ids[is.na(symbols)]
  symbols
}

# ── Align samples between counts and metadata ─────────────────────────────────
align_samples <- function(counts_mat, meta_sub, exclude_samples_str = NULL) {
  common <- intersect(colnames(counts_mat), rownames(meta_sub))
  cat("Common samples between counts and metadata:", length(common), "\n")
  if (length(common) == 0) {
    cat("Count matrix columns (first 10):",
        paste(head(colnames(counts_mat), 10), collapse = ", "), "\n")
    cat("Metadata Sample_Names (first 10):",
        paste(head(rownames(meta_sub), 10), collapse = ", "), "\n")
    stop("No common samples found! Check that Sample_Name matches count matrix columns.")
  }

  # Apply user-confirmed exclusions
  if (!is.null(exclude_samples_str) && nchar(trimws(exclude_samples_str)) > 0) {
    user_excl <- trimws(strsplit(exclude_samples_str, ",")[[1]])
    user_excl <- user_excl[nchar(user_excl) > 0]
    if (length(user_excl) > 0) {
      cat("Excluding samples:", paste(user_excl, collapse = ", "), "\n")
      common <- setdiff(common, user_excl)
      cat("Remaining after exclusion:", length(common), "\n")
    }
  }

  # Ensure consistent ordering between counts and metadata
  counts_out <- counts_mat[, common, drop = FALSE]
  meta_out   <- meta_sub[common, , drop = FALSE]
  stopifnot("Sample order mismatch between counts and metadata" =
              identical(colnames(counts_out), rownames(meta_out)))

  list(counts = counts_out, meta = meta_out)
}

# ── Low-read and CPM filtering ────────────────────────────────────────────────
filter_samples_and_genes <- function(counts_sub, meta_sub, read_thresh = 1000000,
                                      group_col = "group", cpm_thresh = 1,
                                      min_prop = 0.75) {
  # Low-read filter
  lib_size <- colSums(counts_sub)
  keep_samp <- lib_size >= read_thresh
  cat("Samples passing read filter (>=", read_thresh, "):",
      sum(keep_samp), "/", length(keep_samp), "\n")
  counts_sub <- counts_sub[, keep_samp]
  meta_sub <- meta_sub[keep_samp, ]

  # CPM gene filter per group
  cpm_mat <- edgeR::cpm(counts_sub)
  keep_genes <- rep(FALSE, nrow(cpm_mat))
  grp_col <- if (group_col %in% colnames(meta_sub)) group_col else "condition"
  if (grp_col %in% colnames(meta_sub)) {
    meta_sub[[grp_col]] <- as.factor(meta_sub[[grp_col]])
    for (grp in levels(meta_sub[[grp_col]])) {
      samps <- intersect(
        rownames(meta_sub[meta_sub[[grp_col]] == grp, ]),
        colnames(cpm_mat)
      )
      if (!length(samps)) next
      prop <- rowSums(cpm_mat[, samps, drop = FALSE] >= cpm_thresh) / length(samps)
      keep_genes[prop >= min_prop] <- TRUE
    }
  } else {
    # Fallback: filter across all samples
    prop <- rowSums(cpm_mat >= cpm_thresh) / ncol(cpm_mat)
    keep_genes[prop >= min_prop] <- TRUE
  }
  cat("Genes after CPM filter:", sum(keep_genes), "/", nrow(counts_sub), "\n")
  counts_sub <- counts_sub[keep_genes, ]

  list(counts = counts_sub, meta = meta_sub)
}

# ── Batch correction ──────────────────────────────────────────────────────────
apply_batch_correction <- function(counts, metadata, group_col = "Type") {
  if (!"Batch" %in% colnames(metadata)) {
    cat("No Batch column — skipping batch correction\n")
    return(list(counts = counts, corrected = FALSE))
  }
  batch <- as.factor(metadata$Batch)
  if (nlevels(batch) < 2) {
    cat("Only 1 batch — skipping\n")
    return(list(counts = counts, corrected = FALSE))
  }

  grp <- as.factor(metadata[[group_col]])
  tryCatch({
    if (exists("ComBat_seq", mode = "function")) {
      corrected <- ComBat_seq(counts = counts, batch = batch,
                               group = grp, shrink = FALSE)
      # ComBat_seq can return non-integers — round them
      corrected <- round(corrected)
      storage.mode(corrected) <- "integer"
      cat("ComBat_seq batch correction applied\n")
    } else {
      cpm_log <- edgeR::cpm(counts, log = TRUE, prior.count = 1)
      corr_log <- sva::ComBat(dat = cpm_log, batch = batch,
                               mod = model.matrix(~grp), par.prior = TRUE)
      corrected <- round(pmax(2^corr_log - 1, 0))
      storage.mode(corrected) <- "integer"
      cat("sva::ComBat batch correction applied\n")
    }
    list(counts = corrected, corrected = TRUE)
  }, error = function(e) {
    cat("Batch correction failed:", e$message, "— using original counts\n")
    list(counts = counts, corrected = FALSE)
  })
}

# ── PCA outlier detection (writes JSON + plot for app) ────────────────────────
detect_pca_outliers <- function(counts, meta, outdir,
                                 group_col = "condition",
                                 sd_threshold = 2.0) {
  suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })

  dds_qc <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~1)
  vst_qc <- vst(dds_qc, blind = TRUE)
  pca_q <- prcomp(t(assay(vst_qc)))
  pca_data <- as.data.frame(pca_q$x[, 1:min(2, ncol(pca_q$x))])
  pct_var <- round(100 * pca_q$sdev^2 / sum(pca_q$sdev^2))

  pca_outliers <- c()
  if (group_col %in% colnames(meta)) {
    for (grp in unique(meta[[group_col]])) {
      samps <- rownames(meta[meta[[group_col]] == grp, , drop = FALSE])
      if (length(samps) < 3) next
      gdat <- pca_data[samps, , drop = FALSE]
      dists <- sqrt(rowSums(sweep(gdat, 2, colMeans(gdat), "-")^2))
      thr <- mean(dists) + sd_threshold * sd(dists)
      pca_outliers <- c(pca_outliers, names(dists[dists > thr]))
    }
  }

  # Build per-sample report
  pca_report <- lapply(rownames(pca_data), function(s) {
    grp_val <- if (group_col %in% colnames(meta)) meta[s, group_col] else ""
    list(
      sample    = s,
      pc1       = round(pca_data[s, 1], 3),
      pc2       = round(if (ncol(pca_data) >= 2) pca_data[s, 2] else 0, 3),
      condition = as.character(grp_val),
      flagged   = s %in% pca_outliers
    )
  })

  result <- list(
    layer   = "expression",
    status  = if (length(pca_outliers) > 0) "outliers_found" else "all_pass",
    flagged = as.list(pca_outliers),
    samples = pca_report
  )

  write_json(result, file.path(outdir, "pca_outlier_flag.json"), pretty = TRUE,
             auto_unbox = TRUE)

  # ── Generate PCA scatter plot for the QC review UI ──────────────────────────
  plot_df <- pca_data
  plot_df$sample <- rownames(pca_data)
  plot_df$group  <- if (group_col %in% colnames(meta)) meta[rownames(pca_data), group_col] else "All"
  plot_df$status <- ifelse(plot_df$sample %in% pca_outliers, "Outlier", "Pass")

  has_batch_col <- "Batch" %in% colnames(meta)
  if (has_batch_col) {
    plot_df$batch <- as.factor(meta[rownames(pca_data), "Batch"])
  }

  p_qc <- ggplot(plot_df, aes(x = PC1, y = PC2, color = group, label = sample)) +
    geom_point(aes(shape = status), size = 4, alpha = 0.8) +
    scale_shape_manual(values = c("Pass" = 16, "Outlier" = 4),
                       name = "QC Status") +
    geom_text_repel(size = 2.8, max.overlaps = 20, show.legend = FALSE) +
    labs(
      title    = "PCA \u2014 Expression QC (Layer 2)",
      subtitle = paste0(
        length(pca_outliers), " outlier(s) flagged | ",
        nrow(pca_data), " samples | ",
        if (has_batch_col) "batch-corrected for visualisation" else "no batch correction"
      ),
      x = paste0("PC1 (", pct_var[1], "%)"),
      y = paste0("PC2 (", pct_var[2], "%)")
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      legend.position = "bottom"
    )

  # Highlight outliers with red ring
  if (length(pca_outliers) > 0) {
    outlier_df <- plot_df[plot_df$status == "Outlier", , drop = FALSE]
    p_qc <- p_qc +
      geom_point(data = outlier_df, aes(x = PC1, y = PC2),
                 shape = 1, size = 7, color = "red", stroke = 1.5,
                 inherit.aes = FALSE)
  }

  tryCatch({
    ggsave(file.path(outdir, "pca_qc_plot.png"), p_qc,
           width = 10, height = 7, dpi = 200)
    cat("PCA QC plot saved:", file.path(outdir, "pca_qc_plot.png"), "\n")
  }, error = function(e) {
    cat("Warning: Could not save PCA QC plot:", e$message, "\n")
  })

  cat("PCA outlier check: flagged", length(pca_outliers), "sample(s)\n")
  if (length(pca_outliers) > 0) {
    cat("  Flagged:", paste(pca_outliers, collapse = ", "), "\n")
  }

  list(
    outliers = pca_outliers,
    pca_data = pca_data,
    vst_obj  = vst_qc
  )
}
