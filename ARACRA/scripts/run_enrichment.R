#!/usr/bin/env Rscript
# ============================================================
#  run_enrichment.R — Pathway enrichment / over-representation
#  analysis (ORA) on significant DESeq2 DEGs.
#
#  Runs AFTER run_deseq2.R. Gene list = significant DEGs
#  (Custom_Filtered_DEGs.csv); background universe = every gene
#  DESeq2 actually tested (All_Results.csv) — not the whole genome,
#  so the hypergeometric test isn't biased by genes this experiment
#  never had power to detect.
#
#  Databases, each queried through its own peer-reviewed R/Bioconductor
#  interface rather than an ad hoc download:
#    GO BP/MF/CC — clusterProfiler::enrichGO
#      Ashburner et al. 2000, Nat Genet 25:25-29 (Gene Ontology)
#      Wu et al. 2021, Innovation 2:100141 (clusterProfiler v4)
#    KEGG        — clusterProfiler::enrichKEGG
#      Kanehisa & Goto 2000, Nucleic Acids Res 28:27-30
#    Reactome    — ReactomePA::enrichPathway
#      Fabregat et al. 2018, Nucleic Acids Res 46:D649-D655
#      Yu & He 2016, Mol Biosyst 12:477-479
#    MSigDB H/C2 — clusterProfiler::enricher + msigdbr gene sets
#      Liberzon et al. 2015, Cell Syst 1:417-425
#
#  Every database is optional at the package level: if its R package
#  isn't installed, that database is skipped (not a hard failure) and
#  the reason is recorded in enrichment_summary.json for the GUI to
#  surface — same "degrade, don't crash" policy as the rest of the
#  pipeline's optional modules.
# ============================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--degs",          type = "character", help = "Custom_Filtered_DEGs.csv from run_deseq2.R"),
  make_option("--all_results",   type = "character", help = "All_Results.csv from run_deseq2.R (background universe)"),
  make_option("--outdir",        type = "character", default = "."),
  make_option("--databases",     type = "character", default = "GO_BP,KEGG,REACTOME,MSIGDB_H",
              help = "Comma-separated subset of GO_BP,GO_MF,GO_CC,KEGG,REACTOME,MSIGDB_H,MSIGDB_C2"),
  make_option("--padj_cutoff",   type = "double",    default = 0.05),
  make_option("--qvalue_cutoff", type = "double",    default = 0.2),
  make_option("--min_geneset",   type = "integer",   default = 10),
  make_option("--max_geneset",   type = "integer",   default = 500)
)

opt <- parse_args(OptionParser(option_list = option_list))

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
utils_path <- file.path(script_dir, "utils.R")
if (!file.exists(utils_path)) utils_path <- file.path(getwd(), "utils.R")
if (!file.exists(utils_path)) stop("Cannot find utils.R")
source(utils_path)   # gives us jsonlite::write_json + org.Hs.eg.db already loaded

suppressPackageStartupMessages(library(AnnotationDbi))

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

cat("=== Pathway Enrichment (ORA) ===\n")

databases <- toupper(trimws(strsplit(opt$databases, ",")[[1]]))
databases <- databases[databases != ""]
cat("Requested databases:", paste(databases, collapse = ", "), "\n\n")

results_summary <- list()

finish <- function(status, reason = NULL, extra = list()) {
  payload <- list(status = status, databases_requested = databases, results = results_summary)
  if (!is.null(reason)) payload$reason <- reason
  payload <- c(payload, extra)
  write_json(payload, file.path(opt$outdir, "enrichment_summary.json"),
             pretty = TRUE, auto_unbox = TRUE)
  cat("\n=== Pathway Enrichment:", status, "===\n")
  quit(status = 0)
}

# ── Load DEGs & background ──────────────────────────────────────────────────
if (!file.exists(opt$degs) || file.info(opt$degs)$size == 0) {
  cat("No DEG file found — nothing to enrich.\n")
  finish("skipped", "no_degs")
}

degs_df <- read.csv(opt$degs, stringsAsFactors = FALSE)
if (nrow(degs_df) == 0) {
  cat("Zero significant DEGs at the configured threshold — nothing to enrich.\n")
  finish("skipped", "zero_degs")
}
bg_df <- read.csv(opt$all_results, stringsAsFactors = FALSE)

deg_ensembl <- unique(degs_df$ensembl_id)
bg_ensembl  <- unique(bg_df$ensembl_id)
cat("DEGs:", length(deg_ensembl), "| Background universe:", length(bg_ensembl), "\n")

# ENSEMBL -> ENTREZID (KEGG and Reactome index by Entrez, not Ensembl)
deg_entrez <- tryCatch(
  as.character(na.omit(unique(mapIds(org.Hs.eg.db, keys = deg_ensembl, column = "ENTREZID",
                                      keytype = "ENSEMBL", multiVals = "first")))),
  error = function(e) character(0))
bg_entrez <- tryCatch(
  as.character(na.omit(unique(mapIds(org.Hs.eg.db, keys = bg_ensembl, column = "ENTREZID",
                                      keytype = "ENSEMBL", multiVals = "first")))),
  error = function(e) character(0))
cat("Entrez-mapped: ", length(deg_entrez), "DEGs /", length(bg_entrez), "background\n\n")

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  cat("clusterProfiler not installed — skipping all ORA.\n")
  cat("  Install with: BiocManager::install(\"clusterProfiler\")\n")
  finish("skipped", "clusterProfiler_missing")
}
suppressPackageStartupMessages(library(clusterProfiler))

# ── Save one database's result: CSV + dotplot + barplot ────────────────────
save_result <- function(ego, tag, label) {
  df <- tryCatch(as.data.frame(ego), error = function(e) NULL)
  if (is.null(ego) || is.null(df) || nrow(df) == 0) {
    cat(" ", label, ": 0 significant terms\n")
    results_summary[[tag]] <<- list(ran = TRUE, n_terms = 0, label = label)
    return(invisible(NULL))
  }
  write.csv(df, file.path(opt$outdir, paste0("enrichment_", tag, ".csv")), row.names = FALSE)
  cat(" ", label, ":", nrow(df), "significant terms\n")

  n_show <- min(20, nrow(df))
  # GO/Reactome/MSigDB term names run long ("maturation of SSU-rRNA from
  # tricistronic rRNA transcript..."). At the old 0.3in/category, 20 terms in
  # a 6in-tall plot left no room for that text, so labels visually merged
  # into each other. Wrap onto multiple lines instead of letting them run
  # into their neighbors, and give each category enough height to hold a
  # wrapped 2-3 line label without crowding the one above/below it.
  .wrap_labels <- ggplot2::scale_y_discrete(labels = scales::label_wrap(40))
  plot_height <- max(6, n_show * 0.55)
  tryCatch({
    p1 <- enrichplot::dotplot(ego, showCategory = n_show) + .wrap_labels + ggplot2::ggtitle(label)
    ggplot2::ggsave(file.path(opt$outdir, paste0("enrichment_", tag, "_dotplot.png")),
                     p1, width = 11, height = plot_height, dpi = 300, limitsize = FALSE)
  }, error = function(e) cat("    (dotplot skipped:", e$message, ")\n"))
  tryCatch({
    p2 <- graphics::barplot(ego, showCategory = n_show) + .wrap_labels + ggplot2::ggtitle(label)
    ggplot2::ggsave(file.path(opt$outdir, paste0("enrichment_", tag, "_barplot.png")),
                     p2, width = 11, height = plot_height, dpi = 300, limitsize = FALSE)
  }, error = function(e) cat("    (barplot skipped:", e$message, ")\n"))

  results_summary[[tag]] <<- list(ran = TRUE, n_terms = nrow(df), label = label,
                                   top_term = df$Description[1])
}

# ── GO: Biological Process / Molecular Function / Cellular Component ───────
for (ont in c("BP", "MF", "CC")) {
  tag <- paste0("GO_", ont)
  if (!(tag %in% databases)) next
  cat("Running GO", ont, "(Ashburner et al. 2000)...\n")
  ego <- tryCatch(
    enrichGO(gene = deg_ensembl, universe = bg_ensembl, OrgDb = org.Hs.eg.db,
             keyType = "ENSEMBL", ont = ont, pAdjustMethod = "BH",
             pvalueCutoff = opt$padj_cutoff, qvalueCutoff = opt$qvalue_cutoff,
             minGSSize = opt$min_geneset, maxGSSize = opt$max_geneset,
             readable = TRUE),
    error = function(e) { cat("  GO", ont, "FAILED:", e$message, "\n"); NULL })
  save_result(ego, tag, paste("GO", ont))
}

# ── KEGG ─────────────────────────────────────────────────────────────────────
if ("KEGG" %in% databases) {
  if (length(deg_entrez) == 0) {
    cat("KEGG: no Entrez-mapped DEGs — skipping\n")
    results_summary[["KEGG"]] <- list(ran = FALSE, reason = "no_entrez_ids")
  } else {
    cat("Running KEGG (Kanehisa & Goto 2000) — needs internet access to the KEGG REST API...\n")
    ekegg <- tryCatch(
      enrichKEGG(gene = deg_entrez, universe = bg_entrez, organism = "hsa",
                 pAdjustMethod = "BH", pvalueCutoff = opt$padj_cutoff,
                 qvalueCutoff = opt$qvalue_cutoff,
                 minGSSize = opt$min_geneset, maxGSSize = opt$max_geneset),
      error = function(e) { cat("  KEGG FAILED (likely no internet access):", e$message, "\n"); NULL })
    if (!is.null(ekegg)) {
      ekegg <- tryCatch(setReadable(ekegg, org.Hs.eg.db, keyType = "ENTREZID"),
                         error = function(e) ekegg)
    }
    save_result(ekegg, "KEGG", "KEGG")
  }
}

# ── Reactome ─────────────────────────────────────────────────────────────────
if ("REACTOME" %in% databases) {
  if (!requireNamespace("ReactomePA", quietly = TRUE)) {
    cat("ReactomePA not installed — skipping Reactome.\n")
    cat("  Install with: BiocManager::install(\"ReactomePA\")\n")
    results_summary[["REACTOME"]] <- list(ran = FALSE, reason = "package_missing")
  } else if (length(deg_entrez) == 0) {
    cat("Reactome: no Entrez-mapped DEGs — skipping\n")
    results_summary[["REACTOME"]] <- list(ran = FALSE, reason = "no_entrez_ids")
  } else {
    suppressPackageStartupMessages(library(ReactomePA))
    cat("Running Reactome (Fabregat et al. 2018; Yu & He 2016)...\n")
    ereact <- tryCatch(
      enrichPathway(gene = deg_entrez, universe = bg_entrez, organism = "human",
                     pAdjustMethod = "BH", pvalueCutoff = opt$padj_cutoff,
                     qvalueCutoff = opt$qvalue_cutoff,
                     minGSSize = opt$min_geneset, maxGSSize = opt$max_geneset,
                     readable = TRUE),
      error = function(e) { cat("  Reactome FAILED:", e$message, "\n"); NULL })
    save_result(ereact, "REACTOME", "Reactome")
  }
}

# ── MSigDB (Hallmark / C2) via msigdbr ────────────────────────────────────────
run_msigdb <- function(category, tag, label) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    cat(label, ": msigdbr not installed — skipping.\n")
    cat("  Install with: install.packages('msigdbr')\n")
    results_summary[[tag]] <<- list(ran = FALSE, reason = "package_missing")
    return(invisible(NULL))
  }
  cat("Running", label, "(Liberzon et al. 2015)...\n")
  sets <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", category = category),
                    error = function(e) { cat("  ", label, "gene sets FAILED to load:", e$message, "\n"); NULL })
  if (is.null(sets)) {
    results_summary[[tag]] <<- list(ran = FALSE, reason = "load_failed")
    return(invisible(NULL))
  }
  t2g <- unique(as.data.frame(sets)[, c("gs_name", "ensembl_gene")])
  em <- tryCatch(
    enricher(gene = deg_ensembl, universe = bg_ensembl, TERM2GENE = t2g,
             pAdjustMethod = "BH", pvalueCutoff = opt$padj_cutoff,
             qvalueCutoff = opt$qvalue_cutoff,
             minGSSize = opt$min_geneset, maxGSSize = opt$max_geneset),
    error = function(e) { cat("  ", label, "FAILED:", e$message, "\n"); NULL })
  save_result(em, tag, label)
}
if ("MSIGDB_H"  %in% databases) run_msigdb("H",  "MSIGDB_H",  "MSigDB Hallmark")
if ("MSIGDB_C2" %in% databases) run_msigdb("C2", "MSIGDB_C2", "MSigDB C2 (curated)")

finish("completed", extra = list(
  n_degs        = length(deg_ensembl),
  n_degs_entrez = length(deg_entrez),
  n_background  = length(bg_ensembl),
  padj_cutoff   = opt$padj_cutoff,
  qvalue_cutoff = opt$qvalue_cutoff
))
