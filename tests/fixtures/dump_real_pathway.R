#!/usr/bin/env Rscript
# ============================================================
#  dump_real_pathway.R — pick a real, well-sized GO:BP term and a
#  deterministic background gene sample, both as real ENSEMBL IDs.
#
#  Why: 04_fixture.sh's main fixture uses fake gene IDs (ENSG_RESP_0000 etc)
#  that never match anything in org.Hs.eg.db, so run_dromics.R's
#  pathway-level tPOD code path (GO/KEGG/MSigDB gene-set mapping) has never
#  actually been exercised by the test suite. This script supplies real IDs
#  so 05_pathway_fixture.sh can build a fixture that path can actually run on.
#
#  Also times the ENSEMBL -> GO mapping over all ~41k human genes — the one
#  part of run_dromics.R's pathway step ("Step 1") that had never been
#  timed, since it always silently no-ops on the fake-ID fixture.
# ============================================================
suppressPackageStartupMessages({
  library(optparse)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GO.db)
})

option_list <- list(
  make_option("--outdir",   type = "character", default = "."),
  make_option("--min_size", type = "integer",    default = 15),
  make_option("--max_size", type = "integer",    default = 25)
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

t0 <- Sys.time()
cat("Fetching all human ENSEMBL IDs...\n")
all_ensembl <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
cat("  n =", length(all_ensembl), "\n")

t1 <- Sys.time()
cat("Mapping ENSEMBL -> GO (same query as run_dromics.R's pathway Step 1)...\n")
raw <- AnnotationDbi::select(org.Hs.eg.db, keys = all_ensembl,
                              columns = c("ENSEMBL", "GO", "ONTOLOGY"),
                              keytype = "ENSEMBL")
cat("  n rows =", nrow(raw), " (",
    round(as.numeric(Sys.time() - t1, units = "secs"), 1), "s )\n")

bp <- unique(raw[!is.na(raw$GO) & !is.na(raw$ONTOLOGY) & raw$ONTOLOGY == "BP",
                  c("ENSEMBL", "GO")])
cat("GO:BP rows:", nrow(bp), "(", length(unique(bp$GO)), "terms )\n")

sizes <- table(bp$GO)
candidates <- names(sizes)[sizes >= opt$min_size & sizes <= opt$max_size]
cat("Candidate GO:BP terms sized", opt$min_size, "-", opt$max_size, ":",
    length(candidates), "\n")
if (length(candidates) == 0) stop("No GO:BP term found in the requested size range")

# Deterministic pick: sorted GO ID, first match. Not random, so re-running
# this script (even on a different machine/Bioconductor version, modulo
# annotation drift) picks the same term as long as it's still in range.
chosen_go <- sort(candidates)[1]
chosen_genes <- sort(unique(bp$ENSEMBL[bp$GO == chosen_go]))
chosen_name <- tryCatch(
  AnnotationDbi::select(GO.db::GO.db, keys = chosen_go,
                         columns = "TERM", keytype = "GOID")$TERM,
  error = function(e) chosen_go
)
cat("\nChosen pathway:", chosen_go, "-", chosen_name,
    "(", length(chosen_genes), "genes )\n")

# 400 other real, GO:BP-annotated genes, excluding the chosen pathway's
# members, picked by sorted ENSEMBL ID (deterministic, no RNG dependency).
bg_pool <- sort(setdiff(unique(bp$ENSEMBL), chosen_genes))
bg_genes <- bg_pool[round(seq(1, length(bg_pool), length.out = 400))]

write.csv(data.frame(ensembl_id = chosen_genes),
          file.path(opt$outdir, "real_pathway_genes.csv"), row.names = FALSE)
write.csv(data.frame(ensembl_id = bg_genes),
          file.path(opt$outdir, "real_background_genes.csv"), row.names = FALSE)
writeLines(c(chosen_go, chosen_name, as.character(length(chosen_genes))),
           file.path(opt$outdir, "real_pathway_id.txt"))

cat("\nTotal time:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")
cat("Wrote real_pathway_genes.csv, real_background_genes.csv, real_pathway_id.txt to",
    opt$outdir, "\n")
