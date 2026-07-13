#!/usr/bin/env Rscript
# =============================================================================
#  run_tcplfit2_epa.R — EPA gene-level HTTr tPOD workflow
#
#  Implements the gene-level (DESeq2 + tcplfit2) branch of:
#    Harrill JA, Everett LJ, Haggard DE, Bundy JL, Willis CM, Shah I,
#    Paul Friedman K, Basili D, Middleton A, Judson RS (2024).
#    "Exploring the effects of experimental parameters and data modeling
#     approaches on in vitro transcriptomic point-of-departure estimates."
#    Toxicology 501:153694.  §2.5 (data processing) and §2.10.1 (gene_05,
#    gene_abs5, gene_min).
#
#  THIS IS NOT the DRomics-matched script. Preprocessing deliberately departs
#  from utils.R because the EPA workflow specifies a different one. Deltas:
#
#    | Step              | EPA (here)                      | ARACRA/DRomics        |
#    |-------------------|---------------------------------|-----------------------|
#    | Modelling unit    | PROBE-level counts              | gene-level (summed)   |
#    | Probe -> gene     | max |L2FC| in either direction  | sum of probe counts   |
#    | Feature filter    | mean count < 5 removed          | CPM>=1 in >=75% group |
#    | Batch/plate       | covariate in DESeq2 design      | ComBat-seq            |
#    | Response fitted   | shrunken L2FC per dose group    | VST/rlog per sample   |
#    | Gene pre-select   | NONE (hitcall decides)          | itemselect trend test |
#    | Criterion         | AIC                             | AICc                  |
#    | Noise band        | cross-chemical, 2 lowest concs  | per-model residual SD |
#
#  PIPELINE (whole panel, in one run — the cutoff needs all chemicals)
#    1. per chemical: DESeq2 on probe counts, jointly with plate-matched
#       vehicle controls, design ~ plate + dose_group
#       - probes with mean count < 5 dropped within each modelling subset
#       - Wald p-values BEFORE shrinkage; BH within each concentration;
#         independentFiltering = FALSE
#       - moderated L2FC with NORMAL shrinkage, per dose-group-vs-vehicle
#    2. probe -> gene: highest-magnitude L2FC of any probe, sign preserved
#    3. gene retention: L2FC present in >=95% of treatments; missing -> 0
#    4. panel noise band, per gene: pool L2FC at the TWO LOWEST concentrations
#       across ALL chemicals; cutoff = bounds of the 95% interval of that
#       distribution; onesd = SD of it.  [see CAVEAT below]
#    5. tcplfit2 per gene: bmed = 0, BMR = 1.349 x onesd, AIC, no pre-selection,
#       force.fit, bmd_low_bnd = 0.1 (BMD >= lowest conc / 10)
#    6. active if hitcall > 0.9  ->  tPODs: gene_05 / gene_abs5 / gene_min
#
#  CAVEAT ON THE CUTOFF (read this before you write the methods section)
#    Harrill 2024 spells the cutoff rule out only for SIGNATURES (§2.7): take
#    the scores at the two lowest concentrations across all 44 chemicals and
#    set the cutoff to the bounds of the 95% CI of that distribution. §2.10.1
#    says only "fit using tcplfit2, starting with L2FC data" — it does not
#    state the gene-level cutoff derivation. Step 4 above applies the signature
#    rule per gene, which is the natural reading but IS AN INFERENCE. Two
#    alternatives are exposed via --cutoff_method:
#      panel_ci  (default) : 95% interval of the cross-chemical null, per gene
#      panel_bmad          : 3 x MAD of the same null, per gene (ToxCast norm)
#      fixed               : --cutoff_value on the log2 scale
#    Report whichever you use, and say it is an inference if it is panel_ci.
#
#  USAGE
#    python3 build_probe_matrix.py --counts_dir idxstats/ --outfile probe_counts.csv
#    Rscript run_tcplfit2_epa.R \
#      --counts probe_counts.csv \
#      --probe_map probe_to_gene.tsv \
#      --metadata metadata.xlsx \
#      --control DMSO \
#      --outdir results/epa_tpod
# =============================================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--counts",        type = "character",
              help = "PROBE-level count matrix (from build_probe_matrix.py)"),
  make_option("--probe_map",     type = "character", default = NULL,
              help = "probe_to_gene.tsv. If omitted, rows are treated as genes already (NOT EPA-faithful — warns)"),
  make_option("--metadata",      type = "character",
              help = "Metadata: Sample_Name, Treatment, Dose, Type[, Plate/Batch]"),
  make_option("--control",       type = "character", default = "DMSO",
              help = "DEFAULT vehicle. Used only for chemicals with no vehicle assignment."),
  make_option("--vehicle_col",   type = "character", default = "Vehicle",
              help = "Metadata column naming each sample's vehicle (e.g. DMSO / MeOH)"),
  make_option("--vehicle_map",   type = "character", default = NULL,
              help = "CSV Treatment,Vehicle — use when metadata has no vehicle column"),
  make_option("--noise_by_vehicle", type = "logical", default = FALSE,
              help = "Compute the panel noise band separately per vehicle (needs enough chemicals per vehicle)"),
  make_option("--treatments",    type = "character", default = "all",
              help = "Comma-separated chemicals, or 'all' (panel cutoff needs the full panel)"),
  make_option("--outdir",        type = "character", default = "."),

  # ---- preprocessing (EPA §2.5) ----
  make_option("--min_mean_count", type = "double",  default = 5,
              help = "Drop probes with mean count < this within each modelling subset"),
  make_option("--plate_col",      type = "character", default = "Plate",
              help = "Plate column for the DESeq2 covariate; falls back to Batch, then dropped"),
  make_option("--gene_presence",  type = "double",  default = 0.95,
              help = "Gene must have an L2FC in >= this fraction of treatments"),
  make_option("--exclude_doses",  type = "character", default = NULL,
              help = "CSV with Treatment,Dose — cytotoxic concentrations to drop (EPA drops >50%% PI+/casp+)"),
  make_option("--shrinkage",      type = "character", default = "normal",
              help = "normal (EPA) | apeglm | ashr | none"),

  # ---- noise band ----
  make_option("--cutoff_method",  type = "character", default = "panel_ci",
              help = "panel_ci | panel_bmad | fixed"),
  make_option("--n_low_conc",     type = "integer",  default = 2L,
              help = "Number of lowest concentrations forming the null (EPA: 2)"),
  make_option("--cutoff_value",   type = "double",   default = 0.585,
              help = "Only for --cutoff_method fixed (0.585 = 1.5-fold)"),
  make_option("--cutoff_floor",   type = "double",   default = 0.05,
              help = "Lower bound on the per-gene cutoff, guards against zero-variance genes"),

  # ---- tcplfit2 (EPA §2.10.1) ----
  make_option("--fitmodels",      type = "character",
              default = "cnst,hill,poly1,poly2,pow,exp2,exp3,exp4,exp5",
              help = "Harrill 2024 lists cnst + Hill + 2 poly + pow + 4 exp. Add gnls for Harrill 2021."),
  make_option("--aicc",           type = "logical",  default = FALSE,
              help = "EPA uses AIC. TRUE switches to AICc (DRomics convention)."),
  make_option("--bmr_scale",      type = "double",   default = 1.349,
              help = "BMR = bmr_scale x onesd. 1.349 is fixed convention (Thomas 2007) — do not tune."),
  make_option("--bmd_low_bnd",    type = "double",   default = 0.1,
              help = "BMD >= bmd_low_bnd x lowest tested conc (EPA: 1 order of magnitude)"),
  make_option("--bmd_up_bnd",     type = "double",   default = 1,
              help = "BMD <= bmd_up_bnd x highest tested conc. NA to disable."),
  make_option("--hitcall_thresh", type = "double",   default = 0.9,
              help = "Active if hitcall > this (EPA: 0.9)"),
  make_option("--tc_min",         type = "double",   default = 0,
              help = "Optional top/cutoff floor (0 = off; signature-level EPA used TC>=1)"),
  make_option("--force_fit",      type = "logical",  default = TRUE),
  make_option("--bidirectional",  type = "logical",  default = TRUE),

  # ---- extras ----
  make_option("--also_dromics_tpods", type = "logical", default = TRUE,
              help = "Additionally report the 5 distribution tPODs used in run_dromics.R (clearly labelled non-EPA)"),
  make_option("--ncores",         type = "integer",  default = 0L)
)

opt <- parse_args(OptionParser(option_list = option_list))

suppressPackageStartupMessages({
  library(DESeq2); library(dplyr); library(readxl); library(stringr)
  library(jsonlite); library(ggplot2); library(parallel); library(tcplfit2)
})

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)
n_cores <- if (opt$ncores > 0) opt$ncores else max(1, parallel::detectCores() - 1)

cat("=== EPA gene-level HTTr tPOD workflow ===\n")
cat("Reference: Harrill et al. 2024, Toxicology 501:153694 (§2.5, §2.10.1)\n\n")

# =============================================================================
#  LOAD
# =============================================================================
read_matrix <- function(path) {
  sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  raw <- read.table(path, header = TRUE, sep = sep, row.names = 1,
                    check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(raw); mode(m) <- "numeric"
  m[is.na(m)] <- 0
  storage.mode(m) <- "integer"
  # tolerate SRR-style column decoration
  cl <- str_extract(colnames(m), "SRR[0-9]+|ERR[0-9]+|DRR[0-9]+")
  colnames(m) <- ifelse(is.na(cl), colnames(m), cl)
  m
}

counts <- read_matrix(opt$counts)
cat("Counts:", nrow(counts), "features x", ncol(counts), "samples\n")

meta <- if (grepl("\\.xlsx?$", opt$metadata, ignore.case = TRUE))
  as.data.frame(readxl::read_excel(opt$metadata)) else
  read.csv(opt$metadata, stringsAsFactors = FALSE)
rownames(meta) <- meta$Sample_Name

common <- intersect(colnames(counts), rownames(meta))
if (length(common) == 0) stop("No overlap between count columns and Sample_Name")
counts <- counts[, common, drop = FALSE]
meta   <- meta[common, , drop = FALSE]
cat("Matched samples:", length(common), "\n")

meta$dose_numeric <- ifelse(meta$Type == "control", 0,
                            suppressWarnings(as.numeric(as.character(meta$Dose))))
if (any(is.na(meta$dose_numeric))) stop("Non-numeric Dose in treated samples")

# Plate covariate
plate_col <- NULL
for (cc in c(opt$plate_col, "Plate", "Batch")) {
  if (!is.null(cc) && cc %in% colnames(meta) &&
      nlevels(as.factor(meta[[cc]])) > 1) { plate_col <- cc; break }
}
if (is.null(plate_col)) {
  cat("No usable plate/batch column — design reduces to ~ dose_group\n")
  cat("  (EPA models counts ~ plate + dose group; note the deviation.)\n")
} else {
  cat("Plate covariate:", plate_col, "(",
      nlevels(as.factor(meta[[plate_col]])), "levels )\n")
}

# Probe -> gene map
probe_map <- NULL
if (!is.null(opt$probe_map)) {
  pm <- read.table(opt$probe_map, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
  pm$gene <- ifelse(nchar(pm$ensembl_id) > 0, pm$ensembl_id, pm$gene_symbol)
  pm <- pm[nchar(pm$gene) > 0, ]
  probe_map <- setNames(pm$gene, pm$probe_name)
  sym_map   <- setNames(pm$gene_symbol, pm$gene)
  sym_map   <- sym_map[!duplicated(names(sym_map))]
  hit <- sum(rownames(counts) %in% names(probe_map))
  cat("Probe map:", hit, "/", nrow(counts), "rows mapped to",
      length(unique(probe_map)), "genes\n")
} else {
  warning("No --probe_map: treating rows as genes. This SKIPS the EPA ",
          "max-|L2FC| probe aggregation and is NOT EPA-faithful.")
  cat("!! WARNING: no probe map — max|L2FC| aggregation skipped (not EPA-faithful)\n")
  sym_map <- setNames(rownames(counts), rownames(counts))
}

# Cytotoxic concentrations to exclude
drop_key <- character(0)
if (!is.null(opt$exclude_doses) && file.exists(opt$exclude_doses)) {
  ex <- read.csv(opt$exclude_doses, stringsAsFactors = FALSE)
  drop_key <- paste(ex$Treatment, as.numeric(ex$Dose), sep = "|")
  cat("Cytotoxic concentrations excluded:", length(drop_key), "\n")
}

# =============================================================================
#  VEHICLE RESOLUTION — each chemical is compared ONLY against its own vehicle
#
#  A panel can mix vehicles (e.g. most chemicals in DMSO, poorly-soluble ones in
#  MeOH). Pooling all control wells into one reference would push a vehicle
#  effect straight into every L2FC. EPA's design compares each chemical against
#  plate-matched wells of ITS OWN vehicle; so does this.
#
#  Resolution order:
#    1. metadata column named by --vehicle_col (applies to treated AND control rows)
#    2. --vehicle_map CSV (Treatment,Vehicle)
#    3. fall back to --control for everything (single-vehicle panel)
# =============================================================================
is_ctrl <- meta$Type == "control"

meta$vehicle <- NA_character_
if (opt$vehicle_col %in% colnames(meta)) {
  meta$vehicle <- as.character(meta[[opt$vehicle_col]])
  cat("Vehicle source: metadata column '", opt$vehicle_col, "'\n", sep = "")
} else if (!is.null(opt$vehicle_map) && file.exists(opt$vehicle_map)) {
  vm  <- read.csv(opt$vehicle_map, stringsAsFactors = FALSE)
  lut <- setNames(as.character(vm$Vehicle), as.character(vm$Treatment))
  meta$vehicle <- unname(lut[meta$Treatment])
  cat("Vehicle source: --vehicle_map (", length(lut), "chemicals mapped )\n")
} else {
  cat("Vehicle source: none — assuming single vehicle '", opt$control, "'\n", sep = "")
}

# Control wells: their vehicle IS their treatment name unless stated otherwise
meta$vehicle[is_ctrl & is.na(meta$vehicle)] <- meta$Treatment[is_ctrl & is.na(meta$vehicle)]
# Treated wells with no assignment default to --control
meta$vehicle[is.na(meta$vehicle) | meta$vehicle == ""] <- opt$control

vehicles <- sort(unique(meta$vehicle[is_ctrl]))
if (length(vehicles) == 0) stop("No control (Type == 'control') samples found")

ctrl_by_vehicle <- lapply(vehicles, function(v)
  rownames(meta)[is_ctrl & meta$vehicle == v])
names(ctrl_by_vehicle) <- vehicles

cat("Vehicles detected:", length(vehicles), "\n")
for (v in vehicles)
  cat("  ", v, ": ", length(ctrl_by_vehicle[[v]]), " control wells\n", sep = "")

# Chemicals = everything that is not a control well. Vehicle names are excluded
# by construction (their rows are Type == 'control').
chems <- if (opt$treatments == "all") {
  sort(unique(meta$Treatment[!is_ctrl]))
} else trimws(strsplit(opt$treatments, ",")[[1]])
chems <- setdiff(chems, vehicles)   # belt and braces

chem_vehicle <- sapply(chems, function(ch) {
  v <- unique(meta$vehicle[meta$Treatment == ch & !is_ctrl])
  if (length(v) != 1) {
    cat("!! ", ch, ": ", length(v), " vehicles assigned (",
        paste(v, collapse = ", "), ") — using the first\n", sep = "")
  }
  v[1]
})
names(chem_vehicle) <- chems

cat("\nChemical -> vehicle:\n")
for (ch in chems) cat("  ", ch, " -> ", chem_vehicle[[ch]], "\n", sep = "")

# Fail loudly rather than silently borrowing the wrong controls
bad <- chems[!(chem_vehicle %in% vehicles)]
if (length(bad))
  stop("No control wells for vehicle(s) of: ", paste(bad, collapse = ", "),
       ". Available vehicles: ", paste(vehicles, collapse = ", "))
thin <- chems[sapply(chems, function(ch)
  length(ctrl_by_vehicle[[chem_vehicle[[ch]]]]) < 2)]
if (length(thin))
  stop("Fewer than 2 control wells for the vehicle of: ",
       paste(thin, collapse = ", "))

if (length(chems) < 3 && opt$cutoff_method == "panel_ci")
  cat("!! WARNING: panel_ci pools the", opt$n_low_conc,
      "lowest concs ACROSS chemicals.\n   With", length(chems),
      "chemical(s) the null is thin — consider --cutoff_method fixed\n")

# =============================================================================
#  STEP 1-2 — DESeq2 per chemical, probe-level, then max|L2FC| to gene
# =============================================================================
l2fc_by_chem <- list()   # chem -> genes x concentrations matrix
padj_by_chem <- list()

for (chem in chems) {
  veh <- chem_vehicle[[chem]]
  cat("\n--- ", chem, "  (vehicle: ", veh, ") ---\n", sep = "")
  trt <- rownames(meta)[meta$Treatment == chem & !is_ctrl]
  if (length(drop_key)) {
    k <- paste(chem, meta[trt, "dose_numeric"], sep = "|")
    n_dropped <- sum(k %in% drop_key)
    if (n_dropped > 0) {
      trt <- trt[!(k %in% drop_key)]
      cat("  cytotox filter: dropped", n_dropped, "samples\n")
    }
  }
  if (length(trt) == 0) { cat("  no samples — skipped\n"); next }

  # controls: same vehicle first, then plate-matched within that vehicle
  ctl_pool <- ctrl_by_vehicle[[veh]]
  ctl <- ctl_pool
  if (!is.null(plate_col)) {
    plates  <- unique(meta[trt, plate_col])
    matched <- ctl_pool[meta[ctl_pool, plate_col] %in% plates]
    if (length(matched) >= 2) {
      ctl <- matched
    } else {
      cat("  <2 plate-matched ", veh, " wells — using all ", length(ctl_pool),
          " ", veh, " wells (plate matching lost)\n", sep = "")
    }
  }
  cat("  controls used:", length(ctl), veh, "wells\n")

  sm  <- c(trt, ctl)
  cd  <- meta[sm, , drop = FALSE]
  cd$dose_group <- factor(make.names(paste0("d", cd$dose_numeric)))
  ref <- make.names("d0")
  if (!(ref %in% levels(cd$dose_group))) { cat("  no dose 0 — skipped\n"); next }
  cd$dose_group <- relevel(cd$dose_group, ref = ref)

  cnt <- counts[, sm, drop = FALSE]
  keep <- rowMeans(cnt) >= opt$min_mean_count          # EPA: mean count < 5 out
  cnt  <- cnt[keep, , drop = FALSE]
  cat("  samples:", ncol(cnt), "| probes after mean-count filter:", nrow(cnt),
      "/", length(keep), "\n")

  use_plate <- !is.null(plate_col) &&
    nlevels(as.factor(cd[[plate_col]])) > 1 &&
    !any(table(cd[[plate_col]], cd$dose_group) == 0 &
           rowSums(table(cd[[plate_col]], cd$dose_group)) > 0 &
           colSums(table(cd[[plate_col]], cd$dose_group)) > 0) # crude confound check
  if (use_plate) {
    cd$plate <- factor(cd[[plate_col]])
    design_f <- ~ plate + dose_group
  } else {
    design_f <- ~ dose_group
    if (!is.null(plate_col))
      cat("  plate confounded with dose or single-level — using ~ dose_group\n")
  }

  dds <- DESeqDataSetFromMatrix(cnt, colData = cd, design = design_f)

  # Wald p-values BEFORE shrinkage (betaPrior = FALSE)
  dds_np <- tryCatch(DESeq(dds, quiet = TRUE),
                     error = function(e) { cat("  DESeq failed:", e$message, "\n"); NULL })
  if (is.null(dds_np)) next

  # Shrunken L2FC (normal prior = what EPA used, DESeq2 v1.24 lfcShrink type='normal')
  dds_bp <- NULL
  if (opt$shrinkage == "normal") {
    dds_bp <- tryCatch(DESeq(dds, betaPrior = TRUE, quiet = TRUE),
                       error = function(e) {
                         cat("  betaPrior refit failed:", e$message,
                             "— falling back to ashr\n"); NULL })
  }

  dose_levels <- setdiff(levels(cd$dose_group), ref)
  dose_vals   <- sapply(dose_levels, function(g)
    unique(cd$dose_numeric[cd$dose_group == g])[1])
  ord <- order(dose_vals)
  dose_levels <- dose_levels[ord]; dose_vals <- dose_vals[ord]

  lfc_mat  <- matrix(NA_real_, nrow(cnt), length(dose_levels),
                     dimnames = list(rownames(cnt), as.character(dose_vals)))
  padj_mat <- lfc_mat

  for (i in seq_along(dose_levels)) {
    g <- dose_levels[i]
    res <- results(dds_np, contrast = c("dose_group", g, ref),
                   independentFiltering = FALSE, pAdjustMethod = "BH")
    # BH within this concentration (independent filtering off)
    padj_mat[, i] <- p.adjust(res$pvalue, method = "BH")

    lfc <- if (!is.null(dds_bp)) {
      results(dds_bp, contrast = c("dose_group", g, ref),
              independentFiltering = FALSE)$log2FoldChange
    } else if (opt$shrinkage == "none") {
      res$log2FoldChange
    } else {
      sh <- tryCatch(lfcShrink(dds_np, contrast = c("dose_group", g, ref),
                               type = if (opt$shrinkage == "normal") "ashr"
                                      else opt$shrinkage, quiet = TRUE),
                     error = function(e) res)
      sh$log2FoldChange
    }
    lfc_mat[, i] <- lfc
  }

  lfc_mat[is.na(lfc_mat)] <- 0     # EPA: missing L2FC set to zero

  # ---- probe -> gene: highest magnitude L2FC in EITHER direction ----
  if (!is.null(probe_map)) {
    gvec <- probe_map[rownames(lfc_mat)]
    ok   <- !is.na(gvec)
    lm2  <- lfc_mat[ok, , drop = FALSE]; pm2 <- padj_mat[ok, , drop = FALSE]
    gv   <- gvec[ok]
    genes <- unique(gv)
    gl <- matrix(0, length(genes), ncol(lm2),
                 dimnames = list(genes, colnames(lm2)))
    gp <- matrix(NA_real_, length(genes), ncol(lm2),
                 dimnames = list(genes, colnames(lm2)))
    idx_by_gene <- split(seq_along(gv), gv)
    for (gname in genes) {
      ii <- idx_by_gene[[gname]]
      sub <- lm2[ii, , drop = FALSE]
      win <- apply(abs(sub), 2, which.max)            # max |L2FC|, sign kept
      gl[gname, ] <- sub[cbind(win, seq_len(ncol(sub)))]
      subp <- pm2[ii, , drop = FALSE]
      gp[gname, ] <- subp[cbind(win, seq_len(ncol(subp)))]
    }
    cat("  probes -> genes:", nrow(lm2), "->", nrow(gl), "(max |L2FC| rule)\n")
  } else {
    gl <- lfc_mat; gp <- padj_mat
  }

  l2fc_by_chem[[chem]] <- gl
  padj_by_chem[[chem]] <- gp
  write.csv(gl, file.path(opt$outdir, paste0("l2fc_gene_", make.names(chem), ".csv")))
}

if (length(l2fc_by_chem) == 0) stop("No chemical produced an L2FC matrix")

# =============================================================================
#  STEP 3 — gene retention: L2FC in >= 95% of treatments
# =============================================================================
all_genes <- sort(unique(unlist(lapply(l2fc_by_chem, rownames))))
n_treat   <- sum(sapply(l2fc_by_chem, ncol))     # chemical x concentration
present   <- sapply(all_genes, function(g)
  sum(sapply(l2fc_by_chem, function(m) if (g %in% rownames(m)) ncol(m) else 0)))
keep_genes <- all_genes[present / n_treat >= opt$gene_presence]
cat("\nGene retention (>=", 100 * opt$gene_presence, "% of", n_treat,
    "treatments):", length(keep_genes), "/", length(all_genes), "\n")

# pad missing genes with 0 (EPA: missing L2FC set to zero)
for (chem in names(l2fc_by_chem)) {
  m <- l2fc_by_chem[[chem]]
  miss <- setdiff(keep_genes, rownames(m))
  if (length(miss)) {
    add <- matrix(0, length(miss), ncol(m), dimnames = list(miss, colnames(m)))
    m <- rbind(m, add)
  }
  l2fc_by_chem[[chem]] <- m[keep_genes, , drop = FALSE]
}

# =============================================================================
#  STEP 4 — PANEL NOISE BAND from the n lowest concentrations, all chemicals
# =============================================================================
cat("\n=== Panel noise band (", opt$cutoff_method, ") ===\n", sep = "")

# NOTE ON MIXED VEHICLES
#   L2FC values are already expressed relative to each chemical's OWN vehicle,
#   so a DMSO chemical and a MeOH chemical are on a common, zero-centred scale
#   and can legitimately share a null. That is the default (--noise_by_vehicle
#   FALSE) and it keeps the null wide enough to estimate a 95% interval.
#   If a vehicle has its own noise characteristics AND enough chemicals to
#   support it, --noise_by_vehicle TRUE computes a separate band per vehicle.
#   With only 1-2 chemicals on a minor vehicle, do NOT turn this on: the null
#   collapses to 2-4 values per gene and the cutoff becomes meaningless.
fitted_chems <- names(l2fc_by_chem)
noise_groups <- if (opt$noise_by_vehicle) {
  split(fitted_chems, unlist(chem_vehicle[fitted_chems]))
} else {
  list(panel = fitted_chems)
}

if (opt$noise_by_vehicle) {
  cat("Noise band computed SEPARATELY per vehicle:\n")
  for (g in names(noise_groups)) {
    n_ch <- length(noise_groups[[g]])
    cat("  ", g, ": ", n_ch, " chemical(s) -> ",
        n_ch * opt$n_low_conc, " null values per gene\n", sep = "")
    if (n_ch < 3)
      cat("     !! too few chemicals for a stable 95% interval — ",
          "consider --noise_by_vehicle FALSE\n", sep = "")
  }
} else {
  cat("Noise band POOLED across all vehicles (L2FC are vehicle-relative)\n")
}

low_concs <- function(m) {
  cn <- as.numeric(colnames(m))
  m[, order(cn)[seq_len(min(opt$n_low_conc, ncol(m)))], drop = FALSE]
}

cutoff_by_group <- list(); onesd_by_group <- list(); noise_meta <- list()

for (g in names(noise_groups)) {
  null_mat <- do.call(cbind, lapply(l2fc_by_chem[noise_groups[[g]]], low_concs))

  onesd_vec <- apply(null_mat, 1, stats::sd)
  med_sd <- stats::median(onesd_vec[is.finite(onesd_vec) & onesd_vec > 0], na.rm = TRUE)
  if (!is.finite(med_sd) || med_sd <= 0) med_sd <- opt$cutoff_floor
  onesd_vec[!is.finite(onesd_vec) | onesd_vec <= 0] <- med_sd

  cutoff_vec <- switch(
    opt$cutoff_method,
    "panel_ci"   = apply(null_mat, 1, function(v)
      max(abs(stats::quantile(v, c(0.025, 0.975), na.rm = TRUE)))),  # 95% interval bounds
    "panel_bmad" = 3 * apply(null_mat, 1, stats::mad),
    "fixed"      = rep(opt$cutoff_value, nrow(null_mat)),
    stop("Unknown --cutoff_method")
  )
  names(cutoff_vec) <- rownames(null_mat)
  cutoff_vec[!is.finite(cutoff_vec)] <- opt$cutoff_floor
  cutoff_vec <- pmax(cutoff_vec, opt$cutoff_floor)

  cutoff_by_group[[g]] <- cutoff_vec
  onesd_by_group[[g]]  <- onesd_vec

  cat("  [", g, "] null =", ncol(null_mat), "values/gene | median cutoff",
      round(stats::median(cutoff_vec), 4), "| median onesd",
      round(stats::median(onesd_vec), 4), "\n")

  noise_meta[[g]] <- list(
    chemicals      = noise_groups[[g]],
    n_null_values  = ncol(null_mat),
    median_cutoff  = round(stats::median(cutoff_vec), 5),
    median_onesd   = round(stats::median(onesd_vec), 5))

  write.csv(data.frame(gene = names(cutoff_vec),
                       gene_symbol = unname(sym_map[names(cutoff_vec)]),
                       cutoff = round(cutoff_vec, 5),
                       onesd  = round(onesd_vec[names(cutoff_vec)], 5)),
            file.path(opt$outdir, paste0("noise_band_", make.names(g), ".csv")),
            row.names = FALSE)
}
cat("BMR = ", opt$bmr_scale, " x onesd\n", sep = "")

# =============================================================================
#  STEP 5-6 — tcplfit2 per gene, per chemical; tPODs
# =============================================================================
fitmodels <- trimws(strsplit(opt$fitmodels, ",")[[1]])
up_bnd <- if (is.na(opt$bmd_up_bnd)) NULL else opt$bmd_up_bnd
cat("\nModels:", paste(fitmodels, collapse = ", "),
    "| AICc:", opt$aicc, "| bmd_low_bnd:", opt$bmd_low_bnd,
    "| bmd_up_bnd:", ifelse(is.null(up_bnd), "none", up_bnd), "\n")

tpod_table <- list()
panel_summary <- list()

for (chem in names(l2fc_by_chem)) {
  m    <- l2fc_by_chem[[chem]]
  conc <- as.numeric(colnames(m))
  ord  <- order(conc); conc <- conc[ord]; m <- m[, ord, drop = FALSE]

  veh   <- chem_vehicle[[chem]]
  grp   <- if (opt$noise_by_vehicle) veh else "panel"
  cutoff_vec <- cutoff_by_group[[grp]]
  onesd_vec  <- onesd_by_group[[grp]]

  cat("\n=== Fitting", chem, "| vehicle", veh, "| noise band:", grp,
      "|", nrow(m), "genes x", length(conc), "concentrations ===\n")

  fit_one <- function(g) {
    row <- list(conc = conc, resp = as.numeric(m[g, ]), bmed = 0,
                cutoff = as.numeric(cutoff_vec[[g]]),
                onesd  = as.numeric(onesd_vec[[g]]),
                name = g, assay = chem)
    tryCatch(
      tcplfit2::concRespCore(
        row, fitmodels = fitmodels, conthits = TRUE, aicc = opt$aicc,
        force.fit = opt$force_fit, bidirectional = opt$bidirectional,
        verbose = FALSE, do.plot = FALSE, bmr_scale = opt$bmr_scale,
        bmd_low_bnd = opt$bmd_low_bnd, bmd_up_bnd = up_bnd),
      error = function(e)
        data.frame(name = g, assay = chem, fit_method = NA_character_,
                   hitcall = NA_real_, bmd = NA_real_, bmdl = NA_real_,
                   bmdu = NA_real_, top = NA_real_, ac50 = NA_real_,
                   cutoff = as.numeric(cutoff_vec[[g]]),
                   fit_error = conditionMessage(e), stringsAsFactors = FALSE))
  }

  t0 <- Sys.time()
  res <- if (.Platform$OS.type == "unix" && n_cores > 1)
    parallel::mclapply(rownames(m), fit_one, mc.cores = n_cores) else
    lapply(rownames(m), fit_one)
  fits <- dplyr::bind_rows(res)
  cat("  fitted in", round(difftime(Sys.time(), t0, units = "mins"), 2), "min |",
      sum(is.na(fits$fit_method)), "failures\n")

  fits$gene        <- fits$name
  fits$gene_symbol <- unname(sym_map[fits$gene])
  fits$gene_symbol[is.na(fits$gene_symbol) | fits$gene_symbol == ""] <-
    fits$gene[is.na(fits$gene_symbol) | fits$gene_symbol == ""]
  if (!"top_over_cutoff" %in% colnames(fits) &&
      all(c("top", "cutoff") %in% colnames(fits)))
    fits$top_over_cutoff <- abs(fits$top) / fits$cutoff

  write.csv(fits, file.path(opt$outdir,
                            paste0("tcplfit2_fits_", make.names(chem), ".csv")),
            row.names = FALSE)

  # ---- ACTIVE genes: hitcall > 0.9 (EPA §2.10.1) ----
  act <- fits[!is.na(fits$hitcall) & fits$hitcall > opt$hitcall_thresh &
                !is.na(fits$bmd) & fits$bmd > 0 &
                !is.na(fits$fit_method) & fits$fit_method != "cnst", ]
  if (opt$tc_min > 0 && "top_over_cutoff" %in% colnames(act))
    act <- act[!is.na(act$top_over_cutoff) & act$top_over_cutoff >= opt$tc_min, ]
  act <- act[order(act$bmd), ]
  cat("  active genes (hitcall >", opt$hitcall_thresh, "):", nrow(act), "\n")

  if (nrow(act) == 0) {
    tpod_table[[chem]] <- data.frame(
      treatment = chem, vehicle = veh, n_active = 0,
      gene_05 = NA, gene_abs5 = NA, gene_min = NA)
    panel_summary[[chem]] <- list(vehicle = veh, noise_band = grp,
                                  n_active = 0, tpods_epa = NULL)
    cat("  0 active genes — tPODs = NA (EPA convention)\n")
    next
  }

  bmds <- sort(act$bmd)
  gene_05   <- unname(stats::quantile(bmds, 0.05))
  gene_abs5 <- if (length(bmds) >= 5) bmds[5] else NA_real_
  gene_min  <- bmds[1]

  cat("  gene_05  :", signif(gene_05, 4), "uM\n")
  cat("  gene_abs5:", signif(gene_abs5, 4), "uM\n")
  cat("  gene_min :", signif(gene_min, 4), "uM\n")

  tp <- list(
    gene_05   = list(value = round(gene_05, 6),
                     label = "5th percentile of active-gene BMDs",
                     ref = "Harrill et al. 2024 §2.10.1"),
    gene_abs5 = list(value = round(gene_abs5, 6),
                     label = "5th lowest active-gene BMD",
                     ref = "Harrill et al. 2024 §2.10.1"),
    gene_min  = list(value = round(gene_min, 6),
                     label = "lowest active-gene BMD",
                     ref = "Harrill et al. 2024 §2.10.1 (over-sensitive; 10-100x low vs ER mPOD)")
  )

  # Non-EPA cross-reference tPODs (same code as run_dromics.R) — clearly labelled
  tp_extra <- NULL
  if (opt$also_dromics_tpods) {
    n <- length(bmds)
    tp_extra <- list()
    if (n >= 25) {
      tp_extra$rank25 <- round(bmds[25], 6)          # Reardon 2021
      tp_extra$perc05 <- round(unname(stats::quantile(bmds, 0.05)), 6)
      lb <- log10(bmds); ec <- seq_len(n) / n
      xn <- (lb - min(lb)) / max(1e-12, max(lb) - min(lb))
      tp_extra$first_mode <- round(bmds[which.max(abs(ec - xn))], 6)  # Kneedle
    }
    if (n >= 10) tp_extra$perc10 <- round(unname(stats::quantile(bmds, 0.10)), 6)
    tp_extra$median_all <- round(stats::median(bmds), 6)
    tp_extra$NOTE <- "NOT part of the EPA workflow; for cross-reference with run_dromics.R only"
  }

  write.csv(act, file.path(opt$outdir,
                           paste0("active_genes_", make.names(chem), ".csv")),
            row.names = FALSE)

  tpod_table[[chem]] <- data.frame(
    treatment = chem, vehicle = veh, n_active = nrow(act),
    gene_05 = gene_05, gene_abs5 = gene_abs5, gene_min = gene_min)
  panel_summary[[chem]] <- list(
    vehicle = veh, noise_band = grp,
    n_genes_fitted = nrow(fits), n_active = nrow(act),
    median_bmd = round(stats::median(bmds), 6),
    tpods_epa = tp, tpods_non_epa = tp_extra,
    model_distribution = as.list(table(act$fit_method)))

  # BMD accumulation plot with tPOD markers
  df <- data.frame(bmd = bmds, frac = seq_along(bmds) / length(bmds))
  mk <- data.frame(
    method = c("gene_05", "gene_abs5", "gene_min"),
    value  = c(gene_05, gene_abs5, gene_min))
  mk <- mk[!is.na(mk$value), ]
  p <- ggplot(df, aes(bmd, frac)) +
    geom_step(colour = "#1a6b3a", linewidth = 0.9) +
    geom_vline(data = mk, aes(xintercept = value, colour = method),
               linetype = "dashed", linewidth = 0.7) +
    scale_x_log10() +
    labs(title = paste("Active-gene BMD accumulation —", chem),
         subtitle = paste0(nrow(act), " active genes (hitcall > ",
                           opt$hitcall_thresh, ") | EPA gene-level tcplfit2"),
         x = expression(paste("BMD (", mu, "M, log"[10], ")")),
         y = "Cumulative fraction of active genes", colour = "tPOD") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(opt$outdir, paste0("bmd_accumulation_", make.names(chem), ".png")),
         p, width = 9, height = 6, dpi = 300)
}

# =============================================================================
#  PANEL OUTPUT
# =============================================================================
tp_df <- dplyr::bind_rows(tpod_table)
write.csv(tp_df, file.path(opt$outdir, "tpods_epa_gene_level.csv"), row.names = FALSE)

cat("\n=== Panel tPODs (uM) ===\n")
print(tp_df, row.names = FALSE)

write_json(list(
  workflow  = "EPA gene-level HTTr tPOD (Harrill et al. 2024, Toxicology 501:153694)",
  tcplfit2_version = as.character(utils::packageVersion("tcplfit2")),
  deseq2_version   = as.character(utils::packageVersion("DESeq2")),
  preprocessing = list(
    unit_modelled     = if (is.null(probe_map)) "gene (NO probe map — deviation)" else "probe",
    probe_to_gene     = "max |L2FC| in either direction",
    min_mean_count    = opt$min_mean_count,
    design            = if (is.null(plate_col)) "~ dose_group (no plate column)"
                        else paste("~", plate_col, "+ dose_group"),
    shrinkage         = opt$shrinkage,
    pvalues           = "Wald before shrinkage; BH within concentration; independentFiltering=FALSE",
    gene_presence_min = opt$gene_presence,
    cytotox_excluded  = length(drop_key)
  ),
  vehicles = list(
    detected          = vehicles,
    control_wells     = lapply(ctrl_by_vehicle, length),
    chemical_vehicle  = as.list(chem_vehicle),
    matching          = "each chemical modelled against plate-matched wells of its OWN vehicle"
  ),
  noise_band = list(
    method          = opt$cutoff_method,
    n_lowest_concs  = opt$n_low_conc,
    per_vehicle     = opt$noise_by_vehicle,
    groups          = noise_meta,
    bmr_scale       = opt$bmr_scale,
    CAVEAT          = paste("Gene-level cutoff derivation is inferred from the",
                            "signature-level rule in Harrill 2024 §2.7; §2.10.1",
                            "does not state it explicitly.")
  ),
  fitting = list(
    fitmodels = fitmodels, aicc = opt$aicc, force_fit = opt$force_fit,
    bidirectional = opt$bidirectional,
    bmd_low_bnd = opt$bmd_low_bnd, bmd_up_bnd = opt$bmd_up_bnd,
    hitcall_thresh = opt$hitcall_thresh, tc_min = opt$tc_min
  ),
  per_chemical = panel_summary
), file.path(opt$outdir, "tpod_epa_summary.json"),
  pretty = TRUE, auto_unbox = TRUE, na = "null")

cat("\n=== Done ===\n")
cat("Outputs in:", opt$outdir, "\n")
cat("  tpods_epa_gene_level.csv       — gene_05 / gene_abs5 / gene_min per chemical\n")
cat("  panel_noise_band.csv           — per-gene cutoff + onesd\n")
cat("  tcplfit2_fits_<chem>.csv       — every gene, every model, hitcalls\n")
cat("  active_genes_<chem>.csv        — hitcall > threshold\n")
cat("  tpod_epa_summary.json          — full provenance incl. cutoff caveat\n")
