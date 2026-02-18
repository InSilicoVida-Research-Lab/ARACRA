#!/usr/bin/env python3
"""
Integrated RNA-seq Processing Pipeline - GUI Module with DEG Analysis
Streamlit-based interface combining RNA-seq pipeline and DEG analysis
"""

import streamlit as st
import pandas as pd
import threading
import time
from pathlib import Path
import os
import json
import traceback
from datetime import datetime
import re
import subprocess
import tempfile
import shutil
import zipfile
import numpy as np
import glob     # <--- Make sure this is here
import sys      # <--- ADD THIS

# Import core RNA-seq functionality
from rnaseq_pipeline_core import (
    PipelineConfig,
    ProcessingResult,
    ThreadSafeLogger,
    ResourceMonitor,
    validate_hisat2_index,
    validate_tools,
    create_directories,
    run_pipeline_thread,
    export_results_to_csv,
    create_processing_summary,
    validate_accession_format
)

# Page configuration
try:
    st.set_page_config(
        page_title="ARACRA",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded"
    )
except:
    pass


def auto_detect_configuration():
    """
    Aggressively detect system resources and tool paths.
    Derives paths from the active Python environment to ensure accuracy.
    """
    config = {
        'cpu_count': 1,
        'default_threads': 1,
        'trimmomatic_jar': '',
        'adapters': '',
        'base_dir': str(Path.home() / "rnaseq_pipeline")
    }
    
    # --- 1. CPU Detection ---
    try:
        total_cpus = os.cpu_count() or 4
        # Requirement: Process at half the number of cores available
        default_threads = max(1, int(total_cpus * 0.5))
        config['cpu_count'] = total_cpus
        config['default_threads'] = default_threads
    except:
        config['cpu_count'] = 4
        config['default_threads'] = 2

    # --- 2. Environment Path Detection ---
    # Try environment variable first, fall back to Python executable path
    # (e.g., /home/user/conda/envs/aracra/bin/python -> /home/user/conda/envs/aracra)
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if not conda_prefix:
        conda_prefix = os.path.dirname(os.path.dirname(sys.executable))
    
    # --- 3. Trimmomatic JAR Search ---
    search_patterns = [
        # Standard Bioconda location
        os.path.join(conda_prefix, "share", "trimmomatic*", "trimmomatic.jar"),
        # Alternate location (sometimes in java folder)
        os.path.join(conda_prefix, "share", "java", "trimmomatic*", "trimmomatic.jar"),
        # Linux system default
        "/usr/share/java/trimmomatic.jar",
        "/usr/share/java/trimmomatic*.jar"
    ]

    for pattern in search_patterns:
        found = glob.glob(pattern)
        if found:
            config['trimmomatic_jar'] = found[0]
            break
            
    # --- 4. Adapter Search ---
    # If JAR was found, look for adapters relative to it
    if config['trimmomatic_jar']:
        jar_dir = os.path.dirname(config['trimmomatic_jar'])
        
        adapter_patterns = [
            # Standard folder structure: .../share/trimmomatic/adapters/
            os.path.join(jar_dir, "adapters", "TruSeq3-PE.fa"),
            # Sometimes adapters are in the root of the share folder
            os.path.join(jar_dir, "TruSeq3-PE.fa"),
            # Check for ANY .fa file in adapters folder as fallback
            os.path.join(jar_dir, "adapters", "*.fa")
        ]
        
        for pattern in adapter_patterns:
            found = glob.glob(pattern)
            if found:
                config['adapters'] = found[0]
                break
    
    # --- 5. Fallbacks if nothing found ---
    if not config['trimmomatic_jar']:
        config['trimmomatic_jar'] = "/path/to/trimmomatic.jar" # Placeholder
        
    if not config['adapters']:
        config['adapters'] = "/path/to/adapters.fa" # Placeholder

    return config

def init_session_state():
    """Initialize session state variables"""
    
    # Run auto-detection
    detected_config = auto_detect_configuration()
    
    defaults = {
        # Auto-detected defaults
        'detected_cpu_max': detected_config['cpu_count'],
        'detected_cpu_default': detected_config['default_threads'],
        'detected_trim_jar': detected_config['trimmomatic_jar'],
        'detected_adapters': detected_config['adapters'],
        'detected_base_dir': detected_config['base_dir'],

        # RNA-seq pipeline states
        'processing_status': "idle",
        'log_messages': [],
        'accession_numbers': [],
        'pipeline_thread': None,
        'progress_value': 0,
        'current_step': "",
        'total_accessions': 0,
        'completed_accessions': 0,
        'pipeline_results': {},
        'last_log_update': time.time(),
        'thread_error': None,
        'pipeline_config': None,
        'current_memory': 0,
        'current_cpu': 0,
        'peak_memory': 0,
        'peak_cpu': 0,
        'resource_monitor': None,
        'cancelled_accessions': set(),
        'processing_start_time': None,
        'estimated_completion': None,
        'pipeline_logger': None,
        
        # DEG analysis states
        'deg_analysis_complete': False,
        'deg_output_dir': None,
        'deg_summary': None,
        'counts_path': None,
        'metadata_path': None,
        'metadata_df': None,
        'chemicals': None,
        'conditions': None,

        # Dose-response analysis states
        'dr_analysis_complete': False,
        'dr_output_dir': None,
    }
    
    for key, default_value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = default_value

def add_custom_css():
    """Add custom CSS styling"""
    st.markdown("""
    <style>
        .main-header {
            font-size: 2.5rem;
            color: #1f77b4;
            text-align: center;
            margin-bottom: 2rem;
            font-weight: 600;
        }
        .status-running { color: #28a745; font-weight: bold; }
        .status-error { color: #dc3545; font-weight: bold; }
        .status-idle { color: #6c757d; font-weight: bold; }
        .deg-header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 1.5rem;
            border-radius: 10px;
            text-align: center;
            margin-bottom: 1rem;
        }
        .metric-card {
            background: white;
            padding: 1rem;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 1rem;
        }
    </style>
    """, unsafe_allow_html=True)

# --- Analysis Script Generation and Execution ---

def generate_deg_r_script(params):
    """Generate comprehensive R-ODAF script for DEG analysis"""
    
    r_script = f"""#!/usr/bin/env Rscript
#================================================================================
# INTEGRATED R-ODAF PIPELINE: QC + DIFFERENTIAL EXPRESSION ANALYSIS
# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
#================================================================================

# Clear workspace and set up environment
rm(list = ls())
cat("=== INTEGRATED R-ODAF ANALYSIS PIPELINE ===\\n")

#================================================================================
# STEP 1: LOAD REQUIRED LIBRARIES
#================================================================================

cat("\\n=== STEP 1: LOADING LIBRARIES ===\\n")

suppressPackageStartupMessages({{
  library(dplyr)          # Data manipulation
  library(readxl)         # Excel file reading
  library(DESeq2)         # RNA-seq analysis
  library(edgeR)          # RNA-seq analysis
  library(sva)            # Batch correction
  library(ggplot2)        # Plotting
  library(ggrepel)        # Plot labels
  library(stringr)        # String manipulation
  library(org.Hs.eg.db)   # Human gene annotation
}})

cat("Libraries loaded successfully!\\n")

#================================================================================
# STEP 2: LOADING INTEGRATED FUNCTIONS
#================================================================================

cat("\\n=== STEP 2: LOADING INTEGRATED FUNCTIONS ===\\n")

# Source ComBat-seq functions if available
tryCatch({{
  source("/home/saurav/Downloads/ComBat_seq.R")
  source("/home/saurav/Downloads/helper_seq.R")
  cat("ComBat-seq functions loaded successfully\\n")
}}, error = function(e) {{
  cat("Warning: ComBat-seq functions not found. Will use sva::ComBat instead.\\n")
}})

#================================================================================
# QC FUNCTIONS
#================================================================================

read_featurecounts_output <- function(file_path) {{
  cat("\\n=== READING FEATURECOUNTS OUTPUT ===\\n")
  
  first_lines <- readLines(file_path, n = 5)
  cat("Examining file structure...\\n")
  
  if (grepl("^#", first_lines[1])) {{
    skip_lines <- 1
  }} else {{
    skip_lines <- 0
  }}
  
  tryCatch({{
    raw_data <- read.table(file_path, header = TRUE, sep = "\\t", 
                           comment.char = "", stringsAsFactors = FALSE, 
                           skip = skip_lines, check.names = FALSE)
  }}, error = function(e1) {{
    raw_data <<- read.table(file_path, header = TRUE, sep = "\\t", 
                            comment.char = "", stringsAsFactors = FALSE, 
                            skip = skip_lines, fill = TRUE, check.names = FALSE)
  }})
  
  cat("Successfully loaded data with dimensions:", nrow(raw_data), "x", ncol(raw_data), "\\n")
  
  if ("Length" %in% colnames(raw_data)) {{
    length_col_idx <- which(colnames(raw_data) == "Length")
  }} else {{
    length_col_idx <- 6
  }}
  
  count_cols <- (length_col_idx + 1):ncol(raw_data)
  count_matrix <- raw_data[, count_cols, drop = FALSE]
  
  if ("Geneid" %in% colnames(raw_data)) {{
    rownames(count_matrix) <- raw_data$Geneid
  }} else {{
    rownames(count_matrix) <- raw_data[, 1]
  }}
  
  sample_names <- colnames(count_matrix)
  clean_names <- str_extract(sample_names, "SRR[0-9]+")
  
  if (any(is.na(clean_names))) {{
    alt_names <- str_extract(sample_names, "[A-Za-z]+[0-9]+")
    clean_names[is.na(clean_names)] <- alt_names[is.na(clean_names)]
    
    still_na <- is.na(clean_names)
    if (any(still_na)) {{
      clean_names[still_na] <- basename(sample_names[still_na])
      clean_names[still_na] <- gsub("\\\\.[^.]*$", "", clean_names[still_na])
    }}
  }}
  
  colnames(count_matrix) <- clean_names
  count_matrix <- as.matrix(count_matrix)
  mode(count_matrix) <- "numeric"
  
  if (any(is.na(count_matrix))) {{
    cat("Warning: Found NA values in count matrix. Converting to 0...\\n")
    count_matrix[is.na(count_matrix)] <- 0
  }}
  
  cat("Successfully parsed featureCounts output:\\n")
  cat("- Genes:", nrow(count_matrix), "\\n")
  cat("- Samples:", ncol(count_matrix), "\\n")
  
  return(list(counts = count_matrix, original_names = sample_names))
}}

validate_metadata_with_batch <- function(count_data, metadata_df) {{
  cat("\\n=== VALIDATING METADATA WITH BATCH INFO ===\\n")
  
  if (!"Run" %in% colnames(metadata_df)) {{
    stop("Error: 'Run' column not found in metadata.")
  }}
  
  rownames(metadata_df) <- metadata_df$Run
  
  count_samples <- colnames(count_data$counts)
  metadata_samples <- rownames(metadata_df)
  common_samples <- intersect(count_samples, metadata_samples)
  
  cat("Count data samples:", length(count_samples), "\\n")
  cat("Metadata samples:", length(metadata_samples), "\\n")
  cat("Common samples:", length(common_samples), "\\n")
  
  if (length(common_samples) == 0) {{
    stop("Error: No common samples found between count data and metadata.")
  }}
  
  batch_cols <- grep("batch|Batch|BATCH", colnames(metadata_df), value = TRUE)
  cat("Potential batch columns found:", paste(batch_cols, collapse = ", "), "\\n")
  
  aligned_counts <- count_data$counts[, common_samples]
  aligned_metadata <- metadata_df[common_samples, ]
  
  return(list(counts = aligned_counts, metadata = aligned_metadata))
}}

apply_combat_batch_correction <- function(counts, metadata) {{
  cat("\\n=== APPLYING BATCH CORRECTION ===\\n")
  
  batch_col <- NULL
  possible_batch_cols <- c("batch", "Batch", "BATCH", "batch_id", "Batch_ID")
  
  for (col in possible_batch_cols) {{
    if (col %in% colnames(metadata)) {{
      batch_col <- col
      break
    }}
  }}
  
  if (is.null(batch_col)) {{
    cat("No batch column found. Batch correction will be skipped.\\n")
    return(list(corrected_counts = counts, batch_corrected = FALSE))
  }}
  
  common_samples <- intersect(colnames(counts), rownames(metadata))
  counts <- counts[, common_samples]
  metadata <- metadata[common_samples, ]
  
  batch_factor <- as.factor(metadata[[batch_col]])
  
  if (length(unique(batch_factor)) < 2) {{
    cat("Only one batch detected. Batch correction not applicable.\\n")
    return(list(corrected_counts = counts, batch_corrected = FALSE))
  }}
  
  if (!"group" %in% colnames(metadata)) {{
    if ("condition" %in% colnames(metadata)) {{
      group_factor <- as.factor(metadata$condition)
    }} else {{
      cat("Error: No group or condition column found.\\n")
      return(list(corrected_counts = counts, batch_corrected = FALSE))
    }}
  }} else {{
    group_factor <- as.factor(metadata$group)
  }}
  
  tryCatch({{
    if (!exists("ComBat_seq")) {{
      library(sva)
      cpm_data <- edgeR::cpm(counts, log = TRUE, prior.count = 1)
      corrected_log_data <- ComBat(
        dat = cpm_data,
        batch = batch_factor,
        mod = model.matrix(~group_factor),
        par.prior = TRUE
      )
      
      corrected_counts_matrix <- 2^corrected_log_data - 1
      corrected_counts_matrix[corrected_counts_matrix < 0] <- 0
      corrected_counts_matrix <- round(corrected_counts_matrix)
      
      cat("ComBat (sva) batch correction completed successfully.\\n")
      return(list(corrected_counts = corrected_counts_matrix, batch_corrected = TRUE))
    }} else {{
      corrected_counts_matrix <- ComBat_seq(
        counts = counts,
        batch = batch_factor,
        group = group_factor,
        shrink = FALSE
      )
      
      cat("ComBat-seq batch correction completed successfully.\\n")
      return(list(corrected_counts = corrected_counts_matrix, batch_corrected = TRUE))
    }}
  }}, error = function(e) {{
    cat("ERROR: Batch correction failed:", e$message, "\\n")
    cat("Returning original counts.\\n")
    return(list(corrected_counts = counts, batch_corrected = FALSE))
  }})
}}

remove_low_read_samples <- function(counts, metadata, threshold = {params.get('read_threshold', 1000000)}) {{
  cat("\\n=== R-ODAF SAMPLE FILTERING (Read Count) ===\\n")
  library_size <- colSums(counts)
  
  low_read_samples <- names(library_size[library_size < threshold])
  
  if (length(low_read_samples) > 0) {{
    cat("Removing samples with <", format(threshold, scientific = FALSE), "reads:", length(low_read_samples), "\\n")
    for (sample in low_read_samples) {{
      reads <- format(library_size[sample], big.mark = ",")
      cat("  -", sample, "(", reads, "reads )\\n")
    }}
    counts <- counts[, !(colnames(counts) %in% low_read_samples)]
    metadata <- metadata[!(rownames(metadata) %in% low_read_samples), ]
  }} else {{
    cat("All samples pass read count filter\\n")
  }}
  
  cat("Samples retained:", ncol(counts), "\\n")
  return(list(counts = counts, metadata = metadata, removed_samples = low_read_samples))
}}

filter_genes_by_cpm <- function(counts, metadata) {{
  cat("\\n=== R-ODAF RELEVANCE FILTER ===\\n")
  cpm_matrix <- edgeR::cpm(counts)
  genes_to_keep_mask <- rep(FALSE, nrow(cpm_matrix))
  metadata$group <- as.factor(metadata$group)
  
  for (group_level in levels(metadata$group)) {{
    group_samples <- rownames(metadata[metadata$group == group_level, ])
    group_samples_in_matrix <- intersect(group_samples, colnames(cpm_matrix))
    if(length(group_samples_in_matrix) == 0) next
    
    replicates_in_group_cpm <- cpm_matrix[, group_samples_in_matrix, drop = FALSE]
    proportion_expressed <- rowSums(replicates_in_group_cpm >= {params.get('cpm_threshold', 1)}) / length(group_samples_in_matrix)
    genes_to_keep_mask[proportion_expressed >= {params.get('min_prop_expressed', 0.75)}] <- TRUE
    
    cat("Group", group_level, ":", sum(proportion_expressed >= {params.get('min_prop_expressed', 0.75)}), "genes pass filter\\n")
  }}
  
  cat("Original gene count:", nrow(counts), "\\n")
  cat("Genes passing relevance filter:", sum(genes_to_keep_mask), "\\n")
  
  return(counts[genes_to_keep_mask, ])
}}

generate_pca_plot <- function(vst_data, metadata, title) {{
  pca <- prcomp(t(assay(vst_data)))
  pca_data <- as.data.frame(pca$x)
  pca_data <- cbind(pca_data, metadata)
  percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2))
  
  cat("PC variance explained: PC1 =", percent_var[1], "%, PC2 =", percent_var[2], "%\\n")
  
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = rownames(pca_data))) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(size = 3, max.overlaps = 15, box.padding = 0.5) +
    labs(
      title = title,
      x = paste0("PC1: ", percent_var[1], "% variance"),
      y = paste0("PC2: ", percent_var[2], "% variance"),
      color = "Group"
    ) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, size = 16), legend.position = "bottom")
  
  batch_cols <- c("batch", "Batch", "BATCH", "batch_id", "Batch_ID")
  batch_col_found <- NULL
  
  for (col in batch_cols) {{
    if (col %in% colnames(metadata)) {{
      batch_col_found <- col
      break
    }}
  }}
  
  if (!is.null(batch_col_found)) {{
    p <- p + aes(shape = as.factor(metadata[[batch_col_found]])) + labs(shape = "Batch")
  }}
  
  return(p)
}}

identify_pca_outliers <- function(vst_data, metadata, pc_dims = 1:2, method = "combined") {{
  cat("\\n=== IMPROVED PCA OUTLIER DETECTION BY TREATMENT GROUP ===\\n")
  
  pca <- prcomp(t(assay(vst_data)))
  pca_data <- as.data.frame(pca$x[, pc_dims])
  
  outliers <- c()
  outlier_info <- list()
  
  if ("condition" %in% colnames(metadata)) {{
    treatment_groups <- unique(metadata$condition)
    cat("Using 'condition' column for grouping:", paste(treatment_groups, collapse = ", "), "\\n")
  }} else if ("group" %in% colnames(metadata)) {{
    metadata$treatment <- ifelse(metadata$group == "Control", "Control", "Treated")
    treatment_groups <- unique(metadata$treatment)
    cat("Created treatment groups from 'group' column:", paste(treatment_groups, collapse = ", "), "\\n")
  }} else {{
    stop("No suitable grouping column found for treatment-based outlier detection")
  }}
  
  for (treatment_level in treatment_groups) {{
    if ("condition" %in% colnames(metadata)) {{
      group_samples <- rownames(metadata[metadata$condition == treatment_level, ])
    }} else {{
      group_samples <- rownames(metadata[metadata$treatment == treatment_level, ])
    }}
    n_samples <- length(group_samples)
    
    cat("\\nAnalyzing treatment group:", treatment_level, "(n =", n_samples, ")\\n")
    
    if (n_samples < 3) {{
      cat("  Skipping - insufficient samples for outlier detection\\n")
      next
    }}
    
    group_pca_data <- pca_data[group_samples, , drop = FALSE]
    centroid <- colMeans(group_pca_data)
    distances <- sqrt(rowSums(sweep(group_pca_data, 2, centroid, '-')^2))
    
    # Method 1: IQR method
    q1 <- quantile(distances, 0.25)
    q3 <- quantile(distances, 0.75)
    iqr <- q3 - q1
    iqr_threshold <- q3 + (1.2 * iqr)
    
    # Method 2: Z-score method
    mean_dist <- mean(distances)
    sd_dist <- sd(distances)
    z_threshold <- 2.0
    z_outlier_threshold <- mean_dist + (z_threshold * sd_dist)
    
    # Method 3: Modified Z-score using MAD
    median_dist <- median(distances)
    mad_dist <- mad(distances)
    modified_z_threshold <- 2.5
    mad_outlier_threshold <- median_dist + (modified_z_threshold * mad_dist)
    
    iqr_outliers <- names(distances[distances > iqr_threshold])
    z_outliers <- names(distances[distances > z_outlier_threshold])
    mad_outliers <- names(distances[distances > mad_outlier_threshold])
    
    if (method == "combined") {{
      all_potential <- unique(c(iqr_outliers, z_outliers, mad_outliers))
      group_outliers <- c()
      for (sample in all_potential) {{
        detection_count <- sum(c(
          sample %in% iqr_outliers,
          sample %in% z_outliers,
          sample %in% mad_outliers
        ))
        if (detection_count >= 2) {{
          group_outliers <- c(group_outliers, sample)
        }}
      }}
    }}
    
    if (length(group_outliers) > 0) {{
      for (outlier in group_outliers) {{
        outlier_specific_group <- metadata[outlier, "group"]
        cat("      ", outlier, "(", outlier_specific_group, ", distance:", round(distances[outlier], 3), ")\\n")
        outlier_info[[outlier]] <- list(
          treatment = treatment_level,
          specific_group = outlier_specific_group,
          distance = distances[outlier],
          methods_detected = c(
            if(outlier %in% iqr_outliers) "IQR",
            if(outlier %in% z_outliers) "Z-score", 
            if(outlier %in% mad_outliers) "MAD"
          )
        )
      }}
      outliers <- c(outliers, group_outliers)
    }}
  }}
  
  cat("\\nTotal outliers identified:", length(outliers), "\\n")
  
  return(list(outliers = outliers, outlier_info = outlier_info))
}}

#================================================================================
# DEG ANALYSIS FUNCTIONS
#================================================================================

extract_rodaf_results <- function(dds, comparison_name, output_prefix, fdr_threshold_1 = {params.get('fdr_strict', 0.01)}, fdr_threshold_2 = {params.get('fdr_relaxed', 0.05)}) {{
  cat("\\n--- Extracting results for:", comparison_name, "---\\n")
  
  res <- results(dds, name = comparison_name)
  res_df <- as.data.frame(res)
  res_df$ensembl_id <- rownames(res_df)
  res_df <- res_df[, c("ensembl_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
  res_df <- res_df[!is.na(res_df$padj), ]
  res_df <- res_df[order(res_df$padj), ]
  
  sig_001 <- res_df[res_df$padj < fdr_threshold_1, ]
  sig_005 <- res_df[res_df$padj < fdr_threshold_2, ]
  
  cat("Total genes with valid p-values:", nrow(res_df), "\\n")
  cat("Significant at FDR <", fdr_threshold_1, ":", nrow(sig_001), "genes\\n")
  cat("Significant at FDR <", fdr_threshold_2, ":", nrow(sig_005), "genes\\n")
  
  write.csv(res_df, paste0(output_prefix, "_all_results.csv"), row.names = FALSE)
  write.csv(sig_001, paste0(output_prefix, "_FDR001.csv"), row.names = FALSE)
  write.csv(sig_005, paste0(output_prefix, "_FDR005.csv"), row.names = FALSE)
  
  cat("Results saved with prefix:", output_prefix, "\\n")
  
  return(list(all = res_df, sig_001 = sig_001, sig_005 = sig_005))
}}

#================================================================================
# MAIN PIPELINE FUNCTION
#================================================================================

run_complete_rodaf_pipeline <- function(counts_path, metadata_path, output_dir, 
                                        treatment_to_analyze, control_to_use,
                                        read_threshold = {params.get('read_threshold', 1000000)}, 
                                        outlier_method = "combined") {{
  
  cat("\\n=== STARTING COMPLETE R-ODAF PIPELINE ===\\n")
  
  if (!dir.exists(output_dir)) {{ 
    dir.create(output_dir, recursive = TRUE) 
  }}
  
  cat("\\n--- PART I: QC PIPELINE ---\\n")
  
  # Load count data
  count_data <- read_featurecounts_output(counts_path)
  
  # Load metadata
  cat("\\n=== LOADING METADATA ===\\n")
  if (grepl("\\\\.xlsx$", metadata_path, ignore.case = TRUE) || grepl("\\\\.xls$", metadata_path, ignore.case = TRUE)) {{
    metadata_full <- as.data.frame(readxl::read_excel(metadata_path))
  }} else if (grepl("\\\\.csv$", metadata_path, ignore.case = TRUE)) {{
    metadata_full <- read.csv(metadata_path, stringsAsFactors = FALSE)
  }} else {{
    tryCatch({{
      metadata_full <- read.csv(metadata_path, stringsAsFactors = FALSE)
    }}, error = function(e) {{
      metadata_full <<- as.data.frame(readxl::read_excel(metadata_path))
    }})
  }}
  
  # Validate and align
  aligned_data <- validate_metadata_with_batch(count_data, metadata_full)
  
  # Subset data
  metadata_subset <- aligned_data$metadata %>%
    filter(chemical %in% c(treatment_to_analyze, control_to_use)) %>%
    mutate(
      group = if_else(
        chemical == control_to_use,
        "Control",
        paste0("Treated_", gsub("\\\\.", "p", as.character(dose)), "uM")
      ),
      condition = if_else(chemical == control_to_use, "Control", "Treated")
    )
  
  common_samples <- intersect(colnames(aligned_data$counts), rownames(metadata_subset))
  counts_subset <- aligned_data$counts[, common_samples]
  metadata_subset <- metadata_subset[common_samples, ]
  
  # Apply QC steps
  read_filter_results <- remove_low_read_samples(counts_subset, metadata_subset, threshold = read_threshold)
  filtered_counts <- filter_genes_by_cpm(read_filter_results$counts, read_filter_results$metadata)
  batch_results <- apply_combat_batch_correction(filtered_counts, read_filter_results$metadata)
  
  # PCA and outlier detection
  dds_qc <- DESeqDataSetFromMatrix(countData = batch_results$corrected_counts, 
                                   colData = read_filter_results$metadata, 
                                   design = ~ 1)
  vst_qc <- vst(dds_qc, blind = TRUE)
  
  pca_plot <- generate_pca_plot(vst_qc, read_filter_results$metadata, "QC: Final PCA")
  ggsave(file.path(output_dir, "Final_PCA_Plot.png"), pca_plot, width = 12, height = 8, dpi = 300)
  
  outlier_results <- identify_pca_outliers(vst_qc, read_filter_results$metadata, method = outlier_method)
  
  # Final clean data
  remaining_samples <- setdiff(rownames(read_filter_results$metadata), outlier_results$outliers)
  final_counts <- batch_results$corrected_counts[, remaining_samples]
  final_metadata <- read_filter_results$metadata[remaining_samples, ]
  
  cat("\\nQC PIPELINE COMPLETE\\n")
  cat("Final samples:", ncol(final_counts), "\\n")
  cat("Final genes:", nrow(final_counts), "\\n")
  
  cat("\\n--- PART II: DIFFERENTIAL EXPRESSION ANALYSIS ---\\n")
  
  # Check for dose-response analysis
  unique_groups <- unique(final_metadata$group)
  treatment_groups <- unique_groups[unique_groups != "Control"]
  
  dose_results <- NULL
  if (length(treatment_groups) > 1) {{
    cat("\\n### R-ODAF ANALYSIS 1: DOSE-RESPONSE (EACH DOSE vs. CONTROL)\\n")
    
    final_metadata$group <- as.factor(final_metadata$group)
    final_metadata$group <- relevel(final_metadata$group, ref = "Control")
    
    dds_dose <- DESeqDataSetFromMatrix(countData = final_counts, colData = final_metadata, design = ~ group)
    dds_dose <- DESeq(dds_dose)
    
    comparisons <- resultsNames(dds_dose)[grepl("^group_Treated", resultsNames(dds_dose))]
    dose_results <- list()
    
    for (comp_name in comparisons) {{
      clean_comp_name <- gsub("_vs_", "-vs-", sub("group_", "", comp_name))
      output_prefix <- file.path(output_dir, paste0("RODAF_DGE_DOSE_", clean_comp_name))
      dose_results[[comp_name]] <- extract_rodaf_results(dds_dose, comp_name, output_prefix)
    }}
  }}
  
  # OVERALL ANALYSIS
  cat("\\n### R-ODAF ANALYSIS 2: OVERALL (ALL TREATED vs. CONTROL)\\n")
  
  final_metadata$condition <- as.factor(final_metadata$condition)
  final_metadata$condition <- relevel(final_metadata$condition, ref = "Control")
  
  dds_overall <- DESeqDataSetFromMatrix(countData = final_counts, colData = final_metadata, design = ~ condition)
  dds_overall <- DESeq(dds_overall)
  
  comp_name_overall <- "condition_Treated_vs_Control"
  output_prefix_overall <- file.path(output_dir, "RODAF_DGE_OVERALL_Treated-vs-Control")
  overall_results <- extract_rodaf_results(dds_overall, comp_name_overall, output_prefix_overall)
  
  # GENE SYMBOL MAPPING
  cat("\\n--- Adding Gene Symbol Mapping and Custom Filtering ---\\n")
  
  tryCatch({{
    gene_symbols <- mapIds(org.Hs.eg.db, 
                           keys = overall_results$all$ensembl_id, 
                           column = "SYMBOL", 
                           keytype = "ENSEMBL", 
                           multiVals = "first")
    
    overall_results_with_symbols <- overall_results$all
    overall_results_with_symbols$gene_symbol <- gene_symbols[overall_results_with_symbols$ensembl_id]
    
    overall_results_with_symbols <- overall_results_with_symbols[, c(
      "ensembl_id", "gene_symbol", "baseMean", 
      "log2FoldChange", "lfcSE", "pvalue", "padj"
    )]
    
    overall_results_with_symbols <- overall_results_with_symbols[order(overall_results_with_symbols$padj), ]
    
    cat("Successfully mapped", sum(!is.na(gene_symbols)), "out of", length(gene_symbols), "genes to symbols\\n")
    
  }}, error = function(e) {{
    cat("Error in gene symbol mapping:", e$message, "\\n")
    overall_results_with_symbols <- overall_results$all
    overall_results_with_symbols$gene_symbol <- NA
  }})
  
  # Custom filtering
  custom_filtered_degs <- overall_results_with_symbols %>%
    filter(!is.na(padj) & padj < {params.get('fdr_relaxed', 0.05)} & abs(log2FoldChange) > {params.get('log2fc_threshold', 1)}) %>%
    arrange(padj)
  
  upregulated_genes <- custom_filtered_degs %>%
    filter(log2FoldChange > {params.get('log2fc_threshold', 1)}) %>%
    arrange(desc(log2FoldChange))
  
  downregulated_genes <- custom_filtered_degs %>%
    filter(log2FoldChange < -{params.get('log2fc_threshold', 1)}) %>%
    arrange(log2FoldChange)
  
  # Save results
  write.csv(overall_results_with_symbols, file.path(output_dir, "Overall_Results_With_Symbols.csv"), row.names = FALSE)
  write.csv(custom_filtered_degs, file.path(output_dir, "Custom_Filtered_DEGs.csv"), row.names = FALSE)
  write.csv(upregulated_genes, file.path(output_dir, "Upregulated_Genes.csv"), row.names = FALSE)
  write.csv(downregulated_genes, file.path(output_dir, "Downregulated_Genes.csv"), row.names = FALSE)
  
  # Generate volcano plot
  create_volcano_plot <- function(results_df, title, output_path) {{
    volcano_data <- results_df %>%
      mutate(
        log10_padj = -log10(padj),
        significance = case_when(
          is.na(padj) ~ "Not significant",
          padj < {params.get('fdr_relaxed', 0.05)} & log2FoldChange > {params.get('log2fc_threshold', 1)} ~ "Upregulated",
          padj < {params.get('fdr_relaxed', 0.05)} & log2FoldChange < -{params.get('log2fc_threshold', 1)} ~ "Downregulated",
          padj < {params.get('fdr_relaxed', 0.05)} ~ "Significant (small effect)",
          TRUE ~ "Not significant"
        )
      )
    
    p <- ggplot(volcano_data, aes(x = log2FoldChange, y = log10_padj, color = significance)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_manual(values = c(
        "Upregulated" = "red",
        "Downregulated" = "blue",
        "Significant (small effect)" = "orange",
        "Not significant" = "grey"
      )) +
      geom_vline(xintercept = c(-{params.get('log2fc_threshold', 1)}, {params.get('log2fc_threshold', 1)}), linetype = "dashed", alpha = 0.7) +
      geom_hline(yintercept = -log10({params.get('fdr_relaxed', 0.05)}), linetype = "dashed", alpha = 0.7) +
      labs(
        title = title,
        x = "log2 Fold Change",
        y = "-log10(Adjusted P-value)",
        color = "Significance"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "bottom"
      )
    
    if (nrow(volcano_data[volcano_data$significance %in% c("Upregulated", "Downregulated"), ]) > 0) {{
      top_genes <- volcano_data %>%
        filter(significance %in% c("Upregulated", "Downregulated")) %>%
        arrange(padj) %>%
        head(10)
      
      p <- p + geom_text_repel(
        data = top_genes,
        aes(label = ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, ensembl_id)),
        size = 3,
        max.overlaps = 10,
        box.padding = 0.5
      )
    }}
    
    ggsave(output_path, p, width = 12, height = 10, dpi = 300)
    return(p)
  }}
  
  volcano_plot <- create_volcano_plot(
    overall_results_with_symbols,
    paste("Volcano Plot:", treatment_to_analyze, "vs", control_to_use),
    file.path(output_dir, "Volcano_Plot_Overall.png")
  )
  cat("Volcano plot saved to:", file.path(output_dir, "Volcano_Plot_Overall.png"), "\\n")
  
  # Generate summary report
  generate_summary_report <- function() {{
    report_file <- file.path(output_dir, "R-ODAF_Analysis_Summary_Report.txt")
    
    sink(report_file)
    
    cat("================================================================================\\n")
    cat("R-ODAF COMPREHENSIVE ANALYSIS REPORT\\n")
    cat("================================================================================\\n\\n")
    
    cat("Analysis Parameters:\\n")
    cat("- Treatment Analyzed:", treatment_to_analyze, "\\n")
    cat("- Control Used:", control_to_use, "\\n")
    cat("- Read Count Threshold:", format(read_threshold, scientific = FALSE), "reads\\n")
    cat("- Outlier Detection Method:", outlier_method, "\\n")
    cat("- Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")
    
    cat("QC PIPELINE RESULTS:\\n")
    cat("- Final samples retained:", ncol(final_counts), "\\n")
    cat("- Final genes analyzed:", nrow(final_counts), "\\n")
    cat("- Samples removed (low reads):", length(read_filter_results$removed_samples), "\\n")
    cat("- Samples removed (outliers):", length(outlier_results$outliers), "\\n")
    cat("- Batch correction applied:", batch_results$batch_corrected, "\\n\\n")
    
    cat("DIFFERENTIAL EXPRESSION ANALYSIS RESULTS:\\n\\n")
    
    if (!is.null(dose_results)) {{
      cat("Dose-Response Analysis:\\n")
      for (comparison in names(dose_results)) {{
        clean_name <- gsub("_vs_", " vs ", sub("group_", "", comparison))
        cat("  ", clean_name, ":\\n")
        cat("    FDR < {params.get('fdr_strict', 0.01)}:", nrow(dose_results[[comparison]]$sig_001), "genes\\n")
        cat("    FDR < {params.get('fdr_relaxed', 0.05)}:", nrow(dose_results[[comparison]]$sig_005), "genes\\n")
      }}
      cat("\\n")
    }}
    
    cat("Overall Analysis (All Treated vs Control):\\n")
    cat("  FDR < {params.get('fdr_strict', 0.01)}:", nrow(overall_results$sig_001), "genes\\n")
    cat("  FDR < {params.get('fdr_relaxed', 0.05)}:", nrow(overall_results$sig_005), "genes\\n")
    cat("  Custom filtered (padj < {params.get('fdr_relaxed', 0.05)} & |log2FC| > {params.get('log2fc_threshold', 1)}):", nrow(custom_filtered_degs), "genes\\n")
    cat("    - Upregulated (log2FC > {params.get('log2fc_threshold', 1)}):", nrow(upregulated_genes), "genes\\n")
    cat("    - Downregulated (log2FC < -{params.get('log2fc_threshold', 1)}):", nrow(downregulated_genes), "genes\\n\\n")
    
    if (nrow(custom_filtered_degs) > 0) {{
      cat("TOP 20 DIFFERENTIALLY EXPRESSED GENES (Custom Filter):\\n")
      cat("Rank\\tGene Symbol\\tEnsembl ID\\tlog2FC\\tAdjusted P-value\\n")
      
      top_genes <- head(custom_filtered_degs, 20)
      for (i in 1:nrow(top_genes)) {{
        gene_name <- ifelse(is.na(top_genes$gene_symbol[i]) || top_genes$gene_symbol[i] == "", 
                            "N/A", top_genes$gene_symbol[i])
        cat(sprintf("%d\\t%s\\t%s\\t%.3f\\t%.2e\\n", 
                    i, gene_name, top_genes$ensembl_id[i], 
                    top_genes$log2FoldChange[i], top_genes$padj[i]))
      }}
    }}
    
    cat("\\n================================================================================\\n")
    cat("FILES GENERATED:\\n")
    cat("- Final_PCA_Plot.png: PCA visualization after QC\\n")
    cat("- Volcano_Plot_Overall.png: Volcano plot for overall comparison\\n")
    cat("- Overall_Results_With_Symbols.csv: All genes with symbols\\n")
    cat("- Custom_Filtered_DEGs.csv: Significant DEGs (padj<{params.get('fdr_relaxed', 0.05)}, |log2FC|>{params.get('log2fc_threshold', 1)})\\n")
    cat("- Upregulated_Genes.csv: Significantly upregulated genes\\n")
    cat("- Downregulated_Genes.csv: Significantly downregulated genes\\n")
    cat("- RODAF_DGE_OVERALL_*: Standard R-ODAF outputs (FDR {params.get('fdr_strict', 0.01)}, {params.get('fdr_relaxed', 0.05)})\\n")
    
    if (!is.null(dose_results)) {{
      cat("- RODAF_DGE_DOSE_*: Dose-response analysis results\\n")
    }}
    
    cat("================================================================================\\n")
    
    sink()
    
    cat("Comprehensive summary report saved to:", report_file, "\\n")
  }}
  
  generate_summary_report()
  
  # Save summary statistics
  summary_stats <- list(
    total_genes = nrow(overall_results_with_symbols),
    sig_strict = nrow(overall_results$sig_001),
    sig_relaxed = nrow(overall_results$sig_005),
    sig_custom = nrow(custom_filtered_degs),
    upregulated = nrow(upregulated_genes),
    downregulated = nrow(downregulated_genes),
    samples_analyzed = ncol(final_counts),
    genes_analyzed = nrow(final_counts),
    treatment = treatment_to_analyze,
    control = control_to_use,
    fdr_strict = {params.get('fdr_strict', 0.01)},
    fdr_relaxed = {params.get('fdr_relaxed', 0.05)},
    log2fc_threshold = {params.get('log2fc_threshold', 1)},
    batch_corrected = batch_results$batch_corrected,
    outliers_removed = length(outlier_results$outliers),
    low_read_samples_removed = length(read_filter_results$removed_samples)
  )
  
  # Save JSON summary for Python to read
  library(jsonlite)
  write_json(summary_stats, file.path(output_dir, "analysis_summary.json"), pretty = TRUE)
  
  cat("\\n=== CUSTOM FILTERING RESULTS (padj < {params.get('fdr_relaxed', 0.05)} & |log2FC| > {params.get('log2fc_threshold', 1)}) ===\\n")
  cat("Total DEGs meeting criteria:", nrow(custom_filtered_degs), "\\n")
  cat("Upregulated genes (log2FC > {params.get('log2fc_threshold', 1)}):", nrow(upregulated_genes), "\\n")
  cat("Downregulated genes (log2FC < -{params.get('log2fc_threshold', 1)}):", nrow(downregulated_genes), "\\n")
  
  if (nrow(custom_filtered_degs) > 0) {{
    cat("\\nTop 10 DEGs (by adjusted p-value):\\n")
    top_degs <- head(custom_filtered_degs, 10)
    for (i in 1:nrow(top_degs)) {{
      gene_name <- ifelse(is.na(top_degs$gene_symbol[i]) || top_degs$gene_symbol[i] == "", 
                          top_degs$ensembl_id[i], 
                          top_degs$gene_symbol[i])
      cat(sprintf("  %d. %s (log2FC: %.2f, padj: %.2e)\\n", 
                  i, gene_name, top_degs$log2FoldChange[i], top_degs$padj[i]))
    }}
  }}
  
  cat("\\n=== R-ODAF PIPELINE COMPLETION SUMMARY ===\\n")
  cat("Output directory:", output_dir, "\\n")
  cat("\\nCOMPLETE R-ODAF PIPELINE FINISHED!\\n")
  
  return(list(success = TRUE))
}}

#================================================================================
# STEP 3: SET YOUR FILE PATHS AND PARAMETERS
#================================================================================

cat("\\n=== STEP 3: CONFIGURING ANALYSIS PARAMETERS ===\\n")

# File paths from Python parameters
counts_path <- "{params['counts_file']}"
metadata_path <- "{params['metadata_file']}"
output_dir <- "{params['output_dir']}"

# Analysis parameters
TREATMENT_TO_ANALYZE <- "{params['treatment']}"
CONTROL_TO_USE <- "{params['control']}"
READ_THRESHOLD <- {params.get('read_threshold', 1000000)}
OUTLIER_METHOD <- "combined"

cat("Analysis parameters set:\\n")
cat("- Treatment:", TREATMENT_TO_ANALYZE, "\\n")
cat("- Control:", CONTROL_TO_USE, "\\n")
cat("- Read threshold:", format(READ_THRESHOLD, scientific = FALSE), "\\n")
cat("- Outlier method:", OUTLIER_METHOD, "\\n")
cat("- Output directory:", output_dir, "\\n")

#================================================================================
# STEP 4: RUN THE COMPLETE INTEGRATED PIPELINE
#================================================================================

cat("\\n=== STEP 4: RUNNING COMPLETE INTEGRATED R-ODAF PIPELINE ===\\n")

tryCatch({{
  pipeline_results <- run_complete_rodaf_pipeline(
    counts_path = counts_path,
    metadata_path = metadata_path,
    output_dir = output_dir,
    treatment_to_analyze = TREATMENT_TO_ANALYZE,
    control_to_use = CONTROL_TO_USE,
    read_threshold = READ_THRESHOLD,
    outlier_method = OUTLIER_METHOD
  )
  
  if (pipeline_results$success) {{
    cat("\\n=== ANALYSIS COMPLETE ===\\n")
    cat("Results saved to:", output_dir, "\\n\\n")
    
    # Read and display summary
    summary_file <- file.path(output_dir, "analysis_summary.json")
    if (file.exists(summary_file)) {{
      summary <- jsonlite::fromJSON(summary_file)
      cat("Summary:\\n")
      cat("- Total genes analyzed:", summary$total_genes, "\\n")
      cat("- Significant (FDR < {params.get('fdr_relaxed', 0.05)}):", summary$sig_relaxed, "\\n")
      cat("- Upregulated:", summary$upregulated, "\\n")
      cat("- Downregulated:", summary$downregulated, "\\n")
      cat("\\nCheck the output directory for detailed results and plots.\\n")
    }}
  }}
  
}}, error = function(e) {{
  cat("\\nERROR in analysis pipeline:\\n")
  cat(toString(e), "\\n")
  quit(status = 1)
}})
"""
    
    return r_script
    
def generate_dromics_r_script(params):
    """Generate DRomics R script with integrated BMD analysis"""
    
    r_script = f"""#!/usr/bin/env Rscript
#================================================================================
# INTEGRATED DROMICS DOSE-RESPONSE ANALYSIS WITH BMD CALCULATION
# Following R-ODAF QC methodology with BMD analysis
# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
#================================================================================

# Clear workspace
rm(list = ls())
cat("=== INTEGRATED DROMICS PIPELINE WITH R-ODAF QC AND BMD ANALYSIS ===\\n")

#================================================================================
# STEP 1: LOAD REQUIRED LIBRARIES
#================================================================================

cat("\\n=== STEP 1: LOADING LIBRARIES ===\\n")

suppressPackageStartupMessages({{
  library(dplyr)          # Data manipulation
  library(readxl)         # Excel file reading
  library(DESeq2)         # For initial QC
  library(edgeR)          # For CPM filtering
  library(sva)            # Batch correction
  library(ggplot2)        # Plotting
  library(ggrepel)        # Plot labels
  library(stringr)        # String manipulation
  library(org.Hs.eg.db)   # Gene annotation
  library(DRomics)        # Dose-response analysis
  library(jsonlite)       # JSON output
  library(parallel)       # For parallel processing in BMD bootstrap
}})

cat("Libraries loaded successfully!\\n")

#================================================================================
# STEP 2: LOADING QC FUNCTIONS FROM DEG PIPELINE
#================================================================================

cat("\\n=== STEP 2: LOADING INTEGRATED QC FUNCTIONS ===\\n")

# Source ComBat-seq functions if available
tryCatch({{
  source("/home/saurav/Downloads/ComBat_seq.R")
  source("/home/saurav/Downloads/helper_seq.R")
  cat("ComBat-seq functions loaded successfully\\n")
}}, error = function(e) {{
  cat("Warning: ComBat-seq functions not found. Will use sva::ComBat instead.\\n")
}})

# Read featureCounts output function (same as DEG)
read_featurecounts_output <- function(file_path) {{
  cat("\\n=== READING FEATURECOUNTS OUTPUT ===\\n")
  
  first_lines <- readLines(file_path, n = 5)
  
  if (grepl("^#", first_lines[1])) {{
    skip_lines <- 1
  }} else {{
    skip_lines <- 0
  }}
  
  tryCatch({{
    raw_data <- read.table(file_path, header = TRUE, sep = "\\t", 
                           comment.char = "", stringsAsFactors = FALSE, 
                           skip = skip_lines, check.names = FALSE)
  }}, error = function(e1) {{
    raw_data <<- read.table(file_path, header = TRUE, sep = "\\t", 
                            comment.char = "", stringsAsFactors = FALSE, 
                            skip = skip_lines, fill = TRUE, check.names = FALSE)
  }})
  
  cat("Successfully loaded data with dimensions:", nrow(raw_data), "x", ncol(raw_data), "\\n")
  
  if ("Length" %in% colnames(raw_data)) {{
    length_col_idx <- which(colnames(raw_data) == "Length")
  }} else {{
    length_col_idx <- 6
  }}
  
  count_cols <- (length_col_idx + 1):ncol(raw_data)
  count_matrix <- raw_data[, count_cols, drop = FALSE]
  
  if ("Geneid" %in% colnames(raw_data)) {{
    rownames(count_matrix) <- raw_data$Geneid
  }} else {{
    rownames(count_matrix) <- raw_data[, 1]
  }}
  
  sample_names <- colnames(count_matrix)
  clean_names <- str_extract(sample_names, "SRR[0-9]+")
  
  if (any(is.na(clean_names))) {{
    alt_names <- str_extract(sample_names, "[A-Za-z]+[0-9]+")
    clean_names[is.na(clean_names)] <- alt_names[is.na(clean_names)]
    
    still_na <- is.na(clean_names)
    if (any(still_na)) {{
      clean_names[still_na] <- basename(sample_names[still_na])
      clean_names[still_na] <- gsub("\\\\.[^.]*$", "", clean_names[still_na])
    }}
  }}
  
  colnames(count_matrix) <- clean_names
  count_matrix <- as.matrix(count_matrix)
  mode(count_matrix) <- "numeric"
  
  if (any(is.na(count_matrix))) {{
    cat("Warning: Found NA values in count matrix. Converting to 0...\\n")
    count_matrix[is.na(count_matrix)] <- 0
  }}
  
  cat("Successfully parsed featureCounts output:\\n")
  cat("- Genes:", nrow(count_matrix), "\\n")
  cat("- Samples:", ncol(count_matrix), "\\n")
  
  return(list(counts = count_matrix, original_names = sample_names))
}}

# Validate metadata with batch
validate_metadata_with_batch <- function(count_data, metadata_df) {{
  cat("\\n=== VALIDATING METADATA WITH BATCH INFO ===\\n")
  
  if (!"Run" %in% colnames(metadata_df)) {{
    stop("Error: 'Run' column not found in metadata.")
  }}
  
  rownames(metadata_df) <- metadata_df$Run
  
  count_samples <- colnames(count_data$counts)
  metadata_samples <- rownames(metadata_df)
  common_samples <- intersect(count_samples, metadata_samples)
  
  cat("Count data samples:", length(count_samples), "\\n")
  cat("Metadata samples:", length(metadata_samples), "\\n")
  cat("Common samples:", length(common_samples), "\\n")
  
  if (length(common_samples) == 0) {{
    stop("Error: No common samples found between count data and metadata.")
  }}
  
  batch_cols <- grep("batch|Batch|BATCH", colnames(metadata_df), value = TRUE)
  cat("Potential batch columns found:", paste(batch_cols, collapse = ", "), "\\n")
  
  aligned_counts <- count_data$counts[, common_samples]
  aligned_metadata <- metadata_df[common_samples, ]
  
  return(list(counts = aligned_counts, metadata = aligned_metadata))
}}

# Apply batch correction (same as DEG)
apply_combat_batch_correction <- function(counts, metadata) {{
  cat("\\n=== APPLYING BATCH CORRECTION ===\\n")
  
  batch_col <- NULL
  possible_batch_cols <- c("batch", "Batch", "BATCH", "batch_id", "Batch_ID")
  
  for (col in possible_batch_cols) {{
    if (col %in% colnames(metadata)) {{
      batch_col <- col
      break
    }}
  }}
  
  if (is.null(batch_col)) {{
    cat("No batch column found. Batch correction will be skipped.\\n")
    return(list(corrected_counts = counts, batch_corrected = FALSE))
  }}
  
  common_samples <- intersect(colnames(counts), rownames(metadata))
  counts <- counts[, common_samples]
  metadata <- metadata[common_samples, ]
  
  batch_factor <- as.factor(metadata[[batch_col]])
  
  if (length(unique(batch_factor)) < 2) {{
    cat("Only one batch detected. Batch correction not applicable.\\n")
    return(list(corrected_counts = counts, batch_corrected = FALSE))
  }}
  
  # For dose-response, we need to consider dose groups
  if ("dose" %in% colnames(metadata)) {{
    dose_factor <- as.factor(metadata$dose)
  }} else {{
    dose_factor <- as.factor(metadata$group)
  }}
  
  tryCatch({{
    if (!exists("ComBat_seq")) {{
      library(sva)
      cpm_data <- edgeR::cpm(counts, log = TRUE, prior.count = 1)
      corrected_log_data <- ComBat(
        dat = cpm_data,
        batch = batch_factor,
        mod = model.matrix(~dose_factor),
        par.prior = TRUE
      )
      
      corrected_counts_matrix <- 2^corrected_log_data - 1
      corrected_counts_matrix[corrected_counts_matrix < 0] <- 0
      corrected_counts_matrix <- round(corrected_counts_matrix)
      
      cat("ComBat (sva) batch correction completed successfully.\\n")
      return(list(corrected_counts = corrected_counts_matrix, batch_corrected = TRUE))
    }} else {{
      corrected_counts_matrix <- ComBat_seq(
        counts = counts,
        batch = batch_factor,
        group = dose_factor,
        shrink = FALSE
      )
      
      cat("ComBat-seq batch correction completed successfully.\\n")
      return(list(corrected_counts = corrected_counts_matrix, batch_corrected = TRUE))
    }}
  }}, error = function(e) {{
    cat("ERROR: Batch correction failed:", e$message, "\\n")
    cat("Returning original counts.\\n")
    return(list(corrected_counts = counts, batch_corrected = FALSE))
  }})
}}

# Remove low read samples (same as DEG)
remove_low_read_samples <- function(counts, metadata, threshold = {params.get('read_threshold', 1000000)}) {{
  cat("\\n=== R-ODAF SAMPLE FILTERING (Read Count) ===\\n")
  library_size <- colSums(counts)
  
  low_read_samples <- names(library_size[library_size < threshold])
  
  if (length(low_read_samples) > 0) {{
    cat("Removing samples with <", format(threshold, scientific = FALSE), "reads:", length(low_read_samples), "\\n")
    for (sample in low_read_samples) {{
      reads <- format(library_size[sample], big.mark = ",")
      cat("  -", sample, "(", reads, "reads )\\n")
    }}
    counts <- counts[, !(colnames(counts) %in% low_read_samples)]
    metadata <- metadata[!(rownames(metadata) %in% low_read_samples), ]
  }} else {{
    cat("All samples pass read count filter\\n")
  }}
  
  cat("Samples retained:", ncol(counts), "\\n")
  return(list(counts = counts, metadata = metadata, removed_samples = low_read_samples))
}}

# Filter genes by CPM (adapted for dose groups)
filter_genes_by_cpm <- function(counts, metadata) {{
  cat("\\n=== R-ODAF RELEVANCE FILTER (DOSE-AWARE) ===\\n")
  cpm_matrix <- edgeR::cpm(counts)
  genes_to_keep_mask <- rep(FALSE, nrow(cpm_matrix))
  
  # Group by dose for filtering
  if ("dose" %in% colnames(metadata)) {{
    metadata$dose_group <- as.factor(metadata$dose)
  }} else {{
    metadata$dose_group <- as.factor(metadata$group)
  }}
  
  for (dose_level in levels(metadata$dose_group)) {{
    dose_samples <- rownames(metadata[metadata$dose_group == dose_level, ])
    dose_samples_in_matrix <- intersect(dose_samples, colnames(cpm_matrix))
    if(length(dose_samples_in_matrix) == 0) next
    
    replicates_in_dose_cpm <- cpm_matrix[, dose_samples_in_matrix, drop = FALSE]
    proportion_expressed <- rowSums(replicates_in_dose_cpm >= {params.get('cpm_threshold', 1)}) / length(dose_samples_in_matrix)
    genes_to_keep_mask[proportion_expressed >= {params.get('min_prop_expressed', 0.75)}] <- TRUE
    
    cat("Dose", dose_level, ":", sum(proportion_expressed >= {params.get('min_prop_expressed', 0.75)}), "genes pass filter\\n")
  }}
  
  cat("Original gene count:", nrow(counts), "\\n")
  cat("Genes passing relevance filter:", sum(genes_to_keep_mask), "\\n")
  
  return(counts[genes_to_keep_mask, ])
}}

# PCA outlier detection (adapted for dose-response)
identify_pca_outliers <- function(vst_data, metadata, pc_dims = 1:2, method = "combined") {{
  cat("\\n=== IMPROVED PCA OUTLIER DETECTION BY DOSE GROUP ===\\n")
  
  pca <- prcomp(t(assay(vst_data)))
  pca_data <- as.data.frame(pca$x[, pc_dims])
  
  outliers <- c()
  outlier_info <- list()
  
  # Group by dose for outlier detection
  if ("dose" %in% colnames(metadata)) {{
    dose_groups <- unique(metadata$dose)
    cat("Using dose levels for grouping:", paste(dose_groups, collapse = ", "), "\\n")
    grouping_col <- "dose"
  }} else {{
    dose_groups <- unique(metadata$group)
    cat("Using group column for grouping:", paste(dose_groups, collapse = ", "), "\\n")
    grouping_col <- "group"
  }}
  
  for (dose_level in dose_groups) {{
    group_samples <- rownames(metadata[metadata[[grouping_col]] == dose_level, ])
    n_samples <- length(group_samples)
    
    cat("\\nAnalyzing dose group:", dose_level, "(n =", n_samples, ")\\n")
    
    if (n_samples < 3) {{
      cat("  Skipping - insufficient samples for outlier detection\\n")
      next
    }}
    
    group_pca_data <- pca_data[group_samples, , drop = FALSE]
    centroid <- colMeans(group_pca_data)
    distances <- sqrt(rowSums(sweep(group_pca_data, 2, centroid, '-')^2))
    
    # Method 1: IQR method
    q1 <- quantile(distances, 0.25)
    q3 <- quantile(distances, 0.75)
    iqr <- q3 - q1
    iqr_threshold <- q3 + (1.5 * iqr)
    
    # Method 2: Z-score method
    mean_dist <- mean(distances)
    sd_dist <- sd(distances)
    z_threshold <- 2.5
    z_outlier_threshold <- mean_dist + (z_threshold * sd_dist)
    
    # Method 3: Modified Z-score using MAD
    median_dist <- median(distances)
    mad_dist <- mad(distances)
    modified_z_threshold <- 3.0
    mad_outlier_threshold <- median_dist + (modified_z_threshold * mad_dist)
    
    iqr_outliers <- names(distances[distances > iqr_threshold])
    z_outliers <- names(distances[distances > z_outlier_threshold])
    mad_outliers <- names(distances[distances > mad_outlier_threshold])
    
    if (method == "combined") {{
      all_potential <- unique(c(iqr_outliers, z_outliers, mad_outliers))
      group_outliers <- c()
      for (sample in all_potential) {{
        detection_count <- sum(c(
          sample %in% iqr_outliers,
          sample %in% z_outliers,
          sample %in% mad_outliers
        ))
        if (detection_count >= 2) {{
          group_outliers <- c(group_outliers, sample)
        }}
      }}
    }}
    
    if (length(group_outliers) > 0) {{
      for (outlier in group_outliers) {{
        cat("      ", outlier, "(distance:", round(distances[outlier], 3), ")\\n")
        outlier_info[[outlier]] <- list(
          dose = dose_level,
          distance = distances[outlier],
          methods_detected = c(
            if(outlier %in% iqr_outliers) "IQR",
            if(outlier %in% z_outliers) "Z-score", 
            if(outlier %in% mad_outliers) "MAD"
          )
        )
      }}
      outliers <- c(outliers, group_outliers)
    }}
  }}
  
  cat("\\nTotal outliers identified:", length(outliers), "\\n")
  
  return(list(outliers = outliers, outlier_info = outlier_info))
}}

#================================================================================
# BMD ANALYSIS HELPER FUNCTIONS
#================================================================================

# Function to calculate BMDs with multiple parameter sets
calculate_multiple_bmds <- function(fit_object) {{
  cat("\\n=== CALCULATING BMDs WITH MULTIPLE PARAMETERS ===\\n")
  
  # Robustly handle parameters from Python
  z_vals_str <- "{params.get('bmd_z_values', '1.0')}"
  x_vals_str <- "{params.get('bmd_x_values', '10')}"
  
  z_values <- as.numeric(unlist(strsplit(z_vals_str, ",")))
  if (length(z_values) == 0 || all(is.na(z_values))) {{
    cat("Warning: Invalid z_values provided. Using default: 1.0\\n")
    z_values <- c(1.0)
  }}
  
  x_values <- as.numeric(unlist(strsplit(x_vals_str, ",")))
  if (length(x_values) == 0 || all(is.na(x_values))) {{
    cat("Warning: Invalid x_values provided. Using default: 10\\n")
    x_values <- c(10)
  }}

  cat("Using z-values:", paste(z_values, collapse=", "), "\\n")
  cat("Using x-values:", paste(x_values, collapse=", "), "\\n")

  bmd_results_list <- list()
  
  for (z in z_values) {{
    for (x in x_values) {{
      param_name <- paste0("z", z, "_x", x)
      cat(paste("Calculating BMD for z =", z, ", x =", x, "...\\n"))
      
      tryCatch({{
        bmd_result <- bmdcalc(f = fit_object, z = z, x = x)
        bmd_results_list[[param_name]] <- bmd_result
        
        # Count valid BMDs
        valid_count <- sum(!is.na(bmd_result$res$BMD.zSD) & bmd_result$res$BMD.zSD > 0)
        cat(paste("  Valid BMDs for z =", z, ":", valid_count, "\\n"))
        
        if (valid_count > 0) {{
          median_bmd <- median(bmd_result$res$BMD.zSD, na.rm = TRUE)
          cat(paste("  Median BMD:", round(median_bmd, 4), "\\n"))
        }}
        
      }}, error = function(e) {{
        cat(paste("  Error calculating BMD for z =", z, ", x =", x, ":", e$message, "\\n"))
      }})
    }}
  }}
  
  return(bmd_results_list)
}}

# Function to select best BMD parameters
select_best_bmd_params <- function(bmd_results_list) {{
  best_params <- NULL
  most_valid <- 0
  
  if (length(bmd_results_list) == 0) {{
      cat("Warning: No BMD results were generated. Cannot select best parameters.\\n")
      return(NULL)
  }}
  
  for (param_name in names(bmd_results_list)) {{
    bmd_result <- bmd_results_list[[param_name]]
    valid_count <- sum(!is.na(bmd_result$res$BMD.zSD) & bmd_result$res$BMD.zSD > 0)
    
    cat(paste("Parameter set", param_name, "- Valid BMDs:", valid_count, "\\n"))
    
    if (valid_count > most_valid) {{
      most_valid <- valid_count
      best_params <- param_name
    }}
  }}
  
  if (!is.null(best_params)) {{
    cat(paste("\\nSelected best parameter set:", best_params, "with", most_valid, "valid BMDs\\n"))
    return(bmd_results_list[[best_params]])
  }} else {{
    cat("Warning: No valid BMD results found across all parameter sets.\\n")
    return(NULL)
  }}
}}

# Function to add gene annotations to BMD results
add_gene_annotations_to_bmd <- function(bmd_results, gene_mapping) {{
  cat("Adding gene symbol annotations to BMD results...\\n")
  
  bmd_results$gene_symbol <- gene_mapping$gene_symbol[match(bmd_results$id, gene_mapping$ensembl_id)]
  
  # Use Ensembl ID if symbol is missing
  bmd_results$gene_symbol[is.na(bmd_results$gene_symbol)] <- 
    bmd_results$id[is.na(bmd_results$gene_symbol)]
  
  cat(paste("Added symbols for", sum(!is.na(bmd_results$gene_symbol)), "genes\\n"))
  
  return(bmd_results)
}}

#================================================================================
# MAIN DOSE-RESPONSE PIPELINE WITH BMD
#================================================================================

run_complete_dromics_pipeline <- function(counts_path, metadata_path, output_dir, 
                                          treatment_to_analyze, control_to_use,
                                          read_threshold = {params.get('read_threshold', 1000000)}, 
                                          outlier_method = "combined",
                                          perform_bmd = {str(params.get('perform_bmd', True)).upper()},
                                          bmd_bootstrap = {str(params.get('bmd_bootstrap', False)).upper()}) {{
  
  cat("\\n=== STARTING COMPLETE DROMICS PIPELINE WITH R-ODAF QC AND BMD ===\\n")
  
  if (!dir.exists(output_dir)) {{ 
    dir.create(output_dir, recursive = TRUE) 
  }}
  
  cat("\\n--- PART I: R-ODAF QC PIPELINE ---\\n")
  
  # Load count data
  count_data <- read_featurecounts_output(counts_path)
  
  # Load metadata
  cat("\\n=== LOADING METADATA ===\\n")
  if (grepl("\\\\.xlsx$", metadata_path, ignore.case = TRUE) || grepl("\\\\.xls$", metadata_path, ignore.case = TRUE)) {{
    metadata_full <- as.data.frame(readxl::read_excel(metadata_path))
  }} else if (grepl("\\\\.csv$", metadata_path, ignore.case = TRUE)) {{
    metadata_full <- read.csv(metadata_path, stringsAsFactors = FALSE)
  }} else {{
    tryCatch({{
      metadata_full <- read.csv(metadata_path, stringsAsFactors = FALSE)
    }}, error = function(e) {{
      metadata_full <<- as.data.frame(readxl::read_excel(metadata_path))
    }})
  }}
  
  # Validate and align
  aligned_data <- validate_metadata_with_batch(count_data, metadata_full)
  
  # Subset data for dose-response
  metadata_subset <- aligned_data$metadata %>%
    filter(chemical %in% c(treatment_to_analyze, control_to_use))
  
  # Add dose information
  metadata_subset$dose_numeric <- ifelse(metadata_subset$chemical == control_to_use, 0, metadata_subset$dose)
  
  common_samples <- intersect(colnames(aligned_data$counts), rownames(metadata_subset))
  counts_subset <- aligned_data$counts[, common_samples]
  metadata_subset <- metadata_subset[common_samples, ]
  
  cat("\\nDose levels in analysis:", paste(sort(unique(metadata_subset$dose_numeric)), collapse = ", "), "\\n")
  
  # Apply QC steps
  read_filter_results <- remove_low_read_samples(counts_subset, metadata_subset, threshold = read_threshold)
  filtered_counts <- filter_genes_by_cpm(read_filter_results$counts, read_filter_results$metadata)
  batch_results <- apply_combat_batch_correction(filtered_counts, read_filter_results$metadata)
  
  # PCA and outlier detection
  dds_qc <- DESeqDataSetFromMatrix(countData = batch_results$corrected_counts, 
                                   colData = read_filter_results$metadata, 
                                   design = ~ 1)
  vst_qc <- vst(dds_qc, blind = TRUE)
  
  # Generate PCA plot
  generate_pca_plot <- function(vst_data, metadata, title) {{
    pca <- prcomp(t(assay(vst_data)))
    pca_data <- as.data.frame(pca$x)
    pca_data <- cbind(pca_data, metadata)
    percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2))
    
    # Color by dose
    pca_data$dose_factor <- as.factor(pca_data$dose_numeric)
    
    p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = dose_factor, label = rownames(pca_data))) +
      geom_point(size = 4, alpha = 0.8) +
      geom_text_repel(size = 3, max.overlaps = 15, box.padding = 0.5) +
      scale_color_viridis_d(name = "Dose") +
      labs(
        title = title,
        x = paste0("PC1: ", percent_var[1], "% variance"),
        y = paste0("PC2: ", percent_var[2], "% variance")
      ) +
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5, size = 16), legend.position = "bottom")
    
    # Add batch shape if exists
    batch_cols <- c("batch", "Batch", "BATCH", "batch_id", "Batch_ID")
    batch_col_found <- NULL
    
    for (col in batch_cols) {{
      if (col %in% colnames(metadata)) {{
        batch_col_found <- col
        break
      }}
    }}
    
    if (!is.null(batch_col_found)) {{
      p <- p + aes(shape = as.factor(metadata[[batch_col_found]])) + labs(shape = "Batch")
    }}
    
    return(p)
  }}
  
  pca_plot <- generate_pca_plot(vst_qc, read_filter_results$metadata, "QC: PCA After Preprocessing")
  ggsave(file.path(output_dir, "QC_PCA_Plot.png"), pca_plot, width = 12, height = 8, dpi = 300)
  
  outlier_results <- identify_pca_outliers(vst_qc, read_filter_results$metadata, method = outlier_method)
  
  # Final clean data
  remaining_samples <- setdiff(rownames(read_filter_results$metadata), outlier_results$outliers)
  final_counts <- batch_results$corrected_counts[, remaining_samples]
  final_metadata <- read_filter_results$metadata[remaining_samples, ]
  
  cat("\\nQC PIPELINE COMPLETE\\n")
  cat("Final samples:", ncol(final_counts), "\\n")
  cat("Final genes:", nrow(final_counts), "\\n")
  
  cat("\\n--- PART II: DROMICS DOSE-RESPONSE ANALYSIS ---\\n")
  
  # Prepare data for DRomics
  dose_vector <- final_metadata$dose_numeric
  samplenames_vector <- rownames(final_metadata)
  
  # Extract batch information if available
  batch_vector <- NULL
  for (col in c("batch", "Batch", "BATCH")) {{
    if (col %in% colnames(final_metadata)) {{
      batch_vector <- final_metadata[[col]]
      break
    }}
  }}
  
  # Gene symbol mapping
  cat("\\n=== CREATING GENE SYMBOL MAPPING ===\\n")
  all_gene_ids <- rownames(final_counts)
  all_gene_symbols <- mapIds(org.Hs.eg.db, 
                             keys = all_gene_ids, 
                             column = "SYMBOL", 
                             keytype = "ENSEMBL", 
                             multiVals = "first")
  
  missing_symbols <- is.na(all_gene_symbols)
  all_gene_symbols[missing_symbols] <- names(all_gene_symbols)[missing_symbols]
  
  gene_mapping <- data.frame(
    ensembl_id = names(all_gene_symbols),
    gene_symbol = all_gene_symbols,
    stringsAsFactors = FALSE
  )
  
  cat("Gene symbols mapped:", sum(!missing_symbols), "out of", length(all_gene_symbols), "\\n")
  
  # Format data for DRomics
  cat("\\n=== FORMATTING DATA FOR DROMICS ===\\n")
  dromics_raw_data <- formatdata4DRomics(
    signalmatrix = final_counts,
    dose = dose_vector,
    samplenames = samplenames_vector
  )
  
  # Import and normalize in DRomics
  num_samples <- ncol(final_counts)
  transformation_method <- if ("{params.get('transformation_method', 'auto')}" == "auto") {{
    if (num_samples < 30) "rlog" else "vst"
  }} else {{
    "{params.get('transformation_method', 'auto')}"
  }}
  cat("Using transformation method:", transformation_method, "for", num_samples, "samples\\n")
  
  o <- RNAseqdata(
    file = dromics_raw_data,
    check = TRUE,
    transfo.method = transformation_method,
    transfo.blind = TRUE,
    round.counts = FALSE
  )
  
  cat("DRomics data import complete\\n")
  
  # Generate DRomics QC plots
  png(file.path(output_dir, "DRomics_QC_plot.png"), width = 800, height = 600)
  plot(o)
  dev.off()
  
  if (!is.null(batch_vector)) {{
    png(file.path(output_dir, "DRomics_PCA_with_batch.png"), width = 800, height = 600)
    PCAdataplot(o, batch = batch_vector, label = TRUE)
    dev.off()
  }} else {{
    png(file.path(output_dir, "DRomics_PCA.png"), width = 800, height = 600)
    PCAdataplot(o, label = TRUE)
    dev.off()
  }}
  
  # Select significantly responding genes
  cat("\\n=== SELECTING DOSE-RESPONSIVE GENES ===\\n")
  s_quad <- itemselect(
    omicdata = o,
    select.method = "{params.get('selection_method', 'quadratic')}",
    FDR = {params.get('fdr_threshold', 0.05)}
  )
  
  n_selected <- length(s_quad$selectindex)
  cat("Genes selected as significantly dose-responsive:", n_selected, "\\n")
  
  # Export selected genes with symbols
  if (n_selected > 0) {{
    selected_genes <- rownames(o$data)[s_quad$selectindex]
    selected_genes_symbols <- gene_mapping$gene_symbol[match(selected_genes, gene_mapping$ensembl_id)]
    
    selected_genes_df <- data.frame(
      ensembl_id = selected_genes,
      gene_symbol = selected_genes_symbols,
      stringsAsFactors = FALSE
    )
    
    write.csv(selected_genes_df, file.path(output_dir, "selected_responding_genes.csv"), row.names = FALSE)
  }}
  
  # Initialize BMD results storage
  bmd_analysis_performed <- FALSE
  bmd_summary <- list()
  
  # Fit dose-response models
  if (n_selected > 0) {{
    cat("\\n=== FITTING DOSE-RESPONSE MODELS ===\\n")
    
    f <- drcfit(
      itemselect = s_quad,
      information.criterion = "{params.get('criterion', 'AIC')}",
      progressbar = FALSE,
      parallel = "no"
    )
    
    n_fitted <- nrow(f$fitres)
    cat("Successfully fitted models:", n_fitted, "\\n")
    
    if (n_fitted > 0) {{
      # Add gene symbols to results
      f$fitres$gene_symbol <- gene_mapping$gene_symbol[match(f$fitres$id, gene_mapping$ensembl_id)]
      
      # Reorder columns for clarity
      f$fitres <- f$fitres[, c("id", "gene_symbol", setdiff(names(f$fitres), c("id", "gene_symbol")))]
      
      # Export fitted models
      write.csv(f$fitres, file.path(output_dir, "dose_response_models.csv"), row.names = FALSE)
      
      # Plot dose-response curves
      n_plots <- min({params.get('n_plots', 12)}, n_fitted)
      
      png(file.path(output_dir, "dose_response_curves.png"), width = 1200, height = 800)
      plot(f, items = n_plots, dose_log_transfo = {params.get('log_transform_dose', 'FALSE')})
      dev.off()
      
      cat("Dose-response curves plotted for top", n_plots, "genes\\n")
      
      #================================================================================
      # PART III: BMD ANALYSIS (INTEGRATED)
      #================================================================================
      
      if (perform_bmd && n_fitted > 0) {{
        cat("\\n--- PART III: BENCHMARK DOSE (BMD) ANALYSIS ---\\n")
        
        # Calculate BMDs with multiple parameter sets
        bmd_results_multiple <- calculate_multiple_bmds(f)
        
        # Select the best parameter set
        r_bmd <- select_best_bmd_params(bmd_results_multiple)
        
        # Diagnostic block
        cat("\\n--- PRE-BOOTSTRAP DIAGNOSTICS ---\\n")
        cat("Value of 'bmd_bootstrap' parameter:", bmd_bootstrap, "\\n")
        cat("Object 'r_bmd' is null:", is.null(r_bmd), "\\n")
        if (!is.null(r_bmd)) {{
            cat("Dimensions of r_bmd$res:", dim(r_bmd$res)[1], "rows,", dim(r_bmd$res)[2], "cols\\n")
        }}
        cat("-------------------------------------\\n")

        if (!is.null(r_bmd)) {{
          # Perform bootstrap if requested
          if (bmd_bootstrap) {{
            cat("\\n=== ROBUST BOOTSTRAP FOR BMD CONFIDENCE INTERVALS ===\\n")
            n_iter <- {params.get('n_bootstrap', 500)}
            n_cores <- min(parallel::detectCores() - 1, 4)
            cat(paste("Target iterations:", n_iter, "| Cores to use:", n_cores, "\\n"))

            b_bmd <- NULL
            
            # Try 'snow' backend first (more compatible, manual setup)
            cat("\\n1. Attempting bootstrap with 'snow' parallel backend (manual setup)...\\n")
            cl <- NULL
            tryCatch({{
              cl <- parallel::makeCluster(n_cores)
              cat("   - Cluster created successfully.\\n")
              parallel::clusterEvalQ(cl, library(DRomics))
              cat("   - DRomics library loaded on workers.\\n")
              parallel::clusterExport(cl, varlist = "r_bmd", envir = environment())
              cat("   - Fit object 'r_bmd' exported to workers.\\n")

              b_bmd <- bmdboot(r = r_bmd, niter = n_iter, progressbar = FALSE, cl = cl)
              cat("Bootstrap with 'snow' backend SUCCEEDED.\\n")
              
            }}, error = function(e_snow) {{
              cat("-> Bootstrap with 'snow' backend FAILED. Reason:", e_snow$message, "\\n")
            }}, finally = {{
              if (!is.null(cl)) {{
                parallel::stopCluster(cl)
                cat("   - Cluster stopped.\\n")
              }}
            }})
            
            # If both parallel methods failed, fall back to sequential
            if (is.null(b_bmd)) {{
              cat("\\n2. All parallel methods failed. Falling back to non-parallel bootstrap. This may be very slow...\\n")
              tryCatch({{
                b_bmd <- bmdboot(r = r_bmd, niter = n_iter, progressbar = TRUE, parallel = "no")
                cat("Non-parallel bootstrap SUCCEEDED.\\n")
              }}, error = function(e_seq) {{
                cat("!!! CRITICAL: Non-parallel bootstrap also FAILED. Reason:", e_seq$message, "\\n")
                cat("!!! Cannot calculate confidence intervals. Proceeding without bootstrap results.\\n")
              }})
            }}

            # Filter results based on whether bootstrap succeeded
            if (!is.null(b_bmd)) {{
              filtered_bmd_results <- bmdfilter(res = b_bmd$res, BMDfilter = "definedCI", BMDtype = "zSD")
            }} else {{
              cat("WARNING: Proceeding with BMD values but without confidence intervals due to bootstrap failure.\\n")
              filtered_bmd_results <- bmdfilter(res = r_bmd$res, BMDfilter = "definedBMD", BMDtype = "zSD")
            }}
            
          }} else {{
            # Filter without bootstrap if not requested
            cat("Bootstrap not requested. Filtering based on defined BMDs.\\n")
            filtered_bmd_results <- bmdfilter(res = r_bmd$res, BMDfilter = "definedBMD", BMDtype = "zSD")
          }}
          
          cat(paste("After filtering:", nrow(filtered_bmd_results), "genes retained\\n"))
          
          # Add gene annotations to BMD results
          if (nrow(filtered_bmd_results) > 0) {{
            final_bmd_results <- add_gene_annotations_to_bmd(filtered_bmd_results, gene_mapping)
            
            # Sort by BMD (most sensitive first)
            if ("BMD.zSD" %in% colnames(final_bmd_results)) {{
              final_bmd_results <- final_bmd_results[order(final_bmd_results$BMD.zSD, na.last = TRUE), ]
            }}
            
            # Calculate BMD statistics
            valid_bmds <- final_bmd_results[!is.na(final_bmd_results$BMD.zSD) & final_bmd_results$BMD.zSD > 0, ]
            
            if (nrow(valid_bmds) > 0) {{
              bmd_analysis_performed <- TRUE
              
              # Calculate quartiles
              quartiles <- quantile(valid_bmds$BMD.zSD, c(0.25, 0.5, 0.75))
              
              # Store BMD summary
              bmd_summary <- list(
                total_bmds = nrow(valid_bmds),
                min_bmd = min(valid_bmds$BMD.zSD),
                max_bmd = max(valid_bmds$BMD.zSD),
                median_bmd = median(valid_bmds$BMD.zSD),
                mean_bmd = mean(valid_bmds$BMD.zSD),
                q1_bmd = quartiles[1],
                q2_bmd = quartiles[2],
                q3_bmd = quartiles[3],
                bootstrap_performed = bmd_bootstrap
              )
              
              cat("\\n=== BMD RESULTS SUMMARY ===\\n")
              cat(paste("Total genes with valid BMDs:", bmd_summary$total_bmds, "\\n"))
              cat(paste("BMD range:", round(bmd_summary$min_bmd, 6), "-", 
                        round(bmd_summary$max_bmd, 6), "μM\\n"))
              cat(paste("Median BMD:", round(bmd_summary$median_bmd, 6), "μM\\n"))
              cat(paste("Mean BMD:", round(bmd_summary$mean_bmd, 6), "μM\\n"))
              
              # Export BMD results
              write.csv(final_bmd_results, file.path(output_dir, "bmd_results.csv"), row.names = FALSE)
              
              # Export top sensitive genes
              most_sensitive <- head(valid_bmds, 50)
              write.csv(most_sensitive, file.path(output_dir, "top50_sensitive_genes_bmd.csv"), row.names = FALSE)
              
              # Create BMD distribution plot
              if (nrow(valid_bmds) >= 5) {{
                bmd_plot <- ggplot(valid_bmds, aes(x = BMD.zSD)) +
                  geom_histogram(bins = min(30, max(10, nrow(valid_bmds)/10)), 
                               fill = "darkblue", alpha = 0.7, color = "white") +
                  labs(
                    title = "Distribution of Benchmark Doses (BMD)",
                    subtitle = paste("Total genes analyzed:", nrow(valid_bmds)),
                    x = "BMD (μM)",
                    y = "Number of Genes"
                  ) +
                  theme_minimal() +
                  theme(plot.title = element_text(size = 14, face = "bold")) +
                  geom_vline(aes(xintercept = median(BMD.zSD)), 
                           color = "red", linetype = "dashed", size = 1) +
                  geom_vline(aes(xintercept = quantile(BMD.zSD, 0.25)), 
                           color = "orange", linetype = "dashed", size = 0.8) +
                  geom_vline(aes(xintercept = quantile(BMD.zSD, 0.75)), 
                           color = "orange", linetype = "dashed", size = 0.8) +
                  annotate("text", x = median(valid_bmds$BMD.zSD) * 1.1, y = Inf, 
                         label = paste("Median =", round(median(valid_bmds$BMD.zSD), 4), "μM"), 
                         vjust = 1.5, color = "red", size = 4)
                
                ggsave(file.path(output_dir, "bmd_distribution.png"), bmd_plot, 
                      width = 10, height = 6, dpi = 300)
                
                # Create sensitivity ranking plot
                if (nrow(valid_bmds) >= 20) {{
                  top_plot_genes <- head(valid_bmds, 25)
                  
                  sensitivity_plot <- ggplot(top_plot_genes, 
                                           aes(x = reorder(gene_symbol, -BMD.zSD), y = BMD.zSD)) +
                    geom_col(fill = "steelblue", alpha = 0.8) +
                    coord_flip() +
                    labs(
                      title = paste("Top 25 Most Sensitive Genes to", treatment_to_analyze),
                      x = "Gene Symbol",
                      y = "BMD (μM)",
                      subtitle = "Lower BMD = Higher Sensitivity"
                    ) +
                    theme_minimal() +
                    theme(axis.text.y = element_text(size = 9))
                  
                  ggsave(file.path(output_dir, "top25_sensitive_genes_bmd.png"), 
                        sensitivity_plot, width = 10, height = 8, dpi = 300)
                }}
              }}
              
              cat("\\nBMD analysis complete!\\n")
            }}
          }}
        }}
      }}
    }}
  }} else {{
    f <- NULL
    n_fitted <- 0
  }}
  
  # Generate comprehensive summary report
  generate_summary_report <- function() {{
    report_file <- file.path(output_dir, "DRomics_Analysis_Report.txt")
    
    sink(report_file)
    
    cat("================================================================================\\n")
    cat("DROMICS DOSE-RESPONSE ANALYSIS REPORT WITH R-ODAF QC AND BMD\\n")
    cat("================================================================================\\n\\n")
    
    cat("Analysis Parameters:\\n")
    cat("- Treatment Analyzed:", treatment_to_analyze, "\\n")
    cat("- Control Used:", control_to_use, "\\n")
    cat("- Read Count Threshold:", format(read_threshold, scientific = FALSE), "reads\\n")
    cat("- Outlier Detection Method:", outlier_method, "\\n")
    cat("- Selection Method:", "{params.get('selection_method', 'quadratic')}", "\\n")
    cat("- FDR Threshold:", {params.get('fdr_threshold', 0.05)}, "\\n")
    cat("- Model Criterion:", "{params.get('criterion', 'AIC')}", "\\n")
    cat("- BMD Analysis:", perform_bmd, "\\n")
    if (perform_bmd) {{
      cat("- BMD Bootstrap:", bmd_bootstrap, "\\n")
    }}
    cat("- Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")
    
    cat("QC PIPELINE RESULTS:\\n")
    cat("- Initial samples:", ncol(counts_subset), "\\n")
    cat("- Final samples retained:", ncol(final_counts), "\\n")
    cat("- Initial genes:", nrow(counts_subset), "\\n")
    cat("- Final genes analyzed:", nrow(final_counts), "\\n")
    cat("- Samples removed (low reads):", length(read_filter_results$removed_samples), "\\n")
    cat("- Samples removed (outliers):", length(outlier_results$outliers), "\\n")
    cat("- Batch correction applied:", batch_results$batch_corrected, "\\n\\n")
    
    cat("DOSE-RESPONSE ANALYSIS RESULTS:\\n")
    cat("- Dose levels analyzed:", length(unique(dose_vector)), "\\n")
    cat("- Samples per dose:\\n")
    print(table(dose_vector))
    cat("\\n")
    cat("- Significantly dose-responsive genes:", n_selected, "\\n")
    cat("- Successfully fitted models:", n_fitted, "\\n\\n")
    
    if (bmd_analysis_performed) {{
      cat("BMD ANALYSIS RESULTS:\\n")
      cat("- Total genes with valid BMDs:", bmd_summary$total_bmds, "\\n")
      cat("- BMD range:", round(bmd_summary$min_bmd, 6), "-", 
          round(bmd_summary$max_bmd, 6), "μM\\n")
      cat("- Median BMD:", round(bmd_summary$median_bmd, 6), "μM\\n")
      cat("- Q1 BMD:", round(bmd_summary$q1_bmd, 6), "μM\\n")
      cat("- Q3 BMD:", round(bmd_summary$q3_bmd, 6), "μM\\n\\n")
    }}
    
    if (n_fitted > 0) {{
      cat("TOP 20 DOSE-RESPONSIVE GENES:\\n")
      cat("Rank\\tGene Symbol\\tEnsembl ID\\tModel\\tAIC")
      if (bmd_analysis_performed) cat("\\tBMD (μM)")
      cat("\\n")
      
      top_genes <- head(f$fitres, 20)
      for (i in 1:nrow(top_genes)) {{
        cat(sprintf("%d\\t%s\\t%s\\t%s\\t%.2f", 
                    i, 
                    top_genes$gene_symbol[i], 
                    top_genes$id[i],
                    top_genes$best.model[i],
                    top_genes$AIC.model[i]))
        
        # Add BMD if available
        if (bmd_analysis_performed && exists("final_bmd_results")) {{
          bmd_val <- final_bmd_results$BMD.zSD[final_bmd_results$id == top_genes$id[i]]
          if (length(bmd_val) > 0 && !is.na(bmd_val[1])) {{
            cat(sprintf("\\t%.6f", bmd_val[1]))
          }} else {{
            cat("\\tNA")
          }}
        }}
        cat("\\n")
      }}
    }}
    
    cat("\\n================================================================================\\n")
    cat("FILES GENERATED:\\n")
    cat("- QC_PCA_Plot.png: PCA visualization after QC\\n")
    cat("- DRomics_QC_plot.png: DRomics quality control visualization\\n")
    cat("- DRomics_PCA.png: DRomics PCA visualization\\n")
    if (!is.null(batch_vector)) {{
      cat("- DRomics_PCA_with_batch.png: PCA with batch information\\n")
    }}
    if (n_selected > 0) {{
      cat("- selected_responding_genes.csv: List of significantly dose-responsive genes\\n")
    }}
    if (n_fitted > 0) {{
      cat("- dose_response_models.csv: Fitted model parameters\\n")
      cat("- dose_response_curves.png: Dose-response curves\\n")
    }}
    if (bmd_analysis_performed) {{
      cat("- bmd_results.csv: Complete BMD results\\n")
      cat("- top50_sensitive_genes_bmd.csv: Top 50 most sensitive genes by BMD\\n")
      cat("- bmd_distribution.png: BMD distribution plot\\n")
      cat("- top25_sensitive_genes_bmd.png: Top 25 sensitive genes visualization\\n")
    }}
    cat("- dromics_summary.json: Analysis summary in JSON format\\n")
    cat("================================================================================\\n")
    
    sink()
    
    cat("Report saved to:", report_file, "\\n")
  }}
  
  generate_summary_report()
  
  # Save summary statistics
  summary_stats <- list(
    # Input data
    total_genes_initial = nrow(counts_subset),
    total_samples_initial = ncol(counts_subset),
    total_genes_final = nrow(final_counts),
    total_samples_final = ncol(final_counts),
    
    # QC metrics
    samples_removed_low_reads = length(read_filter_results$removed_samples),
    samples_removed_outliers = length(outlier_results$outliers),
    batch_corrected = batch_results$batch_corrected,
    
    # Dose-response analysis
    treatment = treatment_to_analyze,
    control = control_to_use,
    dose_levels = length(unique(dose_vector)),
    dose_distribution = as.list(table(dose_vector)),
    transformation_method = transformation_method,
    selection_method = "{params.get('selection_method', 'quadratic')}",
    fdr_threshold = {params.get('fdr_threshold', 0.05)},
    genes_selected = n_selected,
    models_fitted = n_fitted,
    
    # BMD analysis
    bmd_performed = bmd_analysis_performed,
    bmd_summary = if(bmd_analysis_performed) bmd_summary else NULL,
    
    # Parameters
    read_threshold = read_threshold,
    cpm_threshold = {params.get('cpm_threshold', 1)},
    min_prop_expressed = {params.get('min_prop_expressed', 0.75)}
  )
  
  # Save JSON summary
  library(jsonlite)
  write_json(summary_stats, file.path(output_dir, "dromics_summary.json"), pretty = TRUE)
  
  cat("\\n=== DROMICS PIPELINE WITH R-ODAF QC AND BMD COMPLETE ===\\n")
  cat("Output directory:", output_dir, "\\n")
  
  return(list(success = TRUE))
}}

#================================================================================
# STEP 3: SET FILE PATHS AND PARAMETERS (DROMICS)
#================================================================================

cat("\\n=== STEP 3: CONFIGURING DROMICS ANALYSIS PARAMETERS ===\\n")

# File paths from Python parameters
counts_path <- "{params['counts_file']}"
metadata_path <- "{params['metadata_file']}"
output_dir <- "{params['output_dir']}"

# Analysis parameters
TREATMENT_TO_ANALYZE <- "{params['treatment']}"
CONTROL_TO_USE <- "{params['control']}"
READ_THRESHOLD <- {params.get('read_threshold', 1000000)}
OUTLIER_METHOD <- "combined"
PERFORM_BMD <- {str(params.get('perform_bmd', True)).upper()}
BMD_BOOTSTRAP <- {str(params.get('bmd_bootstrap', False)).upper()}

cat("Dromics analysis parameters set:\\n")
cat("- Treatment:", TREATMENT_TO_ANALYZE, "\\n")
cat("- Control:", CONTROL_TO_USE, "\\n")
cat("- Read threshold:", format(READ_THRESHOLD, scientific = FALSE), "\\n")
cat("- Outlier method:", OUTLIER_METHOD, "\\n")
cat("- Perform BMD:", PERFORM_BMD, "\\n")
cat("- BMD Bootstrap:", BMD_BOOTSTRAP, "\\n")
cat("- Output directory:", output_dir, "\\n")

#================================================================================
# STEP 4: RUN THE COMPLETE INTEGRATED DROMICS PIPELINE
#================================================================================

cat("\\n=== STEP 4: RUNNING COMPLETE INTEGRATED DROMICS PIPELINE ===\\n")

tryCatch({{
  pipeline_results <- run_complete_dromics_pipeline(
    counts_path = counts_path,
    metadata_path = metadata_path,
    output_dir = output_dir,
    treatment_to_analyze = TREATMENT_TO_ANALYZE,
    control_to_use = CONTROL_TO_USE,
    read_threshold = READ_THRESHOLD,
    outlier_method = OUTLIER_METHOD,
    perform_bmd = PERFORM_BMD,
    bmd_bootstrap = BMD_BOOTSTRAP
  )
  
  if (pipeline_results$success) {{
    cat("\\n=== DROMICS ANALYSIS COMPLETE ===\\n")
    cat("Results saved to:", output_dir, "\\n\\n")
    
    # Read and display summary
    summary_file <- file.path(output_dir, "dromics_summary.json")
    if (file.exists(summary_file)) {{
      summary <- jsonlite::fromJSON(summary_file)
      cat("Summary:\\n")
      cat("- Genes after QC:", summary$total_genes_final, "\\n")
      cat("- Samples after QC:", summary$total_samples_final, "\\n")
      cat("- Dose-responsive genes:", summary$genes_selected, "\\n")
      cat("- Models fitted:", summary$models_fitted, "\\n")
      
      if (summary$bmd_performed) {{
        cat("\\nBMD Analysis:\\n")
        cat("- Genes with valid BMDs:", summary$bmd_summary$total_bmds, "\\n")
        cat("- Median BMD:", round(summary$bmd_summary$median_bmd, 6), "μM\\n")
        cat("- BMD range:", round(summary$bmd_summary$min_bmd, 6), "-", 
            round(summary$bmd_summary$max_bmd, 6), "μM\\n")
      }}
      
      cat("\\nCheck the output directory for detailed results and plots.\\n")
    }}
  }}
  
}}, error = function(e) {{
  cat("\\nERROR in dromics analysis pipeline:\\n")
  cat(toString(e), "\\n")
  quit(status = 1)
}})
"""
    return r_script
def run_r_analysis(params, progress_callback=None):
    """Execute R analysis for DEG with live progress updates."""
    try:
        r_script = generate_deg_r_script(params)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False) as f:
            f.write(r_script)
            script_path = f.name
        
        script_copy = os.path.join(params['output_dir'], 'analysis_script_deg.R')
        with open(script_copy, 'w') as f:
            f.write(r_script)
        
        process = subprocess.Popen(
            ['Rscript', '--vanilla', script_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            universal_newlines=True
        )
        
        log_lines = []
        for line in process.stdout:
            log_lines.append(line)
            if progress_callback:
                progress_callback(line.strip())

        process.wait()
        stderr = process.stderr.read()
        
        os.remove(script_path)
        
        if process.returncode == 0:
            return True, "DEG analysis completed successfully.", log_lines
        else:
            # Combine stdout and stderr for a complete log on failure
            full_log = "".join(log_lines) + "\n--- STDERR ---\n" + stderr
            return False, f"R script failed with exit code {process.returncode}. See full log for details.", [full_log]
    
    except Exception as e:
        return False, str(e), []
def run_dromics_analysis(params, progress_callback=None):
    """Execute R analysis for DRomics with live progress updates."""
    try:
        r_script = generate_dromics_r_script(params)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False) as f:
            f.write(r_script)
            script_path = f.name

        script_copy = os.path.join(params['output_dir'], 'analysis_script_dromics.R')
        with open(script_copy, 'w') as f:
            f.write(r_script)
        
        process = subprocess.Popen(
            ['Rscript', '--vanilla', script_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            universal_newlines=True
        )

        log_lines = []
        for line in process.stdout:
            log_lines.append(line)
            if progress_callback:
                progress_callback(line.strip())

        process.wait()
        stderr = process.stderr.read()
        
        os.remove(script_path)
        
        if process.returncode == 0:
            return True, "DRomics analysis completed successfully.", log_lines
        else:
            # Combine stdout and stderr for a complete log on failure
            full_log = "".join(log_lines) + "\n--- STDERR ---\n" + stderr
            return False, f"R script failed with exit code {process.returncode}. See full log for details.", [full_log]
    except Exception as e:
        return False, str(e), []

def check_r_dependencies():
    """Check if R and required packages are installed for all analyses."""
    required_packages = [
        'DESeq2', 'edgeR', 'dplyr', 'ggplot2', 'ggrepel',
        'sva', 'org.Hs.eg.db', 'readxl', 'stringr', 'DRomics', 'jsonlite'
    ]
    
    try:
        r_check = subprocess.run(['R', '--version'], capture_output=True, text=True)
        if r_check.returncode != 0:
            return False, "R is not installed", []
        
        check_script = """
        packages <- c({})
        missing <- packages[!packages %in% installed.packages()[,"Package"]]
        if(length(missing) > 0) {{
            cat(paste(missing, collapse=","))
        }}
        """.format(','.join([f'"{p}"' for p in required_packages]))
        
        result = subprocess.run(
            ['R', '--slave', '-e', check_script],
            capture_output=True, text=True
        )
        
        missing = result.stdout.strip().split(',') if result.stdout.strip() else []
        missing = [p for p in missing if p]
        
        if missing:
            return False, f"Missing R packages", missing
        else:
            return True, "All R dependencies installed", []
            
    except FileNotFoundError:
        return False, "R is not installed", []
    except Exception as e:
        return False, str(e), []


def start_processing_pipeline(accession_numbers: list, config: PipelineConfig) -> bool:
    """Start the processing pipeline"""
    try:
        is_valid, errors = config.validate()
        if not is_valid:
            st.error("Configuration validation failed:")
            for error in errors:
                st.error(f"  • {error}")
            return False
        
        directories = [config.output_dir, config.trimmed_dir, config.hisat2_dir]
        created_dirs, failed_dirs = create_directories(directories)
        
        if failed_dirs:
            st.error("Failed to create directories:")
            for failed_dir in failed_dirs:
                st.error(f"  • {failed_dir}")
            return False
        
        missing_tools = validate_tools()
        if missing_tools:
            st.error(f"Missing required tools: {', '.join(missing_tools)}")
            return False
        
        stop_file = os.path.join(config.output_dir, ".stop_pipeline")
        if os.path.exists(stop_file):
            os.remove(stop_file)
        
        st.session_state.processing_status = "initializing"
        st.session_state.total_accessions = len(accession_numbers)
        st.session_state.completed_accessions = 0
        st.session_state.current_step = "Initializing pipeline..."
        st.session_state.thread_error = None
        st.session_state.pipeline_config = config
        st.session_state.cancelled_accessions = set()
        st.session_state.pipeline_results = {
            'total_processed': 0,
            'successful': 0,
            'failed': 0,
            'detailed_results': [],
            'start_time': time.time()
        }
        
        global_logger = ThreadSafeLogger()
        st.session_state.pipeline_logger = global_logger
        
        pipeline_thread = threading.Thread(
            target=run_pipeline_thread,
            args=(accession_numbers, config),
            daemon=True,
            name="RNASeqPipeline"
        )
        
        st.session_state.pipeline_thread = pipeline_thread
        pipeline_thread.start()
        
        time.sleep(1)
        if pipeline_thread.is_alive():
            st.success("Pipeline started successfully!")
            return True
        else:
            st.error("Failed to start pipeline thread")
            return False
            
    except Exception as e:
        st.error(f"Error starting pipeline: {str(e)}")
        st.session_state.processing_status = "error"
        return False
        
        
def get_json_metric(summary_dict, key, default=0):
    """Safely extracts a numeric metric from a JSON object that might be a single-element list."""
    value = summary_dict.get(key, default)
    if isinstance(value, list) and value:
        return value[0]
    return value
def save_parameters_log(params: dict, filename: str = "analysis_parameters.json"):
    """Saves the analysis parameters to a JSON file in the output directory."""
    try:
        output_dir = params.get('output_dir')
        if not output_dir:
            st.warning("Could not save parameters log: output directory not specified in params.")
            return

        # Ensure the directory exists, although it should have been created just before this call
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        log_path = os.path.join(output_dir, filename)
        
        with open(log_path, 'w') as f:
            # We will convert Path objects or other non-serializable types to strings for logging
            serializable_params = {k: str(v) for k, v in params.items()}
            json.dump(serializable_params, f, indent=4)
            
    except Exception as e:
        st.error(f"Failed to write parameters log to {log_path}: {e}")    
    
def main():
    """Main application function with dose-response integration"""
    try:
        init_session_state()
        add_custom_css()
        
        st.markdown("<h1 class='main-header'>ARACRA: RNA-seq Pipeline & Analysis Suite for Chemical Risk Assessment</h1>", unsafe_allow_html=True)
        
        tab1, tab2, tab3, tab4, tab5 = st.tabs([
            "🧬 RNA-seq Pipeline", 
            "📊 DEG Analysis", 
            "📈 Dose-Response",
            "📋 Results", 
            "⚙️ Configuration"
        ])
        
        # RNA-seq Pipeline Tab
        with tab1:
            st.markdown("## RNA-seq Processing Pipeline")
            
            with st.sidebar:
                st.markdown("### Pipeline Configuration")
                
                uploaded_file = st.file_uploader(
                    "Excel file with accession numbers",
                    type=['xlsx', 'xls'],
                    help="Upload an Excel file containing SRA accession numbers"
                )
                
                if uploaded_file is not None:
                    try:
                        df = pd.read_excel(uploaded_file)
                        st.success("Excel file loaded!")
                        
                        column_options = df.columns.tolist()
                        selected_column = st.selectbox("Select accession column:", column_options)
                        
                        if selected_column:
                            raw_accessions = df[selected_column].dropna().astype(str).tolist()
                            clean_accessions = [acc.strip() for acc in raw_accessions if acc.strip()]
                            
                            valid_accessions, invalid_accessions = validate_accession_format(clean_accessions)
                            
                            st.session_state.accession_numbers = valid_accessions
                            
                            if valid_accessions:
                                st.success(f"Loaded {len(valid_accessions)} valid accessions")
                            
                            if invalid_accessions:
                                st.warning(f"Found {len(invalid_accessions)} invalid accessions")
                    
                    except Exception as e:
                        st.error(f"Error reading Excel file: {e}")
                
                # --- AUTO-DETECT MODIFICATION START ---
                st.markdown("#### Tool Paths")
                
                # Use session state variables populated by auto_detect_configuration()
                base_output_dir = st.text_input(
                    "Base Output Directory", 
                    value=st.session_state.detected_base_dir
                )
                
                trimmomatic_jar = st.text_input(
                    "Trimmomatic JAR Path", 
                    value=st.session_state.detected_trim_jar,
                    help="Auto-detected from Conda environment if available"
                )
                
                adapters = st.text_input(
                    "Adapter Sequences", 
                    value=st.session_state.detected_adapters,
                    help="Auto-detected from Conda environment if available"
                )
                
                index_path = st.text_input("HISAT2 Index Path", value="/path/to/hisat2/index/genome")
                reference_gtf = st.text_input("Reference GTF", value="/path/to/genome.gtf")
                
                st.markdown("#### Processing Options")
                
                # Use session state variables for CPU limits
                max_cpu = st.session_state.detected_cpu_max
                default_threads = st.session_state.detected_cpu_default
                
                st.caption(f"System detected: {max_cpu} cores. Defaulting to {default_threads} (50%).")
                
                num_processes = st.slider(
                    "Parallel processes", 
                    min_value=1, 
                    max_value=max_cpu, 
                    value=default_threads
                )
                # --- AUTO-DETECT MODIFICATION END ---

                featurecounts_threads = st.slider("featureCounts threads", 1, 32, 8)
                cleanup_intermediates = st.checkbox("Cleanup intermediate files", value=True)
                max_retries = st.slider("Max retries per step", 1, 5, 3)
            
            col1, col2 = st.columns([3, 1])
            
            with col1:
                with st.expander("Pipeline Steps", expanded=False):
                    st.markdown("""
                    **Step 1: Download** - Fetch raw sequences from SRA database using `prefetch`
                    **Step 2: Convert** - Convert SRA files to FASTQ format using `fasterq-dump`
                    **Step 3: Trim** - Remove adapter sequences and low-quality bases using `Trimmomatic`
                    **Step 4: Align** - Map reads to reference genome using `HISAT2`
                    **Step 5: Count (Combined)** - Generate gene expression count matrix for all samples using `featureCounts`
                    """)
                
                val_col1, val_col2 = st.columns(2)
                
                with val_col1:
                    if st.button("Validate Tools"):
                        with st.spinner("Checking required tools..."):
                            missing_tools = validate_tools()
                            if missing_tools:
                                st.error(f"Missing tools: {', '.join(missing_tools)}")
                            else:
                                st.success("All required tools are available!")
                
                with val_col2:
                    if st.button("Validate Configuration"):
                        base_path = Path(base_output_dir)
                        temp_config = PipelineConfig(
                            output_dir=str(base_path / "fastq"),
                            trimmed_dir=str(base_path / "trimmed"),
                            hisat2_dir=str(base_path / "alignment"),
                            trimmomatic_jar=trimmomatic_jar,
                            adapters=adapters,
                            index_path=index_path,
                            reference_gtf=reference_gtf,
                            num_processes=num_processes,
                            featurecounts_threads=featurecounts_threads,
                            max_retries=max_retries,
                            cleanup_intermediates=cleanup_intermediates
                        )
                        
                        is_valid, errors = temp_config.validate()
                        if is_valid:
                            st.success("Configuration is valid!")
                        else:
                            st.error("Configuration errors:")
                            for error in errors:
                                st.error(f"  • {error}")
            
            with col2:
                st.markdown("### Status")
                
                status_color = {
                    "idle": "status-idle",
                    "running": "status-running", 
                    "completed": "status-running",
                    "error": "status-error",
                    "stopped": "status-error"
                }.get(st.session_state.processing_status, "status-idle")
                
                st.markdown(f"<div class='{status_color}'>Status: {st.session_state.processing_status.title()}</div>", 
                           unsafe_allow_html=True)
                
                st.metric("Loaded Accessions", len(st.session_state.accession_numbers))
                
                if st.session_state.processing_status == "running":
                    if st.session_state.total_accessions > 0:
                        progress = st.session_state.completed_accessions / st.session_state.total_accessions
                        st.progress(progress)
                        st.metric("Progress", f"{st.session_state.completed_accessions}/{st.session_state.total_accessions}")
            
            st.markdown("### Processing Controls")
            
            ctrl_col1, ctrl_col2, ctrl_col3, ctrl_col4 = st.columns(4)
            
            with ctrl_col1:
                can_start = (len(st.session_state.accession_numbers) > 0 and 
                            st.session_state.processing_status not in ["running", "initializing"])
                
                if st.button("Start Pipeline", disabled=not can_start):
                    base_path = Path(base_output_dir)
                    config = PipelineConfig(
                        output_dir=str(base_path / "fastq"),
                        trimmed_dir=str(base_path / "trimmed"),
                        hisat2_dir=str(base_path / "alignment"),
                        trimmomatic_jar=trimmomatic_jar,
                        adapters=adapters,
                        index_path=index_path,
                        reference_gtf=reference_gtf,
                        num_processes=num_processes,
                        featurecounts_threads=featurecounts_threads,
                        max_retries=max_retries,
                        cleanup_intermediates=cleanup_intermediates
                    )
                    
                    success = start_processing_pipeline(st.session_state.accession_numbers, config)
                    if not success:
                        st.error("Failed to start pipeline - check configuration")
            
            with ctrl_col2:
                can_stop = st.session_state.processing_status == "running"
                
                if st.button("Stop Pipeline", disabled=not can_stop):
                    if st.session_state.get('pipeline_config'):
                        stop_file = os.path.join(st.session_state.pipeline_config.output_dir, ".stop_pipeline")
                        try:
                            Path(stop_file).touch()
                            st.warning("Stop signal sent to pipeline...")
                            st.session_state.processing_status = "stopping"
                        except:
                            st.error("Could not create stop file")
            
            with ctrl_col3:
                if st.button("Reset"):
                    for key in ['processing_status', 'pipeline_results', 'completed_accessions', 
                               'total_accessions', 'current_step', 'thread_error']:
                        if key in st.session_state:
                            if key == 'processing_status':
                                st.session_state[key] = "idle"
                            else:
                                st.session_state[key] = type(st.session_state[key])() if hasattr(st.session_state[key], '__class__') else None
                    st.success("Pipeline state reset")
            
            with ctrl_col4:
                can_test = (len(st.session_state.accession_numbers) > 0 and 
                           st.session_state.processing_status not in ["running", "initializing"])
                
                if st.button("Test Single", disabled=not can_test):
                    if st.session_state.accession_numbers:
                        test_accession = st.session_state.accession_numbers[0]
                        st.info(f"Testing with: {test_accession}")
                        
                        base_path = Path(base_output_dir)
                        config = PipelineConfig(
                            output_dir=str(base_path / "fastq"),
                            trimmed_dir=str(base_path / "trimmed"),
                            hisat2_dir=str(base_path / "alignment"),
                            trimmomatic_jar=trimmomatic_jar,
                            adapters=adapters,
                            index_path=index_path,
                            reference_gtf=reference_gtf,
                            num_processes=1,
                            featurecounts_threads=featurecounts_threads,
                            max_retries=1,
                            cleanup_intermediates=False
                        )
                        
                        start_processing_pipeline([test_accession], config)
            
            if st.session_state.processing_status in ["running", "initializing"]:
                st.markdown("### Live Progress")
                
                try:
                    if st.session_state.get('pipeline_config'):
                        results_file = os.path.join(st.session_state.pipeline_config.output_dir, ".pipeline_results.json")
                        if os.path.exists(results_file):
                            with open(results_file, 'r') as f:
                                file_results = json.load(f)
                            
                            st.session_state.pipeline_results.update(file_results)
                            st.session_state.processing_status = file_results.get('processing_status', 'running')
                            st.session_state.completed_accessions = file_results.get('total_processed', 0)
                except:
                    pass
                
                if st.session_state.current_step:
                    st.info(f"Current: {st.session_state.current_step}")
                
                if st.session_state.total_accessions > 0:
                    progress = st.session_state.completed_accessions / st.session_state.total_accessions
                    st.progress(progress, text=f"Overall progress: {progress:.1%}")
                
                if hasattr(st.session_state, 'pipeline_logger'):
                    recent_logs = st.session_state.pipeline_logger.get_messages()[-15:]
                    if recent_logs:
                        st.markdown("#### Live Log")
                        log_text = "\n".join(recent_logs)
                        st.text_area("", value=log_text, height=250, disabled=True)
                
                time.sleep(2)
                st.rerun()
            
            if (st.session_state.pipeline_results and 
                st.session_state.pipeline_results.get('total_processed', 0) > 0):
                
                st.markdown("### Pipeline Results")
                
                results = st.session_state.pipeline_results
                
                metric_col1, metric_col2, metric_col3, metric_col4 = st.columns(4)
                
                with metric_col1:
                    st.metric("Total Processed", results.get('total_processed', 0))
                
                with metric_col2:
                    st.metric("Successful", results.get('successful', 0))
                
                with metric_col3:
                    st.metric("Failed", results.get('failed', 0))
                
                with metric_col4:
                    if 'total_time' in results:
                        hours = results['total_time'] / 3600
                        if hours >= 1:
                            st.metric("Processing Time", f"{hours:.1f}h")
                        else:
                            minutes = results['total_time'] / 60
                            if minutes >= 1:
                                st.metric("Processing Time", f"{minutes:.1f}m")
                            else:
                                st.metric("Processing Time", f"{results['total_time']:.0f}s")
                
                if (st.session_state.get('pipeline_config') and 
                    results.get('processing_status') == 'completed'):
                    count_matrix_path = os.path.join(st.session_state.pipeline_config.hisat2_dir, "combined_counts.out")
                    if os.path.exists(count_matrix_path):
                        st.success("Count matrix is ready for DEG analysis!")
                        if st.button("Auto-load Count Matrix for DEG Analysis"):
                            st.session_state.counts_path = count_matrix_path
                            st.success("Count matrix loaded! Switch to DEG Analysis tab.")

        # DEG Analysis Tab
        with tab2:
            st.markdown("<div class='deg-header'><h2>Differential Expression Analysis</h2><p>R-ODAF DEG Pipeline</p></div>", 
                       unsafe_allow_html=True)
            
            col1, col2 = st.columns([3, 1])
            
            with col2:
                st.markdown("### System Check")
                if st.button("Check R Dependencies", key="deg_r_check"):
                    with st.spinner("Checking R and packages..."):
                        r_ready, message, missing = check_r_dependencies()
                        if r_ready:
                            st.success(message)
                        else:
                            st.error(message)
                            if missing:
                                packages_str = ', '.join([f'"{p}"' for p in missing])
                                st.code(f"# Install in R:\ninstall.packages(c({packages_str}))")
            
            with col1:
                st.markdown("### Input Data")
                
                deg_col1, deg_col2 = st.columns(2)
                
                with deg_col1:
                    st.markdown("#### Count Matrix")
                    
                    if st.session_state.get('counts_path'):
                        st.success(f"Count matrix loaded: {os.path.basename(st.session_state.counts_path)}")
                        if st.button("Clear Count Matrix"):
                            st.session_state.counts_path = None
                            st.rerun()
                    else:
                        counts_file = st.file_uploader(
                            "Upload count matrix file",
                            type=['txt', 'tsv', 'csv', 'out'],
                            help="FeatureCounts output or CSV matrix",
                            key="deg_counts_upload"
                        )
                        
                        if counts_file:
                            temp_dir = tempfile.mkdtemp()
                            counts_path = os.path.join(temp_dir, counts_file.name)
                            with open(counts_path, 'wb') as f:
                                f.write(counts_file.getvalue())
                            st.session_state.counts_path = counts_path
                            st.success(f"Count matrix uploaded: {counts_file.name}")
                
                with deg_col2:
                    st.markdown("#### Metadata")
                    
                    if st.session_state.get('metadata_path'):
                        st.success(f"Metadata loaded: {os.path.basename(st.session_state.metadata_path)}")
                        if st.button("Clear Metadata"):
                            st.session_state.metadata_path = None
                            st.session_state.metadata_df = None
                            st.session_state.chemicals = None
                            st.session_state.conditions = None
                            st.rerun()
                    else:
                        metadata_file = st.file_uploader(
                            "Upload metadata file",
                            type=['xlsx', 'xls', 'csv'],
                            help="Excel or CSV file with sample information",
                            key="deg_metadata_upload"
                        )
                        
                        if metadata_file:
                            temp_dir = tempfile.mkdtemp()
                            metadata_path = os.path.join(temp_dir, metadata_file.name)
                            with open(metadata_path, 'wb') as f:
                                f.write(metadata_file.getvalue())
                            st.session_state.metadata_path = metadata_path
                            
                            if metadata_file.name.endswith('.csv'):
                                metadata_df = pd.read_csv(metadata_path)
                            else:
                                metadata_df = pd.read_excel(metadata_path)
                            
                            st.session_state.metadata_df = metadata_df
                            
                            if 'chemical' in metadata_df.columns:
                                st.session_state.chemicals = metadata_df['chemical'].unique()
                            elif 'condition' in metadata_df.columns:
                                st.session_state.conditions = metadata_df['condition'].unique()
                            
                            st.success(f"Metadata uploaded: {metadata_file.name}")
            
            if st.session_state.get('metadata_df') is not None:
                with st.expander("View Metadata Preview"):
                    st.dataframe(st.session_state.metadata_df.head(10))
            
            st.markdown("### Analysis Parameters")
            
            param_col1, param_col2 = st.columns(2)
            
            with param_col1:
                st.markdown("#### Experimental Design")
                
                if st.session_state.get('chemicals') is not None:
                    treatment = st.selectbox(
                        "Treatment Group",
                        options=st.session_state.chemicals,
                        index=0 if len(st.session_state.chemicals) > 0 else None,
                        key="deg_treatment"
                    )
                    control = st.selectbox(
                        "Control Group",
                        options=st.session_state.chemicals,
                        index=1 if len(st.session_state.chemicals) > 1 else 0,
                        key="deg_control"
                    )
                elif st.session_state.get('conditions') is not None:
                    treatment = st.selectbox(
                        "Treatment Group",
                        options=st.session_state.conditions,
                        index=0 if len(st.session_state.conditions) > 0 else None,
                        key="deg_treatment"
                    )
                    control = st.selectbox(
                        "Control Group",
                        options=st.session_state.conditions,
                        index=1 if len(st.session_state.conditions) > 1 else 0,
                        key="deg_control"
                    )
                else:
                    treatment = st.text_input("Treatment Group Name", key="deg_treatment")
                    control = st.text_input("Control Group Name", value="Control", key="deg_control")
            
            with param_col2:
                st.markdown("#### Statistical Thresholds")
                
                fdr_strict = st.number_input(
                    "Strict FDR",
                    min_value=0.001,
                    max_value=0.1,
                    value=0.01,
                    step=0.001,
                    format="%.3f",
                    key="deg_fdr_strict"
                )
                
                fdr_relaxed = st.number_input(
                    "Relaxed FDR",
                    min_value=0.01,
                    max_value=0.2,
                    value=0.05,
                    step=0.01,
                    format="%.2f",
                    key="deg_fdr_relaxed"
                )
                
                log2fc = st.number_input(
                    "Log2 FC Threshold",
                    min_value=0.0,
                    max_value=5.0,
                    value=1.0,
                    step=0.1,
                    format="%.1f",
                    key="deg_log2fc"
                )
            
            with st.expander("Advanced Filtering Options"):
                adv_col1, adv_col2, adv_col3 = st.columns(3)
                
                with adv_col1:
                    read_threshold = st.number_input(
                        "Min Reads/Sample",
                        min_value=0,
                        max_value=10000000,
                        value=1000000,
                        step=100000,
                        key="deg_read_thresh"
                    )
                
                with adv_col2:
                    cpm_threshold = st.number_input(
                        "CPM Threshold",
                        min_value=0.0,
                        max_value=10.0,
                        value=1.0,
                        step=0.1,
                        key="deg_cpm_thresh"
                    )
                
                with adv_col3:
                    min_prop = st.slider(
                        "Min Proportion",
                        min_value=0.0,
                        max_value=1.0,
                        value=0.75,
                        step=0.05,
                        key="deg_min_prop"
                    )
                
                batch_correction = st.checkbox(
                    "Apply Batch Correction (if batch column exists)",
                    value=True,
                    key="deg_batch_correction"
                )
            
            deg_output_dir_val = st.text_input(
                "DEG Analysis Output Directory",
                value=str(Path.home() / "deg_results" / datetime.now().strftime("%Y%m%d_%H%M%S"))
            )
            
            st.markdown("---")
            
            can_run_deg = all([
                st.session_state.get('counts_path'),
                st.session_state.get('metadata_path'),
                treatment,
                control,
                deg_output_dir_val
            ])
            
            if st.button("🚀 Run DEG Analysis", disabled=not can_run_deg, type="primary", use_container_width=True):
                params = {
                    'counts_file': st.session_state.counts_path,
                    'metadata_file': st.session_state.metadata_path,
                    'output_dir': deg_output_dir_val,
                    'treatment': treatment,
                    'control': control,
                    'fdr_strict': fdr_strict,
                    'fdr_relaxed': fdr_relaxed,
                    'log2fc_threshold': log2fc,
                    'read_threshold': int(read_threshold),
                    'cpm_threshold': cpm_threshold,
                    'min_prop_expressed': min_prop,
                    'batch_correction': batch_correction
                }
                
                Path(deg_output_dir_val).mkdir(parents=True, exist_ok=True)
                
                # ADDED: Save parameters to a log file
                save_parameters_log(params, filename="deg_parameters.json")
                
                with st.spinner("🧬 Running R-ODAF analysis... This may take several minutes"):
                    status_text = st.empty()
                    
                    def progress_updater(message):
                        status_text.text(message)
                    
                    success, message, output = run_r_analysis(params, progress_callback=progress_updater)
                    
                    status_text.empty()
                    
                    if success:
                        st.success("✅ DEG analysis completed successfully!")
                        st.session_state.deg_analysis_complete = True
                        st.session_state.deg_output_dir = deg_output_dir_val
                        
                        summary_path = os.path.join(deg_output_dir_val, "analysis_summary.json")
                        if os.path.exists(summary_path):
                            with open(summary_path, 'r') as f:
                                summary = json.load(f)
                                st.session_state.deg_summary = summary
                            
                            st.markdown("### 📊 Analysis Summary")
                            col1, col2, col3, col4 = st.columns(4)
                            
                            with col1:
                                total_genes = get_json_metric(summary, 'total_genes')
                                st.metric("Total Genes", f"{total_genes:,}")
                            
                            with col2:
                                sig_relaxed = get_json_metric(summary, 'sig_relaxed')
                                fdr_relaxed_val = get_json_metric(summary, 'fdr_relaxed', 0.05)
                                st.metric(f"FDR < {fdr_relaxed_val}", f"{sig_relaxed:,}")
                            
                            with col3:
                                upreg = get_json_metric(summary, 'upregulated')
                                st.metric("⬆ Upregulated", f"{upreg:,}")
                            
                            with col4:
                                downreg = get_json_metric(summary, 'downregulated')
                                st.metric("⬇ Downregulated", f"{downreg:,}")
                        
                        st.info("📁 Results saved to: " + deg_output_dir_val)
                        st.success("🎯 Switch to the Results tab to view detailed analysis!")
                        
                        log_file = os.path.join(deg_output_dir_val, "R-ODAF_Analysis_Summary_Report.txt")
                        if os.path.exists(log_file):
                            with st.expander("View Analysis Summary Report", expanded=False):
                                with open(log_file, 'r') as f:
                                    st.text(f.read())
                    else:
                        st.error(f"❌ DEG analysis failed")
                        with st.expander("Error Details", expanded=True):
                            st.code(message[:2000])

        # Dose-Response Analysis Tab
        with tab3:
            st.markdown("<div class='deg-header'><h2>Dose-Response Analysis with BMD</h2><p>DRomics Pipeline with Benchmark Dose Analysis</p></div>", 
                       unsafe_allow_html=True)
            
            col1, col2 = st.columns([3, 1])
            
            with col2:
                st.markdown("### System Check")
                if st.button("Check DRomics Dependencies", key="dromics_check"):
                    with st.spinner("Checking R and packages..."):
                        r_ready, message, missing = check_r_dependencies()
                        if r_ready:
                            st.success(message)
                        else:
                            st.error(message)
                            if missing:
                                packages_str = ', '.join([f'"{p}"' for p in missing])
                                st.code(f"# Install in R:\ninstall.packages(c({packages_str}))")
            
            with col1:
                st.markdown("### Input Data")
                st.info("💡 Tip: You can use the same count matrix and metadata from the RNA-seq pipeline or DEG analysis")
                
                dr_col1, dr_col2 = st.columns(2)
                
                with dr_col1:
                    st.markdown("#### Count Matrix")
                    
                    use_existing_counts = False
                    if st.session_state.get('counts_path'):
                        use_existing_counts = st.checkbox(
                            f"Use loaded count matrix ({os.path.basename(st.session_state.counts_path)})",
                            value=True,
                            key="use_existing_counts"
                        )
                    
                    if use_existing_counts:
                        dr_counts_path = st.session_state.counts_path
                        st.success("Using existing count matrix")
                    else:
                        dr_counts_file = st.file_uploader(
                            "Upload count matrix",
                            type=['txt', 'tsv', 'csv', 'out'],
                            key="dr_counts_upload"
                        )
                        
                        if dr_counts_file:
                            temp_dir = tempfile.mkdtemp()
                            dr_counts_path = os.path.join(temp_dir, dr_counts_file.name)
                            with open(dr_counts_path, 'wb') as f:
                                f.write(dr_counts_file.getvalue())
                            st.success(f"Count matrix uploaded: {dr_counts_file.name}")
                        else:
                            dr_counts_path = None
                
                with dr_col2:
                    st.markdown("#### Metadata with Dose Info")
                    
                    use_existing_metadata = False
                    if st.session_state.get('metadata_path'):
                        use_existing_metadata = st.checkbox(
                            f"Use loaded metadata ({os.path.basename(st.session_state.metadata_path)})",
                            value=True,
                            key="use_existing_metadata"
                        )
                    
                    dr_metadata_df = None
                    if use_existing_metadata:
                        dr_metadata_path = st.session_state.metadata_path
                        dr_metadata_df = st.session_state.metadata_df
                        st.success("Using existing metadata")
                    else:
                        dr_metadata_file = st.file_uploader(
                            "Upload metadata with dose column",
                            type=['xlsx', 'xls', 'csv'],
                            key="dr_metadata_upload"
                        )
                        
                        if dr_metadata_file:
                            temp_dir = tempfile.mkdtemp()
                            dr_metadata_path = os.path.join(temp_dir, dr_metadata_file.name)
                            with open(dr_metadata_path, 'wb') as f:
                                f.write(dr_metadata_file.getvalue())
                            
                            if dr_metadata_file.name.endswith('.csv'):
                                dr_metadata_df = pd.read_csv(dr_metadata_path)
                            else:
                                dr_metadata_df = pd.read_excel(dr_metadata_path)
                            
                            st.success(f"Metadata uploaded: {dr_metadata_file.name}")
                        else:
                            dr_metadata_path = None
                            dr_metadata_df = None
            
            if dr_metadata_df is not None:
                with st.expander("View Metadata with Dose Information"):
                    if 'dose' in dr_metadata_df.columns:
                        st.success("✓ Dose column found")
                        dose_info = dr_metadata_df.groupby(['chemical', 'dose']).size().reset_index(name='n_samples')
                        st.dataframe(dose_info, use_container_width=True, hide_index=True)
                    else:
                        st.warning("⚠️ No 'dose' column found in metadata. Please ensure your metadata contains dose information.")
                    
                    st.dataframe(dr_metadata_df.head(10), use_container_width=True)
            
            st.markdown("### Analysis Parameters")
            
            dr_param_col1, dr_param_col2 = st.columns(2)
            
            with dr_param_col1:
                st.markdown("#### Experimental Design")
                dr_treatment, dr_control = None, None
                if dr_metadata_df is not None and 'chemical' in dr_metadata_df.columns:
                    chemicals = dr_metadata_df['chemical'].unique()
                    dr_treatment = st.selectbox(
                        "Treatment Chemical",
                        options=chemicals,
                        index=0 if len(chemicals) > 0 else None,
                        key="dr_treatment",
                        help="Select the chemical treatment to analyze"
                    )
                    dr_control = st.selectbox(
                        "Control Chemical",
                        options=chemicals,
                        index=1 if len(chemicals) > 1 else 0,
                        key="dr_control",
                        help="Select the control/vehicle treatment"
                    )
                    
                    if dr_treatment and 'dose' in dr_metadata_df.columns:
                        treatment_doses = dr_metadata_df[dr_metadata_df['chemical'] == dr_treatment]['dose'].unique()
                        st.info(f"Doses for {dr_treatment}: {sorted(treatment_doses)}")
                else:
                    dr_treatment = st.text_input("Treatment Chemical", key="dr_treatment_text")
                    dr_control = st.text_input("Control Chemical", value="Control", key="dr_control_text")
                
                selection_method = st.selectbox(
                    "Gene Selection Method",
                    options=["quadratic", "linear", "ANOVA"],
                    index=0,
                    key="dr_selection_method",
                    help="Method for selecting dose-responsive genes"
                )
            
            with dr_param_col2:
                st.markdown("#### Statistical Parameters")
                
                dr_fdr = st.number_input(
                    "FDR Threshold",
                    min_value=0.01,
                    max_value=0.2,
                    value=0.05,
                    step=0.01,
                    format="%.2f",
                    key="dr_fdr",
                    help="False Discovery Rate threshold for gene selection"
                )
                
                model_criterion = st.selectbox(
                    "Model Selection Criterion",
                    options=["AIC", "BIC"],
                    index=0,
                    key="dr_criterion",
                    help="Information criterion for model selection"
                )
                
                n_plots = st.number_input(
                    "Number of Curves to Plot",
                    min_value=1,
                    max_value=50,
                    value=12,
                    step=1,
                    key="dr_n_plots",
                    help="Number of top dose-response curves to visualize"
                )
            
            st.markdown("### 📊 Benchmark Dose (BMD) Analysis")
            
            bmd_col1, bmd_col2 = st.columns(2)
            
            with bmd_col1:
                perform_bmd = st.checkbox(
                    "Perform BMD Analysis",
                    value=True,
                    key="perform_bmd",
                    help="Calculate Benchmark Doses for dose-responsive genes"
                )
                
                if perform_bmd:
                    bmd_bootstrap = st.checkbox(
                        "Bootstrap for Confidence Intervals",
                        value=False,
                        key="bmd_bootstrap",
                        help="Perform bootstrap to calculate BMD confidence intervals (slower but more accurate)"
                    )
                    
                    z_values = st.text_input(
                        "BMD z-values",
                        value="0.1, 0.5, 1.0",
                        key="bmd_z_values",
                        help="Comma-separated z-values for BMD calculation (lower = more sensitive)"
                    )
            
            with bmd_col2:
                if perform_bmd:
                    x_values = st.text_input(
                        "BMD x-values (fold change)",
                        value="10",
                        key="bmd_x_values",
                        help="Comma-separated x-values for BMD calculation"
                    )
                    
                    if bmd_bootstrap:
                        n_bootstrap = st.slider(
                            "Bootstrap Iterations",
                            min_value=100,
                            max_value=5000,
                            value=500,
                            step=100,
                            key="n_bootstrap",
                            help="Number of bootstrap iterations (more = slower but more accurate)"
                        )
            
            with st.expander("Advanced Options"):
                adv_col1, adv_col2 = st.columns(2)
                
                with adv_col1:
                    log_transform_dose = st.checkbox(
                        "Log-transform dose for plotting",
                        value=False,
                        key="dr_log_dose",
                        help="Apply log transformation to dose axis in plots"
                    )
                    
                    dr_read_threshold = st.number_input(
                        "Min Reads/Sample (DR)",
                        min_value=0,
                        max_value=10000000,
                        value=1000000,
                        step=100000,
                        key="dr_read_thresh"
                    )
                
                with adv_col2:
                    transformation_method = st.selectbox(
                        "Normalization Method",
                        options=["auto", "rlog", "vst"],
                        index=0,
                        key="dr_transform",
                        help="auto: rlog for <30 samples, vst for ≥30 samples"
                    )
                    
                    dr_cpm_threshold = st.number_input(
                        "CPM Threshold (DR)",
                        min_value=0.0,
                        max_value=10.0,
                        value=1.0,
                        step=0.1,
                        key="dr_cpm_thresh"
                    )
            
            dr_output_dir_val = st.text_input(
                "Output Directory for Dose-Response Results",
                value=str(Path.home() / "dromics_results" / datetime.now().strftime("%Y%m%d_%H%M%S")),
                key="dr_output_dir_input"
            )
            
            st.markdown("---")
            
            can_run_dr = all([
                'dr_counts_path' in locals() and dr_counts_path,
                'dr_metadata_path' in locals() and dr_metadata_path,
                dr_treatment,
                dr_control,
                dr_output_dir_val
            ])
            
            if not can_run_dr:
                missing = []
                if not ('dr_counts_path' in locals() and dr_counts_path): missing.append("count matrix")
                if not ('dr_metadata_path' in locals() and dr_metadata_path): missing.append("metadata")
                if not dr_treatment: missing.append("treatment selection")
                if not dr_control: missing.append("control selection")
                if not dr_output_dir_val: missing.append("output directory")
                st.warning(f"Missing required inputs: {', '.join(missing)}")
            
            if st.button("🚀 Run Dose-Response Analysis with BMD", 
                        disabled=not can_run_dr, 
                        type="primary", 
                        use_container_width=True,
                        key="run_dromics"):
                
                params = {
                    'counts_file': dr_counts_path,
                    'metadata_file': dr_metadata_path,
                    'output_dir': dr_output_dir_val,
                    'treatment': dr_treatment,
                    'control': dr_control,
                    'selection_method': selection_method,
                    'fdr_threshold': dr_fdr,
                    'criterion': model_criterion,
                    'n_plots': n_plots,
                    'log_transform_dose': str(log_transform_dose).upper(),
                    'perform_bmd': perform_bmd,
                    'bmd_bootstrap': bmd_bootstrap if perform_bmd else False,
                    'read_threshold': dr_read_threshold,
                    'cpm_threshold': dr_cpm_threshold,
                    'min_prop_expressed': 0.75
                }
                
                if perform_bmd:
                    params['bmd_z_values'] = z_values
                    params['bmd_x_values'] = x_values
                    if bmd_bootstrap:
                        params['n_bootstrap'] = n_bootstrap
                
                if transformation_method != "auto":
                    params['transformation_method'] = transformation_method
                
                Path(dr_output_dir_val).mkdir(parents=True, exist_ok=True)

                # ADDED: Save parameters to a log file
                save_parameters_log(params, filename="dromics_parameters.json")
                
                with st.spinner("🧬 Running DRomics dose-response analysis with BMD... This may take several minutes"):
                    status_text = st.empty()
                    
                    def progress_updater(message):
                        status_text.text(message)

                    success, message, output = run_dromics_analysis(params, progress_callback=progress_updater)
                    
                    if success:
                        st.success("✅ Dose-response analysis with BMD completed successfully!")
                        st.session_state.dr_analysis_complete = True
                        st.session_state.dr_output_dir = dr_output_dir_val
                        
                        summary_path = os.path.join(dr_output_dir_val, "dromics_summary.json")
                        if os.path.exists(summary_path):
                            with open(summary_path, 'r') as f:
                                dr_summary = json.load(f)
                            
                            st.markdown("### 📊 Analysis Summary")
                            col1, col2, col3, col4 = st.columns(4)
                            
                            with col1:
                                genes_final = get_json_metric(dr_summary, 'total_genes_final')
                                st.metric("Genes Analyzed", f"{genes_final:,}")
                            with col2:
                                dose_levels = get_json_metric(dr_summary, 'dose_levels')
                                st.metric("Dose Levels", dose_levels)
                            with col3:
                                genes_selected = get_json_metric(dr_summary, 'genes_selected')
                                st.metric("Responding Genes", 
                                         genes_selected,
                                         delta=f"FDR < {dr_fdr}")
                            with col4:
                                models_fitted = get_json_metric(dr_summary, 'models_fitted')
                                st.metric("Models Fitted", models_fitted)
                            
                            if dr_summary.get('bmd_performed') and dr_summary.get('bmd_summary'):
                                st.markdown("### 🎯 BMD Analysis Results")
                                bmd_summary = dr_summary['bmd_summary']
                                
                                bmd_col1, bmd_col2, bmd_col3, bmd_col4 = st.columns(4)
                                
                                with bmd_col1:
                                    total_bmds = get_json_metric(bmd_summary, 'total_bmds')
                                    st.metric("Genes with BMD", total_bmds)
                                with bmd_col2:
                                    median_bmd = get_json_metric(bmd_summary, 'median_bmd')
                                    st.metric("Median BMD", f"{median_bmd:.3f} μM")
                                with bmd_col3:
                                    min_bmd = get_json_metric(bmd_summary, 'min_bmd')
                                    st.metric("Min BMD", f"{min_bmd:.3f} μM")
                                with bmd_col4:
                                    max_bmd = get_json_metric(bmd_summary, 'max_bmd')
                                    st.metric("Max BMD", f"{max_bmd:.3f} μM")
                        
                        st.markdown("### 📈 Visualizations")
                        
                        viz_tabs = st.tabs(["Quality Control", "PCA Analysis", "Dose-Response Curves", "BMD Analysis"])
                        
                        with viz_tabs[0]:
                            qc_plot = os.path.join(dr_output_dir_val, "DRomics_QC_plot.png")
                            if os.path.exists(qc_plot): st.image(qc_plot, caption="Quality Control Visualization", use_column_width=True)
                        with viz_tabs[1]:
                            pca_batch_plot = os.path.join(dr_output_dir_val, "DRomics_PCA_with_batch.png")
                            pca_plot = os.path.join(dr_output_dir_val, "DRomics_PCA.png")
                            if os.path.exists(pca_batch_plot): st.image(pca_batch_plot, caption="PCA with Batch Information", use_column_width=True)
                            elif os.path.exists(pca_plot): st.image(pca_plot, caption="PCA Plot", use_column_width=True)
                        with viz_tabs[2]:
                            curves_plot = os.path.join(dr_output_dir_val, "dose_response_curves.png")
                            if os.path.exists(curves_plot): st.image(curves_plot, caption=f"Top {n_plots} Dose-Response Curves", use_column_width=True)
                            else: st.info("No dose-response curves generated.")
                        with viz_tabs[3]:
                            if perform_bmd:
                                bmd_dist_plot = os.path.join(dr_output_dir_val, "bmd_distribution.png")
                                bmd_sens_plot = os.path.join(dr_output_dir_val, "top25_sensitive_genes_bmd.png")
                                if os.path.exists(bmd_dist_plot): st.image(bmd_dist_plot, caption="BMD Distribution", use_column_width=True)
                                if os.path.exists(bmd_sens_plot): st.image(bmd_sens_plot, caption="Top 25 Most Sensitive Genes by BMD", use_column_width=True)
                                bmd_results_file = os.path.join(dr_output_dir_val, "top50_sensitive_genes_bmd.csv")
                                if os.path.exists(bmd_results_file):
                                    st.markdown("#### Top Sensitive Genes by BMD")
                                    try:
                                        bmd_df = pd.read_csv(bmd_results_file)
                                        if not bmd_df.empty:
                                            display_cols = ['gene_symbol', 'BMD.zSD', 'model']
                                            display_cols = [col for col in display_cols if col in bmd_df.columns]
                                            st.dataframe(bmd_df[display_cols].head(20), use_container_width=True, hide_index=True)
                                    except: pass
                            else:
                                st.info("BMD analysis was not performed.")
                        
                        st.markdown("### 📥 Download Results")
                        st.info(f"📁 All results saved to: {dr_output_dir_val}")
                        with st.expander("View Analysis Report", expanded=False):
                            report_file = os.path.join(dr_output_dir_val, "DRomics_Analysis_Report.txt")
                            if os.path.exists(report_file):
                                with open(report_file, 'r') as f:
                                    st.text(f.read())
                    else:
                        st.error(f"❌ Dose-response analysis failed")
                        with st.expander("Error Details", expanded=True):
                            st.code(message[:2000])
        
        # Results Tab
        with tab4:
            st.markdown("## Analysis Results")
            
            result_tabs = st.tabs(["RNA-seq Pipeline Results", "DEG Analysis Results", "Dose-Response Results"])

            with result_tabs[0]:
                if (st.session_state.pipeline_results and 
                    st.session_state.pipeline_results.get('total_processed', 0) > 0):
                    results = st.session_state.pipeline_results
                    col1, col2 = st.columns(2)
                    with col1:
                        csv_data = export_results_to_csv(results)
                        st.download_button(
                            label="Download CSV Results",
                            data=csv_data,
                            file_name=f"rnaseq_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                            mime="text/csv"
                        )
                    with col2:
                        summary_text = create_processing_summary(results)
                        st.download_button(
                            label="Download Summary",
                            data=summary_text,
                            file_name=f"rnaseq_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
                            mime="text/plain"
                        )
                else:
                    st.info("No pipeline results available. Run the RNA-seq pipeline first.")
            
            with result_tabs[1]:
                if st.session_state.get('deg_analysis_complete') and st.session_state.get('deg_output_dir'):
                    output_dir = st.session_state.deg_output_dir
                    
                    if st.session_state.get('deg_summary'):
                        summary = st.session_state.deg_summary
                        st.markdown("#### Summary Statistics")
                        
                        metric_col1, metric_col2, metric_col3, metric_col4 = st.columns(4)
                        with metric_col1:
                            total_genes = get_json_metric(summary, 'total_genes')
                            st.metric("Total Genes", f"{total_genes:,}")
                        with metric_col2:
                            sig_strict = get_json_metric(summary, 'sig_strict')
                            fdr_strict_val = get_json_metric(summary, 'fdr_strict', 0.01)
                            st.metric(f"FDR < {fdr_strict_val}", f"{sig_strict:,}")
                        with metric_col3:
                            upreg = get_json_metric(summary, 'upregulated')
                            log2fc_val = get_json_metric(summary, 'log2fc_threshold', 1.0)
                            st.metric("Upregulated", 
                                     f"{upreg:,}",
                                     delta=f"FC > {2**log2fc_val:.1f}")
                        with metric_col4:
                            downreg = get_json_metric(summary, 'downregulated')
                            log2fc_val = get_json_metric(summary, 'log2fc_threshold', 1.0)
                            st.metric("Downregulated", 
                                     f"{downreg:,}",
                                     delta=f"FC < {1/(2**log2fc_val):.1f}")
                    
                    st.markdown("#### Visualizations")
                    viz_col1, viz_col2 = st.columns(2)
                    with viz_col1:
                        pca_file = os.path.join(output_dir, "Final_PCA_Plot.png")
                        if os.path.exists(pca_file):
                            st.image(pca_file, caption="PCA Plot", use_column_width=True)
                    with viz_col2:
                        volcano_file = os.path.join(output_dir, "Volcano_Plot_Overall.png")
                        if os.path.exists(volcano_file):
                            st.image(volcano_file, caption="Volcano Plot", use_column_width=True)
                    
                    st.markdown("#### Top Differentially Expressed Genes")
                    sig_file = os.path.join(output_dir, "Custom_Filtered_DEGs.csv")
                    if os.path.exists(sig_file):
                        try:
                            sig_df = pd.read_csv(sig_file)
                            if not sig_df.empty:
                                st.dataframe(sig_df.head(20), use_container_width=True, hide_index=True)
                            else:
                                st.info("No significant genes found with current thresholds")
                        except pd.errors.EmptyDataError:
                             st.info("No significant genes found with current thresholds")
                    
                    st.markdown("#### Download DEG Results")
                    dl_col1, dl_col2, dl_col3 = st.columns(3)
                    with dl_col1:
                        all_genes_file = os.path.join(output_dir, "Overall_Results_With_Symbols.csv")
                        if os.path.exists(all_genes_file):
                            with open(all_genes_file, 'rb') as f:
                                st.download_button("All Genes CSV", data=f.read(), file_name="all_genes.csv", mime="text/csv")
                    with dl_col2:
                        sig_custom_file = os.path.join(output_dir, "Custom_Filtered_DEGs.csv")
                        if os.path.exists(sig_custom_file):
                            with open(sig_custom_file, 'rb') as f:
                                st.download_button("Significant DEGs CSV", data=f.read(), file_name="significant_degs.csv", mime="text/csv")
                    with dl_col3:
                        zip_path = os.path.join(output_dir, "deg_results.zip")
                        with zipfile.ZipFile(zip_path, 'w') as zipf:
                            for root, _, files in os.walk(output_dir):
                                for file in files:
                                    if not file.endswith('.zip'):
                                        zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), output_dir))
                        with open(zip_path, 'rb') as f:
                            st.download_button("All DEG Results (ZIP)", data=f.read(), file_name=f"deg_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip", mime="application/zip")
                else:
                    st.info("No DEG analysis results available. Run the DEG analysis first.")

            with result_tabs[2]:
                if st.session_state.get('dr_analysis_complete') and st.session_state.get('dr_output_dir'):
                    dr_output_dir_res = st.session_state.dr_output_dir
                    
                    summary_path = os.path.join(dr_output_dir_res, "dromics_summary.json")
                    if os.path.exists(summary_path):
                        with open(summary_path, 'r') as f:
                            dr_summary = json.load(f)
                        st.markdown("#### Dose-Response Analysis Summary")
                        metric_col1, metric_col2, metric_col3, metric_col4 = st.columns(4)
                        
                        with metric_col1:
                            genes_final = get_json_metric(dr_summary, 'total_genes_final')
                            st.metric("Genes Analyzed", f"{genes_final:,}")
                        with metric_col2:
                            dose_levels = get_json_metric(dr_summary, 'dose_levels')
                            st.metric("Dose Levels", dose_levels)
                        with metric_col3:
                            genes_selected = get_json_metric(dr_summary, 'genes_selected')
                            st.metric("Responding Genes", genes_selected)
                        with metric_col4:
                            models_fitted = get_json_metric(dr_summary, 'models_fitted')
                            st.metric("Models Fitted", models_fitted)
                    
                    st.markdown("#### Top Dose-Responsive Genes")
                    models_file = os.path.join(dr_output_dir_res, "dose_response_models.csv")
                    if os.path.exists(models_file):
                        try:
                            models_df = pd.read_csv(models_file)
                            if not models_df.empty:
                                display_cols = [col for col in ['gene_symbol', 'id', 'best.model', 'AIC.model', 'BMD.zSD'] if col in models_df.columns]
                                top_models = models_df[display_cols].head(20).copy()
                                for col in top_models.select_dtypes(include=np.number).columns:
                                    top_models[col] = top_models[col].round(3)
                                st.dataframe(top_models, use_container_width=True, hide_index=True)
                            else:
                                st.info("No models were successfully fitted")
                        except pd.errors.EmptyDataError:
                            st.info("No models were successfully fitted")
                    
                    st.markdown("#### Visualizations")
                    viz_col1, viz_col2 = st.columns(2)
                    with viz_col1:
                        qc_plot = os.path.join(dr_output_dir_res, "DRomics_QC_plot.png")
                        if os.path.exists(qc_plot):
                            st.image(qc_plot, caption="Quality Control", use_column_width=True)
                    with viz_col2:
                        pca_plot = os.path.join(dr_output_dir_res, "DRomics_PCA.png")
                        if os.path.exists(pca_plot):
                            st.image(pca_plot, caption="PCA Analysis", use_column_width=True)
                    curves_plot = os.path.join(dr_output_dir_res, "dose_response_curves.png")
                    if os.path.exists(curves_plot):
                        st.markdown("#### Dose-Response Curves")
                        st.image(curves_plot, use_column_width=True)
                    
                    st.markdown("#### Download Dose-Response Results")
                    dl_col1, dl_col2, dl_col3 = st.columns(3)
                    with dl_col1:
                        selected_file = os.path.join(dr_output_dir_res, "selected_responding_genes.csv")
                        if os.path.exists(selected_file):
                            with open(selected_file, 'rb') as f:
                                st.download_button("Responding Genes CSV", data=f.read(), file_name="responding_genes.csv", mime="text/csv")
                    with dl_col2:
                        models_file = os.path.join(dr_output_dir_res, "dose_response_models.csv")
                        if os.path.exists(models_file):
                            with open(models_file, 'rb') as f:
                                st.download_button("Fitted Models CSV", data=f.read(), file_name="dose_response_models.csv", mime="text/csv")
                    with dl_col3:
                        zip_path = os.path.join(dr_output_dir_res, "dromics_results.zip")
                        with zipfile.ZipFile(zip_path, 'w') as zipf:
                            for root, _, files in os.walk(dr_output_dir_res):
                                for file in files:
                                    if not file.endswith('.zip'):
                                        zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), dr_output_dir_res))
                        with open(zip_path, 'rb') as f:
                            st.download_button("All DR Results (ZIP)", data=f.read(), file_name=f"dromics_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip", mime="application/zip")
                else:
                    st.info("No dose-response analysis results available. Run the dose-response analysis first.")

        # Configuration Tab
        with tab5:
            st.markdown("## System Configuration")
            
            config_col1, config_col2 = st.columns(2)
            
            with config_col1:
                st.markdown("### RNA-seq Pipeline Tools")
                
                if st.button("Check Pipeline Tools", key="config_pipeline_tools"):
                    with st.spinner("Checking tools..."):
                        missing_tools = validate_tools()
                        if missing_tools:
                            st.error("Missing tools:")
                            for tool in missing_tools:
                                st.write(f"❌ {tool}")
                            st.info("Install missing tools using conda or your package manager")
                        else:
                            st.success("All pipeline tools are available!")
                            for tool in ['prefetch', 'fasterq-dump', 'hisat2', 'samtools', 'featureCounts', 'java']:
                                st.write(f"✅ {tool}")
            
            with config_col2:
                st.markdown("### Analysis (R)")
                
                if st.button("Check R Environment", key="config_r_env"):
                    with st.spinner("Checking R and packages..."):
                        r_ready, message, missing = check_r_dependencies()
                        if r_ready:
                            st.success("R environment is ready!")
                            st.write("✅ R installed")
                            st.write("✅ All required packages available")
                        else:
                            st.error(f"R environment issue: {message}")
                            if missing:
                                st.markdown("**Missing packages:**")
                                for pkg in missing:
                                    st.write(f"❌ {pkg}")
                                
                                cran_packages = [p for p in missing if p not in ["DESeq2", "edgeR", "org.Hs.eg.db", "DRomics"]]
                                install_cran_str = ""
                                if cran_packages:
                                    packages_str = ', '.join([f'"{p}"' for p in cran_packages])
                                    install_cran_str = f"install.packages(c({packages_str}))"

                                st.code(f"""
# Install missing packages in R:
{install_cran_str}

# For Bioconductor packages:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "org.Hs.eg.db"))

# For DRomics (if missing):
# install.packages("DRomics")
""")
            
            st.markdown("### About")
            st.info("""
            **ARACRA: Automated RNA-seq and Comprehensive 'omics Analysis for Chemical Risk Assessment**
            
            This application combines:
            - **RNA-seq Processing Pipeline**: Download, convert, trim, align, and count RNA-seq data from SRA.
            - **DEG Analysis (R-ODAF)**: Differential expression analysis using DESeq2 and edgeR.
            - **Dose-Response Analysis (DRomics)**: Modeling dose-dependent effects and calculating Benchmark Doses (BMD).
            
            **Workflow:**
            1. Process raw RNA-seq data to generate a count matrix.
            2. Upload a metadata file with sample information.  
            3. Run differential expression and/or dose-response analysis.
            4. View and download results and visualizations.
            
            **Version:** 1.2  
            **Dependencies:** Python 3.8+, R 4.0+, Bioconductor packages
            """)

    except Exception as e:
        st.error(f"Critical application error: {e}")
        st.code(traceback.format_exc())
        
# This block is needed to run the app
if __name__ == "__main__":
    main()
