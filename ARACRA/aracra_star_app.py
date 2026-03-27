#!/usr/bin/env python3
"""
ARACRA Pipeline — Streamlit App (v2)
Two-phase design:
  Phase 1: SRA → fastp → STAR/HISAT2 → QC → featureCounts/Salmon → MultiQC
  Phase 2: (after user reviews QC) → DESeq2 → DRomics/BMD

Key change from v1: preprocessing and analysis are decoupled.
After Phase 1, the user sees QC outlier results and chooses which samples
to exclude BEFORE analysis runs. No more silent auto-exclusion.
"""

import streamlit as st
import subprocess, os, signal, json, time, re, shutil
import pandas as pd
from pathlib import Path
from datetime import datetime

# ─────────────────────────────────────────────────────────────────────────────
#  Configuration & Environment
# ─────────────────────────────────────────────────────────────────────────────

PIPELINE_DIR = Path(__file__).parent.absolute()
MAIN_NF      = PIPELINE_DIR / "main.nf"
ENV_FILE     = PIPELINE_DIR / ".env"
STATE_FILE   = PIPELINE_DIR / ".pipeline_state.json"
PID_FILE     = PIPELINE_DIR / ".pipeline.pid"


def _load_env() -> dict:
    """Parse the .env file written by setup.sh."""
    env = {}
    if ENV_FILE.exists():
        for line in ENV_FILE.read_text().splitlines():
            line = line.strip()
            if line and "=" in line and not line.startswith("#"):
                k, v = line.split("=", 1)
                env[k.strip()] = v.strip()
    return env


ENV = _load_env()


def _env(key: str, fallback: str = "") -> str:
    return ENV.get(key, fallback)


CONDA_BASE       = _env("CONDA_BASE")
ENV_NAME         = _env("ENV_NAME", "test_ARACRA")
ENV_PATH         = _env("ENV_PATH")
ENV_BIN          = _env("ENV_BIN") or (str(Path(ENV_PATH) / "bin") if ENV_PATH else "")
DEFAULT_WORK     = _env("WORK_DIR",  str(Path.home() / "aracra_star_work/work"))
DEFAULT_OUT      = _env("OUT_DIR",   str(Path.home() / "aracra_star_work/results"))
LOG_FILE         = Path(_env("LOG_FILE", str(Path.home() / "aracra_star_work/pipeline.log")))
STAR_INDEX       = _env("STAR_INDEX",       str(Path.home() / "databases/hg38_reference/star_index_hg38"))
HISAT2_INDEX     = _env("HISAT2_INDEX",     str(Path.home() / "databases/hg38_reference/hisat2_index_hg38/genome_tran"))
SALMON_INDEX     = _env("SALMON_INDEX",     str(Path.home() / "databases/hg38_reference/salmon_index_hg38"))
GTF_PATH         = _env("GTF_PATH",         str(Path.home() / "databases/hg38_reference/gencode.v44.primary_assembly.annotation.gtf"))
BED_PATH         = _env("BED_PATH",         str(Path.home() / "databases/annotations/hg38_RefSeq.bed"))
HK_BED_PATH      = _env("HK_BED_PATH",     str(Path.home() / "databases/annotations/hg38_housekeeping.bed"))
REFFLAT_PATH     = _env("REFFLAT_PATH",     str(Path.home() / "databases/annotations/hg38_refFlat.txt"))
SCREEN_CONF_PATH = _env("SCREEN_CONF_PATH", str(Path.home() / "databases/fastq_screen_genomes/FastQ_Screen_Genomes/fastq_screen.conf"))
RECOMMENDED_ALIGNER = _env("RECOMMENDED_ALIGNER", "star")
SETUP_DONE       = ENV_FILE.exists() and bool(ENV_BIN)


def _get_ram_gb() -> float:
    try:
        text = open("/proc/meminfo").read()
        return int(re.search(r"MemTotal:\s+(\d+)", text).group(1)) / 1_048_576
    except Exception:
        return 0.0


SYSTEM_RAM_GB = _get_ram_gb()

# ─────────────────────────────────────────────────────────────────────────────
#  Process Environment & Launcher
# ─────────────────────────────────────────────────────────────────────────────


def _process_env() -> dict:
    env = os.environ.copy()
    if ENV_BIN:
        env["PATH"] = f"{ENV_BIN}:{env.get('PATH', '')}"
    return env


def launch_pipeline(cmd: list, work_dir: str) -> int:
    """Launch a Nextflow command as a detached background process."""
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(LOG_FILE, "a") as f:
        f.write(f"\n{'=' * 60}\nStarted : {datetime.now()}\nCommand : {' '.join(cmd)}\n{'=' * 60}\n\n")
    log_fd = os.open(str(LOG_FILE), os.O_WRONLY | os.O_APPEND | os.O_CREAT)
    proc = subprocess.Popen(
        cmd, cwd=work_dir,
        stdout=log_fd, stderr=subprocess.STDOUT,
        start_new_session=True,
        env=_process_env(),
    )
    os.close(log_fd)
    return proc.pid


def stop_pipeline():
    pid = _load_pid()
    if not pid:
        return
    try:
        os.killpg(os.getpgid(pid), signal.SIGTERM)
    except Exception:
        try:
            os.kill(pid, signal.SIGTERM)
        except Exception:
            pass
    _clear_pid()


# ─────────────────────────────────────────────────────────────────────────────
#  State & PID Management
# ─────────────────────────────────────────────────────────────────────────────


def _save_state(s: dict):
    STATE_FILE.write_text(json.dumps(s, default=str))


def _load_state() -> dict:
    try:
        return json.loads(STATE_FILE.read_text()) if STATE_FILE.exists() else {}
    except Exception:
        return {}


def _save_pid(pid: int):
    PID_FILE.write_text(str(pid))


def _load_pid() -> int | None:
    try:
        return int(PID_FILE.read_text().strip()) if PID_FILE.exists() else None
    except Exception:
        return None


def _clear_pid():
    if PID_FILE.exists():
        PID_FILE.unlink()


def _is_running() -> bool:
    pid = _load_pid()
    if not pid:
        return False
    try:
        os.kill(pid, 0)  # check if process exists
        try:
            cmd = open(f"/proc/{pid}/cmdline", "rb").read().decode(errors="replace")
            # Check for nextflow OR java (nextflow runs as java process)
            return "nextflow" in cmd.lower() or ("java" in cmd.lower() and "nextflow" in cmd.lower())
        except FileNotFoundError:
            # /proc not available (non-Linux): be conservative, check start time
            start = state.get("start_time")
            if start:
                try:
                    elapsed = (datetime.now() - datetime.fromisoformat(start)).total_seconds()
                    return elapsed < 86400  # assume it's ours if < 24h old
                except Exception:
                    pass
            return False
    except (ProcessLookupError, PermissionError):
        return False


def _elapsed(state: dict) -> str:
    t0 = state.get("start_time")
    if not t0:
        return ""
    try:
        s = int((datetime.now() - datetime.fromisoformat(t0)).total_seconds())
        return f"{s // 3600:02d}:{(s % 3600) // 60:02d}:{s % 60:02d}"
    except Exception:
        return ""


# ─────────────────────────────────────────────────────────────────────────────
#  Log Reading & Progress Parsing
# ─────────────────────────────────────────────────────────────────────────────


def log_tail(n: int = 80) -> str:
    if not LOG_FILE.exists():
        return "Log will appear here when pipeline starts\u2026"
    try:
        # Read from end of file — much faster for large logs
        with open(LOG_FILE, "rb") as f:
            f.seek(0, 2)  # seek to end
            size = f.tell()
            # Read last ~20KB which is plenty for 80 lines
            read_size = min(size, 20480)
            f.seek(max(0, size - read_size))
            chunk = f.read().decode(errors="replace")
        lines = chunk.splitlines()
        return "\n".join(lines[-n:])
    except Exception as e:
        return f"Error reading log: {e}"


def parse_progress() -> dict:
    progress = {}
    if not LOG_FILE.exists():
        return progress
    try:
        for line in LOG_FILE.read_text(errors="replace").splitlines()[-300:]:
            m = re.search(
                r"process\s*>\s*(\w+)[^\[]*\[(\d+)%\]\s+(\d+)\s+of\s+(\d+)", line
            )
            if m:
                progress[m.group(1)] = {
                    "pct": int(m.group(2)),
                    "done": int(m.group(3)),
                    "total": int(m.group(4)),
                }
    except Exception:
        pass
    return progress


PROCESS_LABELS = {
    "PARSE_METADATA":           "0. Parse metadata",
    "DOWNLOAD_SRA":             "1. Download SRA",
    "FASTP":                    "2. Fastp QC / Trim",
    "FASTQ_SCREEN":             "3. FastQ Screen",
    "STAR_LOAD_GENOME":         "4a. Load genome",
    "STAR_ALIGN":               "4b. STAR alignment",
    "STAR_REMOVE_GENOME":       "4c. Unload genome",
    "QC_OUTLIER_CHECK":         "4d. Alignment QC",
    "HISAT2_ALIGN":             "4b. HISAT2 alignment",
    "SAMTOOLS_INDEX":           "5a. Index BAMs",
    "SAMTOOLS_FLAGSTAT":        "5b. flagstat",
    "SAMTOOLS_IDXSTATS":        "5c. idxstats",
    "INFER_STRANDEDNESS":       "6. Strandedness",
    "GENEBODY_COVERAGE_HK":     "7a. Coverage (HK)",
    "GENEBODY_COVERAGE_FULL":   "7b. Coverage (full)",
    "READ_DISTRIBUTION":        "8. Read distribution",
    "PICARD_RNA_METRICS":       "9a. Picard RNA",
    "QUALIMAP_BAMQC":           "9b. Qualimap BAM QC",
    "QUALIMAP_RNASEQ":          "9c. Qualimap RNA-seq",
    "FEATURECOUNTS":            "10. featureCounts",
    "SALMON_QUANT":             "11a. Salmon (aligned)",
    "SALMON_PSEUDO":            "11b. Salmon pseudomap",
    "MERGE_SALMON_PSEUDO":      "11c. Merge pseudo counts",
    "MAKE_COUNT_MATRIX":        "12a. Count matrix",
    "MERGE_SALMON":             "12b. Merge Salmon",
    "MULTIQC":                  "13. MultiQC",
    "DESEQ2_ANALYSIS":          "14. DESeq2",
    "DROMICS_ANALYSIS":         "15. DRomics / BMD",
}


def detect_completion() -> str | None:
    if not LOG_FILE.exists():
        return None
    try:
        tail = "\n".join(LOG_FILE.read_text(errors="replace").splitlines()[-20:])
        if "SUCCESS" in tail or "\u2705" in tail:
            return "completed"
        if "FAILED" in tail or "\u274c" in tail or "ERROR ~" in tail:
            return "failed"
    except Exception:
        pass
    return None


# ─────────────────────────────────────────────────────────────────────────────
#  System & Tool Checks
# ─────────────────────────────────────────────────────────────────────────────

TOOLS = {
    "nextflow":      ["nextflow", "-version"],
    "prefetch":      ["prefetch", "--version"],
    "fasterq-dump":  ["fasterq-dump", "--version"],
    "fastp":         ["fastp", "--version"],
    "STAR":          ["STAR", "--version"],
    "hisat2":        ["hisat2", "--version"],
    "samtools":      ["samtools", "--version"],
    "featureCounts": ["featureCounts", "-v"],
    "salmon":        ["salmon", "--version"],
    "fastq_screen":  ["fastq_screen", "--version"],
    "Rscript":       ["Rscript", "--version"],
    "multiqc":       ["multiqc", "--version"],
}


@st.cache_data(ttl=60)
def check_tools() -> dict:
    result = {}
    env = _process_env()
    for name, cmd in TOOLS.items():
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=10, env=env)
            ver = (r.stdout + r.stderr).strip().split("\n")[0][:60]
            result[name] = {"ok": True, "ver": ver or "OK"}
        except FileNotFoundError:
            result[name] = {"ok": False, "ver": "NOT FOUND"}
        except Exception as e:
            result[name] = {"ok": False, "ver": str(e)[:60]}
    return result


@st.cache_data(ttl=300)
def check_r_packages() -> tuple[bool, list]:
    script = """
pkgs <- c("DESeq2","edgeR","sva","org.Hs.eg.db","GO.db","KEGGREST","AnnotationDbi",
          "dplyr","readxl","ggplot2","ggrepel","stringr","jsonlite","DRomics","optparse")
miss <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
if (length(miss)) cat("MISSING:", paste(miss, collapse=",")) else cat("OK")
"""
    try:
        r = subprocess.run(
            ["Rscript", "-e", script],
            capture_output=True, text=True, timeout=30, env=_process_env(),
        )
        out = (r.stdout + r.stderr).strip()
        if "MISSING:" in out:
            return False, out.replace("MISSING:", "").strip().split(",")
        return True, []
    except Exception as e:
        return False, [str(e)]


def get_sysinfo() -> dict:
    info = {}
    try:
        mem = open("/proc/meminfo").read()
        total_kb = int(re.search(r"MemTotal:\s+(\d+)", mem).group(1))
        avail_kb = int(re.search(r"MemAvailable:\s+(\d+)", mem).group(1))
        info["ram_total"] = f"{total_kb / 1_048_576:.1f} GB"
        info["ram_avail"] = f"{avail_kb / 1_048_576:.1f} GB"
    except Exception:
        info["ram_total"] = info["ram_avail"] = "?"
    info["cpus"] = str(os.cpu_count() or "?")
    try:
        d = os.statvfs(DEFAULT_OUT if Path(DEFAULT_OUT).exists() else str(Path.home()))
        info["disk"] = f"{(d.f_bavail * d.f_frsize) / 1e9:.0f} GB free"
    except Exception:
        info["disk"] = "?"
    return info


# ─────────────────────────────────────────────────────────────────────────────
#  Helpers
# ─────────────────────────────────────────────────────────────────────────────


def _make_run_name(treatment: str, control: str) -> str:
    def _short(s):
        s = re.sub(r"[^A-Za-z0-9]+", "_", s.strip())[:20].strip("_")
        return s or "X"
    return f"{_short(treatment)}_vs_{_short(control)}"


# ── Column mapping helpers ───────────────────────────────────────────────────

STANDARD_COLS = ["Sample_Name", "Treatment", "Dose", "Batch", "Type"]

def _guess_column(df_columns: list, standard: str) -> str:
    """Try to auto-detect which user column maps to a standard name."""
    candidates = {
        "Sample_Name": ["sample_name", "sample", "samplename", "sample_id", "sampleid",
                        "srr", "accession", "run", "geo_accession"],
        "Treatment":   ["treatment", "chemical", "compound", "substance", "drug",
                        "exposure", "agent", "chemical_name", "chem"],
        "Dose":        ["dose", "concentration", "conc", "dose_um", "dose_µm",
                        "dose_mg", "level", "amount"],
        "Batch":       ["batch", "plate", "run_batch", "experiment", "replicate_group",
                        "block", "experiment_batch"],
        "Type":        ["type", "sample_type", "group", "condition", "role",
                        "control_treatment", "trt_type"],
    }
    # Exact match first
    for col in df_columns:
        if col == standard:
            return col
    # Case-insensitive match
    lower_map = {c.lower().strip(): c for c in df_columns}
    for alias in candidates.get(standard, []):
        if alias in lower_map:
            return lower_map[alias]
    return ""


def _show_column_mapper(df: pd.DataFrame, key_prefix: str = "colmap") -> dict:
    """Show column mapping dropdowns. Returns {standard_name: user_column}."""
    cols = list(df.columns)
    none_option = "— not available —"
    col_map = {}

    st.markdown("**Map your columns** (auto-detected where possible)")
    mapping_cols = st.columns(len(STANDARD_COLS))
    for i, std in enumerate(STANDARD_COLS):
        with mapping_cols[i]:
            guessed = _guess_column(cols, std)
            options = [none_option] + cols
            default_idx = options.index(guessed) if guessed in options else 0
            required = std in ("Sample_Name", "Treatment", "Dose", "Type")
            label = f"**{std}**" if required else std
            chosen = st.selectbox(
                label, options, index=default_idx,
                key=f"{key_prefix}_{std}",
                help=f"{'Required' if required else 'Optional'}: select the column for {std}",
            )
            if chosen != none_option:
                col_map[std] = chosen
    return col_map


def _apply_column_mapping(df: pd.DataFrame, col_map: dict) -> pd.DataFrame:
    """Rename user columns to standard names. Returns a new DataFrame."""
    df_out = df.copy()
    rename_map = {user_col: std_name for std_name, user_col in col_map.items()
                  if user_col in df_out.columns and user_col != std_name}
    if rename_map:
        df_out = df_out.rename(columns=rename_map)
    # Add missing optional columns with defaults
    if "Batch" not in df_out.columns:
        df_out["Batch"] = 1
    if "Type" not in df_out.columns and "Treatment" in df_out.columns:
        # Try to infer: if a value contains "control" (case-insensitive), mark as control
        df_out["Type"] = df_out["Treatment"].apply(
            lambda x: "control" if "control" in str(x).lower() else "treatment"
        )
    return df_out


def _subset_for_comparison(
    count_matrix_path: str,
    meta_df: pd.DataFrame,
    treatment: str,
    control: str,
    outdir: str,
    run_name: str,
) -> tuple:
    """
    Subset count matrix and metadata to only samples matching treatment/control.
    Returns (subsetted_counts_path, subsetted_metadata_path, n_samples) or (None, None, 0) on error.
    """
    if meta_df is None or "Treatment" not in meta_df.columns or "Sample_Name" not in meta_df.columns:
        return None, None, 0

    # Filter metadata to selected treatment + control
    meta_sub = meta_df[meta_df["Treatment"].isin([treatment, control])].copy()
    if meta_sub.empty:
        return None, None, 0

    sample_names = set(meta_sub["Sample_Name"].astype(str).tolist())

    # Read count matrix and subset columns
    cm_path = Path(count_matrix_path)
    first_line = cm_path.read_text().split("\n", 1)[0]
    sep = "\t" if "\t" in first_line else ","
    df_counts = pd.read_csv(cm_path, sep=sep)

    # Find the gene ID column (first column or named gene_id/Geneid)
    gene_col = df_counts.columns[0]
    sample_cols = [c for c in df_counts.columns if c in sample_names]

    if not sample_cols:
        return None, None, 0

    df_counts_sub = df_counts[[gene_col] + sample_cols]

    # Save subsetted files
    run_dir = Path(outdir) / "runs" / run_name
    run_dir.mkdir(parents=True, exist_ok=True)
    sub_counts_path = run_dir / "subset_counts.csv"
    sub_meta_path = run_dir / "subset_metadata.csv"

    df_counts_sub.to_csv(sub_counts_path, index=False)
    # Only keep samples that are in both metadata and count matrix
    meta_sub = meta_sub[meta_sub["Sample_Name"].astype(str).isin(sample_cols)]
    meta_sub.to_csv(sub_meta_path, index=False)

    return str(sub_counts_path), str(sub_meta_path), len(meta_sub)


def list_runs(outdir: str) -> list[dict]:
    runs_dir = Path(outdir) / "runs"
    if not runs_dir.exists():
        return []
    runs = []
    for d in sorted(runs_dir.iterdir(), reverse=True):
        if d.is_dir():
            runs.append({
                "name":        d.name,
                "path":        d,
                "has_deg":     (d / "deg_results" / "analysis_summary.json").exists(),
                "has_dromics": (d / "dromics_results" / "dromics_summary.json").exists(),
                "mtime":       datetime.fromtimestamp(d.stat().st_mtime).strftime("%Y-%m-%d %H:%M"),
            })
    return runs


# ─────────────────────────────────────────────────────────────────────────────
#  Outlier QC — load & display
# ─────────────────────────────────────────────────────────────────────────────


def _unbox(v, default=None):
    """Unwrap R jsonlite list-wrapped scalars: [1.23] → 1.23."""
    if isinstance(v, list):
        return v[0] if v else default
    return v if v is not None else default


def _flatten_qc_json(d: dict) -> dict:
    """Normalise R jsonlite quirks (scalars wrapped in arrays) to plain Python."""
    if "flagged" in d:
        flat = []
        for x in d["flagged"]:
            if isinstance(x, list):
                flat.extend(str(i) for i in x if i)
            elif x:
                flat.append(str(x))
        d["flagged"] = flat

    for s in d.get("samples", []):
        for k in ("sample", "condition", "dose", "p"):
            if k in s:
                s[k] = str(_unbox(s[k], ""))
        for k in ("flagged", "confirmed_outlier", "iqr_flag", "grubbs_flag"):
            if k in s:
                s[k] = bool(_unbox(s[k], False))
        for k in ("pc1", "pc2", "pct_unique", "total_reads_M",
                   "pct_multi", "pct_too_short", "G", "G_crit", "fence", "q1", "q3", "value"):
            if k in s:
                try:
                    s[k] = float(_unbox(s[k], 0))
                except (ValueError, TypeError):
                    s[k] = 0.0
        for dk in ("iqr_details", "grubbs_details"):
            if dk in s and isinstance(s[dk], dict):
                for k2, v2 in s[dk].items():
                    if isinstance(v2, list) and len(v2) == 1:
                        s[dk][k2] = v2[0]
    return d


def load_outlier_flags(outdir: str, run_path=None) -> dict:
    """Load Layer 1 (alignment) and Layer 2 (PCA) outlier JSONs."""
    flags = {}
    align_path = Path(outdir) / "qc" / "outlier" / "outlier_flag.json"
    if align_path.exists():
        try:
            flags["alignment"] = _flatten_qc_json(json.loads(align_path.read_text()))
        except Exception:
            pass
    if run_path:
        for sub in ("deg_results", "dromics_results"):
            pca_path = Path(run_path) / sub / "pca_outlier_flag.json"
            if pca_path.exists():
                try:
                    flags["pca"] = _flatten_qc_json(json.loads(pca_path.read_text()))
                    break
                except Exception:
                    pass
    return flags


def _show_qc_review(flags: dict, all_samples: list[str], key_prefix: str = "full", run_name: str = "") -> list[str]:
    """
    Display QC tables and sample-exclusion checkboxes.
    Returns list of samples the user chose to exclude.
    key_prefix differentiates widgets when called from multiple tabs.
    run_name is used to locate PCA plots on disk.
    """
    align = flags.get("alignment", {})
    pca   = flags.get("pca", {})

    align_outliers = [
        s["sample"] for s in align.get("samples", []) if s.get("confirmed_outlier")
    ]
    pca_outliers = [
        s["sample"] for s in pca.get("samples", []) if s.get("flagged")
    ]
    auto_exclude = set(align_outliers) | set(pca_outliers)
    any_flagged  = bool(auto_exclude)

    if any_flagged:
        st.warning(
            "**QC Warning** — outlier samples detected (pre-ticked below). "
            "Review the tables and adjust before continuing.",
            icon="\u26a0\ufe0f",
        )
    else:
        st.success(
            "**Sample QC** — no outliers detected. "
            "You can still manually exclude samples if needed.",
            icon="\u2714\ufe0f",
        )

    # Layer 1 table
    if align.get("samples"):
        st.markdown("##### Layer 1 — Alignment QC (Grubbs + IQR)")
        rows = []
        for s in align["samples"]:
            rows.append({
                "Sample":      s.get("sample", ""),
                "% Mapped":    s.get("pct_unique", ""),
                "Reads (M)":   s.get("total_reads_M", ""),
                "% Multi":     s.get("pct_multi", ""),
                "% Too Short": s.get("pct_too_short", ""),
                "IQR":         "\u2718" if s.get("iqr_flag")    else "\u2714",
                "Grubbs":      "\u2718" if s.get("grubbs_flag") else "\u2714",
                "Status":      "\u26a0 OUTLIER" if s.get("confirmed_outlier") else "\u2714 OK",
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
        st.caption(
            f"Mean mapping: {align.get('mean_mapping', '?')}% "
            f"\u00b1 {align.get('std_mapping', '?')}%  |  "
            f"n={align.get('n_samples', '?')}  |  "
            "Flagged = both IQR AND Grubbs agree"
        )

    # Layer 2 table
    if pca.get("samples"):
        st.markdown("##### Layer 2 — Expression PCA")

        # ── Display PCA scatter plot if available ──
        # Search comprehensively — the plot could be in qc_results, deg_results, or dromics_results
        _pca_plot_found = False
        _run_name = run_name or state.get("run_name", "")
        _pca_search_paths = [
            # QC-only run output
            Path(outdir) / "runs" / _run_name / "qc_results" / "pca_qc_plot.png",
            # DESeq2 run output (detect_pca_outliers runs inside DESeq2)
            Path(outdir) / "runs" / _run_name / "deg_results" / "pca_qc_plot.png",
            # DRomics run output
            Path(outdir) / "runs" / _run_name / "dromics_results" / "pca_qc_plot.png",
            # Alignment-level QC directory
            Path(outdir) / "qc" / "outlier" / "pca_qc_plot.png",
        ]
        for _pca_candidate in _pca_search_paths:
            if _pca_candidate.exists():
                st.image(str(_pca_candidate),
                         caption="PCA \u2014 Expression QC (outliers marked with \u00d7 and red ring)",
                         use_container_width=True)
                _pca_plot_found = True
                break

        # Fallback: only look in the current run's subdirectories, not across all runs
        if not _pca_plot_found and _run_name:
            for _subdir in ("qc_results", "deg_results", "dromics_results"):
                _fb = Path(outdir) / "runs" / _run_name / _subdir / "pca_qc_plot.png"
                if _fb.exists():
                    st.image(str(_fb),
                             caption="PCA \u2014 Expression QC (outliers marked with \u00d7 and red ring)",
                             use_container_width=True)
                    _pca_plot_found = True
                    break

        if not _pca_plot_found:
            st.caption(
                "\u2139 PCA plot not found. Run **Step 1: QC Check** in the "
                "DESeq2/DRomics tab to generate it."
            )

        # PCA data table
        rows2 = []
        for s in pca["samples"]:
            rows2.append({
                "Sample": s.get("sample", ""),
                "PC1":    round(float(s.get("pc1", 0)), 2),
                "PC2":    round(float(s.get("pc2", 0)), 2),
                "Group":  s.get("condition", s.get("dose", "")),
                "Status": "\u26a0 OUTLIER" if s.get("flagged") else "\u2714 OK",
            })
        with st.expander(f"PCA coordinates ({len(rows2)} samples)", expanded=False):
            st.dataframe(pd.DataFrame(rows2), use_container_width=True, hide_index=True)

    # Exclusion checkboxes
    st.markdown("---")
    st.markdown(
        "**Select samples to exclude** "
        "*(outliers are pre-ticked \u2014 untick to keep)*"
    )
    # Only use previously-selected exclusions if they came from the CURRENT QC context.
    # This prevents stale exclusions from a previous run from persisting.
    prev = set(st.session_state.get("exclude_samples", []))
    _qc_fresh = st.session_state.get("_qc_exclusions_set_for", "") == key_prefix
    defaults = prev if (prev and _qc_fresh) else auto_exclude

    to_exclude = []
    n_cols = min(4, max(1, len(all_samples)))
    cols = st.columns(n_cols)
    for i, samp in enumerate(all_samples):
        is_outlier = samp in auto_exclude
        label = f"{'⚠ ' if is_outlier else ''}{samp}"
        with cols[i % n_cols]:
            if st.checkbox(label, value=(samp in defaults), key=f"excl_{key_prefix}_{samp}"):
                to_exclude.append(samp)

    st.session_state["exclude_samples"] = to_exclude
    st.session_state["_qc_exclusions_set_for"] = key_prefix
    return to_exclude


# ─────────────────────────────────────────────────────────────────────────────
#  Results Display
# ─────────────────────────────────────────────────────────────────────────────

def _jv(d: dict, k: str, default=0):
    """Extract a value, unwrapping R jsonlite list wrappers."""
    v = d.get(k, default)
    return v[0] if isinstance(v, list) and v else v


def display_results(run_dir: Path):
    deg_dir = run_dir / "deg_results"
    dr_dir  = run_dir / "dromics_results"

    # --- DESeq2 ---
    st.markdown("#### DESeq2")
    summary_path = deg_dir / "analysis_summary.json"
    if summary_path.exists():
        s = json.loads(summary_path.read_text())
        n_samples = _jv(s, "samples_analyzed")
        n_genes = _jv(s, "total_genes")
        n_degs = _jv(s, "sig_relaxed")
        n_up = _jv(s, "upregulated")
        n_down = _jv(s, "downregulated")

        # Custom metric cards
        mc = st.columns(5)
        for col, val, label in zip(mc,
            [n_samples, f"{n_genes:,}", f"{n_degs:,}", f"{n_up:,}", f"{n_down:,}"],
            ["Samples", "Genes Tested", "DEGs (FDR<0.05)", "\u2b06 Upregulated", "\u2b07 Downregulated"]):
            col.markdown(
                f'<div class="metric-card"><div class="val">{val}</div>'
                f'<div class="lab">{label}</div></div>',
                unsafe_allow_html=True)

        # Info line
        pca_out = _jv(s, "pca_outliers_flagged", 0)
        batch_cov = "Yes (design covariate)" if _jv(s, "batch_as_covariate") else "No"
        batch_pca = "Yes" if _jv(s, "batch_corrected_pca") else "No"
        st.caption(
            f"Batch in design: {batch_cov} \u00b7 "
            f"Batch-corrected PCA: {batch_pca} \u00b7 "
            f"PCA outliers flagged: {pca_out} \u00b7 "
            f"Doses analyzed: {_jv(s, 'doses_analyzed', 'N/A')}"
        )

        pc, vc = st.columns(2)
        pca_img = deg_dir / "Final_PCA_Plot.png"
        vol_img = deg_dir / "Volcano_Plot.png"
        if pca_img.exists():
            pca_cap = "PCA \u2014 batch-corrected VST" if _jv(s, "batch_corrected_pca") else "PCA \u2014 VST counts"
            pc.image(str(pca_img), caption=pca_cap, use_container_width=True)
        if vol_img.exists():
            vc.image(str(vol_img), caption="Volcano Plot", use_container_width=True)

        sig_csv = deg_dir / "Custom_Filtered_DEGs.csv"
        if sig_csv.exists():
            with st.expander(f"\U0001f4ca Top DEGs ({min(30, n_degs)} shown)", expanded=True):
                st.dataframe(pd.read_csv(sig_csv).head(30), use_container_width=True, hide_index=True)
            dl_cols = st.columns(4)
            for col, fname in zip(dl_cols, ["All_Results.csv", "Custom_Filtered_DEGs.csv",
                                             "Upregulated_Genes.csv", "Downregulated_Genes.csv"]):
                fp = deg_dir / fname
                if fp.exists():
                    col.download_button(
                        f"\u2b07 {fname.replace('_', ' ').replace('.csv', '')}",
                        fp.read_bytes(), file_name=fname, mime="text/csv",
                        key=f"dl_{run_dir.name}_{fname}",
                    )
    else:
        st.info("No DESeq2 results in this run.")

    st.divider()

    # --- DRomics ---
    st.markdown("#### DRomics / BMD")
    dp = dr_dir / "dromics_summary.json"
    if dp.exists():
        ds = json.loads(dp.read_text())
        n_genes_dr   = _jv(ds, "total_genes_final")
        n_responsive = _jv(ds, "genes_selected")
        n_models     = _jv(ds, "models_fitted")

        # ── Top metric cards ──────────────────────────────────────────────────
        metrics = [
            (f"{n_genes_dr:,}",   "Genes Analyzed"),
            (f"{n_responsive:,}", "Responsive Genes"),
            (f"{n_models:,}",     "Models Fitted"),
        ]
        if ds.get("bmd_performed") and ds.get("bmd_summary"):
            bms = ds["bmd_summary"]
            median_bmd = round(float(_jv(bms, "median_bmd", 0)), 4)
            metrics.append((f"{median_bmd} \u03bcM", "Median Gene BMD"))
        if ds.get("tpod_performed") and ds.get("tpod_summary"):
            tpod_val = round(float(_jv(ds["tpod_summary"], "tpod_value", 0)), 4)
            metrics.append((f"{tpod_val} \u03bcM", "tPOD (NTP 2018)"))

        mc = st.columns(len(metrics))
        for col, (val, label) in zip(mc, metrics):
            col.markdown(
                f'<div class="metric-card"><div class="val">{val}</div>'
                f'<div class="lab">{label}</div></div>',
                unsafe_allow_html=True)

        # ── Batch / PCA info line ─────────────────────────────────────────────
        st.caption(
            f"Samples: {_jv(ds, 'total_samples_final')} \u00b7 "
            f"Dose levels: {_jv(ds, 'dose_levels')} \u00b7 "
            f"PCA outliers: {_jv(ds, 'pca_outliers_flagged', 0)}"
        )

        # ── Model quality panel ───────────────────────────────────────────────
        mq = ds.get("model_quality")
        if mq:
            with st.expander("\U0001f4cf Model Quality Metrics", expanded=True):
                q1, q2, q3, q4 = st.columns(4)
                q1.metric("Median pseudo-R²",  f"{_jv(mq, 'median_pseudo_R2', '?')}")
                q2.metric("Median SNR",         f"{_jv(mq, 'median_SNR', '?')}")
                q3.metric("R² ≥ 0.5",           f"{_jv(mq, 'pct_R2_above_0.5', '?')}%")
                q4.metric("SNR ≥ 2",            f"{_jv(mq, 'pct_SNR_above_2', '?')}%")

                qa, qb = st.columns(2)
                with qa:
                    st.caption(
                        f"Mean pseudo-R²: {_jv(mq, 'mean_pseudo_R2', '?')} \u00b7 "
                        f"R² ≥ 0.3: {_jv(mq, 'pct_R2_above_0.3', '?')}% \u00b7 "
                        f"SNR ≥ 1: {_jv(mq, 'pct_SNR_above_1', '?')}% \u00b7 "
                        f"Median SDres: {_jv(mq, 'median_SDres', '?')}"
                    )
                with qb:
                    model_dist = mq.get("model_distribution") or {}
                    trend_dist = mq.get("trend_distribution") or {}
                    if model_dist:
                        model_str = " \u00b7 ".join(f"{k}: {v}" for k, v in model_dist.items())
                        st.caption(f"Models: {model_str}")
                    if trend_dist:
                        trend_str = " \u00b7 ".join(f"{k}: {v}" for k, v in trend_dist.items())
                        st.caption(f"Trends: {trend_str}")

        # ── NTP quality filter summary ────────────────────────────────────────
        if ds.get("bmd_performed") and ds.get("bmd_summary"):
            bms = ds["bmd_summary"]
            ntp = bms.get("ntp_quality_filters") or {}
            if ntp:
                with st.expander("\U0001f6e1\ufe0f NTP 2018 / EFSA Quality Filters Applied", expanded=False):
                    f1, f2, f3, f4 = st.columns(4)
                    f1.metric("Removed: BMD > max dose",
                               _jv(ntp, "n_removed_maxdose", 0),
                               help=f"max dose = {round(float(_jv(ntp, 'max_dose', 0)), 4)} µM")
                    f2.metric("Flagged: BMD extrapolated",
                               _jv(ntp, "n_extrapolated", 0),
                               help=f"threshold = lowest_dose / {_jv(ntp, 'extrap_factor', 10)}")
                    f3.metric("Removed: BMDU/BMDL ratio",
                               _jv(ntp, "n_removed_ratio", 0),
                               help=f"ratio threshold = {_jv(ntp, 'ratio_threshold', 40)}")
                    f4.metric("Removed: fold-change filter",
                               _jv(ntp, "n_fold_change_removed", 0),
                               help=f"min |log2FC| = {_jv(ntp, 'fold_change_min', 0)}")
                    st.caption(
                        f"Bootstrap CI available: {'Yes' if ntp.get('has_bootstrap_ci') else 'No'} \u00b7 "
                        f"Gene selection method: {_jv(ntp, 'select_method', '?')} \u00b7 "
                        f"Total BMDs before filters: {_jv(bms, 'total_bmds', '?')} \u00b7 "
                        f"After filters: {_jv(bms, 'total_bmds_filtered', '?')}"
                    )

        # ── tPOD highlight box ────────────────────────────────────────────────
        if ds.get("tpod_performed") and ds.get("tpod_summary"):
            tpod = ds["tpod_summary"]
            st.markdown(
                f'<div class="box-info" style="border-color:rgba(255,87,34,.4);background:rgba(255,87,34,.04)">'
                f'\U0001f3af <b>Transcriptomic Point of Departure (tPOD): '
                f'{round(float(_jv(tpod, "tpod_value", 0)), 4)} \u03bcM</b><br>'
                f'Pathway: <b>{_jv(tpod, "tpod_pathway", "?")}</b> '
                f'({_jv(tpod, "tpod_pathway_id", "")}) &nbsp;\u00b7&nbsp; '
                f'Category: {_jv(tpod, "tpod_category", "?")} &nbsp;\u00b7&nbsp; '
                f'{_jv(tpod, "tpod_n_genes", "?")} genes &nbsp;\u00b7&nbsp; '
                f'{_jv(tpod, "tpod_coverage", "?")}% coverage &nbsp;\u00b7&nbsp; '
                f'Method: NTP 2018 gene set &nbsp;\u00b7&nbsp; '
                f'{_jv(tpod, "n_pathways_tested", "?")} pathways tested<br>'
                f'<span style="color:#7a90a4;font-size:.75rem">'
                f'Gene-level median BMD: {round(float(_jv(tpod, "gene_level_median", 0)), 4)} \u03bcM '
                f'&nbsp;\u00b7&nbsp; Gene-level Q25: {round(float(_jv(tpod, "gene_level_q25", 0)), 4)} \u03bcM'
                f'</span></div>',
                unsafe_allow_html=True)

            # ── All tPOD methods comparison table ─────────────────────────────
            tpod_methods = tpod.get("tpod_methods") or {}
            if tpod_methods:
                with st.expander("\U0001f4cb All tPOD Methods Comparison", expanded=True):
                    rows = []
                    for key, m in tpod_methods.items():
                        if isinstance(m, dict):
                            rows.append({
                                "Method":       m.get("label", key),
                                "tPOD (µM)":    round(float(m.get("value", 0)), 6),
                                "Reference":    m.get("ref", ""),
                                "Genes Used":   m.get("n_genes_used", ""),
                                "Pathway":      m.get("pathway_name", ""),
                                "Coverage (%)": m.get("coverage_pct", ""),
                            })
                    if rows:
                        st.dataframe(
                            pd.DataFrame(rows).sort_values("tPOD (µM)"),
                            use_container_width=True, hide_index=True
                        )

            # ── Per-collection tPOD breakdown ──────────────────────────────────
            coll_tpod = tpod.get("collection_tpod") or {}
            if coll_tpod:
                with st.expander("\U0001f5c2\ufe0f tPOD by Gene Set Collection", expanded=False):
                    crows = []
                    for cat, info in coll_tpod.items():
                        if isinstance(info, dict):
                            crows.append({
                                "Collection":   cat,
                                "tPOD (µM)":    round(float(info.get("tpod_value", 0)), 6),
                                "Pathway":      info.get("pathway_name", ""),
                                "Pathway ID":   info.get("pathway_id", ""),
                                "Genes":        info.get("n_genes", ""),
                                "Coverage (%)": info.get("coverage_pct", ""),
                                "Pathways Tested": info.get("n_pathways", ""),
                            })
                    if crows:
                        st.dataframe(
                            pd.DataFrame(crows).sort_values("tPOD (µM)"),
                            use_container_width=True, hide_index=True
                        )

            # ── Top 5 pathways ─────────────────────────────────────────────────
            top5 = tpod.get("top5_pathways") or []
            if top5:
                with st.expander("\U0001f3c6 Top 5 Most Sensitive Pathways", expanded=False):
                    t5rows = []
                    for pw in top5:
                        if isinstance(pw, dict):
                            t5rows.append({
                                "Pathway":      pw.get("name", ""),
                                "ID":           pw.get("id", ""),
                                "Category":     pw.get("category", ""),
                                "Median BMD (µM)": round(float(pw.get("median_bmd", 0)), 6),
                                "Genes":        pw.get("n_genes", ""),
                                "Coverage (%)": pw.get("coverage", ""),
                            })
                    if t5rows:
                        st.dataframe(pd.DataFrame(t5rows), use_container_width=True, hide_index=True)

        # ── Plots row 1: dose-response + BMD distribution ─────────────────────
        c1, c2 = st.columns(2)
        for ci, fname, cap in [(c1, "dose_response_curves.png", "Dose-Response Curves"),
                                (c2, "bmd_distribution.png",    "BMD Distribution")]:
            fp = dr_dir / fname
            if fp.exists():
                ci.image(str(fp), caption=cap, use_container_width=True)

        # ── Plots row 2: pathway sensitivity + top 25 genes ───────────────────
        pw_plot = dr_dir / "pathway_sensitivity.png"
        top25_img = dr_dir / "top25_sensitive_genes.png"
        if pw_plot.exists() or top25_img.exists():
            pc1, pc2 = st.columns(2)
            if pw_plot.exists():
                pc1.image(str(pw_plot), caption="Pathway Sensitivity (NTP 2018 tPOD)",
                          use_container_width=True)
            if top25_img.exists():
                pc2.image(str(top25_img), caption="Top 25 Most Sensitive Genes (by BMD)",
                          use_container_width=True)

        # ── Pathway table ─────────────────────────────────────────────────────
        pw_csv = dr_dir / "pathway_bmd_summary.csv"
        if pw_csv.exists():
            with st.expander("\U0001f4ca Pathway BMD Summary (NTP 2018)", expanded=False):
                st.dataframe(pd.read_csv(pw_csv).head(50), use_container_width=True, hide_index=True)

        # ── Model quality table ───────────────────────────────────────────────
        mq_csv = dr_dir / "model_quality_metrics.csv"
        if mq_csv.exists():
            with st.expander("\U0001f4c8 Model Fit Quality (pseudo-R² / SNR per gene)", expanded=False):
                st.dataframe(pd.read_csv(mq_csv).head(50), use_container_width=True, hide_index=True)

        # ── Downloads ─────────────────────────────────────────────────────────
        dl_files = [
            "bmd_results.csv",
            "top50_sensitive_genes.csv",
            "dose_response_models.csv",
            "model_quality_metrics.csv",
            "model_fit_summary.csv",
            "pathway_bmd_summary.csv",
            "selected_responding_genes.csv",
        ]
        available = [f for f in dl_files if (dr_dir / f).exists()]
        if available:
            dl_cols = st.columns(min(len(available), 4))
            for i, fname in enumerate(available):
                fp = dr_dir / fname
                dl_cols[i % 4].download_button(
                    f"\u2b07 {fname.replace('_', ' ').replace('.csv', '')}",
                    fp.read_bytes(), file_name=fname, mime="text/csv",
                    key=f"dl_{run_dir.name}_{fname}",
                )
    else:
        st.info("No DRomics results in this run.")


# ═════════════════════════════════════════════════════════════════════════════
#  PAGE CONFIG + CSS
# ═════════════════════════════════════════════════════════════════════════════

st.set_page_config(
    page_title="ARACRA Pipeline",
    page_icon="\U0001f9ec",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;500;600&family=DM+Sans:wght@300;400;500;600;700&display=swap');
html,body,[class*="css"]{font-family:'DM Sans',sans-serif;-webkit-font-smoothing:antialiased}
/* ── Header ── */
.header{background:linear-gradient(135deg,#0a1118 0%,#0f1d2b 60%,#0d2234 100%);border-left:3px solid #00d4aa;border-radius:0 10px 10px 0;padding:1.3rem 2rem;margin-bottom:1.6rem;box-shadow:0 4px 20px rgba(0,0,0,.25)}
.header h1{font-family:'JetBrains Mono',monospace;color:#00d4aa;font-size:1.5rem;margin:0;letter-spacing:-.3px;font-weight:600}
.header p{color:#7a90a4;font-size:.82rem;margin:5px 0 0 0;letter-spacing:.2px}
/* ── Section labels ── */
.sec{font-family:'JetBrains Mono',monospace;font-size:.65rem;color:#00b894;text-transform:uppercase;letter-spacing:2.5px;margin:1.2rem 0 .4rem 0;border-bottom:1px solid rgba(0,212,170,.15);padding-bottom:4px}
/* ── Token chips ── */
.tok{display:inline-block;background:rgba(0,212,170,.06);border:1px solid rgba(0,212,170,.2);color:#00d4aa;font-family:'JetBrains Mono',monospace;font-size:.65rem;padding:2px 8px;border-radius:4px;margin:2px;transition:background .15s}
.tok:hover{background:rgba(0,212,170,.12)}
/* ── Status indicators ── */
.s-idle{font-family:'JetBrains Mono',monospace;color:#6c7a89;font-size:.95rem;font-weight:600}
.s-run{font-family:'JetBrains Mono',monospace;color:#4da6ff;font-size:.95rem;font-weight:600;animation:pulse 2s ease-in-out infinite}
.s-ok{font-family:'JetBrains Mono',monospace;color:#00d4aa;font-size:.95rem;font-weight:600}
.s-err{font-family:'JetBrains Mono',monospace;color:#ff6b6b;font-size:.95rem;font-weight:600}
.s-stp{font-family:'JetBrains Mono',monospace;color:#ffa726;font-size:.95rem;font-weight:600}
.s-review{font-family:'JetBrains Mono',monospace;color:#b388ff;font-size:.95rem;font-weight:600}
@keyframes pulse{0%,100%{opacity:1}50%{opacity:.6}}
/* ── Tool check ── */
.tok-ok{color:#00d4aa;font-family:'JetBrains Mono',monospace;font-size:.75rem;display:block;padding:1px 0}
.tok-err{color:#ff6b6b;font-family:'JetBrains Mono',monospace;font-size:.75rem;display:block;padding:1px 0}
/* ── Info/warn boxes ── */
.box-info{background:rgba(0,212,170,.04);border:1px solid rgba(0,212,170,.18);border-radius:8px;padding:.55rem .9rem;color:#00d4aa;font-size:.8rem;margin:4px 0;backdrop-filter:blur(4px)}
.box-warn{background:rgba(255,167,38,.04);border:1px solid rgba(255,167,38,.2);border-radius:8px;padding:.55rem .9rem;color:#ffa726;font-size:.8rem;margin:4px 0;backdrop-filter:blur(4px)}
.box-err{background:rgba(255,107,107,.04);border:1px solid rgba(255,107,107,.2);border-radius:8px;padding:.55rem .9rem;color:#ff6b6b;font-size:.8rem;margin:4px 0}
.box-sys{background:rgba(15,25,35,.6);border:1px solid rgba(30,45,61,.6);border-radius:8px;padding:.6rem .9rem;font-family:'JetBrains Mono',monospace;font-size:.7rem;color:#7a90a4;margin:4px 0}
.tool-group{background:rgba(15,25,35,.5);border:1px solid rgba(30,45,61,.5);border-radius:8px;padding:.7rem .9rem;margin:4px 0}
/* ── Metric cards ── */
.metric-card{background:rgba(0,212,170,.03);border:1px solid rgba(0,212,170,.12);border-radius:10px;padding:.8rem 1rem;text-align:center;margin:2px 0}
.metric-card .val{font-family:'JetBrains Mono',monospace;font-size:1.4rem;font-weight:600;color:#00d4aa}
.metric-card .lab{font-size:.72rem;color:#7a90a4;margin-top:2px;text-transform:uppercase;letter-spacing:1px}
/* ── Phase indicator ── */
.phase-badge{display:inline-block;font-family:'JetBrains Mono',monospace;font-size:.65rem;padding:3px 10px;border-radius:12px;letter-spacing:1px;font-weight:500}
.phase-1{background:rgba(77,166,255,.12);color:#4da6ff;border:1px solid rgba(77,166,255,.25)}
.phase-2{background:rgba(179,136,255,.12);color:#b388ff;border:1px solid rgba(179,136,255,.25)}
.phase-qc{background:rgba(255,167,38,.12);color:#ffa726;border:1px solid rgba(255,167,38,.25)}
/* ── Better scrollbar ── */
::-webkit-scrollbar{width:6px}
::-webkit-scrollbar-track{background:transparent}
::-webkit-scrollbar-thumb{background:rgba(0,212,170,.2);border-radius:3px}
::-webkit-scrollbar-thumb:hover{background:rgba(0,212,170,.35)}
/* ── Tab styling ── */
.stTabs [data-baseweb="tab-list"]{gap:8px}
.stTabs [data-baseweb="tab"]{font-family:'DM Sans',sans-serif;font-weight:500;font-size:.85rem}
</style>
""", unsafe_allow_html=True)


# ═════════════════════════════════════════════════════════════════════════════
#  SESSION STATE
# ═════════════════════════════════════════════════════════════════════════════

_defaults = {
    "accessions":       [],
    "meta_df":          None,
    "meta_path":        None,
    "show_results":     False,
    "chemicals":        [],
    "tools":            None,
    "fdr_strict":       0.01,
    "fdr_relaxed":      0.05,
    "log2fc":           1.0,
    "read_thresh":      1_000_000,
    "dr_fdr":           0.05,
    "dr_criterion":     "AICc",
    "perform_bmd":      True,
    "bmd_bootstrap":    False,
    "dr_select_method": "quadratic",
    "dr_transfo":       "auto",
    "dr_bmd_z":         1.0,
    "dr_bmd_x":         10.0,
    "dr_niter":         1000,
    "bmdu_bmdl_ratio":  40.0,
    "bmd_max_dose_filter": True,
    "bmd_extrap_factor": 10.0,
    "fold_change_min":  0.0,
    "use_msigdb":       True,
    "count_matrix_path": None,
    "count_matrix_df":  None,
    "exclude_samples":  [],
    "col_map":          {},     # user column name → standard name mapping
    "_phase1_transferred": False,
    "_current_run_id":  None,   # tracks which run the current session state belongs to
    "_analysis_qc_run_id": None,  # tracks which run_name the QC was run for
}
for k, v in _defaults.items():
    if k not in st.session_state:
        st.session_state[k] = v

# Auto-populate chemicals list from metadata if available
if not st.session_state.chemicals and st.session_state.meta_df is not None:
    _df = st.session_state.meta_df
    if hasattr(_df, 'columns') and "Treatment" in _df.columns:
        st.session_state.chemicals = sorted(_df["Treatment"].unique().tolist())

# Reconcile running state
state   = _load_state()
running = _is_running()

# If state says running but process is dead → detect final status
if state.get("status") == "running" and not running:
    final = detect_completion()
    state["status"] = final or "failed"
    _save_state(state)
    _clear_pid()
    # Force a clean rerun so the full page re-evaluates with updated state.
    # Without this, the QC gate and inline results don't see the new output files
    # because they rendered in the same pass that discovered completion.
    st.rerun()

# On fresh page load (no prior interaction), don't auto-show old results
# Reset UI flags on genuine fresh page load only (not hot-reload)
if "_app_initialized" not in st.session_state:
    st.session_state._app_initialized = True
    st.session_state.show_results = False
    st.session_state._show_prev_qc = False
    # Clear stale analysis state from previous session to prevent
    # old QC/results from bleeding into the DESeq2/DRomics tab
    st.session_state.exclude_samples = []
    st.session_state._qc_valid = False
    st.session_state._analysis_qc_run_id = None
    st.session_state._current_run_id = None
    st.session_state._phase1_transferred = False

# Also reset stale state when pipeline isn't running and status is from a previous session
if not running and state.get("status") in ("completed", "failed", "stopped"):
    # Check if the state is stale (older than current session)
    state_time = state.get("start_time", "")
    if state_time:
        try:
            from datetime import timedelta
            elapsed = datetime.now() - datetime.fromisoformat(state_time)
            # If the last run finished more than 1 hour ago, consider it stale
            if elapsed > timedelta(hours=1):
                state["_stale"] = True
        except Exception:
            pass

STATUS_MAP = {
    "idle":      ("\u25cb  IDLE",      "s-idle"),
    "running":   ("\u25cf  RUNNING",   "s-run"),
    "completed": ("\u2714  COMPLETED", "s-ok"),
    "failed":    ("\u2718  FAILED",    "s-err"),
    "stopped":   ("\u23f9  STOPPED",   "s-stp"),
    "review":    ("\u2714  REVIEW QC", "s-review"),
}

# ═════════════════════════════════════════════════════════════════════════════
#  HEADER
# ═════════════════════════════════════════════════════════════════════════════

st.markdown("""
<div class="header">
  <h1>\U0001f9ec ARACRA Pipeline</h1>
  <p>RNA-seq &amp; TempO-Seq Analysis &nbsp;\u2502&nbsp;
     SRA \u2192 fastp \u2192 STAR/HISAT2 \u2192 QC \u2192 featureCounts/Salmon \u2192 DESeq2 \u2192 DRomics/BMD
     &nbsp;\u2502&nbsp; <span style="opacity:.6">v2.0</span></p>
</div>""", unsafe_allow_html=True)


# ═════════════════════════════════════════════════════════════════════════════
#  SIDEBAR
# ═════════════════════════════════════════════════════════════════════════════

with st.sidebar:
    sysinfo = get_sysinfo()
    st.markdown(
        f'<div class="box-sys">'
        f'\U0001f5a5 <b>System</b><br>'
        f'&nbsp;&nbsp;RAM: {sysinfo["ram_total"]} ({sysinfo["ram_avail"]} free)<br>'
        f'&nbsp;&nbsp;CPUs: {sysinfo["cpus"]}<br>'
        f'&nbsp;&nbsp;Disk: {sysinfo["disk"]}<br>'
        f'&nbsp;&nbsp;Aligner: {RECOMMENDED_ALIGNER.upper()}'
        f'</div>',
        unsafe_allow_html=True,
    )

    if not SETUP_DONE:
        st.markdown(
            '<div class="box-warn">\u26a0 .env not found \u2014 run setup.sh first.</div>',
            unsafe_allow_html=True,
        )

    st.markdown('<p class="sec">\U0001f4c2 Reference Files</p>', unsafe_allow_html=True)
    star_index   = st.text_input("STAR index",       value=STAR_INDEX,       key="star_idx")
    hisat2_index = st.text_input("HISAT2 index",     value=HISAT2_INDEX,     key="h2_idx")
    salmon_index = st.text_input("Salmon index",     value=SALMON_INDEX,     key="sal_idx")
    gtf          = st.text_input("GTF annotation",   value=GTF_PATH,         key="gtf_in")
    bed          = st.text_input("RefSeq BED",       value=BED_PATH,         key="bed_in")
    hk_bed       = st.text_input("Housekeeping BED", value=HK_BED_PATH,     key="hk_bed_in")
    refflat      = st.text_input("refFlat",          value=REFFLAT_PATH,     key="refflat_in")
    screen_conf  = st.text_input("FastQ Screen conf", value=SCREEN_CONF_PATH, key="screen_in")

    st.markdown('<p class="sec">\U0001f4c1 Directories</p>', unsafe_allow_html=True)
    outdir  = st.text_input("Output directory", value=DEFAULT_OUT,  key="out_in")
    workdir = st.text_input("Work directory",   value=DEFAULT_WORK, key="wk_in")

    st.markdown('<p class="sec">\U0001f5a5 Executor</p>', unsafe_allow_html=True)
    profile = st.selectbox("Profile", ["standard", "slurm"])
    layout  = st.selectbox("Library layout", ["SE", "PE"])
    resume  = st.checkbox("Resume previous run", value=True)

    st.markdown('<p class="sec">\U0001f3af Platform</p>', unsafe_allow_html=True)
    platform = st.selectbox(
        "Sequencing platform",
        ["RNA-seq (whole transcriptome)", "TempO-Seq (targeted)"],
        help="RNA-seq: align to whole genome with STAR/HISAT2 + featureCounts.\n"
             "TempO-Seq: align to BioSpyder probe reference + probe-level counting."
    )
    is_temposeq = platform.startswith("TempO")
    temposeq_manifest_path = None
    if is_temposeq:
        st.info("🎯 TempO-Seq mode: aligns to probe reference, not whole genome. "
                "Upload BioSpyder manifest CSV below.")
        temposeq_manifest = st.file_uploader(
            "BioSpyder manifest CSV", type=["csv"], key="temposeq_manifest"
        )
        if temposeq_manifest:
            mp = Path(outdir) / "temposeq_manifest.csv"
            mp.parent.mkdir(parents=True, exist_ok=True)
            mp.write_bytes(temposeq_manifest.getvalue())
            temposeq_manifest_path = str(mp)
            st.success(f"✔ Manifest saved: {mp.name}")

    st.markdown('<p class="sec">\U0001f52c Downstream</p>', unsafe_allow_html=True)
    run_deg     = st.checkbox("DESeq2", value=True)
    run_dromics = st.checkbox("DRomics / BMD",    value=True)

    st.markdown('<p class="sec">\U0001f527 System Check</p>', unsafe_allow_html=True)
    if st.button("\u25b6 Check Tools", use_container_width=True):
        st.cache_data.clear()
        st.session_state.tools = check_tools()
    if st.session_state.tools:
        for name, info in st.session_state.tools.items():
            cls  = "tok-ok" if info["ok"] else "tok-err"
            icon = "\u2714" if info["ok"] else "\u2718"
            st.markdown(
                f'<span class="{cls}">{icon} {name}: {info["ver"][:40]}</span>',
                unsafe_allow_html=True,
            )
        r_ok, r_miss = check_r_packages()
        if r_ok:
            st.markdown('<span class="tok-ok">\u2714 R packages (12/12)</span>', unsafe_allow_html=True)
        else:
            st.markdown(
                f'<span class="tok-err">\u2718 R missing: {", ".join(r_miss)}</span>',
                unsafe_allow_html=True,
            )

    if running:
        st.divider()
        st.markdown(f'<span class="tok-ok">\u2699 PID: {_load_pid()}</span>', unsafe_allow_html=True)
        st.caption("Detached \u2014 safe to close browser")


# ═════════════════════════════════════════════════════════════════════════════
#  TABS
# ═════════════════════════════════════════════════════════════════════════════

tab_full, tab_analysis, tab_params, tab_results = st.tabs([
    "\U0001f680 Full Pipeline", "\U0001f52c DESeq2 / DRomics", "\u2699 Parameters", "\U0001f4ca Results",
])


# ─────────────────────────────────────────────────────────────────────────────
#  Shared: build Nextflow params dict
# ─────────────────────────────────────────────────────────────────────────────


def _build_params(
    mode: str,
    *,
    aligner: str = "star",
    fastq_source: str = "Download from SRA",
    fastq_dir_val=None,
    trimmed_dir_val=None,
    treatment: str = "Treatment",
    control: str = "Control",
    run_name: str = "",
    run_screen: bool = True,
    run_read_dist: bool = False,
    run_qualimap: bool = False,
    use_fc: bool = True,
    use_salmon: bool = False,
    use_pseudo: bool = False,
    use_strand: bool = True,
    use_picard: bool = False,
    use_cov_hk: bool = True,
    use_cov_full: bool = False,
    no_aligner: bool = False,
    phase: str = "full",
    count_matrix: str = None,
    exclude_samples: list = None,
    platform_temposeq: bool = False,
    temposeq_manifest_path: str = None,
) -> dict:
    """Build params dict for Nextflow -params-file."""
    is_phase1 = (phase == "preprocess")
    is_phase2 = (phase == "analysis")

    params = {
        "metadata":           st.session_state.meta_path or "",
        "layout":             layout,
        "outdir":             outdir,
        "aligner":            "none" if no_aligner else aligner.lower(),
        "star_index":         star_index,
        "hisat2_index":       hisat2_index,
        "salmon_index":       salmon_index,
        "gtf":                gtf,
        "bed":                bed,
        "housekeeping_bed":   hk_bed,
        "refflat":            refflat,
        "fastq_screen_conf":  screen_conf,
        "fastq_dir":          fastq_dir_val,
        "trimmed_dir":        trimmed_dir_val,
        "run_name":           run_name,
        "treatment":          treatment,
        "control":            control,
        # Phase 1: preprocess only
        "run_download":       (fastq_source == "Download from SRA") and not is_phase2,
        "run_fastp":          (fastq_source != "Local FASTQs (already trimmed)") and not is_phase2,
        "run_fastq_screen":   run_screen and not is_phase2,
        "run_star":           (aligner.lower() == "star") and not is_phase2,
        "run_hisat2":         (aligner.lower() == "hisat2") and not is_phase2,
        "run_samtools_stats": (not no_aligner) and not is_phase2,
        "run_strandedness":   use_strand and not is_phase2 and not platform_temposeq,
        "run_coverage_hk":    use_cov_hk and not is_phase2 and not platform_temposeq,
        "run_coverage_full":  use_cov_full and not is_phase2 and not platform_temposeq,
        "run_read_distribution": run_read_dist and (not no_aligner) and not is_phase2 and not platform_temposeq,
        "run_picard":         use_picard and not is_phase2 and not platform_temposeq,
        "run_qualimap":       run_qualimap and (not no_aligner) and not is_phase2 and not platform_temposeq,
        "run_featurecounts":  use_fc and not is_phase2 and not platform_temposeq,
        "run_salmon":         use_salmon and not is_phase2 and not platform_temposeq,
        "run_salmon_pseudo":  use_pseudo and not is_phase2 and not platform_temposeq,
        "run_merge_counts":   not is_phase2,
        "run_multiqc":        not is_phase2,
        # Phase 2: analysis only (disabled during Phase 1)
        "run_deg":            run_deg and not is_phase1,
        "run_dromics":        run_dromics and not is_phase1,
        "fdr_strict":         st.session_state.fdr_strict,
        "fdr_relaxed":        st.session_state.fdr_relaxed,
        "log2fc_threshold":   st.session_state.log2fc,
        "read_threshold":     st.session_state.read_thresh,
        "dr_fdr":             st.session_state.dr_fdr,
        "dr_criterion":       st.session_state.dr_criterion,
        "perform_bmd":        st.session_state.perform_bmd,
        "bmd_bootstrap":      st.session_state.bmd_bootstrap,
        "dr_select_method":   st.session_state.dr_select_method,
        "dr_transfo_method":  st.session_state.dr_transfo,
        "dr_bmd_z":           st.session_state.dr_bmd_z,
        "dr_bmd_x":           st.session_state.dr_bmd_x,
        "dr_niter":           st.session_state.dr_niter,
        "bmdu_bmdl_ratio":    st.session_state.bmdu_bmdl_ratio,
        "bmd_max_dose_filter": st.session_state.bmd_max_dose_filter,
        "bmd_extrap_factor":  st.session_state.bmd_extrap_factor,
        "fold_change_min":    st.session_state.fold_change_min,
        "exclude_samples":    ",".join(exclude_samples or []),
        # TempO-Seq platform
        "platform":           "temposeq" if platform_temposeq else "rnaseq",
        "temposeq_manifest":  temposeq_manifest_path or "",
    }

    # Direct / Phase 2: pass count matrix
    if count_matrix:
        params["count_matrix"] = count_matrix

    return params


def _launch_with_params(params_dict: dict, run_name: str, mode: str = "full", phase: str = "full"):
    """Write params JSON and launch Nextflow."""
    Path(outdir).mkdir(parents=True, exist_ok=True)
    params_file = Path(outdir) / "nextflow_params.json"
    params_file.write_text(json.dumps(params_dict, indent=2))

    cmd = [
        "nextflow", "run", str(MAIN_NF),
        "-profile", profile,
        "-w", workdir,
        "-ansi-log", "false",
        "-params-file", str(params_file),
    ]
    if resume:
        cmd.append("-resume")

    pid = launch_pipeline(cmd, str(PIPELINE_DIR))
    _save_pid(pid)
    _save_state({
        "status":     "running",
        "start_time": datetime.now().isoformat(),
        "pid":        pid,
        "mode":       mode,
        "phase":      phase,
        "run_name":   run_name,
    })
    # Track the active run in session state so other tabs can detect staleness
    st.session_state._current_run_id = f"{run_name}_{datetime.now().isoformat()}"


# ─────────────────────────────────────────────────────────────────────────────
#  Shared: status panel + log viewer
# ─────────────────────────────────────────────────────────────────────────────


def _show_status_panel(key_prefix: str):
    """Show status, progress bars, controls, and log."""
    cur_running = _is_running()  # re-check live, don't use stale global
    cur_state = _load_state()    # re-read from disk for latest
    status = "running" if cur_running else cur_state.get("status", "idle")

    # If process just finished, update state
    if cur_state.get("status") == "running" and not cur_running:
        final = detect_completion()
        status = final or "failed"
        cur_state["status"] = status
        _save_state(cur_state)
        _clear_pid()
        # Force a clean rerun so QC gate / inline results see the new output files.
        # This is safe because on the next pass, status is no longer "running",
        # so this block won't fire again.
        st.rerun()

    # For display: if not running, show idle (clean slate) unless just completed
    is_fresh_completion = False
    if status in ("completed", "failed", "stopped") and not cur_running:
        start_time = cur_state.get("start_time", "")
        if start_time:
            try:
                elapsed = datetime.now() - datetime.fromisoformat(start_time)
                # Show completion status for 10 minutes, then revert to idle
                is_fresh_completion = elapsed.total_seconds() < 600
            except Exception:
                pass

    display_status = status if (cur_running or is_fresh_completion) else "idle"
    lbl, css = STATUS_MAP.get(display_status, ("\u25cb  IDLE", "s-idle"))
    st.markdown(f'<p class="{css}">{lbl}</p>', unsafe_allow_html=True)

    if cur_running:
        elapsed_str = _elapsed(cur_state)
        if cur_state.get("run_name"):
            st.caption(f"Run: **{cur_state['run_name']}** \u00b7 Elapsed: {elapsed_str}")
        else:
            st.caption(f"Elapsed: {elapsed_str}")
        phase = cur_state.get("phase", "full")
        phase_badges = {
            "preprocess": ('<span class="phase-badge phase-1">PHASE 1 \u2014 PREPROCESSING</span>', "Alignment + QC + Counts"),
            "analysis":   ('<span class="phase-badge phase-2">PHASE 2 \u2014 ANALYSIS</span>', "DESeq2 / DRomics"),
            "qc_check":   ('<span class="phase-badge phase-qc">QC CHECK</span>', "PCA outlier detection"),
        }
        badge, desc = phase_badges.get(phase, ("", ""))
        if badge:
            st.markdown(f'{badge} &nbsp;<span style="color:#7a90a4;font-size:.78rem">{desc}</span>',
                        unsafe_allow_html=True)
    elif is_fresh_completion and cur_state.get("run_name"):
        elapsed_str = _elapsed(cur_state)
        st.caption(f"Run: **{cur_state['run_name']}** \u00b7 Duration: {elapsed_str}")
    elif not cur_running and status in ("completed", "failed") and cur_state.get("run_name"):
        # Stale result — show as subtle info
        with st.expander(f"Last run: {cur_state['run_name']}", expanded=False):
            elapsed_str = _elapsed(cur_state)
            st.caption(f"Status: {status} \u00b7 Duration: {elapsed_str}")
            if st.button("Clear old status", key=f"clear_status_{key_prefix}"):
                _save_state({"status": "idle"})
                st.rerun()

    # Progress bars
    prog = parse_progress()
    if prog:
        st.markdown('<p class="sec">\u26a1 Progress</p>', unsafe_allow_html=True)
        for proc, info in prog.items():
            pct = min(info["pct"], 100) / 100
            label = PROCESS_LABELS.get(proc, proc)
            st.progress(pct, text=f"{label} \u2014 {info['done']}/{info['total']}")

    # Controls
    st.markdown('<p class="sec">\U0001f3ae Controls</p>', unsafe_allow_html=True)
    c1, c2 = st.columns(2)
    with c1:
        if st.button("\u23f9 Stop", disabled=not cur_running, use_container_width=True, key=f"stop_{key_prefix}"):
            stop_pipeline()
            _save_state({**cur_state, "status": "stopped"})
            st.rerun()
    with c2:
        if st.button("\u21ba Reset", disabled=cur_running, use_container_width=True, key=f"reset_{key_prefix}"):
            _save_state({"status": "idle"})
            _clear_pid()
            # Clear all session state to prevent stale data from previous runs
            for k in list(st.session_state.keys()):
                if k not in ("tools", "_app_initialized"):  # keep tool check cache + init flag
                    del st.session_state[k]
            st.rerun()

    # Log — only show when running or freshly completed
    if cur_running or is_fresh_completion:
        with st.expander("\U0001f4dc Pipeline Log" + (" (live)" if cur_running else ""), expanded=cur_running):
            st.text_area("", value=log_tail(80), height=250, disabled=True,
                         key=f"log_{key_prefix}", label_visibility="collapsed")
            log_c1, log_c2 = st.columns([1, 3])
            if LOG_FILE.exists():
                log_c1.download_button(
                    "\u2b07 Full log", LOG_FILE.read_bytes(),
                    file_name="pipeline.log", key=f"dl_log_{key_prefix}",
                )

        if is_fresh_completion and status == "completed":
            _completed_phase = cur_state.get("phase", "")
            if _completed_phase == "qc_check":
                st.markdown(
                    '<div class="box-info">\u2714 QC check complete. '
                    'Review the outlier results below.</div>',
                    unsafe_allow_html=True)
            elif _completed_phase == "analysis":
                st.markdown(
                    '<div class="box-info">\u2714 Analysis complete. '
                    'Results are shown below.</div>',
                    unsafe_allow_html=True)
            else:
                st.markdown(
                    '<div class="box-info">\u2714 Pipeline finished successfully. '
                    'Check the <b>Results</b> tab or proceed to <b>Phase 2</b>.</div>',
                    unsafe_allow_html=True)
        elif is_fresh_completion and status == "failed":
            st.markdown(
                '<div class="box-err">\u2718 Pipeline failed. '
                'Check the log for errors. You can <b>Resume</b> after fixing the issue.</div>',
                unsafe_allow_html=True)

        if cur_running:
            refresh_sec = st.select_slider(
                "Auto-refresh interval",
                options=[2, 3, 5, 10, 30],
                value=5, key=f"refresh_{key_prefix}",
            )
            time.sleep(refresh_sec)
            st.rerun()


# ─────────────────────────────────────────────────────────────────────────────
#  Shared: validation
# ─────────────────────────────────────────────────────────────────────────────


def _validate_full(
    aligner, no_aligner, pseudo_only, use_fc, use_salmon, use_pseudo,
    use_strand, use_picard, run_qualimap, use_cov_hk, use_cov_full,
    run_read_dist, treatment, control,
    fastq_source, fastq_dir_val, trimmed_dir_val,
    is_temposeq=False, temposeq_manifest_path=None,
) -> tuple[list, list]:
    """Return (issues, warnings) for the full pipeline."""
    issues, warnings = [], []

    qc_needs_bam = (
        use_strand or use_picard or (run_qualimap and not no_aligner)
        or use_cov_hk or use_cov_full or (run_read_dist and not no_aligner)
    )

    if not MAIN_NF.exists():
        issues.append("main.nf not found")
    if not st.session_state.accessions:
        issues.append("No SRR accessions loaded")
    if not st.session_state.meta_path:
        issues.append("No metadata uploaded")
    if no_aligner and not use_pseudo:
        issues.append("No aligner \u2014 enable Salmon pseudomapping")

    # TempO-Seq: needs manifest, does NOT need genome index/GTF
    if is_temposeq:
        if not temposeq_manifest_path or not Path(temposeq_manifest_path).exists():
            issues.append("TempO-Seq mode requires a BioSpyder manifest CSV \u2014 upload in sidebar")
    else:
        if not Path(gtf).exists():
            issues.append(f"GTF not found: {gtf}")
        if not pseudo_only or qc_needs_bam:
            if aligner == "STAR" and not Path(star_index).is_dir():
                issues.append(f"STAR index not found: {star_index}")
            if aligner == "HISAT2" and not Path(hisat2_index + ".1.ht2").exists():
                issues.append(f"HISAT2 index not found: {hisat2_index}")
            if aligner == "STAR" and 0 < SYSTEM_RAM_GB < 32:
                issues.append(f"STAR needs \u226532 GB RAM (you have {SYSTEM_RAM_GB:.0f} GB)")
    if use_salmon and not Path(salmon_index).is_dir():
        issues.append(f"Salmon index not found: {salmon_index}")
    if use_pseudo and not Path(salmon_index).is_dir():
        issues.append(f"Salmon index not found (pseudomapping): {salmon_index}")
    if treatment and control and treatment == control:
        issues.append("Treatment = Control")
    if fastq_source == "Local FASTQs (raw)" and fastq_dir_val and not Path(fastq_dir_val).is_dir():
        issues.append(f"FASTQ dir not found: {fastq_dir_val}")
    if fastq_source == "Local FASTQs (already trimmed)" and trimmed_dir_val and not Path(trimmed_dir_val).is_dir():
        issues.append(f"Trimmed FASTQ dir not found: {trimmed_dir_val}")

    if no_aligner:
        warnings.append("\u26a1 No-alignment mode \u2014 post-alignment QC disabled.")

    if is_temposeq:
        warnings.append("\U0001f3af TempO-Seq mode \u2014 aligning to probe reference. "
                        "Strandedness, gene body coverage, Picard, Qualimap, and featureCounts are skipped.")

    return issues, warnings


# ═════════════════════════════════════════════════════════════════════════════
#  TAB 1 — FULL PIPELINE (two-phase)
# ═════════════════════════════════════════════════════════════════════════════

with tab_full:
    L, R = st.columns([3, 2], gap="large")

    with L:
        # ── Input section ──
        st.markdown('<p class="sec">\U0001f4e5 SRA Accessions / Metadata</p>', unsafe_allow_html=True)
        method = st.radio(
            "Input method",
            ["Paste accessions", "Upload metadata (xlsx)"],
            horizontal=True, label_visibility="collapsed",
        )

        if method == "Paste accessions":
            raw = st.text_area("", height=120, placeholder="SRR29684827\nSRR29684828\n...")
            if raw.strip():
                valid = [a.strip() for a in raw.splitlines()
                         if a.strip() and re.match(r"^[SED]RR\d+$", a.strip())]
                invalid = [a.strip() for a in raw.splitlines()
                           if a.strip() and not re.match(r"^[SED]RR\d+$", a.strip())]
                st.session_state.accessions = valid
                if invalid:
                    st.warning(f"Ignored: {', '.join(invalid)}")
                if valid and not st.session_state.meta_path:
                    Path(outdir).mkdir(parents=True, exist_ok=True)
                    mp = Path(outdir) / "metadata_input.csv"
                    stub = pd.DataFrame({
                        "Sample_Name": valid,
                        "Treatment":   ["Treatment"] * len(valid),
                        "Dose":        [1.0] * len(valid),
                        "Batch":       [1] * len(valid),
                        "Type":        ["treatment"] * len(valid),
                    })
                    stub.to_csv(mp, index=False)
                    st.session_state.meta_path = str(mp)
                    st.info("Stub metadata created \u2014 upload a proper metadata file before running.")
        else:
            mf = st.file_uploader("Upload metadata", type=["xlsx", "xls", "csv"], key="meta_full")
            if mf:
                try:
                    df = pd.read_csv(mf) if mf.name.endswith(".csv") else pd.read_excel(mf)
                    srr_col = None
                    for c in df.columns:
                        if df[c].astype(str).str.match(r"^[SED]RR\d+$").any():
                            srr_col = c
                            break
                    if srr_col:
                        st.session_state.accessions = df[srr_col].dropna().astype(str).str.strip().tolist()
                    Path(outdir).mkdir(parents=True, exist_ok=True)
                    # Save as-is — column mapping happens in DESeq2/DRomics tab
                    if mf.name.endswith(".csv"):
                        mp = Path(outdir) / "metadata_input.csv"
                        df.to_csv(str(mp), index=False)
                    else:
                        mp = Path(outdir) / "metadata_input.xlsx"
                        mf.seek(0)
                        mp.write_bytes(mf.read())
                    st.session_state.meta_path = str(mp)
                    st.session_state.meta_df = df
                    st.success(f"\u2714 {len(df)} rows, {len(st.session_state.accessions)} SRR IDs found")
                    with st.expander("Preview"):
                        st.dataframe(df.head(10), use_container_width=True, hide_index=True)
                except Exception as e:
                    st.error(f"Error: {e}")

        if st.session_state.accessions:
            n = len(st.session_state.accessions)
            tags = "".join(f'<span class="tok">{a}</span>' for a in st.session_state.accessions[:20])
            if n > 20:
                tags += f'<span class="tok">+{n - 20} more</span>'
            st.markdown(tags, unsafe_allow_html=True)

        # ── FASTQ source ──
        st.markdown('<p class="sec">\U0001f4c1 FASTQ Source</p>', unsafe_allow_html=True)
        fastq_source = st.radio(
            "Where are your reads?",
            ["Download from SRA", "Local FASTQs (raw)", "Local FASTQs (already trimmed)"],
            horizontal=True, key="fastq_source",
        )
        fastq_dir_val = trimmed_dir_val = None
        if fastq_source == "Local FASTQs (raw)":
            fastq_dir_val = st.text_input(
                "Path to raw FASTQ directory",
                value=str(Path(outdir) / "data"), key="fastq_dir_input",
            )
        elif fastq_source == "Local FASTQs (already trimmed)":
            trimmed_dir_val = st.text_input(
                "Path to trimmed FASTQ directory",
                value=str(Path(outdir) / "trimmed"), key="trimmed_dir_input",
            )

        # ── Tool choices ──
        st.markdown('<p class="sec">\U0001f527 Pipeline Tool Choices</p>', unsafe_allow_html=True)
        al1, al2 = st.columns([2, 3])
        with al1:
            default_idx = 1 if RECOMMENDED_ALIGNER == "hisat2" else 0
            aligner = st.radio(
                "Aligner", ["STAR", "HISAT2", "None (pseudo only)"],
                index=default_idx, horizontal=True, key="aligner_choice",
            )
        with al2:
            if aligner == "STAR" and 0 < SYSTEM_RAM_GB < 32:
                st.markdown(
                    f'<div class="box-warn">\u26a0 Low RAM: {SYSTEM_RAM_GB:.0f} GB \u2014 '
                    "STAR needs ~32 GB. <b>HISAT2 recommended.</b></div>",
                    unsafe_allow_html=True,
                )

        tc1, tc2 = st.columns(2)
        with tc1:
            quant_tool = st.radio(
                "Quantification",
                ["featureCounts", "Salmon (aligned)", "Salmon pseudomap", "Both (fC + pseudo)"],
                horizontal=True, key="quant_tool",
            )
        with tc2:
            coverage_mode = st.radio(
                "Coverage",
                ["Housekeeping genes", "Full genome", "Both", "None"],
                horizontal=True, key="coverage_mode",
            )

        tc3, tc4 = st.columns(2)
        with tc3:
            rna_metrics = st.radio(
                "RNA metrics",
                ["RSeQC (strandedness)", "Picard", "Both", "None"],
                horizontal=True, key="rna_metrics",
            )
        with tc4:
            run_screen    = st.checkbox("FastQ Screen", value=True, key="run_screen")
            run_read_dist = st.checkbox("Read distribution", value=False, key="run_read_dist")
            run_qualimap  = st.checkbox("Qualimap", value=False, key="run_qualimap")

        # ── Comparison ──
        if run_deg or run_dromics:
            st.markdown('<p class="sec">\U0001f52c Comparison</p>', unsafe_allow_html=True)
            c1, c2 = st.columns(2)
            chems = st.session_state.chemicals
            with c1:
                treatment = (st.selectbox("Treatment", chems, key="treat_full") if chems
                             else st.text_input("Treatment", value="Bisphenol A", key="treat_full_txt"))
            with c2:
                control = (st.selectbox("Control", chems, index=min(1, len(chems) - 1), key="ctrl_full") if chems
                           else st.text_input("Control", value="Control", key="ctrl_full_txt"))
        else:
            treatment = "Treatment"
            control   = "Control"

        st.markdown('<p class="sec">\U0001f3f7 Run Name</p>', unsafe_allow_html=True)
        run_name_full = st.text_input("Run name", value=_make_run_name(treatment, control), key="run_name_full")

    # ── Derived flags ──
    no_aligner   = aligner == "None (pseudo only)"
    use_fc       = quant_tool in ("featureCounts", "Both (fC + pseudo)") and not no_aligner
    use_salmon   = quant_tool == "Salmon (aligned)" and not no_aligner
    use_pseudo   = quant_tool in ("Salmon pseudomap", "Both (fC + pseudo)") or no_aligner
    use_cov_hk   = coverage_mode in ("Housekeeping genes", "Both") and not no_aligner
    use_cov_full = coverage_mode in ("Full genome", "Both") and not no_aligner
    use_strand   = rna_metrics in ("RSeQC (strandedness)", "Both") and not no_aligner
    use_picard   = rna_metrics in ("Picard", "Both") and not no_aligner
    pseudo_only  = no_aligner or (quant_tool == "Salmon pseudomap" and not use_fc and not use_salmon)

    with R:
        st.markdown('<p class="sec">\U0001f4e1 Status</p>', unsafe_allow_html=True)
        _show_status_panel("full")

        # ── Validation ──
        st.markdown('<p class="sec">\u2705 Validation</p>', unsafe_allow_html=True)
        issues, warnings = _validate_full(
            aligner, no_aligner, pseudo_only, use_fc, use_salmon, use_pseudo,
            use_strand, use_picard, run_qualimap, use_cov_hk, use_cov_full,
            run_read_dist, treatment, control, fastq_source, fastq_dir_val, trimmed_dir_val,
            is_temposeq=is_temposeq, temposeq_manifest_path=temposeq_manifest_path,
        )
        if not issues:
            st.markdown('<div class="box-info">\u2714 All checks passed</div>', unsafe_allow_html=True)
        for i in issues:
            st.markdown(f'<div class="box-warn">\u26a0 {i}</div>', unsafe_allow_html=True)
        for w in warnings:
            st.markdown(f'<div class="box-info">\u2139 {w}</div>', unsafe_allow_html=True)

        # ── Phase 1 Launch Button ──
        st.markdown("---")
        if st.button(
            "\U0001f680 Phase 1: Preprocess (align + QC + counts)",
            disabled=(bool(issues) or running),
            use_container_width=True,
            type="primary",
            key="launch_phase1",
        ):
            Path(outdir).mkdir(parents=True, exist_ok=True)
            acc_file = Path(outdir) / "accessions.txt"
            acc_file.write_text("\n".join(st.session_state.accessions) + "\n")

            params = _build_params(
                "full",
                aligner=aligner, fastq_source=fastq_source,
                fastq_dir_val=fastq_dir_val, trimmed_dir_val=trimmed_dir_val,
                treatment=treatment, control=control, run_name=run_name_full,
                run_screen=run_screen, run_read_dist=run_read_dist,
                run_qualimap=run_qualimap, use_fc=use_fc, use_salmon=use_salmon,
                use_pseudo=use_pseudo, use_strand=use_strand, use_picard=use_picard,
                use_cov_hk=use_cov_hk, use_cov_full=use_cov_full,
                no_aligner=no_aligner, phase="preprocess",
                platform_temposeq=is_temposeq,
                temposeq_manifest_path=temposeq_manifest_path,
            )
            _launch_with_params(params, run_name_full, mode="full", phase="preprocess")
            st.success("\u2714 Phase 1 launched \u2014 preprocessing only (no DESeq2/DRomics yet)")
            time.sleep(1)
            st.rerun()

    # ═══════════════════════════════════════════════════════════════════════
    #  QC REVIEW GATE — full width, between Phase 1 and Phase 2
    # ═══════════════════════════════════════════════════════════════════════

    current_phase  = state.get("phase", "")
    current_status = "running" if running else state.get("status", "idle")

    # Detect count matrix — check all possible locations
    _count_matrix_path = None
    for _candidate in [
        Path(outdir) / "counts" / "gene_counts.csv",
        Path(outdir) / "counts" / "gene_count_matrix.csv",
        Path(outdir) / "counts" / "salmon_pseudo_gene_counts.csv",
        Path(outdir) / "counts" / "salmon_counts_matrix.csv",
    ]:
        if _candidate.exists():
            _count_matrix_path = _candidate
            break

    phase1_done = (
        current_status in ("completed", "review")
        and current_phase == "preprocess"
    )
    phase1_outputs_exist = _count_matrix_path is not None

    # Also detect metadata from Phase 1 output
    _phase1_metadata = Path(outdir) / "metadata" / "metadata.csv"
    if not _phase1_metadata.exists():
        _phase1_metadata = None

    # Show QC review only when Phase 1 actually completed in this session,
    # OR if user explicitly opted in by clicking "Show QC from previous run"
    _show_qc_gate = phase1_done  # strict: only current pipeline run

    if not _show_qc_gate and phase1_outputs_exist and not running:
        # Previous outputs exist but pipeline wasn't just run — offer to load
        if st.button("\U0001f50d Load QC from previous preprocessing run",
                     key="load_prev_qc", use_container_width=False):
            st.session_state._show_prev_qc = True
            st.rerun()
        if st.session_state.get("_show_prev_qc"):
            _show_qc_gate = True

    if _show_qc_gate and not running:
        st.markdown("---")
        st.markdown("### \U0001f50d QC Review \u2014 Review samples before analysis")
        st.info(
            "Preprocessing complete. Review the QC results below and choose which samples "
            "to **exclude**. Then click **Phase 2** to run DESeq2 / DRomics with your choices."
        )

        # Load outlier flags
        outlier_flags = load_outlier_flags(outdir)

        # Collect all sample names — try alignment QC first, then count matrix columns
        all_samples = []
        for s in outlier_flags.get("alignment", {}).get("samples", []):
            name = s.get("sample", "")
            if name and name not in all_samples:
                all_samples.append(name)

        if not all_samples and _count_matrix_path:
            try:
                # Auto-detect separator
                first_line = _count_matrix_path.read_text().split("\n", 1)[0]
                sep = "\t" if "\t" in first_line else ","
                df_cm = pd.read_csv(_count_matrix_path, nrows=0, sep=sep)
                all_samples = [
                    c for c in df_cm.columns
                    if c.lower() not in ("gene_id", "geneid", "transcript_id")
                ]
            except Exception:
                pass

        if all_samples:
            exclude_list = _show_qc_review(outlier_flags, all_samples, key_prefix="full", run_name=run_name_full)
            if exclude_list:
                st.warning(f"Will exclude: **{', '.join(exclude_list)}**")
        else:
            exclude_list = []
            st.caption("No per-sample QC data available \u2014 proceeding with all samples.")

        # MultiQC download
        mqc_report = Path(outdir) / "multiqc" / "multiqc_report.html"
        col_mqc, col_cm, _ = st.columns([1, 1, 2])
        if mqc_report.exists():
            col_mqc.download_button(
                "\u2b07 MultiQC Report", mqc_report.read_bytes(),
                file_name="multiqc_report.html", mime="text/html",
                key="dl_mqc_gate",
            )
        if _count_matrix_path:
            col_cm.download_button(
                "\u2b07 Count Matrix", _count_matrix_path.read_bytes(),
                file_name=_count_matrix_path.name, mime="text/csv",
                key="dl_cm_gate",
            )

        st.markdown("---")

        # ── Transfer to Analysis tab instead of launching Phase 2 here ──
        phase2_issues = []
        if not _count_matrix_path:
            phase2_issues.append(
                "Count matrix not found \u2014 Phase 1 may not have completed successfully."
            )
        meta_for_phase2 = str(_phase1_metadata) if _phase1_metadata else st.session_state.meta_path
        if not meta_for_phase2:
            phase2_issues.append("No metadata file found")

        for i2 in phase2_issues:
            st.markdown(f'<div class="box-warn">\u26a0 {i2}</div>', unsafe_allow_html=True)

        if st.button(
            "\u27a1 Continue to Analysis (DESeq2 / DRomics)",
            disabled=(bool(phase2_issues) or running),
            use_container_width=True,
            type="primary",
            key="transfer_to_analysis",
        ):
            # Auto-populate the analysis tab with Phase 1 outputs
            st.session_state.count_matrix_path = str(_count_matrix_path)
            st.session_state.meta_path = meta_for_phase2
            st.session_state.exclude_samples = exclude_list
            # Clear stale QC state from any previous analysis run
            st.session_state._analysis_qc_run_id = None
            st.session_state._qc_valid = False
            st.session_state._qc_exclusions_set_for = ""
            # Load metadata for chemical selection
            try:
                _meta_df = pd.read_csv(meta_for_phase2)
                st.session_state.meta_df = _meta_df
                if "Treatment" in _meta_df.columns:
                    st.session_state.chemicals = sorted(_meta_df["Treatment"].unique().tolist())
            except Exception:
                pass
            st.session_state._phase1_transferred = True
            st.success(
                "\u2714 Count matrix and metadata transferred to **DESeq2 / DRomics** tab. "
                "Switch to that tab to run Layer 2 QC and analysis."
            )


# ═════════════════════════════════════════════════════════════════════════════
#  TAB 2 — DESeq2 / DRomics ANALYSIS
#  Three-step flow: (1) Load data → (2) QC check + review → (3) Run analysis
#  Supports both uploaded count matrices AND Phase 1 transfers.
#  User can run DESeq2 only, DRomics only, or both independently.
# ═════════════════════════════════════════════════════════════════════════════

with tab_analysis:
    # Show transfer notice if coming from Phase 1
    if st.session_state.get("_phase1_transferred"):
        st.success(
            "\u2714 **Data loaded from Full Pipeline.** Count matrix and metadata from Phase 1 "
            "are ready below. Run QC check, review outliers, then launch analysis."
        )

    st.info(
        "\U0001f52c **DESeq2 / DRomics Analysis** \u2014 "
        "Provide a count matrix (uploaded or from Phase 1). "
        "Step 1: QC check \u2192 Step 2: Review & exclude \u2192 Step 3: Run DESeq2 and/or DRomics."
    )

    dL, dR = st.columns([3, 2], gap="large")

    with dL:
        st.markdown('<p class="sec">\U0001f4ca Count Matrix</p>', unsafe_allow_html=True)

        # Show current source if already loaded
        if st.session_state.count_matrix_path:
            _src = "Phase 1 pipeline" if st.session_state.get("_phase1_transferred") else "uploaded file"
            st.markdown(
                f'<div class="box-info">\u2714 Count matrix loaded from {_src}: '
                f'<code>{Path(st.session_state.count_matrix_path).name}</code></div>',
                unsafe_allow_html=True)
            if st.button("\u21ba Load different count matrix", key="reset_counts"):
                st.session_state.count_matrix_path = None
                st.session_state.count_matrix_df = None
                st.session_state._phase1_transferred = False
                # Clear stale QC state so old results don't appear for new data
                st.session_state._analysis_qc_run_id = None
                st.session_state._qc_valid = False
                st.session_state.exclude_samples = []
                st.session_state._qc_exclusions_set_for = ""
                st.rerun()
        else:
            st.caption("Accepts: featureCounts output, gene_count_matrix.csv, .tsv, .xlsx")
            count_file = st.file_uploader(
                "Upload count matrix",
                type=["out", "csv", "tsv", "txt", "xlsx", "xls"],
                label_visibility="collapsed", key="count_upload",
            )

            if count_file:
                try:
                    fname = count_file.name
                    if fname.endswith((".xlsx", ".xls")):
                        df_counts = pd.read_excel(count_file)
                    elif fname.endswith((".out", ".tsv")):
                        raw_text = count_file.read().decode(errors="replace")
                        count_file.seek(0)
                        lines = raw_text.splitlines()
                        skip = sum(1 for l in lines if l.startswith("#"))
                        from io import StringIO
                        df_counts = pd.read_csv(StringIO("\n".join(lines[skip:])), sep="\t")
                    elif fname.endswith(".csv"):
                        df_counts = pd.read_csv(count_file)
                    else:
                        raw_text = count_file.read().decode(errors="replace")
                        count_file.seek(0)
                        from io import StringIO
                        sep = "\t" if "\t" in raw_text[:500] else ","
                        df_counts = pd.read_csv(StringIO(raw_text), sep=sep, comment="#")

                    st.session_state.count_matrix_df = df_counts
                    Path(outdir).mkdir(parents=True, exist_ok=True)
                    saved_path = Path(outdir) / "uploaded_counts.csv"
                    df_counts.to_csv(saved_path, index=False)
                    st.session_state.count_matrix_path = str(saved_path)
                    st.success(f"\u2714 {df_counts.shape[0]} genes \u00d7 {df_counts.shape[1]} columns")
                    with st.expander("Preview"):
                        st.dataframe(df_counts.head(20), use_container_width=True, hide_index=True)
                except Exception as e:
                    st.error(f"Error: {e}")

        st.markdown('<p class="sec">\U0001f4cb Metadata</p>', unsafe_allow_html=True)
        if st.session_state.meta_path and st.session_state.meta_df is not None:
            _meta_mapped = bool(st.session_state.chemicals) and "Treatment" in st.session_state.meta_df.columns
            st.markdown(
                f'<div class="box-info">\u2714 Metadata loaded: '
                f'{len(st.session_state.meta_df)} samples'
                f'{" (columns mapped)" if _meta_mapped else " — map columns below"}</div>',
                unsafe_allow_html=True)

            # If metadata loaded but not yet mapped (e.g. from Full Pipeline), show mapper
            if not _meta_mapped:
                st.markdown("**Map your metadata columns before selecting treatment/control:**")
                col_map = _show_column_mapper(st.session_state.meta_df, key_prefix="colmap_remap")
                req_mapped = ["Sample_Name", "Treatment", "Dose"]
                missing = [r for r in req_mapped if r not in col_map]
                if missing:
                    st.warning(f"Please map required columns: {', '.join(missing)}")
                elif st.button("\u2714 Apply mapping", key="apply_colmap_remap"):
                    df_mapped = _apply_column_mapping(st.session_state.meta_df, col_map)
                    Path(outdir).mkdir(parents=True, exist_ok=True)
                    mp = Path(outdir) / "metadata_input.csv"
                    df_mapped.to_csv(mp, index=False)
                    st.session_state.meta_df = df_mapped
                    st.session_state.meta_path = str(mp)
                    st.session_state.col_map = col_map
                    st.session_state.chemicals = sorted(df_mapped["Treatment"].unique().tolist())
                    st.success(f"\u2714 Mapped! {df_mapped['Treatment'].nunique()} treatments found")
                    st.rerun()

            if st.button("\u21ba Load different metadata", key="reset_meta_analysis"):
                st.session_state.meta_path = None
                st.session_state.meta_df = None
                st.session_state.chemicals = []
                st.session_state.col_map = {}
                # Clear stale QC — new metadata means old QC is invalid
                st.session_state._analysis_qc_run_id = None
                st.session_state._qc_valid = False
                st.session_state.exclude_samples = []
                st.session_state._qc_exclusions_set_for = ""
                st.rerun()
        else:
            st.caption("Upload your metadata file — you'll map columns in the next step")
            meta_upload = st.file_uploader(
                "Metadata CSV or XLSX",
                type=["csv", "xlsx", "xls"], label_visibility="collapsed", key="meta_analysis",
            )

            if meta_upload:
                try:
                    df_m_raw = (pd.read_csv(meta_upload) if meta_upload.name.endswith(".csv")
                                else pd.read_excel(meta_upload))
                    st.session_state._meta_raw = df_m_raw

                    with st.expander("Raw preview", expanded=False):
                        st.dataframe(df_m_raw.head(10), use_container_width=True, hide_index=True)

                    # Column mapping UI
                    col_map = _show_column_mapper(df_m_raw, key_prefix="colmap_analysis")

                    # Validate required columns are mapped
                    req_mapped = ["Sample_Name", "Treatment", "Dose"]
                    missing = [r for r in req_mapped if r not in col_map]

                    if missing:
                        st.warning(f"Please map required columns: {', '.join(missing)}")
                    else:
                        if st.button("\u2714 Apply mapping & load", key="apply_colmap_analysis"):
                            df_m = _apply_column_mapping(df_m_raw, col_map)
                            Path(outdir).mkdir(parents=True, exist_ok=True)
                            mp = Path(outdir) / "metadata_input.csv"
                            df_m.to_csv(mp, index=False)
                            st.session_state.meta_df = df_m
                            st.session_state.meta_path = str(mp)
                            st.session_state.col_map = col_map
                            st.session_state.chemicals = sorted(df_m["Treatment"].unique().tolist())
                            st.success(f"\u2714 {len(df_m)} samples | {df_m['Treatment'].nunique()} treatments")
                            st.rerun()
                except Exception as e:
                    st.error(f"Error: {e}")

        # ── Comparison settings ──
        st.markdown('<p class="sec">\U0001f52c Comparison</p>', unsafe_allow_html=True)
        dc1, dc2 = st.columns(2)
        chems_a = st.session_state.chemicals

        # Auto-populate chemicals from metadata if empty but metadata is loaded
        if not chems_a and st.session_state.meta_df is not None and "Treatment" in st.session_state.meta_df.columns:
            chems_a = sorted(st.session_state.meta_df["Treatment"].unique().tolist())
            st.session_state.chemicals = chems_a

        # Detect if chemicals changed and reset stale selectbox values
        if chems_a:
            # If the currently selected treatment isn't in the new chemicals list, reset
            if st.session_state.get("treat_analysis") and \
               st.session_state.get("treat_analysis") not in chems_a:
                del st.session_state["treat_analysis"]
            if st.session_state.get("ctrl_analysis") and \
               st.session_state.get("ctrl_analysis") not in chems_a:
                del st.session_state["ctrl_analysis"]

        with dc1:
            if chems_a:
                treatment_a = st.selectbox("Treatment", chems_a, key="treat_analysis")
            else:
                treatment_a = st.text_input("Treatment", value="Bisphenol A", key="treat_analysis_txt")
        with dc2:
            if chems_a:
                # Default to second item (often "Control") if available
                ctrl_idx = min(1, len(chems_a) - 1)
                control_a = st.selectbox("Control", chems_a, index=ctrl_idx, key="ctrl_analysis")
            else:
                control_a = st.text_input("Control", value="Control", key="ctrl_analysis_txt")

        # ── Analysis selection — user chooses which to run ──
        st.markdown('<p class="sec">\U0001f527 Analysis Selection</p>', unsafe_allow_html=True)
        ac1, ac2 = st.columns(2)
        with ac1:
            run_deseq2_here = st.checkbox("Run DESeq2", value=run_deg, key="run_deseq2_analysis")
        with ac2:
            run_dromics_here = st.checkbox("Run DRomics / BMD", value=run_dromics, key="run_dromics_analysis")

        st.markdown('<p class="sec">\U0001f3f7 Run Name</p>', unsafe_allow_html=True)
        run_name_analysis = _make_run_name(treatment_a, control_a)
        st.code(run_name_analysis, language=None)

    with dR:
        st.markdown('<p class="sec">\U0001f4e1 Status</p>', unsafe_allow_html=True)
        _show_status_panel("analysis")

        st.markdown('<p class="sec">\u2705 Validation</p>', unsafe_allow_html=True)
        a_issues = []
        if not MAIN_NF.exists():
            a_issues.append("main.nf not found")
        if not st.session_state.count_matrix_path:
            a_issues.append("No count matrix loaded")
        if not st.session_state.meta_path:
            a_issues.append("No metadata loaded")
        if not run_deseq2_here and not run_dromics_here:
            a_issues.append("Select at least one analysis (DESeq2 or DRomics)")
        if treatment_a and control_a and treatment_a == control_a:
            a_issues.append("Treatment = Control")

        if not a_issues:
            st.markdown('<div class="box-info">\u2714 Ready</div>', unsafe_allow_html=True)
        for i in a_issues:
            st.markdown(f'<div class="box-warn">\u26a0 {i}</div>', unsafe_allow_html=True)

        # ── Step 1: QC Check (Layer 2 — PCA) ──
        st.markdown("---")
        if st.button(
            "\U0001f50d Step 1: Run QC Check (PCA outlier detection)",
            disabled=(bool(a_issues) or running),
            use_container_width=True, type="primary",
            key="analysis_qc",
        ):
            # Subset count matrix and metadata to selected treatment/control
            sub_counts, sub_meta, n_sub = _subset_for_comparison(
                st.session_state.count_matrix_path,
                st.session_state.meta_df,
                treatment_a, control_a,
                outdir, run_name_analysis,
            )
            if not sub_counts or n_sub == 0:
                st.error("Could not subset data for the selected comparison. "
                         "Check that metadata Sample_Name values match count matrix column names.")
            else:
                # Run QC SYNCHRONOUSLY — this is a single Rscript call (seconds, not hours).
                # Running it synchronously means results appear immediately after completion
                # without any async polling / tab-switching nonsense.
                _qc_outdir = Path(outdir) / "runs" / run_name_analysis / "qc_results"
                _qc_outdir.mkdir(parents=True, exist_ok=True)

                _qc_cmd = [
                    "Rscript", str(PIPELINE_DIR / "scripts" / "run_qc_check.R"),
                    "--counts",    sub_counts,
                    "--metadata",  sub_meta,
                    "--treatment", treatment_a,
                    "--control",   control_a,
                    "--outdir",    str(_qc_outdir),
                    "--read_thresh", str(st.session_state.read_thresh),
                    "--mode",      "deseq2",
                ]

                with st.spinner(f"Running PCA outlier detection on {n_sub} samples..."):
                    try:
                        _qc_result = subprocess.run(
                            _qc_cmd,
                            capture_output=True, text=True,
                            timeout=300,  # 5 min max — QC should be <30s
                            env=_process_env(),
                            cwd=str(_qc_outdir),
                        )
                        if _qc_result.returncode != 0:
                            st.error(f"QC check failed (exit code {_qc_result.returncode})")
                            with st.expander("Error details"):
                                st.code(_qc_result.stderr[-2000:] if _qc_result.stderr else "No error output")
                                if _qc_result.stdout:
                                    st.code(_qc_result.stdout[-2000:])
                        else:
                            # Mark QC as done in session state
                            st.session_state._qc_ran_for_treatment = treatment_a
                            st.session_state._qc_ran_for_control = control_a
                            st.session_state._qc_valid = True
                            st.session_state._analysis_qc_run_id = run_name_analysis
                            # Clear stale exclusions from a previous QC run
                            st.session_state.exclude_samples = []
                            st.session_state._qc_exclusions_set_for = ""
                            st.success(f"\u2714 QC check complete for {n_sub} samples")
                            st.rerun()
                    except subprocess.TimeoutExpired:
                        st.error("QC check timed out after 5 minutes. Try running via the Full Pipeline tab.")
                    except FileNotFoundError:
                        st.error("Rscript not found. Check that R is installed and in PATH.")
                    except Exception as e:
                        st.error(f"QC check error: {e}")

    # ── QC Review Gate (full width, below the two columns) ──
    # Since QC runs synchronously, the output file exists immediately after Step 1 completes.
    # We just check: file exists + session state says QC was run for this comparison.

    _analysis_qc_flag = None
    _qc_dir = Path(outdir) / "runs" / run_name_analysis / "qc_results"
    for _fname in ["pca_outlier_flag.json", "qc_summary.json"]:
        _qc_cand = _qc_dir / _fname
        if _qc_cand.exists():
            _analysis_qc_flag = _qc_cand
            break
    if not _analysis_qc_flag and _qc_dir.exists():
        for _f in _qc_dir.glob("*.json"):
            _analysis_qc_flag = _f
            break

    analysis_qc_done = (
        _analysis_qc_flag is not None
        and st.session_state.get("_analysis_qc_run_id") == run_name_analysis
    )

    # Detect stale QC from a previous run/session
    _stale_qc_exists = (
        not analysis_qc_done
        and _analysis_qc_flag is not None
    )

    if _stale_qc_exists and st.session_state.count_matrix_path:
        st.warning(
            "\u26a0 QC results found from a **previous run**. Click **Step 1: Run QC Check** "
            "to generate fresh QC for the current comparison, or click below to load the old results."
        )
        if st.button(
            "\U0001f50d Load previous QC results for review",
            key="load_stale_qc_analysis",
        ):
            st.session_state._analysis_qc_run_id = run_name_analysis
            st.session_state.exclude_samples = []
            st.session_state._qc_exclusions_set_for = ""
            st.rerun()
    elif not analysis_qc_done and st.session_state.count_matrix_path:
        st.caption("Click **Step 1: Run QC Check** to detect PCA outliers before analysis.")

    if analysis_qc_done:
        st.markdown("---")
        st.markdown("### \U0001f50d Step 2: QC Review \u2014 Review PCA outliers before analysis")
        st.info(
            "PCA outlier detection complete. Review the results below, select samples to exclude, "
            "then choose which analysis to run."
        )

        # Load PCA outlier flags
        analysis_flags = {}
        if _analysis_qc_flag and _analysis_qc_flag.exists():
            try:
                raw_qc = json.loads(_analysis_qc_flag.read_text())
                # Handle qc_summary.json format (has all_samples + flagged_samples)
                if "all_samples" in raw_qc and "samples" not in raw_qc:
                    flagged = set(raw_qc.get("flagged_samples", []))
                    raw_qc["samples"] = [
                        {"sample": s, "flagged": s in flagged}
                        for s in raw_qc["all_samples"]
                    ]
                    raw_qc["flagged"] = list(flagged)
                analysis_flags["pca"] = _flatten_qc_json(raw_qc)
            except Exception:
                pass
        # Also try loading the other QC file if first one didn't have sample data
        if not analysis_flags.get("pca", {}).get("samples"):
            for _alt in ["pca_outlier_flag.json", "qc_summary.json"]:
                _alt_path = _qc_dir / _alt
                if _alt_path.exists() and _alt_path != _analysis_qc_flag:
                    try:
                        raw2 = json.loads(_alt_path.read_text())
                        if "all_samples" in raw2 and "samples" not in raw2:
                            flagged2 = set(raw2.get("flagged_samples", []))
                            raw2["samples"] = [
                                {"sample": s, "flagged": s in flagged2}
                                for s in raw2["all_samples"]
                            ]
                            raw2["flagged"] = list(flagged2)
                        analysis_flags["pca"] = _flatten_qc_json(raw2)
                        break
                    except Exception:
                        pass

        # Also load alignment-level flags if available (from Phase 1)
        # But ONLY if the samples match the current count matrix — prevents cross-study contamination
        align_path = Path(outdir) / "qc" / "outlier" / "outlier_flag.json"
        if align_path.exists() and st.session_state.count_matrix_path:
            try:
                _align_data = _flatten_qc_json(json.loads(align_path.read_text()))
                _align_samples = {s.get("sample", "") for s in _align_data.get("samples", [])}
                # Check overlap with count matrix columns
                _cm_path = Path(st.session_state.count_matrix_path)
                _first_line = _cm_path.read_text().split("\n", 1)[0]
                _sep = "\t" if "\t" in _first_line else ","
                _cm_cols = set(pd.read_csv(_cm_path, nrows=0, sep=_sep).columns)
                # Only use alignment flags if >50% of alignment samples are in the count matrix
                _overlap = _align_samples & _cm_cols
                if len(_overlap) > len(_align_samples) * 0.5:
                    analysis_flags["alignment"] = _align_data
                else:
                    pass  # alignment flags are from a different study — skip
            except Exception:
                pass

        # Get sample names — FILTERED to selected treatment/control only
        analysis_all_samples = []

        # First, get the relevant samples from metadata for the selected comparison
        _relevant_samples = set()
        if st.session_state.meta_df is not None and "Treatment" in st.session_state.meta_df.columns:
            _meta = st.session_state.meta_df
            _relevant = _meta[_meta["Treatment"].isin([treatment_a, control_a])]
            if "Sample_Name" in _relevant.columns:
                _relevant_samples = set(_relevant["Sample_Name"].astype(str).tolist())

        # Build sample list from QC flags, filtered to relevant samples
        for layer_key in ("alignment", "pca"):
            for s in analysis_flags.get(layer_key, {}).get("samples", []):
                name = s.get("sample", "")
                if name and name not in analysis_all_samples:
                    if not _relevant_samples or name in _relevant_samples:
                        analysis_all_samples.append(name)

        # Fallback: from count matrix columns, filtered to relevant samples
        if not analysis_all_samples and st.session_state.count_matrix_path:
            try:
                _cm_path = Path(st.session_state.count_matrix_path)
                first_line = _cm_path.read_text().split("\n", 1)[0]
                sep = "\t" if "\t" in first_line else ","
                df_cm = pd.read_csv(_cm_path, nrows=0, sep=sep)
                all_cols = [
                    c for c in df_cm.columns
                    if c.lower() not in ("gene_id", "geneid", "transcript_id", "")
                ]
                if _relevant_samples:
                    analysis_all_samples = [c for c in all_cols if c in _relevant_samples]
                else:
                    analysis_all_samples = all_cols
            except Exception:
                pass

        if _relevant_samples and analysis_all_samples:
            st.caption(f"Showing {len(analysis_all_samples)} samples for "
                       f"**{treatment_a}** vs **{control_a}**")

        # Also filter the QC flags to only show relevant samples
        for layer_key in ("alignment", "pca"):
            if layer_key in analysis_flags and "samples" in analysis_flags[layer_key]:
                if _relevant_samples:
                    analysis_flags[layer_key]["samples"] = [
                        s for s in analysis_flags[layer_key]["samples"]
                        if s.get("sample", "") in _relevant_samples
                    ]

        if analysis_all_samples:
            analysis_exclude = _show_qc_review(analysis_flags, analysis_all_samples, key_prefix="analysis", run_name=run_name_analysis)
            if analysis_exclude:
                st.warning(f"Will exclude: **{', '.join(analysis_exclude)}**")
        else:
            analysis_exclude = []
            st.caption("No per-sample QC data \u2014 proceeding with all samples.")

        st.markdown("---")

        # ── Step 3: Run Analysis — DESeq2 and DRomics independently ──
        st.markdown("### \U0001f52c Step 3: Run Analysis")

        a2_issues = []
        if not st.session_state.count_matrix_path:
            a2_issues.append("No count matrix")
        if not st.session_state.meta_path:
            a2_issues.append("No metadata")

        # Check for subsetted files from Step 1
        _sub_dir = Path(outdir) / "runs" / run_name_analysis
        _sub_counts_file = _sub_dir / "subset_counts.csv"
        _sub_meta_file = _sub_dir / "subset_metadata.csv"
        if _sub_counts_file.exists() and _sub_meta_file.exists():
            _analysis_counts = str(_sub_counts_file)
            _analysis_meta = str(_sub_meta_file)
        else:
            # Will be created when user clicks a run button
            _analysis_counts = st.session_state.count_matrix_path
            _analysis_meta = st.session_state.meta_path

        # Independent buttons for DESeq2 and DRomics
        btn_c1, btn_c2, btn_c3 = st.columns(3)

        with btn_c1:
            if st.button(
                "\U0001f4ca Run DESeq2 only",
                disabled=(bool(a2_issues) or running or not run_deseq2_here),
                use_container_width=True, type="primary",
                key="run_deseq2_only",
            ):
                params = _build_params(
                    "direct",
                    treatment=treatment_a, control=control_a, run_name=run_name_analysis,
                    phase="analysis",
                    count_matrix=_analysis_counts,
                    exclude_samples=analysis_exclude,
                )
                params["metadata"] = _analysis_meta
                params["run_qc_only"] = False
                params["run_deg"] = True
                params["run_dromics"] = False
                _launch_with_params(params, run_name_analysis, mode="direct", phase="analysis")
                st.session_state._qc_valid = False  # reset after launching analysis
                st.success(
                    f"\u2714 DESeq2 launched "
                    f"({'excluding ' + ', '.join(analysis_exclude) if analysis_exclude else 'all samples'})"
                )
                time.sleep(1)
                st.rerun()

        with btn_c2:
            if st.button(
                "\U0001f4c8 Run DRomics only",
                disabled=(bool(a2_issues) or running or not run_dromics_here),
                use_container_width=True, type="primary",
                key="run_dromics_only",
            ):
                params = _build_params(
                    "direct",
                    treatment=treatment_a, control=control_a, run_name=run_name_analysis + "_dromics",
                    phase="analysis",
                    count_matrix=_analysis_counts,
                    exclude_samples=analysis_exclude,
                )
                params["metadata"] = _analysis_meta
                params["run_qc_only"] = False
                params["run_deg"] = False
                params["run_dromics"] = True
                _launch_with_params(params, run_name_analysis + "_dromics", mode="direct", phase="analysis")
                st.session_state._qc_valid = False
                st.success(
                    f"\u2714 DRomics launched "
                    f"({'excluding ' + ', '.join(analysis_exclude) if analysis_exclude else 'all samples'})"
                )
                time.sleep(1)
                st.rerun()

        with btn_c3:
            if st.button(
                "\U0001f52c Run Both",
                disabled=(bool(a2_issues) or running or (not run_deseq2_here and not run_dromics_here)),
                use_container_width=True,
                key="run_both_analysis",
            ):
                params = _build_params(
                    "direct",
                    treatment=treatment_a, control=control_a, run_name=run_name_analysis,
                    phase="analysis",
                    count_matrix=_analysis_counts,
                    exclude_samples=analysis_exclude,
                )
                params["metadata"] = _analysis_meta
                params["run_qc_only"] = False
                params["run_deg"] = run_deseq2_here
                params["run_dromics"] = run_dromics_here
                _launch_with_params(params, run_name_analysis, mode="direct", phase="analysis")
                st.session_state._qc_valid = False
                st.success(
                    f"\u2714 Analysis launched "
                    f"({'excluding ' + ', '.join(analysis_exclude) if analysis_exclude else 'all samples'})"
                )
                time.sleep(1)
                st.rerun()

        if not run_deseq2_here:
            st.caption("\u2139 DESeq2 disabled \u2014 enable it in the Analysis Selection above.")
        if not run_dromics_here:
            st.caption("\u2139 DRomics disabled \u2014 enable it in the Analysis Selection above.")

    # ── Step 4: Inline Results — show analysis results right here ──────────
    # Detect completed analysis run for the current comparison and show results
    # inline without requiring the user to switch to the Results tab.
    _analysis_run_dir = Path(outdir) / "runs" / run_name_analysis
    _analysis_has_deg = (_analysis_run_dir / "deg_results" / "analysis_summary.json").exists()
    _analysis_has_dromics = (_analysis_run_dir / "dromics_results" / "dromics_summary.json").exists()
    _analysis_has_results = _analysis_has_deg or _analysis_has_dromics

    # Also check the _dromics variant run name (DRomics-only runs append "_dromics")
    _dromics_run_dir = Path(outdir) / "runs" / (run_name_analysis + "_dromics")
    if not _analysis_has_dromics and (_dromics_run_dir / "dromics_results" / "dromics_summary.json").exists():
        _analysis_has_dromics = True
        _dromics_results_dir = _dromics_run_dir
    else:
        _dromics_results_dir = _analysis_run_dir

    if _analysis_has_results and not running:
        st.markdown("---")
        st.markdown("### \U0001f4ca Step 4: Results")

        # Show DESeq2 results inline
        if _analysis_has_deg:
            _inline_deg_dir = _analysis_run_dir / "deg_results"
            _inline_summary = _inline_deg_dir / "analysis_summary.json"
            if _inline_summary.exists():
                st.markdown("#### DESeq2")
                _s = json.loads(_inline_summary.read_text())
                _n_samples = _jv(_s, "samples_analyzed")
                _n_genes = _jv(_s, "total_genes")
                _n_degs = _jv(_s, "sig_relaxed")
                _n_up = _jv(_s, "upregulated")
                _n_down = _jv(_s, "downregulated")

                _mc = st.columns(5)
                for _col, _val, _label in zip(_mc,
                    [_n_samples, f"{_n_genes:,}", f"{_n_degs:,}", f"{_n_up:,}", f"{_n_down:,}"],
                    ["Samples", "Genes Tested", "DEGs (FDR<0.05)", "\u2b06 Up", "\u2b07 Down"]):
                    _col.markdown(
                        f'<div class="metric-card"><div class="val">{_val}</div>'
                        f'<div class="lab">{_label}</div></div>',
                        unsafe_allow_html=True)

                _pca_out = _jv(_s, "pca_outliers_flagged", 0)
                _batch_cov = "Yes" if _jv(_s, "batch_as_covariate") else "No"
                _batch_pca = "Yes" if _jv(_s, "batch_corrected_pca") else "No"
                st.caption(
                    f"Batch in design: {_batch_cov} \u00b7 "
                    f"Batch-corrected PCA: {_batch_pca} \u00b7 "
                    f"PCA outliers: {_pca_out}"
                )

                _pc, _vc = st.columns(2)
                _pca_img = _inline_deg_dir / "Final_PCA_Plot.png"
                _vol_img = _inline_deg_dir / "Volcano_Plot.png"
                if _pca_img.exists():
                    _pca_cap = "PCA \u2014 batch-corrected VST" if _jv(_s, "batch_corrected_pca") else "PCA \u2014 VST counts"
                    _pc.image(str(_pca_img), caption=_pca_cap, use_container_width=True)
                if _vol_img.exists():
                    _vc.image(str(_vol_img), caption="Volcano Plot", use_container_width=True)

                _sig_csv = _inline_deg_dir / "Custom_Filtered_DEGs.csv"
                if _sig_csv.exists():
                    with st.expander(f"\U0001f4ca Top DEGs ({min(30, _n_degs)} shown)", expanded=True):
                        st.dataframe(pd.read_csv(_sig_csv).head(30), use_container_width=True, hide_index=True)
                    _dl_cols = st.columns(4)
                    for _dlc, _fname in zip(_dl_cols, ["All_Results.csv", "Custom_Filtered_DEGs.csv",
                                                        "Upregulated_Genes.csv", "Downregulated_Genes.csv"]):
                        _fp = _inline_deg_dir / _fname
                        if _fp.exists():
                            _dlc.download_button(
                                f"\u2b07 {_fname.replace('_', ' ').replace('.csv', '')}",
                                _fp.read_bytes(), file_name=_fname, mime="text/csv",
                                key=f"inline_dl_{run_name_analysis}_{_fname}",
                            )
                st.divider()

        # Show DRomics results inline
        if _analysis_has_dromics:
            _inline_dr_dir = _dromics_results_dir / "dromics_results"
            _inline_dp = _inline_dr_dir / "dromics_summary.json"
            if _inline_dp.exists():
                st.markdown("#### DRomics / BMD")
                _ds = json.loads(_inline_dp.read_text())
                _n_genes_dr   = _jv(_ds, "total_genes_final")
                _n_responsive = _jv(_ds, "genes_selected")
                _n_models     = _jv(_ds, "models_fitted")

                _metrics = [
                    (f"{_n_genes_dr:,}",   "Genes Analyzed"),
                    (f"{_n_responsive:,}", "Responsive Genes"),
                    (f"{_n_models:,}",     "Models Fitted"),
                ]
                if _ds.get("bmd_performed") and _ds.get("bmd_summary"):
                    _bms = _ds["bmd_summary"]
                    _median_bmd = round(float(_jv(_bms, "median_bmd", 0)), 4)
                    _metrics.append((f"{_median_bmd} \u03bcM", "Median Gene BMD"))
                if _ds.get("tpod_performed") and _ds.get("tpod_summary"):
                    _tpod_val = round(float(_jv(_ds["tpod_summary"], "tpod_value", 0)), 4)
                    _metrics.append((f"{_tpod_val} \u03bcM", "tPOD (NTP 2018)"))

                _mc2 = st.columns(len(_metrics))
                for _col2, (_val2, _label2) in zip(_mc2, _metrics):
                    _col2.markdown(
                        f'<div class="metric-card"><div class="val">{_val2}</div>'
                        f'<div class="lab">{_label2}</div></div>',
                        unsafe_allow_html=True)

                st.caption(
                    f"Samples: {_jv(_ds, 'total_samples_final')} \u00b7 "
                    f"Dose levels: {_jv(_ds, 'dose_levels')} \u00b7 "
                    f"PCA outliers: {_jv(_ds, 'pca_outliers_flagged', 0)}"
                )

                # tPOD highlight
                if _ds.get("tpod_performed") and _ds.get("tpod_summary"):
                    _tpod = _ds["tpod_summary"]
                    st.markdown(
                        f'<div class="box-info" style="border-color:rgba(255,87,34,.4);background:rgba(255,87,34,.04)">'
                        f'\U0001f3af <b>tPOD: '
                        f'{round(float(_jv(_tpod, "tpod_value", 0)), 4)} \u03bcM</b> &nbsp;\u00b7&nbsp; '
                        f'Pathway: <b>{_jv(_tpod, "tpod_pathway", "?")}</b> '
                        f'({_jv(_tpod, "tpod_pathway_id", "")}) &nbsp;\u00b7&nbsp; '
                        f'{_jv(_tpod, "tpod_n_genes", "?")} genes &nbsp;\u00b7&nbsp; '
                        f'{_jv(_tpod, "tpod_coverage", "?")}% coverage</div>',
                        unsafe_allow_html=True)

                # Plots
                _c1, _c2 = st.columns(2)
                for _ci, _fname, _cap in [(_c1, "dose_response_curves.png", "Dose-Response Curves"),
                                            (_c2, "bmd_distribution.png",    "BMD Distribution")]:
                    _fp2 = _inline_dr_dir / _fname
                    if _fp2.exists():
                        _ci.image(str(_fp2), caption=_cap, use_container_width=True)

                _pw_plot = _inline_dr_dir / "pathway_sensitivity.png"
                _top25_img = _inline_dr_dir / "top25_sensitive_genes.png"
                if _pw_plot.exists() or _top25_img.exists():
                    _pc3, _pc4 = st.columns(2)
                    if _pw_plot.exists():
                        _pc3.image(str(_pw_plot), caption="Pathway Sensitivity",
                                   use_container_width=True)
                    if _top25_img.exists():
                        _pc4.image(str(_top25_img), caption="Top 25 Sensitive Genes",
                                   use_container_width=True)

                # Downloads
                _dr_dl_files = [
                    "bmd_results.csv", "top50_sensitive_genes.csv",
                    "dose_response_models.csv", "pathway_bmd_summary.csv",
                ]
                _dr_available = [f for f in _dr_dl_files if (_inline_dr_dir / f).exists()]
                if _dr_available:
                    _dr_dl_cols = st.columns(min(len(_dr_available), 4))
                    for _di, _dfname in enumerate(_dr_available):
                        _dfp = _inline_dr_dir / _dfname
                        _dr_dl_cols[_di % 4].download_button(
                            f"\u2b07 {_dfname.replace('_', ' ').replace('.csv', '')}",
                            _dfp.read_bytes(), file_name=_dfname, mime="text/csv",
                            key=f"inline_dr_dl_{run_name_analysis}_{_dfname}",
                        )

        st.caption(
            "\U0001f4cb Full results with additional tables and plots are available in the **Results** tab."
        )


# ═════════════════════════════════════════════════════════════════════════════
#  TAB 3 — PARAMETERS
# ═════════════════════════════════════════════════════════════════════════════

with tab_params:
    st.markdown("#### DESeq2")
    d1, d2, d3 = st.columns(3)
    with d1:
        st.number_input("Strict FDR",  min_value=0.001, max_value=0.1, step=0.001, format="%.3f", key="fdr_strict")
        st.number_input("Relaxed FDR", min_value=0.01,  max_value=0.2, step=0.01, key="fdr_relaxed")
    with d2:
        st.number_input("log2FC threshold",   min_value=0.0, max_value=5.0, step=0.1, key="log2fc")
        st.number_input("Min reads / sample", min_value=0,   step=100000, key="read_thresh")
    with d3:
        st.number_input("CPM threshold", value=1.0, min_value=0.0, step=0.1, key="cpm_thresh")
        st.slider("Min proportion expressed", 0.0, 1.0, 0.75, 0.05, key="min_prop")

    st.divider()
    st.markdown("#### DRomics / BMD")
    b1, b2, b3, b4 = st.columns(4)
    with b1:
        st.number_input("DRomics FDR", min_value=0.001, max_value=0.2, step=0.01, key="dr_fdr")
        st.selectbox(
            "Model criterion", ["AICc", "AIC", "BIC"],
            help="AICc recommended (corrected AIC, prevents overfitting with small samples)",
            key="dr_criterion",
        )
        st.selectbox(
            "Gene selection method", ["quadratic", "linear", "ANOVA"],
            key="dr_select_method",
        )
    with b2:
        st.checkbox("Perform BMD", key="perform_bmd")
        st.checkbox("Bootstrap CIs (slow)", key="bmd_bootstrap")
        st.selectbox(
            "Transformation", ["auto", "vst", "rlog"],
            help="auto = vst for n\u226530, rlog for n<30",
            key="dr_transfo",
        )
    with b3:
        st.number_input("BMD z-value", min_value=0.1, max_value=5.0, step=0.1, help="Signal-to-noise ratio threshold", key="dr_bmd_z")
        st.number_input("BMD x-value (%)", min_value=1.0, max_value=50.0, step=1.0, help="Percent change threshold", key="dr_bmd_x")
        st.number_input("Bootstrap iterations", min_value=100, max_value=20000, step=100, key="dr_niter")
    with b4:
        df_meta = st.session_state.meta_df
        if df_meta is not None and "Dose" in df_meta.columns:
            doses = sorted(df_meta["Dose"].unique())
            st.markdown("**Dose levels:**")
            st.write(", ".join(f"{d} \u03bcM" for d in doses))
            st.dataframe(df_meta.groupby("Dose").size().rename("n"), use_container_width=True)

    st.divider()
    st.markdown("#### NTP 2018 / EFSA Quality Filters")
    st.caption("Applied after BMD calculation to filter unreliable estimates before pathway tPOD.")
    ntp1, ntp2 = st.columns(2)
    with ntp1:
        st.number_input(
            "Max BMDU/BMDL ratio",
            min_value=1.0, max_value=200.0, step=5.0,
            help="Genes with wider CI ratio are removed. NTP default: 40, EFSA: 50.",
            key="bmdu_bmdl_ratio",
        )
        st.checkbox(
            "Exclude BMD > highest dose",
            help="Remove genes with BMD above the highest tested concentration.",
            key="bmd_max_dose_filter",
        )
        st.number_input(
            "Extrapolation factor",
            min_value=1.0, max_value=100.0, step=1.0,
            help="Flag genes with BMD < lowest_dose / factor as extrapolated.",
            key="bmd_extrap_factor",
        )
    with ntp2:
        st.number_input(
            "Min fold-change",
            min_value=0.0, max_value=10.0, step=0.1,
            help="0 = disabled. Post-fit filter: remove genes with absolute fold-change below threshold.",
            key="fold_change_min",
        )
        st.checkbox(
            "Include MSigDB Hallmark + C2 gene sets",
            help="Auto-detected: included when msigdbr R package is installed. ~50 Hallmark and ~6,300 C2 curated gene sets.",
            key="use_msigdb",
        )

    st.divider()
    st.markdown("#### Resource Allocation")
    cpus = os.cpu_count() or 8
    ram  = f"{SYSTEM_RAM_GB:.0f}" if SYSTEM_RAM_GB > 0 else "?"
    st.markdown(f"""
<div class="box-sys">
System: {cpus} CPUs | {ram} GB RAM | Recommended: {RECOMMENDED_ALIGNER.upper()}<br><br>
Download: 1 core \u00d7 3 parallel &nbsp;|&nbsp;
Fastp: 1 core \u00d7 {max(1, cpus - 1)} parallel &nbsp;|&nbsp;
STAR: {max(2, cpus // 2)} cores \u00d7 2 parallel &nbsp;|&nbsp;
Post-QC: 1 core \u00d7 {max(1, cpus - 2)} parallel &nbsp;|&nbsp;
Final: {cpus} cores
</div>""", unsafe_allow_html=True)


# ═════════════════════════════════════════════════════════════════════════════
#  TAB 4 — RESULTS
# ═════════════════════════════════════════════════════════════════════════════

with tab_results:
    st.markdown("#### Analysis Runs")

    # Pipeline-level outputs
    multiqc_report = Path(outdir) / "multiqc" / "multiqc_report.html"
    count_matrix_file = Path(outdir) / "counts" / "gene_count_matrix.csv"

    if multiqc_report.exists() or count_matrix_file.exists():
        qc_cols = st.columns(3)
        if multiqc_report.exists():
            qc_cols[0].download_button(
                "\u2b07 MultiQC Report", multiqc_report.read_bytes(),
                file_name="multiqc_report.html", mime="text/html",
            )
        if count_matrix_file.exists():
            qc_cols[1].download_button(
                "\u2b07 Count Matrix", count_matrix_file.read_bytes(),
                file_name="gene_count_matrix.csv", mime="text/csv",
            )
        timeline = Path(outdir) / "timeline.html"
        if timeline.exists():
            qc_cols[2].download_button(
                "\u2b07 Timeline", timeline.read_bytes(),
                file_name="timeline.html", mime="text/html",
            )
        st.divider()

    runs = list_runs(outdir)

    if not runs:
        st.info("No analysis runs found yet. Run the pipeline first.")
    else:
        rc1, rc2, rc3, rc4 = st.columns([3, 1, 1, 1])
        with rc1:
            run_labels = []
            for r in runs:
                tags = []
                if r["has_deg"]:
                    tags.append("DESeq2")
                if r["has_dromics"]:
                    tags.append("DRomics")
                tag_str = f" [{', '.join(tags)}]" if tags else ""
                run_labels.append(f"{r['name']}  \u2014  {r['mtime']}{tag_str}")
            selected_idx = st.selectbox(
                "Select run", range(len(run_labels)),
                format_func=lambda i: run_labels[i],
                key="run_selector",
            )
        with rc2:
            if st.button("\U0001f4c2 Load", use_container_width=True, key="load_results"):
                st.session_state.show_results = True
        with rc3:
            if st.button("\U0001f504 Refresh", use_container_width=True, key="refresh_results"):
                st.session_state.show_results = False
                st.rerun()
        with rc4:
            if st.button("\U0001f5d1 Delete", use_container_width=True, key="delete_run"):
                shutil.rmtree(runs[selected_idx]["path"])
                st.session_state.show_results = False
                st.success(f"Deleted: {runs[selected_idx]['name']}")
                time.sleep(1)
                st.rerun()

        if st.session_state.show_results:
            selected_run = runs[selected_idx]
            st.markdown(
                f'<div class="box-info">\U0001f4c1 Viewing: <b>{selected_run["name"]}</b> '
                f'| {selected_run["mtime"]}</div>',
                unsafe_allow_html=True,
            )
            display_results(selected_run["path"])
        else:
            st.info(f"Found {len(runs)} run(s). Select one and click **Load** to view.")
