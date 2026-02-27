#!/usr/bin/env python3
"""
ARACRA-STAR Pipeline — Streamlit App
- Full Pipeline: SRA → fastp → STAR → QC → featureCounts/Salmon → DESeq2 → DRomics
- Direct Analysis: Upload count matrix → DESeq2/DRomics
- Tool choices: featureCounts vs Salmon, RSeQC vs Picard, HK vs Full coverage
- Per-run output folders with results browser
"""

import streamlit as st
import subprocess, os, signal, json, time, re, shutil
import pandas as pd
from pathlib import Path
from datetime import datetime

# ══════════════════════════════════════════════════════════════════════════════
# PATHS & CONFIG
# ══════════════════════════════════════════════════════════════════════════════
PIPELINE_DIR = Path(__file__).parent.absolute()
MAIN_NF      = PIPELINE_DIR / "main.nf"
ENV_FILE     = PIPELINE_DIR / ".env"
STATE_FILE   = PIPELINE_DIR / ".pipeline_state.json"
PID_FILE     = PIPELINE_DIR / ".pipeline.pid"

def load_env() -> dict:
    env = {}
    if ENV_FILE.exists():
        for line in ENV_FILE.read_text().splitlines():
            line = line.strip()
            if line and "=" in line and not line.startswith("#"):
                k, v = line.split("=", 1)
                env[k.strip()] = v.strip()
    return env

ENV = load_env()

CONDA_BASE   = ENV.get("CONDA_BASE", "")
ENV_NAME     = ENV.get("ENV_NAME", "rnaseq_pipeline2")
ENV_PATH     = ENV.get("ENV_PATH", "")
ENV_BIN      = ENV.get("ENV_BIN", "")
DEFAULT_WORK = ENV.get("WORK_DIR", str(Path.home() / "aracra_star_work/work"))
DEFAULT_OUT  = ENV.get("OUT_DIR",  str(Path.home() / "aracra_star_work/results"))
LOG_FILE     = Path(ENV.get("LOG_FILE", str(Path.home() / "aracra_star_work/pipeline.log")))

# Reference defaults
STAR_INDEX       = ENV.get("STAR_INDEX",       str(Path.home() / "databases/hg38_reference/star_index_hg38"))
HISAT2_INDEX     = ENV.get("HISAT2_INDEX",     str(Path.home() / "databases/hg38_reference/hisat2_index_hg38/genome_tran"))
SALMON_INDEX     = ENV.get("SALMON_INDEX",     str(Path.home() / "databases/hg38_reference/salmon_index_hg38"))
GTF_PATH         = ENV.get("GTF_PATH",         str(Path.home() / "databases/hg38_reference/gencode.v44.primary_assembly.annotation.gtf"))
BED_PATH         = ENV.get("BED_PATH",         str(Path.home() / "databases/annotations/hg38_RefSeq.bed"))
HK_BED_PATH      = ENV.get("HK_BED_PATH",      str(Path.home() / "databases/annotations/hg38_housekeeping.bed"))
REFFLAT_PATH     = ENV.get("REFFLAT_PATH",     str(Path.home() / "databases/annotations/hg38_refFlat.txt"))
SCREEN_CONF_PATH = ENV.get("SCREEN_CONF_PATH", str(Path.home() / "databases/fastq_screen_genomes/FastQ_Screen_Genomes/fastq_screen.conf"))
RECOMMENDED_ALIGNER = ENV.get("RECOMMENDED_ALIGNER", "star")

# System RAM detection for warnings
def get_ram_gb():
    try:
        mem = open("/proc/meminfo").read()
        return int(re.search(r'MemTotal:\s+(\d+)', mem).group(1)) / 1048576
    except Exception:
        return 0
SYSTEM_RAM_GB = get_ram_gb()

if not ENV_BIN and ENV_PATH:
    ENV_BIN = str(Path(ENV_PATH) / "bin")
elif not ENV_BIN and CONDA_BASE:
    ENV_BIN = str(Path(CONDA_BASE) / "envs" / ENV_NAME / "bin")

SETUP_DONE = ENV_FILE.exists() and bool(ENV_BIN)

# ══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════════════
def generate_run_name(treatment: str, control: str) -> str:
    def short(s):
        s = re.sub(r'[^A-Za-z0-9]+', '_', s.strip())[:20].strip('_')
        return s or "X"
    ts = datetime.now().strftime("%Y%m%d_%H%M")
    return f"{short(treatment)}_vs_{short(control)}_{ts}"

def list_runs(outdir: str) -> list:
    runs_dir = Path(outdir) / "runs"
    if not runs_dir.exists():
        return []
    runs = []
    for d in sorted(runs_dir.iterdir(), reverse=True):
        if d.is_dir():
            has_deg = (d / "deg_results" / "analysis_summary.json").exists()
            has_dr  = (d / "dromics_results" / "dromics_summary.json").exists()
            runs.append({
                "name": d.name, "path": d,
                "has_deg": has_deg, "has_dromics": has_dr,
                "mtime": datetime.fromtimestamp(d.stat().st_mtime).strftime("%Y-%m-%d %H:%M"),
            })
    return runs

def get_sysinfo() -> dict:
    info = {}
    try:
        mem = open("/proc/meminfo").read()
        total_kb = int(re.search(r'MemTotal:\s+(\d+)', mem).group(1))
        avail_kb = int(re.search(r'MemAvailable:\s+(\d+)', mem).group(1))
        info["ram_total"] = f"{total_kb/1048576:.1f} GB"
        info["ram_avail"] = f"{avail_kb/1048576:.1f} GB"
    except Exception:
        info["ram_total"] = info["ram_avail"] = "?"
    info["cpus"] = str(os.cpu_count() or "?")
    try:
        d = os.statvfs(DEFAULT_OUT if Path(DEFAULT_OUT).exists() else Path.home())
        info["disk"] = f"{(d.f_bavail * d.f_frsize) / 1e9:.0f} GB free"
    except Exception:
        info["disk"] = "?"
    return info

# ══════════════════════════════════════════════════════════════════════════════
# STATE MANAGEMENT
# ══════════════════════════════════════════════════════════════════════════════
def save_state(s): STATE_FILE.write_text(json.dumps(s, default=str))
def load_state():
    try: return json.loads(STATE_FILE.read_text()) if STATE_FILE.exists() else {}
    except Exception: return {}

def save_pid(pid): PID_FILE.write_text(str(pid))
def load_pid():
    try: return int(PID_FILE.read_text().strip()) if PID_FILE.exists() else None
    except Exception: return None
def clear_pid():
    if PID_FILE.exists(): PID_FILE.unlink()

def is_running() -> bool:
    pid = load_pid()
    if not pid: return False
    try:
        os.kill(pid, 0)
        cmd = open(f"/proc/{pid}/cmdline", "rb").read().decode(errors="replace")
        return "nextflow" in cmd.lower()
    except Exception: return False

def elapsed(state) -> str:
    if not state.get("start_time"): return ""
    try:
        s = int((datetime.now() - datetime.fromisoformat(state["start_time"])).total_seconds())
        return f"{s//3600:02d}:{(s%3600)//60:02d}:{s%60:02d}"
    except Exception: return ""

# ══════════════════════════════════════════════════════════════════════════════
# PROCESS ENVIRONMENT & LAUNCHER
# ══════════════════════════════════════════════════════════════════════════════
def get_process_env() -> dict:
    env = os.environ.copy()
    if ENV_BIN:
        env["PATH"] = f"{ENV_BIN}:{env.get('PATH', '')}"
    return env

def launch(cmd: list, work_dir: str) -> int:
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(LOG_FILE, "a") as f:
        f.write(f"\n{'='*60}\nRun started : {datetime.now()}\nCommand     : {' '.join(cmd)}\n{'='*60}\n\n")
    fh = open(LOG_FILE, "a")
    proc = subprocess.Popen(
        cmd, cwd=work_dir,
        stdout=fh, stderr=subprocess.STDOUT,
        start_new_session=True, close_fds=True,
        env=get_process_env()
    )
    fh.close()
    return proc.pid

def stop_pipeline():
    pid = load_pid()
    if not pid: return
    try:    os.killpg(os.getpgid(pid), signal.SIGTERM)
    except Exception:
        try: os.kill(pid, signal.SIGTERM)
        except Exception: pass
    clear_pid()

# ══════════════════════════════════════════════════════════════════════════════
# LOG READER & PROGRESS
# ══════════════════════════════════════════════════════════════════════════════
def log_tail(n=80) -> str:
    if not LOG_FILE.exists(): return "Log will appear here when pipeline starts..."
    try:
        lines = LOG_FILE.read_text(errors="replace").splitlines()
        return "\n".join(lines[-n:])
    except Exception as e: return f"Error reading log: {e}"

def parse_progress() -> dict:
    progress = {}
    if not LOG_FILE.exists(): return progress
    try:
        for line in LOG_FILE.read_text(errors="replace").splitlines()[-300:]:
            m = re.search(r'process\s*>\s*(\w+)[^\[]*\[(\d+)%\]\s+(\d+)\s+of\s+(\d+)', line)
            if m: progress[m.group(1)] = {"pct": int(m.group(2)), "done": int(m.group(3)), "total": int(m.group(4))}
    except Exception: pass
    return progress

def detect_final() -> str:
    if not LOG_FILE.exists(): return None
    try:
        tail = "\n".join(LOG_FILE.read_text(errors="replace").splitlines()[-20:])
        if "SUCCESS" in tail or "✅" in tail: return "completed"
        if "FAILED" in tail or "❌" in tail or "ERROR ~" in tail: return "failed"
    except Exception: pass
    return None

# ══════════════════════════════════════════════════════════════════════════════
# TOOL CHECKER
# ══════════════════════════════════════════════════════════════════════════════
TOOLS = {
    "nextflow":      ["nextflow", "-version"],
    "prefetch":      ["prefetch", "--version"],
    "fasterq-dump":  ["fasterq-dump", "--version"],
    "fastp":         ["fastp", "--version"],
    "STAR":          ["STAR", "--version"],
    "samtools":      ["samtools", "--version"],
    "featureCounts": ["featureCounts", "-v"],
    "salmon":        ["salmon", "--version"],
    "fastq_screen":  ["fastq_screen", "--version"],
    "Rscript":       ["Rscript", "--version"],
    "multiqc":       ["multiqc", "--version"],
}

@st.cache_data(ttl=60)
def check_tools():
    out = {}
    proc_env = get_process_env()
    for name, cmd in TOOLS.items():
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=10, env=proc_env)
            ver = (r.stdout + r.stderr).strip().split("\n")[0][:60]
            out[name] = {"ok": True, "ver": ver or "OK"}
        except FileNotFoundError:
            out[name] = {"ok": False, "ver": "NOT FOUND"}
        except Exception as e:
            out[name] = {"ok": False, "ver": str(e)[:60]}
    return out

@st.cache_data(ttl=300)
def check_r_packages():
    script = """
pkgs <- c("DESeq2","edgeR","sva","org.Hs.eg.db","dplyr","readxl",
          "ggplot2","ggrepel","stringr","jsonlite","DRomics","optparse")
miss <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
if (length(miss)) cat("MISSING:", paste(miss, collapse=",")) else cat("OK")
"""
    try:
        r = subprocess.run(["Rscript", "-e", script], capture_output=True, text=True,
                           timeout=30, env=get_process_env())
        out = (r.stdout + r.stderr).strip()
        if "MISSING:" in out:
            return False, out.replace("MISSING:", "").strip().split(",")
        return True, []
    except Exception as e:
        return False, [str(e)]

# ══════════════════════════════════════════════════════════════════════════════
# RESULTS DISPLAY
# ══════════════════════════════════════════════════════════════════════════════
def display_results(run_dir: Path):
    deg_dir = run_dir / "deg_results"
    dr_dir  = run_dir / "dromics_results"

    # ── DESeq2 ──
    st.markdown('<p class="sec">DESeq2 / R-ODAF Results</p>', unsafe_allow_html=True)
    sp = deg_dir / "analysis_summary.json"
    if sp.exists():
        s = json.loads(sp.read_text())
        def jv(d, k, df=0):
            v = d.get(k, df); return v[0] if isinstance(v, list) and v else v
        m1, m2, m3, m4, m5 = st.columns(5)
        m1.metric("Samples",         jv(s, "samples_analyzed"))
        m2.metric("Genes tested",    f"{jv(s, 'total_genes'):,}")
        m3.metric("DEGs (FDR<0.05)", f"{jv(s, 'sig_relaxed'):,}")
        m4.metric("⬆ Up",            f"{jv(s, 'upregulated'):,}")
        m5.metric("⬇ Down",          f"{jv(s, 'downregulated'):,}")
        st.caption(f"Batch corrected: {'Yes' if jv(s, 'batch_corrected') else 'No'} | "
                   f"Outliers removed: {jv(s, 'outliers_removed')}")
        pc, vc = st.columns(2)
        pca_img = deg_dir / "Final_PCA_Plot.png"
        vol_img = deg_dir / "Volcano_Plot.png"
        if pca_img.exists(): pc.image(str(pca_img), caption="PCA", use_container_width=True)
        if vol_img.exists(): vc.image(str(vol_img), caption="Volcano", use_container_width=True)
        sig_csv = deg_dir / "Custom_Filtered_DEGs.csv"
        if sig_csv.exists():
            st.dataframe(pd.read_csv(sig_csv).head(30), use_container_width=True, hide_index=True)
            dl1, dl2, dl3 = st.columns(3)
            for col_btn, fname in [(dl1, "All_Results.csv"), (dl2, "Custom_Filtered_DEGs.csv"),
                                   (dl3, "Upregulated_Genes.csv")]:
                fp = deg_dir / fname
                if fp.exists():
                    with open(fp, "rb") as f:
                        col_btn.download_button(
                            f"⬇ {fname.replace('_', ' ').replace('.csv', '')}",
                            f.read(), file_name=fname, mime="text/csv",
                            key=f"dl_{run_dir.name}_{fname}")
    else:
        st.info("No DESeq2 results in this run.")

    st.divider()

    # ── DRomics ──
    st.markdown('<p class="sec">DRomics / BMD Results</p>', unsafe_allow_html=True)
    dp = dr_dir / "dromics_summary.json"
    if dp.exists():
        ds = json.loads(dp.read_text())
        def jv2(d, k, df=0):
            v = d.get(k, df); return v[0] if isinstance(v, list) and v else v
        b1, b2, b3, b4 = st.columns(4)
        b1.metric("Genes analyzed",   f"{jv2(ds, 'total_genes_final'):,}")
        b2.metric("Responsive genes", f"{jv2(ds, 'genes_selected'):,}")
        b3.metric("Models fitted",    f"{jv2(ds, 'models_fitted'):,}")
        if ds.get("bmd_performed") and ds.get("bmd_summary"):
            bms = ds["bmd_summary"]
            def bv(k):
                v = bms.get(k, 0); return round(v[0] if isinstance(v, list) and v else v, 4)
            b4.metric("Median BMD", f"{bv('median_bmd')} μM")
        c1, c2 = st.columns(2)
        for ci, fname, cap in [(c1, "dose_response_curves.png", "Dose-Response"),
                                (c2, "bmd_distribution.png", "BMD Distribution")]:
            fp = dr_dir / fname
            if fp.exists(): ci.image(str(fp), caption=cap, use_container_width=True)
        fp = dr_dir / "top25_sensitive_genes.png"
        if fp.exists(): st.image(str(fp), caption="Top 25 Sensitive Genes", use_container_width=True)
        dl1, dl2, dl3 = st.columns(3)
        for col_btn, fname in [(dl1, "bmd_results.csv"), (dl2, "top50_sensitive_genes.csv"),
                                (dl3, "dose_response_models.csv")]:
            fp = dr_dir / fname
            if fp.exists():
                with open(fp, "rb") as f:
                    col_btn.download_button(
                        f"⬇ {fname.replace('_', ' ').replace('.csv', '')}",
                        f.read(), file_name=fname, mime="text/csv",
                        key=f"dl_{run_dir.name}_{fname}")
    else:
        st.info("No DRomics results in this run.")

# ══════════════════════════════════════════════════════════════════════════════
# PAGE CONFIG + CSS
# ══════════════════════════════════════════════════════════════════════════════
st.set_page_config(page_title="ARACRA-STAR Pipeline", page_icon="🧬",
                   layout="wide", initial_sidebar_state="expanded")
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;500;600&display=swap');
html,body,[class*="css"]{font-family:'IBM Plex Sans',sans-serif}
.header{background:#0f1923;border-left:4px solid #00d4aa;border-radius:0 8px 8px 0;padding:1.2rem 1.8rem;margin-bottom:1.5rem}
.header h1{font-family:'IBM Plex Mono',monospace;color:#00d4aa;font-size:1.6rem;margin:0;letter-spacing:-0.5px}
.header p{color:#8899aa;font-size:.85rem;margin:4px 0 0 0}
.sec{font-family:'IBM Plex Mono',monospace;font-size:.68rem;color:#00d4aa;text-transform:uppercase;letter-spacing:2px;margin:1rem 0 .3rem 0;border-bottom:1px solid #1e2d3d;padding-bottom:3px}
.tok{display:inline-block;background:#0f1923;border:1px solid #00d4aa44;color:#00d4aa;font-family:'IBM Plex Mono',monospace;font-size:.68rem;padding:1px 6px;border-radius:3px;margin:2px}
.s-idle{font-family:'IBM Plex Mono',monospace;color:#8899aa;font-size:.95rem;font-weight:600}
.s-run {font-family:'IBM Plex Mono',monospace;color:#4da6ff;font-size:.95rem;font-weight:600}
.s-ok  {font-family:'IBM Plex Mono',monospace;color:#00d4aa;font-size:.95rem;font-weight:600}
.s-err {font-family:'IBM Plex Mono',monospace;color:#ff4b4b;font-size:.95rem;font-weight:600}
.s-stp {font-family:'IBM Plex Mono',monospace;color:#ffaa00;font-size:.95rem;font-weight:600}
.tok-ok{color:#00d4aa;font-family:'IBM Plex Mono',monospace;font-size:.78rem;display:block}
.tok-err{color:#ff4b4b;font-family:'IBM Plex Mono',monospace;font-size:.78rem;display:block}
.box-info{background:#001a14;border:1px solid #00d4aa44;border-radius:6px;padding:.5rem .8rem;color:#00d4aa;font-size:.8rem;margin:3px 0}
.box-warn{background:#1a1200;border:1px solid #ffaa0066;border-radius:6px;padding:.5rem .8rem;color:#ffaa00;font-size:.8rem;margin:3px 0}
.box-det{background:#001428;border:1px solid #4da6ff55;border-radius:6px;padding:.5rem .8rem;color:#4da6ff;font-size:.8rem;margin:3px 0}
.box-sys{background:#0f1923;border:1px solid #1e2d3d;border-radius:6px;padding:.5rem .8rem;font-family:'IBM Plex Mono',monospace;font-size:.72rem;color:#8899aa;margin:3px 0}
.box-direct{background:#0a001a;border:1px solid #8855ff66;border-radius:6px;padding:.5rem .8rem;color:#aa88ff;font-size:.8rem;margin:3px 0}
.tool-group{background:#0f1923;border:1px solid #1e2d3d;border-radius:6px;padding:.6rem .8rem;margin:4px 0}
</style>
""", unsafe_allow_html=True)

# ══════════════════════════════════════════════════════════════════════════════
# SESSION STATE
# ══════════════════════════════════════════════════════════════════════════════
defaults = [
    ("accessions", []), ("meta_df", None), ("meta_path", None),
    ("show_results", False), ("chemicals", []), ("tools", None),
    ("fdr_strict", 0.01), ("fdr_relaxed", 0.05), ("log2fc", 1.0),
    ("read_thresh", 1000000), ("dr_fdr", 0.05), ("dr_criterion", "AIC"),
    ("perform_bmd", True), ("bmd_bootstrap", False),
    ("count_matrix_path", None), ("count_matrix_df", None),
]
for k, v in defaults:
    if k not in st.session_state: st.session_state[k] = v

state   = load_state()
running = is_running()
if state.get("status") == "running" and not running:
    final = detect_final()
    state["status"] = final or "failed"
    save_state(state)
    clear_pid()

# ══════════════════════════════════════════════════════════════════════════════
# HEADER
# ══════════════════════════════════════════════════════════════════════════════
st.markdown("""
<div class="header">
  <h1>🧬 ARACRA-STAR Pipeline</h1>
  <p>RNA-seq &amp; TempO-Seq Analysis &nbsp;|&nbsp;
     SRA → fastp → STAR → QC → featureCounts/Salmon → DESeq2 → DRomics/BMD</p>
</div>""", unsafe_allow_html=True)

# ══════════════════════════════════════════════════════════════════════════════
# SIDEBAR
# ══════════════════════════════════════════════════════════════════════════════
with st.sidebar:
    sysinfo = get_sysinfo()
    st.markdown(
        f'<div class="box-sys">🖥 RAM: {sysinfo["ram_total"]} ({sysinfo["ram_avail"]} free)<br>'
        f'⚙ CPUs: {sysinfo["cpus"]}<br>💾 Disk: {sysinfo["disk"]}</div>',
        unsafe_allow_html=True)

    if not SETUP_DONE:
        st.markdown(
            '<div class="box-warn">⚠ .env not found — run setup.sh first or configure paths below.</div>',
            unsafe_allow_html=True)

    # ── Reference Files ──
    st.markdown('<p class="sec">📂 Reference Files</p>', unsafe_allow_html=True)
    star_index   = st.text_input("STAR index",       value=STAR_INDEX,       key="star_idx")
    hisat2_index = st.text_input("HISAT2 index",     value=HISAT2_INDEX,     key="h2_idx")
    salmon_index = st.text_input("Salmon index",     value=SALMON_INDEX,     key="sal_idx")
    gtf          = st.text_input("GTF annotation",   value=GTF_PATH,         key="gtf_in")
    bed          = st.text_input("RefSeq BED",       value=BED_PATH,         key="bed_in")
    hk_bed       = st.text_input("Housekeeping BED", value=HK_BED_PATH,      key="hk_bed_in")
    refflat      = st.text_input("refFlat",          value=REFFLAT_PATH,     key="refflat_in")
    screen_conf  = st.text_input("FastQ Screen conf", value=SCREEN_CONF_PATH, key="screen_in")

    # ── Directories ──
    st.markdown('<p class="sec">📁 Directories</p>', unsafe_allow_html=True)
    outdir  = st.text_input("Output directory", value=DEFAULT_OUT,  key="out_in")
    workdir = st.text_input("Work directory",   value=DEFAULT_WORK, key="wk_in")

    # ── Executor ──
    st.markdown('<p class="sec">🖥 Executor</p>', unsafe_allow_html=True)
    profile = st.selectbox("Profile", ["standard", "slurm"])
    layout  = st.selectbox("Library layout", ["SE", "PE"])
    resume  = st.checkbox("Resume previous run", value=True)

    # ── Analysis ──
    st.markdown('<p class="sec">🔬 Downstream Analysis</p>', unsafe_allow_html=True)
    run_deg     = st.checkbox("DESeq2 / R-ODAF", value=True)
    run_dromics = st.checkbox("DRomics / BMD",    value=True)

    # ── Tool Check ──
    st.markdown('<p class="sec">🔧 System Check</p>', unsafe_allow_html=True)
    if st.button("▶ Check Tools", use_container_width=True):
        st.cache_data.clear()
        st.session_state.tools = check_tools()
    if st.session_state.tools:
        for name, info in st.session_state.tools.items():
            cls  = "tok-ok" if info["ok"] else "tok-err"
            icon = "✔" if info["ok"] else "✘"
            st.markdown(f'<span class="{cls}">{icon} {name}: {info["ver"][:40]}</span>',
                        unsafe_allow_html=True)
        r_ok, r_miss = check_r_packages()
        if r_ok:
            st.markdown('<span class="tok-ok">✔ R packages (12/12)</span>', unsafe_allow_html=True)
        else:
            st.markdown(f'<span class="tok-err">✘ R missing: {", ".join(r_miss)}</span>',
                        unsafe_allow_html=True)

    if running:
        st.divider()
        st.markdown(f'<span class="tok-ok">⚙ PID: {load_pid()}</span>', unsafe_allow_html=True)
        st.markdown('<div class="box-det">🔒 Detached — safe to close browser</div>',
                    unsafe_allow_html=True)

# ══════════════════════════════════════════════════════════════════════════════
# TABS
# ══════════════════════════════════════════════════════════════════════════════
tab1, tab_direct, tab2, tab3 = st.tabs([
    "🚀 Full Pipeline", "⚡ Direct Analysis", "⚙ Parameters", "📊 Results"])

smap = {
    "idle": ("○  IDLE", "s-idle"), "running": ("●  RUNNING", "s-run"),
    "completed": ("✔  COMPLETED", "s-ok"), "failed": ("✘  FAILED", "s-err"),
    "stopped": ("⏹  STOPPED", "s-stp")
}

# ─────────────────────────────────────────────────────────────────────────────
# TAB 1 — FULL PIPELINE
# ─────────────────────────────────────────────────────────────────────────────
with tab1:
    L, R = st.columns([3, 2], gap="large")

    with L:
        # ── SRA Accessions ──
        st.markdown('<p class="sec">📥 SRA Accessions / Metadata</p>', unsafe_allow_html=True)
        method = st.radio("Input method", ["Paste accessions", "Upload metadata (xlsx)"],
                          horizontal=True, label_visibility="collapsed")
        if method == "Paste accessions":
            raw = st.text_area("", height=120,
                               placeholder="SRR29684827\nSRR29684828\n...")
            if raw.strip():
                valid   = [a.strip() for a in raw.splitlines()
                           if a.strip() and re.match(r'^[SED]RR\d+$', a.strip())]
                invalid = [a.strip() for a in raw.splitlines()
                           if a.strip() and not re.match(r'^[SED]RR\d+$', a.strip())]
                st.session_state.accessions = valid
                if invalid: st.warning(f"Ignored: {', '.join(invalid)}")
        else:
            mf = st.file_uploader("Upload metadata Excel/CSV",
                                   type=["xlsx", "xls", "csv"], key="meta_full")
            if mf:
                try:
                    df = pd.read_csv(mf) if mf.name.endswith(".csv") else pd.read_excel(mf)
                    # Auto-detect SRR column
                    srr_col = None
                    for c in df.columns:
                        if df[c].astype(str).str.match(r'^[SED]RR\d+$').any():
                            srr_col = c; break
                    if srr_col:
                        st.session_state.accessions = df[srr_col].dropna().astype(str).str.strip().tolist()
                    # Save metadata
                    Path(outdir).mkdir(parents=True, exist_ok=True)
                    mp = Path(outdir) / "metadata_input.xlsx"
                    if mf.name.endswith(".csv"):
                        df.to_csv(str(mp).replace(".xlsx", ".csv"), index=False)
                        st.session_state.meta_path = str(mp).replace(".xlsx", ".csv")
                    else:
                        mf.seek(0)
                        mp.write_bytes(mf.read())
                        st.session_state.meta_path = str(mp)
                    st.session_state.meta_df = df
                    if "Treatment" in df.columns:
                        st.session_state.chemicals = sorted(df["Treatment"].unique().tolist())
                    st.success(f"✔ {len(df)} rows, {len(st.session_state.accessions)} SRR IDs found")
                    with st.expander("Preview"): st.dataframe(df, use_container_width=True, hide_index=True)
                except Exception as e:
                    st.error(f"Error: {e}")

        if st.session_state.accessions:
            n = len(st.session_state.accessions)
            tags = "".join(f'<span class="tok">{a}</span>' for a in st.session_state.accessions[:20])
            if n > 20: tags += f'<span class="tok">+{n-20} more</span>'
            st.markdown(tags, unsafe_allow_html=True)

        # ── Tool Choices ──
        st.markdown('<p class="sec">🔧 Pipeline Tool Choices</p>', unsafe_allow_html=True)

        # Aligner selection with RAM warning
        al1, al2 = st.columns([2, 3])
        with al1:
            st.markdown('<div class="tool-group"><b>Aligner</b></div>', unsafe_allow_html=True)
            default_idx = 1 if RECOMMENDED_ALIGNER == "hisat2" else 0
            aligner = st.radio("Alignment tool", ["STAR", "HISAT2"],
                               index=default_idx, horizontal=True, key="aligner_choice")
        with al2:
            if aligner == "STAR" and SYSTEM_RAM_GB < 32 and SYSTEM_RAM_GB > 0:
                st.markdown(
                    f'<div class="box-warn">⚠ <b>Low RAM detected: {SYSTEM_RAM_GB:.0f} GB</b><br>'
                    'STAR needs ~32 GB RAM for human genome. '
                    '<b>HISAT2 recommended</b> (needs only ~8 GB).</div>',
                    unsafe_allow_html=True)
            elif aligner == "STAR":
                st.markdown(
                    '<div class="box-info">✔ STAR — fast splice-aware aligner, needs ~32 GB RAM</div>',
                    unsafe_allow_html=True)
            else:
                st.markdown(
                    '<div class="box-info">✔ HISAT2 — memory-efficient (~8 GB RAM), good accuracy</div>',
                    unsafe_allow_html=True)

        tc1, tc2 = st.columns(2)
        with tc1:
            st.markdown('<div class="tool-group"><b>Quantification</b></div>', unsafe_allow_html=True)
            quant_tool = st.radio("Count method", ["featureCounts", "Salmon", "Both"],
                                  horizontal=True, key="quant_tool")

        with tc2:
            st.markdown('<div class="tool-group"><b>Gene Body Coverage</b></div>', unsafe_allow_html=True)
            coverage_mode = st.radio("Coverage", ["Housekeeping genes", "Full genome", "Both", "None"],
                                     horizontal=True, key="coverage_mode")

        tc3, tc4 = st.columns(2)
        with tc3:
            st.markdown('<div class="tool-group"><b>RNA Metrics</b></div>', unsafe_allow_html=True)
            rna_metrics = st.radio("RNA metrics tool", ["RSeQC (strandedness)", "Picard", "Both"],
                                   horizontal=True, key="rna_metrics")

        with tc4:
            st.markdown('<div class="tool-group"><b>Additional QC</b></div>', unsafe_allow_html=True)
            run_screen   = st.checkbox("FastQ Screen (contamination)", value=True, key="run_screen")
            run_read_dist = st.checkbox("Read distribution (RSeQC)", value=False, key="run_read_dist")

        # ── Comparison ──
        st.markdown('<p class="sec">🔬 Comparison</p>', unsafe_allow_html=True)
        c1, c2 = st.columns(2)
        chems = st.session_state.chemicals
        with c1:
            treatment = st.selectbox("Treatment", chems, key="treat_full") if chems else \
                        st.text_input("Treatment", value="Bisphenol A", key="treat_full_txt")
        with c2:
            control = st.selectbox("Control", chems, index=min(1, len(chems)-1), key="ctrl_full") if chems else \
                      st.text_input("Control", value="Control", key="ctrl_full_txt")

        st.markdown('<p class="sec">🏷 Run Name</p>', unsafe_allow_html=True)
        auto_name = generate_run_name(treatment, control)
        run_name_full = st.text_input("Run name", value=auto_name, key="run_name_full")

    with R:
        st.markdown('<p class="sec">📡 Status</p>', unsafe_allow_html=True)
        status = "running" if running else state.get("status", "idle")
        lbl, css = smap.get(status, ("○  IDLE", "s-idle"))
        st.markdown(f'<p class="{css}">{lbl}</p>', unsafe_allow_html=True)
        if status == "running":
            st.caption(f"Elapsed: {elapsed(state)}")
            if state.get("run_name"):
                st.caption(f"Run: {state['run_name']}")

        prog = parse_progress()
        if prog:
            st.markdown('<p class="sec">⚡ Progress</p>', unsafe_allow_html=True)
            labels = {
                "PARSE_METADATA": "0. Parse metadata",
                "DOWNLOAD_SRA": "1. Download SRA",
                "FASTP": "2. Fastp QC/Trim",
                "FASTQ_SCREEN": "3. FastQ Screen",
                "STAR_LOAD_GENOME": "4a. Load genome",
                "STAR_ALIGN": "4b. STAR alignment",
                "STAR_REMOVE_GENOME": "4c. Unload genome",
                "SAMTOOLS_INDEX": "5a. Index BAMs",
                "SAMTOOLS_FLAGSTAT": "5b. Flagstat",
                "SAMTOOLS_IDXSTATS": "5c. Idxstats",
                "INFER_STRANDEDNESS": "6. Strandedness",
                "GENEBODY_COVERAGE_HK": "7a. Coverage (HK)",
                "GENEBODY_COVERAGE_FULL": "7b. Coverage (full)",
                "READ_DISTRIBUTION": "8. Read distribution",
                "PICARD_RNA_METRICS": "9. Picard metrics",
                "FEATURECOUNTS": "10. featureCounts",
                "SALMON_QUANT": "11. Salmon",
                "MAKE_COUNT_MATRIX": "12. Count matrix",
                "MERGE_SALMON": "12b. Merge Salmon",
                "MULTIQC": "13. MultiQC",
                "DESEQ2_ANALYSIS": "14. DESeq2",
                "DROMICS_ANALYSIS": "15. DRomics/BMD",
            }
            for proc, info in prog.items():
                st.progress(info["pct"] / 100,
                            text=f"{labels.get(proc, proc)} — {info['done']}/{info['total']}")

        st.markdown('<p class="sec">✅ Validation</p>', unsafe_allow_html=True)
        issues = []
        if not MAIN_NF.exists():             issues.append("main.nf not found")
        if not st.session_state.accessions:  issues.append("No SRR accessions loaded")
        if not st.session_state.meta_path:   issues.append("No metadata uploaded")
        if not Path(gtf).exists():           issues.append(f"GTF not found: {gtf}")
        if aligner == "STAR" and not Path(star_index).is_dir():
            issues.append(f"STAR index not found: {star_index}")
        if aligner == "HISAT2" and not Path(hisat2_index + ".1.ht2").exists():
            issues.append(f"HISAT2 index not found: {hisat2_index}")
        if aligner == "STAR" and SYSTEM_RAM_GB < 32 and SYSTEM_RAM_GB > 0:
            issues.append(f"STAR needs ≥32 GB RAM (you have {SYSTEM_RAM_GB:.0f} GB) — switch to HISAT2")
        if quant_tool in ("Salmon", "Both") and not Path(salmon_index).is_dir():
            issues.append(f"Salmon index not found: {salmon_index}")
        if treatment and control and treatment == control:
            issues.append("Treatment = Control")

        if not issues:
            st.markdown('<div class="box-info">✔ All checks passed — ready to run</div>',
                        unsafe_allow_html=True)
        else:
            for i in issues:
                st.markdown(f'<div class="box-warn">⚠ {i}</div>', unsafe_allow_html=True)

        # ── Controls ──
        st.markdown('<p class="sec">🎮 Controls</p>', unsafe_allow_html=True)

        if st.button("🚀 Start Full Pipeline",
                     disabled=(bool(issues) or running),
                     use_container_width=True, type="primary"):

            Path(outdir).mkdir(parents=True, exist_ok=True)
            acc_file = Path(outdir) / "accessions.txt"
            acc_file.write_text("\n".join(st.session_state.accessions) + "\n")

            meta_path_to_use = st.session_state.meta_path or str(acc_file)

            use_fc       = quant_tool in ("featureCounts", "Both")
            use_salmon   = quant_tool in ("Salmon", "Both")
            use_cov_hk   = coverage_mode in ("Housekeeping genes", "Both")
            use_cov_full = coverage_mode in ("Full genome", "Both")
            use_strand   = rna_metrics in ("RSeQC (strandedness)", "Both")
            use_picard   = rna_metrics in ("Picard", "Both")

            # Write params to JSON file — avoids CLI quoting issues
            params_dict = {
                "metadata":           meta_path_to_use,
                "layout":             layout,
                "outdir":             outdir,
                "aligner":            aligner.lower(),
                "star_index":         star_index,
                "hisat2_index":       hisat2_index,
                "salmon_index":       salmon_index,
                "gtf":                gtf,
                "bed":                bed,
                "housekeeping_bed":   hk_bed,
                "refflat":            refflat,
                "fastq_screen_conf":  screen_conf,
                "run_download":       True,
                "run_fastp":          True,
                "run_fastq_screen":   run_screen,
                "run_star":           aligner == "STAR",
                "run_hisat2":         aligner == "HISAT2",
                "run_samtools_stats": True,
                "run_strandedness":   use_strand,
                "run_coverage_hk":    use_cov_hk,
                "run_coverage_full":  use_cov_full,
                "run_read_distribution": run_read_dist,
                "run_picard":         use_picard,
                "run_featurecounts":  use_fc,
                "run_salmon":         use_salmon,
                "run_merge_counts":   True,
                "run_multiqc":        True,
                "run_name":           run_name_full,
                "treatment":          treatment,
                "control":            control,
                "run_deg":            run_deg,
                "run_dromics":        run_dromics,
                "fdr_strict":         st.session_state.fdr_strict,
                "fdr_relaxed":        st.session_state.fdr_relaxed,
                "log2fc_threshold":   st.session_state.log2fc,
                "read_threshold":     st.session_state.read_thresh,
                "dr_fdr":             st.session_state.dr_fdr,
                "dr_criterion":       st.session_state.dr_criterion,
                "perform_bmd":        st.session_state.perform_bmd,
                "bmd_bootstrap":      st.session_state.bmd_bootstrap,
            }
            params_file = Path(outdir) / "nextflow_params.json"
            params_file.write_text(json.dumps(params_dict, indent=2))

            cmd = ["nextflow", "run", str(MAIN_NF),
                   "-profile", profile, "-w", workdir, "-ansi-log", "false",
                   "-params-file", str(params_file)]
            if resume: cmd.append("-resume")

            pid = launch(cmd, str(PIPELINE_DIR))
            save_pid(pid)
            save_state({"status": "running", "start_time": datetime.now().isoformat(),
                        "pid": pid, "run_name": run_name_full})
            st.success(f"✔ Pipeline launched → runs/{run_name_full}/")
            time.sleep(1)
            st.rerun()

        b1, b2 = st.columns(2)
        with b1:
            if st.button("⏹ Stop", disabled=not running, use_container_width=True, key="stop_full"):
                stop_pipeline()
                save_state({**state, "status": "stopped"})
                st.rerun()
        with b2:
            if st.button("↺ Reset", disabled=running, use_container_width=True, key="reset_full"):
                save_state({"status": "idle"})
                clear_pid()
                st.rerun()

    # ── Live Log ──
    if status in ("running", "completed", "failed", "stopped"):
        st.markdown('<p class="sec">📜 Live Log</p>', unsafe_allow_html=True)
        st.text_area("", value=log_tail(80), height=320, disabled=True, key="log_box_full")
        c1, c2, _ = st.columns([1, 2, 4])
        with c1:
            if LOG_FILE.exists():
                with open(LOG_FILE, "rb") as f:
                    st.download_button("⬇ Full log", f.read(), file_name="pipeline.log")
        if running:
            time.sleep(3)
            st.rerun()

# ─────────────────────────────────────────────────────────────────────────────
# TAB DIRECT — DIRECT ANALYSIS (skip to DESeq2/DRomics)
# ─────────────────────────────────────────────────────────────────────────────
with tab_direct:
    st.markdown(
        '<div class="box-direct">⚡ <b>Direct Analysis</b> — Upload a count matrix '
        'and metadata to skip straight to DESeq2 / DRomics.</div>',
        unsafe_allow_html=True)

    dL, dR = st.columns([3, 2], gap="large")

    with dL:
        st.markdown('<p class="sec">📊 Count Matrix</p>', unsafe_allow_html=True)
        st.caption("Accepts: featureCounts output, gene_count_matrix.csv, .tsv, .xlsx")
        count_file = st.file_uploader("Upload count matrix",
            type=["out", "csv", "tsv", "txt", "xlsx", "xls"],
            label_visibility="collapsed", key="count_upload")

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
                saved_path = Path(outdir) / "uploaded_counts.out"
                df_counts.to_csv(saved_path, sep="\t", index=False)
                st.session_state.count_matrix_path = str(saved_path)
                st.success(f"✔ {df_counts.shape[0]} genes × {df_counts.shape[1]} columns")
                with st.expander("Preview"):
                    st.dataframe(df_counts.head(20), use_container_width=True, hide_index=True)
            except Exception as e:
                st.error(f"Error: {e}")

        st.markdown('<p class="sec">📋 Metadata</p>', unsafe_allow_html=True)
        st.caption("Required columns: Sample_Name · Treatment · Dose · Batch · Type")
        meta_direct = st.file_uploader("Metadata CSV or XLSX",
            type=["csv", "xlsx", "xls"], label_visibility="collapsed", key="meta_direct")

        if meta_direct:
            try:
                df_m = pd.read_csv(meta_direct) if meta_direct.name.endswith(".csv") else pd.read_excel(meta_direct)
                req = ["Sample_Name", "Treatment", "Dose", "Batch", "Type"]
                miss = [c for c in req if c not in df_m.columns]
                if miss:
                    st.error(f"Missing columns: {', '.join(miss)}")
                else:
                    Path(outdir).mkdir(parents=True, exist_ok=True)
                    mp = Path(outdir) / "metadata_input.csv"
                    df_m.to_csv(mp, index=False)
                    st.session_state.meta_df   = df_m
                    st.session_state.meta_path = str(mp)
                    st.session_state.chemicals = sorted(df_m["Treatment"].unique().tolist())
                    st.success(f"✔ {len(df_m)} samples | {df_m['Treatment'].nunique()} treatments")
                    with st.expander("Preview"):
                        st.dataframe(df_m, use_container_width=True, hide_index=True)
            except Exception as e:
                st.error(f"Error: {e}")

        st.markdown('<p class="sec">🔬 Comparison</p>', unsafe_allow_html=True)
        dc1, dc2 = st.columns(2)
        chems_d = st.session_state.chemicals
        with dc1:
            treatment_d = st.selectbox("Treatment", chems_d, key="treat_direct") if chems_d else \
                          st.text_input("Treatment", value="Bisphenol A", key="treat_direct_txt")
        with dc2:
            control_d = st.selectbox("Control", chems_d, index=min(1, len(chems_d)-1), key="ctrl_direct") if chems_d else \
                        st.text_input("Control", value="Control", key="ctrl_direct_txt")

        st.markdown('<p class="sec">🏷 Run Name</p>', unsafe_allow_html=True)
        auto_name_d = generate_run_name(treatment_d, control_d)
        run_name_direct = st.text_input("Run name", value=auto_name_d, key="run_name_direct")

    with dR:
        st.markdown('<p class="sec">📡 Status</p>', unsafe_allow_html=True)
        status_d = "running" if running else state.get("status", "idle")
        lbl_d, css_d = smap.get(status_d, ("○  IDLE", "s-idle"))
        st.markdown(f'<p class="{css_d}">{lbl_d}</p>', unsafe_allow_html=True)
        if status_d == "running":
            st.caption(f"Elapsed: {elapsed(state)}")

        st.markdown('<p class="sec">✅ Validation</p>', unsafe_allow_html=True)
        d_issues = []
        if not MAIN_NF.exists():                    d_issues.append("main.nf not found")
        if not st.session_state.count_matrix_path:  d_issues.append("No count matrix uploaded")
        if not st.session_state.meta_path:          d_issues.append("No metadata uploaded")
        if not run_deg and not run_dromics:          d_issues.append("Enable DESeq2 and/or DRomics")
        if treatment_d and control_d and treatment_d == control_d:
            d_issues.append("Treatment = Control")

        if not d_issues:
            st.markdown('<div class="box-info">✔ Ready — click Start</div>',
                        unsafe_allow_html=True)
        else:
            for i in d_issues:
                st.markdown(f'<div class="box-warn">⚠ {i}</div>', unsafe_allow_html=True)

        st.markdown('<p class="sec">🎮 Controls</p>', unsafe_allow_html=True)

        if st.button("⚡ Start Direct Analysis",
                     disabled=(bool(d_issues) or running),
                     use_container_width=True, type="primary"):

            # Write params to JSON — avoids CLI quoting issues
            params_d = {
                "count_matrix":     st.session_state.count_matrix_path,
                "metadata":         st.session_state.meta_path,
                "outdir":           outdir,
                "run_name":         run_name_direct,
                "treatment":        treatment_d,
                "control":          control_d,
                "run_deg":          run_deg,
                "run_dromics":      run_dromics,
                "fdr_strict":       st.session_state.fdr_strict,
                "fdr_relaxed":      st.session_state.fdr_relaxed,
                "log2fc_threshold": st.session_state.log2fc,
                "read_threshold":   st.session_state.read_thresh,
                "dr_fdr":           st.session_state.dr_fdr,
                "dr_criterion":     st.session_state.dr_criterion,
                "perform_bmd":      st.session_state.perform_bmd,
                "bmd_bootstrap":    st.session_state.bmd_bootstrap,
            }
            Path(outdir).mkdir(parents=True, exist_ok=True)
            params_file = Path(outdir) / "nextflow_params.json"
            params_file.write_text(json.dumps(params_d, indent=2))

            cmd = ["nextflow", "run", str(MAIN_NF),
                   "-profile", profile, "-w", workdir, "-ansi-log", "false",
                   "-params-file", str(params_file)]

            pid = launch(cmd, str(PIPELINE_DIR))
            save_pid(pid)
            save_state({"status": "running", "start_time": datetime.now().isoformat(),
                        "pid": pid, "mode": "direct", "run_name": run_name_direct})
            st.success(f"✔ Direct analysis launched → runs/{run_name_direct}/")
            time.sleep(1)
            st.rerun()

        db1, db2 = st.columns(2)
        with db1:
            if st.button("⏹ Stop ", disabled=not running, use_container_width=True, key="stop_direct"):
                stop_pipeline()
                save_state({**state, "status": "stopped"})
                st.rerun()
        with db2:
            if st.button("↺ Reset ", disabled=running, use_container_width=True, key="reset_direct"):
                save_state({"status": "idle"})
                clear_pid()
                st.rerun()

    if status_d in ("running", "completed", "failed", "stopped"):
        st.markdown('<p class="sec">📜 Live Log</p>', unsafe_allow_html=True)
        st.text_area("", value=log_tail(80), height=320, disabled=True, key="log_box_direct")
        if running:
            time.sleep(3)
            st.rerun()

# ─────────────────────────────────────────────────────────────────────────────
# TAB 2 — PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────
with tab2:
    st.markdown('<p class="sec">DESeq2 / R-ODAF</p>', unsafe_allow_html=True)
    d1, d2, d3 = st.columns(3)
    with d1:
        st.session_state.fdr_strict  = st.number_input("Strict FDR",  value=st.session_state.fdr_strict,  min_value=0.001, max_value=0.1, step=0.001, format="%.3f")
        st.session_state.fdr_relaxed = st.number_input("Relaxed FDR", value=st.session_state.fdr_relaxed, min_value=0.01,  max_value=0.2, step=0.01)
    with d2:
        st.session_state.log2fc      = st.number_input("log2FC threshold",   value=st.session_state.log2fc,      min_value=0.0, max_value=5.0, step=0.1)
        st.session_state.read_thresh = st.number_input("Min reads / sample", value=st.session_state.read_thresh, min_value=0, step=100000)
    with d3:
        st.number_input("CPM threshold", value=1.0, min_value=0.0, step=0.1, key="cpm_thresh")
        st.slider("Min proportion expressed", 0.0, 1.0, 0.75, 0.05, key="min_prop")

    st.divider()
    st.markdown('<p class="sec">DRomics / BMD</p>', unsafe_allow_html=True)
    b1, b2, b3 = st.columns(3)
    with b1:
        st.session_state.dr_fdr       = st.number_input("DRomics FDR",    value=st.session_state.dr_fdr, min_value=0.001, max_value=0.2, step=0.01)
        st.session_state.dr_criterion = st.selectbox("Model criterion", ["AIC", "BIC"],
                                                      index=0 if st.session_state.dr_criterion == "AIC" else 1)
    with b2:
        st.session_state.perform_bmd   = st.checkbox("Perform BMD",          value=st.session_state.perform_bmd)
        st.session_state.bmd_bootstrap = st.checkbox("Bootstrap CIs (slow)", value=st.session_state.bmd_bootstrap)
    with b3:
        df_meta = st.session_state.meta_df
        if df_meta is not None and "Dose" in df_meta.columns:
            doses = sorted(df_meta["Dose"].unique())
            st.markdown("**Dose levels:**")
            st.write(", ".join(f"{d} μM" for d in doses))
            st.dataframe(df_meta.groupby("Dose").size().rename("n"), use_container_width=True)

    st.divider()
    st.markdown('<p class="sec">Pipeline Architecture</p>', unsafe_allow_html=True)
    cpus = os.cpu_count() or 8
    ram  = f"{SYSTEM_RAM_GB:.0f}" if SYSTEM_RAM_GB > 0 else "?"
    st.markdown(f"""
    <div class="box-sys">
    System: {cpus} CPUs | {ram} GB RAM | Recommended: {RECOMMENDED_ALIGNER.upper()}<br><br>
    Core strategy ({cpus} cores):<br>
    &nbsp;&nbsp;Download: 1 core × 3 parallel (network-bound)<br>
    &nbsp;&nbsp;Fastp: 1 core × {max(1,cpus-1)} parallel<br>
    &nbsp;&nbsp;FastQ Screen: 2 cores (independent)<br>
    &nbsp;&nbsp;STAR: {max(2,cpus//2)} cores × 2 parallel (shared genome memory, ~32 GB RAM)<br>
    &nbsp;&nbsp;HISAT2: {max(2,cpus//2)} cores × 2 parallel (~8 GB RAM)<br>
    &nbsp;&nbsp;Post-align QC: 1 core × {max(1,cpus-2)} parallel<br>
    &nbsp;&nbsp;featureCounts/Salmon: {cpus} cores (all BAMs together)
    </div>
    """, unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────────────────────
# TAB 3 — RESULTS
# ─────────────────────────────────────────────────────────────────────────────
with tab3:
    st.markdown('<p class="sec">📊 Analysis Runs</p>', unsafe_allow_html=True)

    # Also look for pipeline QC results
    multiqc_report = Path(outdir) / "multiqc" / "multiqc_report.html"
    count_matrix   = Path(outdir) / "counts" / "gene_count_matrix.csv"

    if multiqc_report.exists() or count_matrix.exists():
        st.markdown('<p class="sec">🧪 Pipeline QC Outputs</p>', unsafe_allow_html=True)
        qc1, qc2, qc3 = st.columns(3)
        if multiqc_report.exists():
            with open(multiqc_report, "rb") as f:
                qc1.download_button("⬇ MultiQC Report", f.read(),
                                    file_name="multiqc_report.html", mime="text/html")
        if count_matrix.exists():
            with open(count_matrix, "rb") as f:
                qc2.download_button("⬇ Count Matrix", f.read(),
                                    file_name="gene_count_matrix.csv", mime="text/csv")
            df_preview = pd.read_csv(count_matrix, nrows=10)
            st.caption(f"Count matrix: {df_preview.shape[1]-1} samples")

        timeline_file = Path(outdir) / "timeline.html"
        if timeline_file.exists():
            with open(timeline_file, "rb") as f:
                qc3.download_button("⬇ Timeline", f.read(),
                                    file_name="timeline.html", mime="text/html")
        st.divider()

    # DESeq2/DRomics runs
    runs = list_runs(outdir)

    if not runs:
        st.info("No DESeq2/DRomics runs found yet. Run the pipeline or direct analysis first.")
    else:
        rc1, rc2, rc3, rc4 = st.columns([3, 1, 1, 1])
        with rc1:
            run_labels = []
            for r in runs:
                tags = []
                if r["has_deg"]:     tags.append("DESeq2")
                if r["has_dromics"]: tags.append("DRomics")
                tag_str = f" [{', '.join(tags)}]" if tags else ""
                run_labels.append(f"{r['name']}  —  {r['mtime']}{tag_str}")

            selected_idx = st.selectbox("Select run", range(len(run_labels)),
                                        format_func=lambda i: run_labels[i],
                                        key="run_selector")
        with rc2:
            if st.button("📂 Load", use_container_width=True, key="load_results"):
                st.session_state.show_results = True
        with rc3:
            if st.button("🔄 Refresh", use_container_width=True, key="refresh_results"):
                st.session_state.show_results = False
                st.rerun()
        with rc4:
            if st.button("🗑 Delete", use_container_width=True, key="delete_run"):
                shutil.rmtree(runs[selected_idx]["path"])
                st.session_state.show_results = False
                st.success(f"Deleted: {runs[selected_idx]['name']}")
                time.sleep(1)
                st.rerun()

        if st.session_state.show_results:
            selected_run = runs[selected_idx]
            st.markdown(
                f'<div class="box-info">📁 Viewing: <b>{selected_run["name"]}</b> &nbsp;|&nbsp; '
                f'{selected_run["mtime"]}</div>',
                unsafe_allow_html=True)
            display_results(selected_run["path"])
        else:
            st.info(f"Found {len(runs)} run(s). Select one and click **📂 Load** to view.")
