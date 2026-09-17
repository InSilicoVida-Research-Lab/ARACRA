#!/usr/bin/env python3
"""
complexity_report.py — where the complexity in ARACRA actually lives.

Not a code-quality score. The questions it answers are practical:
  * which files would a new maintainer have to read?
  * how many branch points does a run pass through before producing a number?
  * how many knobs can a user turn, and are they wired end to end?
  * what does the install cost in time, disk and external dependencies?

Run:  python3 tests/complexity_report.py [--repo PATH] [--json]
"""

import argparse
import json
import re
import sys
from pathlib import Path


def strip_comments(text, style="hash"):
    if style == "hash":
        text = re.sub(r"(?m)^\s*#.*$", "", text)
    elif style == "slash":
        text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
        text = re.sub(r"(?m)^\s*//.*$", "", text)
    return text


def measure(path, style):
    raw = path.read_text(errors="replace")
    lines = raw.splitlines()
    code = [l for l in strip_comments(raw, style).splitlines() if l.strip()]
    comment = len([l for l in lines if l.strip().startswith(("#", "//", "*"))])
    # Branch points: a rough cyclomatic proxy, comparable across languages.
    branches = len(re.findall(
        r"\b(if|elif|else if|for|while|case|when|catch|except|&&|\|\|)\b",
        strip_comments(raw, style)))
    return {
        "total": len(lines),
        "code": len(code),
        "comment": comment,
        "branches": branches,
    }


FILES = [
    ("ARACRA/aracra_star_app.py",        "python", "Streamlit GUI"),
    ("ARACRA/main.nf",                   "slash",  "Nextflow workflow"),
    ("ARACRA/setup.sh",                  "hash",   "Installer"),
    ("ARACRA/run_app.sh",                "hash",   "Launcher"),
    ("ARACRA/lib/aracra_common.sh",      "hash",   "Shared defaults"),
    ("ARACRA/nextflow.config",           "slash",  "Config template"),
    ("ARACRA/scripts/run_dromics.R",     "hash",   "DRomics / BMD / tPOD"),
    ("ARACRA/scripts/run_deseq2.R",      "hash",   "Differential expression"),
    ("ARACRA/scripts/utils.R",           "hash",   "Shared R helpers"),
    ("ARACRA/scripts/run_qc_check.R",    "hash",   "PCA QC"),
    ("ARACRA/scripts/qc_outlier_check.py", "python", "Grubbs + IQR outliers"),
    ("ARACRA/scripts/temposeq_counts.py",  "python", "TempO-Seq aggregation"),
    ("ARACRA/scripts/parse_temposeq_manifest.py", "python", "TempO-Seq manifest"),
]

# External tools setup.sh installs, with approximate install footprint.
EXTERNAL = [
    ("nextflow", "workflow engine"), ("STAR", "aligner"), ("hisat2", "aligner"),
    ("bowtie2", "used by FastQ Screen"), ("samtools", "BAM handling"),
    ("fastp", "trimming"), ("fastq-screen", "contamination"),
    ("fastqc", "read QC"), ("rseqc", "post-align QC"), ("picard", "RNA metrics"),
    ("qualimap", "BAM QC"), ("salmon", "quantification"),
    ("subread", "featureCounts"), ("multiqc", "QC aggregation"),
    ("sra-tools", "SRA download"), ("pigz", "compression"),
]
R_PKGS = ["DESeq2", "edgeR", "sva", "DRomics", "rtracklayer", "org.Hs.eg.db",
          "GO.db", "KEGGREST", "AnnotationDbi", "optparse", "dplyr", "readxl",
          "ggplot2", "ggrepel", "stringr", "jsonlite"]

DOWNLOADS = [
    ("hg38 genome + GTF",   15, "always"),
    ("STAR index",          28, "unless --skip-index or RAM < 32 GB"),
    ("HISAT2 index",         5, "prebuilt, downloaded"),
    ("Salmon index",         3, "unless --skip-index"),
    ("FastQ Screen genomes", 14, "unless --skip-screen"),
    ("conda environment",    6, "always"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=str(Path(__file__).resolve().parent.parent))
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    root = Path(a.repo)

    rows, totals = [], {"total": 0, "code": 0, "comment": 0, "branches": 0}
    for rel, style, role in FILES:
        p = root / rel
        if not p.exists():
            continue
        m = measure(p, style)
        for k in totals:
            totals[k] += m[k]
        rows.append({"file": rel, "role": role, **m})

    nf = (root / "ARACRA" / "main.nf").read_text()
    app = (root / "ARACRA" / "aracra_star_app.py").read_text()
    dro = (root / "ARACRA" / "scripts" / "run_dromics.R").read_text()

    pipeline = {
        "nextflow_processes": len(re.findall(r"^process\s+\w+", nf, re.M)),
        "nextflow_params": len(set(re.findall(r"params\.(\w+)", nf))),
        "step_toggles": len(set(re.findall(r"params\.(run_\w+)", nf))),
        "entry_modes": 3,          # full pipeline / direct analysis / QC-only
        "platforms": 2,            # RNA-seq / TempO-Seq
        "aligners": 2,             # STAR / HISAT2
        "gui_widgets": len(re.findall(r"st\.(number_input|selectbox|checkbox|radio|text_input|file_uploader)\(", app)),
        "dromics_options": len(re.findall(r'make_option\("--(\w+)"', dro)),
    }

    unwired = sorted(o for o in re.findall(r'make_option\("--(\w+)"', dro)
                     if f"--{o}" not in nf)

    disk_full = sum(d[1] for d in DOWNLOADS)
    disk_lean = sum(d[1] for d in DOWNLOADS if "skip" not in d[2])

    report = {
        "code": {"files": rows, "totals": totals},
        "pipeline": pipeline,
        "unwired_dromics_options": unwired,
        "dependencies": {"external_tools": len(EXTERNAL), "r_packages": len(R_PKGS)},
        "install": {"disk_full_gb": disk_full, "disk_lean_gb": disk_lean,
                    "downloads": DOWNLOADS},
    }

    if a.json:
        print(json.dumps(report, indent=2))
        return

    print("\n\033[1mARACRA — complexity report\033[0m\n")

    print("  \033[1mCode size\033[0m")
    print(f"    {'file':<44} {'code':>6} {'cmt':>5} {'branch':>7}  role")
    for r in sorted(rows, key=lambda x: -x["code"]):
        print(f"    {r['file']:<44} {r['code']:>6} {r['comment']:>5} "
              f"{r['branches']:>7}  {r['role']}")
    print(f"    {'TOTAL':<44} {totals['code']:>6} {totals['comment']:>5} "
          f"{totals['branches']:>7}")
    ratio = totals["comment"] / max(totals["code"], 1)
    print(f"\n    comment-to-code ratio: {ratio:.2f}")

    top = sorted(rows, key=lambda x: -x["code"])[:3]
    share = sum(t["code"] for t in top) / max(totals["code"], 1)
    print(f"    top 3 files hold {share:.0%} of the code "
          f"({', '.join(Path(t['file']).name for t in top)})")

    print("\n  \033[1mPipeline surface\033[0m")
    for k, v in pipeline.items():
        print(f"    {k.replace('_', ' '):<28} {v}")
    paths = (pipeline["platforms"] * pipeline["aligners"]) + pipeline["entry_modes"] - 1
    print(f"    {'distinct run paths':<28} {paths}   "
          "(RNA-seq/TempO-Seq x STAR/HISAT2, plus direct & QC-only)")

    if unwired:
        print(f"\n    \033[1;33mnot passed by main.nf\033[0m (CLI-only): {', '.join(unwired)}")

    print("\n  \033[1mDependencies\033[0m")
    print(f"    external tools               {len(EXTERNAL)}")
    print(f"    R / Bioconductor packages    {len(R_PKGS)}")
    print(f"    every one is a version that can drift; "
          f"see aracra_versions.txt after setup")

    print("\n  \033[1mInstall footprint\033[0m")
    for name, gb, note in DOWNLOADS:
        print(f"    {name:<26} {gb:>3} GB   {note}")
    print(f"    {'full install':<26} {disk_full:>3} GB")
    print(f"    {'leanest viable':<26} {disk_lean:>3} GB   "
          "(--skip-index --skip-screen)")

    print("\n  \033[1mWhere the risk concentrates\033[0m")
    hot = max(rows, key=lambda r: r["branches"])
    print(f"    most branch points : {hot['file']} ({hot['branches']})")
    print(f"    longest file       : {max(rows, key=lambda r: r['code'])['file']}")
    print( "    least testable     : anything requiring hg38 + an aligner —")
    print( "                         which is why 04_fixture.sh enters via Direct Mode")
    print()


if __name__ == "__main__":
    sys.exit(main())
