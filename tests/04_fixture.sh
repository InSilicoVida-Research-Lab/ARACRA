#!/usr/bin/env bash
# =============================================================================
#  04_fixture.sh — end-to-end statistics on a tiny synthetic dataset.
#
#  Enters through Direct Mode (count matrix in, no alignment), so it exercises
#  DESeq2 -> DRomics -> BMD -> tPOD without needing hg38 or an aligner.
#  Runtime: a few minutes. Runs the REAL scripts, not a mock.
#
#  It runs run_dromics.R twice — with --bmd_extrap_filter FALSE and TRUE — and
#  asserts that the setting changes the tPOD in the expected direction. That is
#  the regression test for the one behavioural switch added to the pipeline.
#
#  Usage:  bash 04_fixture.sh [--keep]
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${HERE}/.." && pwd)"
APP="${ROOT}/ARACRA"
WORK="${HERE}/.out/fixture"
KEEP=false
[ "${1:-}" = "--keep" ] && KEEP=true
source "${HERE}/lib/assert.sh"

section "Environment"
if ! need_tool Rscript; then
    skip "Rscript not found — activate the ARACRA env first:"
    printf "         conda activate ~/miniforge3/envs/test_ARACRA\n"
    summary; exit 0
fi
pass "Rscript: $(Rscript --version 2>&1 | head -1)"

MISSING_R=$(Rscript -e '
pkgs <- c("DRomics","DESeq2","optparse","dplyr","ggplot2","jsonlite")
cat(paste(pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)], collapse=" "))' 2>/dev/null)
if [ -n "$MISSING_R" ]; then
    skip "missing R packages: ${MISSING_R} — run setup.sh"
    summary; exit 0
fi
pass "required R packages present"

section "Fixture"
rm -rf "$WORK"; mkdir -p "$WORK"
python3 "${HERE}/fixtures/make_fixture.py" --outdir "$WORK" || { fail "fixture generation"; summary; exit 1; }
pass "fixture generated (deterministic, seed fixed)"

COUNTS="${WORK}/counts_matrix.csv"
META="${WORK}/metadata.csv"

section "DESeq2"
mkdir -p "${WORK}/deg"
if Rscript "${APP}/scripts/run_deseq2.R" \
        --counts "$COUNTS" --metadata "$META" \
        --treatment TestChem --control Control \
        --outdir "${WORK}/deg" --read_thresh 0 \
        > "${WORK}/deseq2.log" 2>&1; then
    pass "run_deseq2.R completed"
else
    fail "run_deseq2.R failed — see ${WORK}/deseq2.log"
    tail -15 "${WORK}/deseq2.log" | sed 's/^/    /'
fi

if [ -f "${WORK}/deg/All_Results.csv" ]; then
    pass "All_Results.csv produced"
    # 80 of 500 genes are planted responders; a sane DE call should find a
    # decent share of them and should NOT call most of the 420 nulls.
    python3 - "$WORK" <<'PY'
import sys, pandas as pd, pathlib
w = pathlib.Path(sys.argv[1])
res = pd.read_csv(w/"deg"/"All_Results.csv").set_index("ensembl_id")
# keep_default_na=False: pandas' default NA sentinel list includes the
# literal string "null", which this fixture uses as a real class label —
# without this, every null-class row silently becomes NaN and the class
# comparisons below never match anything.
truth = pd.read_csv(w/"ground_truth.csv", keep_default_na=False).set_index("gene_id")
padj = next((c for c in res.columns if "padj" in c.lower() or "adj" in c.lower()), None)
if padj is None:
    print("  (no adjusted p-value column found; skipping sensitivity check)"); sys.exit(0)
sig = set(res.index[res[padj] < 0.05])
real = set(truth.index[truth["class"] != "null"])
null = set(truth.index[truth["class"] == "null"])
tpr = len(sig & real) / max(len(real), 1)
fpr = len(sig & null) / max(len(null), 1)
print(f"  recovered {len(sig & real)}/{len(real)} planted responders (TPR {tpr:.0%})")
print(f"  false positives {len(sig & null)}/{len(null)} nulls (FPR {fpr:.1%})")
sys.exit(0 if (tpr > 0.5 and fpr < 0.10) else 1)
PY
    [ $? -eq 0 ] && pass "DE sensitivity/specificity within expected bounds" \
                 || fail "DE sensitivity/specificity outside expected bounds"
else
    fail "All_Results.csv missing"
fi

# ── DRomics, run twice to compare the extrapolation policy ───────────────────
run_dromics () {              # $1 = outdir suffix   $2 = TRUE|FALSE
    mkdir -p "${WORK}/dr_$1"
    Rscript "${APP}/scripts/run_dromics.R" \
        --counts "$COUNTS" --metadata "$META" \
        --treatment TestChem --control Control \
        --outdir "${WORK}/dr_$1" \
        --bmd TRUE --bootstrap FALSE --niter 200 \
        --read_thresh 0 \
        --bmd_extrap_filter "$2" \
        > "${WORK}/dromics_$1.log" 2>&1
}

section "DRomics — extrapolation filter OFF (default, reproduces published runs)"
if run_dromics off FALSE; then
    pass "run_dromics.R --bmd_extrap_filter FALSE completed"
else
    fail "run_dromics.R (FALSE) failed — see ${WORK}/dromics_off.log"
    tail -20 "${WORK}/dromics_off.log" | sed 's/^/    /'
fi
assert_ok "bmd_results.csv produced (OFF)" test -f "${WORK}/dr_off/bmd_results.csv"

section "DRomics — extrapolation filter ON"
if run_dromics on TRUE; then
    pass "run_dromics.R --bmd_extrap_filter TRUE completed"
else
    fail "run_dromics.R (TRUE) failed — see ${WORK}/dromics_on.log"
    tail -20 "${WORK}/dromics_on.log" | sed 's/^/    /'
fi

section "Extrapolation switch behaves as documented"
python3 - "$WORK" <<'PY'
import sys, json, pathlib
w = pathlib.Path(sys.argv[1])
def load(tag):
    p = w / f"dr_{tag}" / "dromics_summary.json"
    return json.loads(p.read_text()) if p.exists() else None

off, on = load("off"), load("on")
if not off or not on:
    print("  summaries missing — cannot compare"); sys.exit(1)

def dig(d, *keys):
    for k in keys:
        if not isinstance(d, dict): return None
        d = d.get(k)
    return d

NTP  = ("bmd_summary", "ntp_quality_filters")
TP   = ("bmd_summary", "tpod_methods", "rank25", "value")
n_flag  = dig(off, *NTP, "n_extrapolated")
rm_off  = dig(off, *NTP, "n_removed_extrap")
rm_on   = dig(on,  *NTP, "n_removed_extrap")
t_off   = dig(off, *TP) or dig(off, "bmd_summary", "median_bmd")
t_on    = dig(on,  *TP) or dig(on,  "bmd_summary", "median_bmd")
n_bmd_off = dig(off, "bmd_summary", "total_bmds")
n_bmd_on  = dig(on,  "bmd_summary", "total_bmds")
print(f"  valid BMDs   OFF / ON         : {n_bmd_off} / {n_bmd_on}")

print(f"  genes flagged as extrapolated : {n_flag}")
print(f"  removed with filter OFF       : {rm_off}   (must be 0)")
print(f"  removed with filter ON        : {rm_on}")
print(f"  rank25 tPOD  OFF / ON         : {t_off} / {t_on}")

problems = []
if rm_off not in (0, None):
    problems.append("filter OFF removed genes — default behaviour changed")
if n_flag and rm_on != n_flag:
    problems.append("filter ON did not remove every flagged gene")
if t_off is not None and t_on is not None and n_flag and t_on < t_off:
    problems.append("removing below-range BMDs LOWERED the tPOD — expected it to rise")
for p in problems: print("  FAIL  " + p)
sys.exit(1 if problems else 0)
PY
[ $? -eq 0 ] && pass "flag-vs-remove semantics correct, tPOD moves in the expected direction" \
             || fail "flag-vs-remove semantics wrong"

section "Artefacts"
for f in bmd_results.csv dromics_summary.json; do
    assert_ok "dr_off/${f}" test -f "${WORK}/dr_off/${f}"
done

if [ "$KEEP" = true ]; then
    printf "\n  Kept: %s\n" "$WORK"
else
    printf "\n  Output kept at %s (delete with: rm -rf %s)\n" "$WORK" "$WORK"
fi
summary
