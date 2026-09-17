#!/usr/bin/env bash
# =============================================================================
#  03_solve.sh — does the version-constrained solve actually resolve?
#
#  This is the check that matters most, because it is the one thing about
#  setup.sh that cannot be verified by reading it. It asks mamba to SOLVE the
#  pinned spec without downloading or installing anything (--dry-run), which
#  takes a couple of minutes instead of an hour.
#
#  It also solves the UNPINNED spec, so if the pinned one fails you immediately
#  know whether the floors are the cause or the channels are simply broken.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${HERE}/.." && pwd)"
APP="${ROOT}/ARACRA"
OUTDIR="${HERE}/.out"
mkdir -p "$OUTDIR"
source "${HERE}/lib/assert.sh"
source "${APP}/lib/aracra_common.sh"

SOLVER=""
need_tool mamba && SOLVER=mamba
[ -z "$SOLVER" ] && need_tool conda && SOLVER=conda

section "Solver"
if [ -z "$SOLVER" ]; then
    skip "no conda/mamba on PATH — cannot test the solve here"
    skip "pinned solve"
    skip "unpinned solve"
    summary
    exit 0
fi
pass "using ${SOLVER} ($(${SOLVER} --version 2>&1 | head -1))"

TMP_PREFIX="${OUTDIR}/dryrun_env"

section "Pinned solve (the spec setup.sh tries first)"
printf "  spec: nextflow%s salmon%s %s\n\n" \
    "$ARACRA_NEXTFLOW_SPEC" "$ARACRA_SALMON_SPEC" "${ARACRA_TOOL_SPECS[*]}"

if "$SOLVER" create --dry-run -p "$TMP_PREFIX" \
        -c conda-forge -c bioconda -y \
        "nextflow${ARACRA_NEXTFLOW_SPEC}" \
        "salmon${ARACRA_SALMON_SPEC}" \
        "${ARACRA_TOOL_SPECS[@]}" \
        python=3.12 openpyxl fastq-screen fastqc picard qualimap pigz \
        > "${OUTDIR}/solve_pinned.log" 2>&1; then
    pass "pinned spec resolves"
    # Report what the floors would actually give you.
    printf "\n  Resolved versions:\n"
    grep -oE '\b(nextflow|salmon|star|hisat2|samtools|fastp|subread|multiqc|rseqc)-[0-9][^ -]*' \
        "${OUTDIR}/solve_pinned.log" | sort -u | sed 's/^/    /' || true
else
    fail "pinned spec does NOT resolve — setup.sh will fall back to unpinned"
    printf "\n  Solver said:\n"
    grep -iE 'unsatisfiable|conflict|nothing provides|cannot' \
        "${OUTDIR}/solve_pinned.log" | head -12 | sed 's/^/    /' || \
        tail -12 "${OUTDIR}/solve_pinned.log" | sed 's/^/    /'
    printf "\n  Raise or drop the offending floor in ARACRA/lib/aracra_common.sh\n"
fi

section "Unpinned solve (the fallback path)"
if "$SOLVER" create --dry-run -p "$TMP_PREFIX" \
        -c conda-forge -c bioconda -y \
        openjdk=17 nextflow python=3.12 openpyxl sra-tools fastp fastq-screen \
        fastqc star hisat2 bowtie2 samtools rseqc picard qualimap \
        "salmon${ARACRA_SALMON_SPEC}" subread multiqc pigz \
        > "${OUTDIR}/solve_unpinned.log" 2>&1; then
    pass "fallback spec resolves — setup.sh can always complete"
else
    fail "fallback spec does NOT resolve either — channel or base-env problem, not the floors"
    tail -12 "${OUTDIR}/solve_unpinned.log" | sed 's/^/    /'
fi

section "Nextflow version guard"
if need_tool nextflow; then
    NF_VER=$(nextflow -version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
    NF_MAJOR=${NF_VER%%.*}
    if [ "$NF_MAJOR" = "24" ]; then
        pass "installed Nextflow ${NF_VER} is inside the declared range"
    else
        fail "installed Nextflow ${NF_VER} is OUTSIDE '>=24.04.0, <25.0.0' — nextflow.config will refuse to run"
    fi
else
    skip "nextflow not on PATH (expected before setup.sh has run)"
fi

printf "\n  Full solver output: %s\n" "${OUTDIR}/solve_*.log"
summary
