#!/usr/bin/env bash
# =============================================================================
#  02_preflight.sh — can THIS machine run ARACRA, and what will setup cost?
#
#  Run this BEFORE setup.sh. It downloads nothing and installs nothing. The
#  point is to find out in 5 seconds what would otherwise surface 90 minutes
#  into a failed install.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${HERE}/.." && pwd)"
APP="${ROOT}/ARACRA"
source "${HERE}/lib/assert.sh"
source "${APP}/lib/aracra_common.sh"

CORES="$(aracra_cpu_cores)"
RAM="$(aracra_ram_gb_int)"
ALIGNER="$(aracra_recommended_aligner "$RAM")"
DISK_HOME=$(df -BG "${HOME}" 2>/dev/null | awk 'NR==2{gsub("G","");print $4}' || echo 0)

section "Machine"
printf "  CPU cores        : %s\n" "$CORES"
printf "  RAM              : %s GB\n" "$RAM"
printf "  Free disk (HOME) : %s GB\n" "$DISK_HOME"
printf "  Architecture     : %s\n" "$(uname -m)"
printf "  Aligner chosen   : %s\n" "$ALIGNER"

section "Requirements"
[ "$CORES" -ge 4 ]  && pass "CPU cores >= 4"            || fail "CPU cores >= 4 (have ${CORES})"
[ "$RAM"   -ge 16 ] && pass "RAM >= 16 GB (minimum)"    || fail "RAM >= 16 GB (have ${RAM})"
if [ "$RAM" -ge 32 ]; then
    pass "RAM >= 32 GB — STAR usable"
else
    skip "RAM < 32 GB — STAR unavailable, HISAT2 will be used (not a failure)"
fi
[ "${DISK_HOME:-0}" -ge 100 ] && pass "Free disk >= 100 GB" \
    || fail "Free disk >= 100 GB (have ${DISK_HOME}) — use --db-dir for a larger drive"

case "$(uname -m)" in
    x86_64)          pass "Architecture x86_64 (main branch)" ;;
    aarch64|arm64)   skip "Architecture arm64 — use dgx_install_branch, not main" ;;
    *)               fail "Unsupported architecture: $(uname -m)" ;;
esac

section "Prerequisites"
for t in curl wget git tar; do
    need_tool "$t" && pass "$t present" || fail "$t missing (apt install $t)"
done
if need_tool conda || need_tool mamba; then
    pass "conda/mamba present — setup will reuse it"
else
    skip "no conda found — setup.sh will install Miniforge3 (adds ~10 min, ~500 MB)"
fi
need_tool docker && skip "docker present (only needed for the container route)" || true

section "Install cost estimate"
# Figures are the sizes setup.sh actually pulls, not guesses at runtime.
EST_INDEX_MIN=0
[ "$ALIGNER" = "star" ] && EST_INDEX_MIN=40   # STAR index build; HISAT2 is prebuilt
printf "  conda environment       ~  6 GB    ~20-40 min\n"
printf "  hg38 genome + GTF       ~ 15 GB    ~15-30 min (network bound)\n"
printf "  STAR index build        ~ 28 GB    ~%s min\n" "$EST_INDEX_MIN"
printf "  FastQ Screen genomes    ~ 14 GB    ~15 min   (skip: --skip-screen)\n"
printf "  %s\n" "----------------------------------------------"
TOTAL_GB=$(( 6 + 15 + 14 + (EST_INDEX_MIN > 0 ? 28 : 3) ))
printf "  total                   ~ %s GB\n" "$TOTAL_GB"
if [ "${DISK_HOME:-0}" -lt "$TOTAL_GB" ]; then
    fail "not enough free disk for a full setup (${DISK_HOME} GB free, ~${TOTAL_GB} GB needed)"
else
    pass "disk sufficient for a full setup (~${TOTAL_GB} GB)"
fi
printf "\n  Leaner options:\n"
printf "    bash setup.sh --skip-screen                  saves ~14 GB\n"
printf "    bash setup.sh --skip-index                   saves ~28 GB and ~40 min\n"
printf "    bash setup.sh --db-dir=/mnt/bigdrive/dbs     puts references elsewhere\n"

section "Network reachability"
for host in https://ftp.ebi.ac.uk https://bioconductor.org https://conda.anaconda.org; do
    if curl -sSf -I --max-time 8 "$host" >/dev/null 2>&1; then
        pass "reachable: $host"
    else
        skip "unreachable here: $host (proxy/firewall, or offline sandbox)"
    fi
done

summary
