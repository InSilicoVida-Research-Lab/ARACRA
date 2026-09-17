#!/usr/bin/env bash
# =============================================================================
#  run_tests.sh — ARACRA test entrypoint
#
#    bash tests/run_tests.sh            static + preflight   (~5 s, no install)
#    bash tests/run_tests.sh --solve    + conda solve dry-run (~2-5 min)
#    bash tests/run_tests.sh --fixture  + end-to-end statistics (~5 min)
#    bash tests/run_tests.sh --pathway  + pathway-level tPOD, real gene IDs (~1-2 min)
#    bash tests/run_tests.sh --all      everything
#    bash tests/run_tests.sh --complexity   report only
#
#  The tiers are ordered by cost. Nothing here installs anything or touches
#  ~/databases; the fixture writes only inside tests/.out/.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${HERE}/.." && pwd)"

RUN_SOLVE=false
RUN_FIXTURE=false
RUN_PATHWAY=false
ONLY_COMPLEXITY=false
case "${1:-}" in
    --solve)      RUN_SOLVE=true ;;
    --fixture)    RUN_FIXTURE=true ;;
    --pathway)    RUN_PATHWAY=true ;;
    --all)        RUN_SOLVE=true; RUN_FIXTURE=true; RUN_PATHWAY=true ;;
    --complexity) ONLY_COMPLEXITY=true ;;
    "")           ;;
    *) echo "unknown option: $1"; sed -n '4,11p' "$0"; exit 2 ;;
esac

RC=0
banner() { printf "\n\033[1m═══ %s ═══\033[0m\n" "$*"; }

if [ "$ONLY_COMPLEXITY" = true ]; then
    python3 "${HERE}/complexity_report.py" --repo "$ROOT"
    exit 0
fi

banner "1. Static checks"
bash "${HERE}/01_static.sh" || RC=1

banner "2. Preflight (this machine)"
bash "${HERE}/02_preflight.sh" || RC=1

if [ "$RUN_SOLVE" = true ]; then
    banner "3. Dependency solve (dry-run)"
    bash "${HERE}/03_solve.sh" || RC=1
else
    printf "\n\033[1;33m  skipped: dependency solve\033[0m  (add --solve)\n"
fi

if [ "$RUN_FIXTURE" = true ]; then
    banner "4. End-to-end fixture"
    bash "${HERE}/04_fixture.sh" || RC=1
else
    printf "\033[1;33m  skipped: end-to-end fixture\033[0m  (add --fixture)\n"
fi

if [ "$RUN_PATHWAY" = true ]; then
    banner "5. Pathway-level tPOD (real gene IDs)"
    bash "${HERE}/05_pathway_fixture.sh" || RC=1
else
    printf "\033[1;33m  skipped: pathway-level tPOD\033[0m  (add --pathway)\n"
fi

banner "Complexity"
python3 "${HERE}/complexity_report.py" --repo "$ROOT"

if [ "$RC" -eq 0 ]; then
    printf "\033[0;32m  All executed checks passed.\033[0m\n"
else
    printf "\033[0;31m  Some checks failed — see above.\033[0m\n"
fi
printf "  A 'skip' is not a pass; it means the check could not run here.\n\n"
exit $RC
