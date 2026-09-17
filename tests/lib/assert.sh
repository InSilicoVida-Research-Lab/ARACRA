#!/usr/bin/env bash
# =============================================================================
#  tests/lib/assert.sh — tiny test helpers. No framework, no dependencies.
#
#  Three outcomes, and the distinction matters:
#    pass  — checked, correct
#    fail  — checked, wrong            → exit code 1
#    skip  — could NOT be checked here (tool absent, no network)
#
#  A skip is never silently counted as a pass. The point of this folder is to
#  be honest about what was actually verified on THIS machine.
# =============================================================================

T_PASS=0; T_FAIL=0; T_SKIP=0; T_FAILED_NAMES=()

if [ -t 1 ]; then
    C_G='\033[0;32m'; C_R='\033[0;31m'; C_Y='\033[1;33m'; C_B='\033[1m'; C_N='\033[0m'
else
    C_G=''; C_R=''; C_Y=''; C_B=''; C_N=''
fi

section() { printf "\n${C_B}── %s ${C_N}\n" "$*"; }

pass() { T_PASS=$((T_PASS+1)); printf "  ${C_G}pass${C_N}  %s\n" "$*"; }
fail() { T_FAIL=$((T_FAIL+1)); T_FAILED_NAMES+=("$1"); printf "  ${C_R}FAIL${C_N}  %s\n" "$*"; }
skip() { T_SKIP=$((T_SKIP+1)); printf "  ${C_Y}skip${C_N}  %s\n" "$*"; }

# assert_ok "name" <command...>
assert_ok() {
    local name="$1"; shift
    if "$@" >/dev/null 2>&1; then pass "$name"; else fail "$name"; fi
}

# assert_contains "name" <file> <string>
assert_contains() {
    local name="$1" file="$2" needle="$3"
    if [ ! -f "$file" ]; then fail "$name (missing file: $file)"; return; fi
    if grep -qF -- "$needle" "$file"; then pass "$name"; else fail "$name"; fi
}

# assert_absent "name" <file> <string>   — string must NOT appear
assert_absent() {
    local name="$1" file="$2" needle="$3"
    if [ ! -f "$file" ]; then fail "$name (missing file: $file)"; return; fi
    if grep -qF -- "$needle" "$file"; then fail "$name"; else pass "$name"; fi
}

# need_tool <tool> — returns 1 if absent, so callers can skip cleanly
need_tool() { command -v "$1" >/dev/null 2>&1; }

summary() {
    printf "\n${C_B}%s${C_N}\n" "────────────────────────────────────────────"
    printf "  ${C_G}%d passed${C_N}   ${C_R}%d failed${C_N}   ${C_Y}%d skipped${C_N}\n" \
        "$T_PASS" "$T_FAIL" "$T_SKIP"
    if [ "$T_FAIL" -gt 0 ]; then
        printf "\n  Failed:\n"
        for n in "${T_FAILED_NAMES[@]}"; do printf "    - %s\n" "$n"; done
    fi
    if [ "$T_SKIP" -gt 0 ]; then
        printf "\n  %d check(s) could not run here — they are NOT passes.\n" "$T_SKIP"
    fi
    printf "\n"
    [ "$T_FAIL" -eq 0 ]
}
