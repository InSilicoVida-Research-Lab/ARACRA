#!/usr/bin/env bash
# =============================================================================
#  01_static.sh — syntax and cross-file consistency. No install, no network.
#  Runs in about a second. This is the one to put in CI.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${HERE}/.." && pwd)"
APP="${ROOT}/ARACRA"
source "${HERE}/lib/assert.sh"

section "Shell syntax"
for f in "${APP}/setup.sh" "${APP}/run_app.sh" "${APP}/lib/aracra_common.sh"; do
    assert_ok "bash -n $(basename "$f")" bash -n "$f"
done

section "Python / R / Groovy parse"
assert_ok "python parse aracra_star_app.py" \
    python3 -c "import ast,sys;ast.parse(open(sys.argv[1]).read())" "${APP}/aracra_star_app.py"
for pyf in "${APP}"/scripts/*.py; do
    assert_ok "python parse $(basename "$pyf")" \
        python3 -c "import ast,sys;ast.parse(open(sys.argv[1]).read())" "$pyf"
done

if need_tool Rscript; then
    for rf in "${APP}"/scripts/*.R; do
        assert_ok "R parse $(basename "$rf")" \
            Rscript -e 'invisible(parse(commandArgs(TRUE)[1]))' "$rf"
    done
else
    skip "R parse of scripts/*.R (Rscript not installed)"
fi

# Brace/quote balance is a weak proxy for a Groovy parser, but it catches the
# class of damage that was actually in nextflow.config (a stray extra quote).
section "Config structure"
python3 - "$APP" <<'PY'
import sys, pathlib
app = pathlib.Path(sys.argv[1])
ok = True
for name in ["nextflow.config", "main.nf"]:
    s = (app / name).read_text()
    for o, c in [("{", "}"), ("(", ")"), ("[", "]")]:
        if s.count(o) != s.count(c):
            print(f"  FAIL  {name}: unbalanced {o}{c}"); ok = False
    if s.count('"') % 2:
        print(f"  FAIL  {name}: odd number of double quotes"); ok = False
sys.exit(0 if ok else 1)
PY
[ $? -eq 0 ] && pass "nextflow.config / main.nf balanced" || fail "nextflow.config / main.nf balanced"

section "Regression guards (defects previously found in main)"
assert_ok      "no stray shell-redirect file (=*)" \
    bash -c "! compgen -G '${APP}/=*' > /dev/null"
assert_absent  "nextflow.config: no doubled opening quote" "${APP}/nextflow.config" '=  ""'
assert_absent  "nextflow.config: no literal ~ path"        "${APP}/nextflow.config" '"~/'
assert_contains "setup.sh: salmon spec is quoted"          "${APP}/setup.sh"        '"salmon${ARACRA_SALMON_SPEC}"'
assert_absent  "setup.sh: rtracklayer version not hardcoded" "${APP}/setup.sh"      'rtracklayer_1.66.0'
assert_contains "nextflow.config: version guard present"   "${APP}/nextflow.config" 'nextflowVersion'

section "Cross-file agreement"

# The app reads keys from .env; the shared writer must emit every one of them.
python3 - "$APP" <<'PY'
import sys, re, pathlib
app = pathlib.Path(sys.argv[1])
needed = set(re.findall(r'_env\("([A-Z_0-9]+)"', (app/"aracra_star_app.py").read_text()))
written = set(re.findall(r'^([A-Z_0-9]+)=', (app/"lib"/"aracra_common.sh").read_text(), re.M))
missing = sorted(needed - written)
print("  app reads %d key(s); writer emits all of them" % len(needed) if not missing
      else "  MISSING from .env writer: %s" % ", ".join(missing))
sys.exit(1 if missing else 0)
PY
[ $? -eq 0 ] && pass ".env key parity (app vs writer)" || fail ".env key parity (app vs writer)"

# Port must agree between the launcher, the shared default and the manual.
PORT_LIB=$(grep -oP 'ARACRA_DEFAULT_PORT:-\K[0-9]+' "${APP}/lib/aracra_common.sh" | head -1)
if [ -n "$PORT_LIB" ] && grep -q "localhost:${PORT_LIB}" "${ROOT}/README.md" \
   && ! grep -qE '\b850[0-9]\b' "${APP}/run_app.sh"; then
    pass "port ${PORT_LIB} consistent (lib / run_app.sh / README)"
else
    fail "port consistent (lib / run_app.sh / README)"
fi

# Every --option run_dromics.R defines should be reachable from main.nf, or it
# is dead config. (Informational: reports, does not fail.)
python3 - "$APP" <<'PY'
import sys, re, pathlib
app = pathlib.Path(sys.argv[1])
opts = set(re.findall(r'make_option\("--([a-z_0-9]+)"', (app/"scripts"/"run_dromics.R").read_text()))
nf   = (app/"main.nf").read_text()
unused = sorted(o for o in opts if f"--{o}" not in nf)
print(f"  run_dromics.R defines {len(opts)} options; {len(opts)-len(unused)} passed by main.nf")
if unused:
    print("  not passed by main.nf (CLI-only): " + ", ".join(unused))
PY
pass "DRomics option surface reported"

section "Documented defaults match code defaults"
python3 - "$ROOT" <<'PY'
import sys, re, pathlib
root = pathlib.Path(sys.argv[1])
r  = (root/"ARACRA"/"scripts"/"run_dromics.R").read_text()
app = (root/"ARACRA"/"aracra_star_app.py").read_text()

def r_default(opt):
    m = re.search(r'make_option\("--%s".*?default\s*=\s*([^,\)]+)' % opt, r, re.S)
    return m.group(1).strip() if m else None

bad = []
# niter: was 20000L in the script but 1000 in the app/README.
if r_default("niter") not in ("1000L", "1000"):
    bad.append("niter: run_dromics.R=%s, app=1000" % r_default("niter"))
if '"dr_niter":         1000' not in app and '"dr_niter": 1000' not in app:
    bad.append("app dr_niter default is not 1000")
for b in bad: print("  FAIL  " + b)
sys.exit(1 if bad else 0)
PY
[ $? -eq 0 ] && pass "bootstrap iteration default agrees (script/app)" \
             || fail "bootstrap iteration default agrees (script/app)"

summary
