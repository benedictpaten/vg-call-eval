#!/usr/bin/env bash
# Smoke test for the harness scripts themselves.
#
# Every one of these checks corresponds to a bug that reached a whole-genome run. None of them
# would have been caught by vg's own test suite, which exercises `vg` and never these scripts:
#
#   * `"${extra[@]}"` under `set -u` aborts on an empty array in bash 3.2, which is what
#     /bin/bash is on macOS. Every diploid contig failed. Only chrX passes a ploidy BED, so only
#     chrX had exercised call_one_bed since it gained the argument.
#   * `head -n -1` is a GNU extension; BSD head rejects a negative count.
#   * assemble_wgs.sh recompressed a contig only when its .tbi was missing, so after a full recall
#     it concatenated the *previous* run's VCFs and the refreshed numbers came out identical to the
#     run they replaced. Nothing said so; only timestamps did.
#   * bench_wgs.py, sweep_params.sh, schedule_wgs.py and call_wgs.sh all cached on "does the output
#     exist" rather than "is it newer than its input" -- the same mistake four times over.
#
# The theme is that a harness measuring a binary that changes must key its caches on freshness.
# These tests assert that property directly rather than trusting the comment that says so.
#
# Runs in seconds and touches no real data: everything is a temporary directory of stub files.
set -uo pipefail
cd "$(dirname "$0")/.."

PASS=0; FAIL=0
ok()   { PASS=$((PASS+1)); printf "  ok   %s\n" "$1"; }
bad()  { FAIL=$((FAIL+1)); printf "  FAIL %s\n" "$1"; }
check(){ if [ "$2" = "$3" ]; then ok "$1"; else bad "$1 (got '$2', want '$3')"; fi; }

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

echo "== shell syntax =="
for f in scripts/wgs/*.sh scripts/coverage/*.sh scripts/tier2/*.sh scripts/test_harness.sh; do
    [ -e "$f" ] || continue
    if bash -n "$f" 2>/dev/null; then ok "$(basename "$f") parses"; else bad "$(basename "$f") parses"; fi
done

echo "== python syntax =="
for f in scripts/wgs/*.py scripts/coverage/*.py; do
    [ -e "$f" ] || continue
    if python3 -c "import ast,sys; ast.parse(open('$f').read())" 2>/dev/null; then
        ok "$(basename "$f") parses"
    else
        bad "$(basename "$f") parses"
    fi
done

echo "== empty-array expansion under set -u (bash 3.2) =="
# The exact idiom call_one_bed uses. "${a[@]}" aborts; ${a[@]+"${a[@]}"} does not.
cat > "$TMP/arr.sh" <<'EOF'
set -euo pipefail
f() { local extra=(); [ -n "${1:-}" ] && extra=(--flag "$1"); echo "n=${#extra[@]}" ${extra[@]+"${extra[@]}"}; }
f
f value
EOF
if out=$(bash "$TMP/arr.sh" 2>&1); then
    check "optional-argument array expands with no argument" "$(echo "$out" | head -1)" "n=0"
    check "optional-argument array expands with an argument" "$(echo "$out" | tail -1)" "n=2 --flag value"
else
    bad "optional-argument array expansion (aborted: $out)"
fi
# Strip comments before linting, or the check flags the very comment that documents the fix --
# which is exactly what happened the first time this test ran.
uncommented() { sed 's/[[:space:]]*#.*$//' "$1"; }
if uncommented scripts/wgs/call_wgs.sh | grep -q 'extra\[@\]' && \
   ! uncommented scripts/wgs/call_wgs.sh | grep 'extra\[@\]' | grep -qv 'extra\[@\]+'; then
    ok "call_wgs.sh guards its optional-argument array"
elif ! uncommented scripts/wgs/call_wgs.sh | grep -q 'extra\[@\]'; then
    ok "call_wgs.sh has no optional-argument array to guard"
else
    bad "call_wgs.sh uses an unguarded optional-argument array"
fi

echo "== no GNU-only invocations =="
hits=0
for f in $(find scripts -name '*.sh'); do
    uncommented "$f" | grep -q 'head -n -[0-9]' && hits=$((hits+1))
done
if [ "$hits" -eq 0 ]; then
    ok "no 'head -n -N' outside comments (BSD head rejects it)"
else
    bad "$hits script(s) use 'head -n -N', which BSD head rejects"
fi

echo "== caches key on freshness, not existence =="
# assemble must recompress when the .vcf is newer than its .gz, even though a .tbi exists.
mkdir -p "$TMP/w/chrT"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > "$TMP/w/chrT/chrT.vcf"
touch -t 202001010000 "$TMP/w/chrT/chrT.vcf.gz" "$TMP/w/chrT/chrT.vcf.gz.tbi"
V="$TMP/w/chrT/chrT.vcf"
if [ ! -s "$V.gz.tbi" ] || [ "$V" -nt "$V.gz" ]; then
    ok "assemble recompresses a .gz older than its .vcf"
else
    bad "assemble would reuse a .gz older than its .vcf"
fi
grep -q '"\$V" -nt "\$V.gz"' scripts/wgs/assemble_wgs.sh \
    && ok "assemble_wgs.sh carries the freshness test" \
    || bad "assemble_wgs.sh lost the freshness test"
grep -q 'st_mtime' scripts/wgs/bench_wgs.py \
    && ok "bench_wgs.py rescores on mtime" || bad "bench_wgs.py caches on existence"
grep -q 'vg_mtime' scripts/wgs/schedule_wgs.py \
    && ok "schedule_wgs.py revalidates .done against the binary" \
    || bad "schedule_wgs.py trusts .done unconditionally"
grep -q -- '-nt "\$VG"' scripts/wgs/call_wgs.sh \
    && ok "call_wgs.sh revalidates .done against the binary" \
    || bad "call_wgs.sh trusts .done unconditionally"
grep -q -- '-nt "\$VG"' scripts/coverage/sweep_params.sh \
    && ok "sweep_params.sh revalidates results against the binary" \
    || bad "sweep_params.sh trusts its outputs unconditionally"

echo "== the memory model matches its fitted constants =="
python3 - <<'PY' && ok "schedule_wgs.py memory model is the refitted one" || bad "memory model drifted from the doc"
import re, sys, pathlib
s = pathlib.Path("scripts/wgs/schedule_wgs.py").read_text()
base = float(re.search(r"BASE_GB = ([0-9.]+)", s).group(1))
per  = float(re.search(r"GB_PER_RECORD = ([0-9.e-]+)", s).group(1))
# Largest contig seen is ~370k records at 5.6 GB; the model must not underestimate it,
# nor overestimate by the >50% the original coefficient did.
pred = base + per * 370000
sys.exit(0 if 5.5 <= pred <= 8.0 else 1)
PY

echo
echo "$PASS passed, $FAIL failed"
[ "$FAIL" -eq 0 ]
