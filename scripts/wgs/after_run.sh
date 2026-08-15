#!/usr/bin/env bash
# Everything after the serial calling run, in the order that protects the deliverable.
#
# Assemble and benchmark first. Those produce the whole-genome VCF and mosaic, which are the point
# of the exercise; the scheduler experiment re-runs contigs and overwrites their per-contig outputs,
# so it must not go first. The assembled files are independent of that churn once written.
set -uo pipefail
cd "$(dirname "$0")/../.."
export PATH="$HOME/CLionProjects/vg/bin:$PATH"
PY=/opt/homebrew/bin/python3.11
step() { echo; echo "=============== $* ==============="; date +%H:%M:%S; }

# --- 0. Fix the resume marker -----------------------------------------------
# `[ -s ]` requires a non-empty file and `touch` makes an empty one, so the skip never fired and a
# resumed run silently recalled every contig. Write the record count instead: non-empty, and it
# says what was done.
step "fix resume marker"
$PY - <<'PY'
import pathlib
p = pathlib.Path("scripts/wgs/call_wgs.sh"); s = p.read_text()
s = s.replace('    if [ -s "$D/$C.done" ]; then',
              '    if [ -f "$D/$C.done" ]; then')
s = s.replace('    touch "$D/$C.done"',
              '    grep -vc "^#" "$D/$C.vcf" > "$D/$C.done"')
p.write_text(s); print("patched call_wgs.sh")
PY
for d in work/wgs/*/; do
    c=$(basename "$d")
    [ -f "$d/$c.done" ] && [ -s "$d/$c.vcf" ] && grep -vc "^#" "$d/$c.vcf" > "$d/$c.done"
done

# --- 1. Assemble -------------------------------------------------------------
step "assemble"
bash scripts/wgs/assemble_wgs.sh 2>&1 | tail -8

# --- 2. Benchmark ------------------------------------------------------------
step "benchmark"
$PY scripts/wgs/bench_wgs.py --work work/wgs --out docs/wgs-results.md 2>&1 | tail -6

step "phasing"
$PY scripts/tier2/phasing_benchmark.py \
    --calls work/wgs/HG002.vcf.gz \
    --truth data/truth/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.vcf.gz \
    --out work/wgs/score/phasing-wgs 2>&1 | tail -14 | tee work/wgs/score/phasing-wgs.txt

# --- 3. Thread A/B -----------------------------------------------------------
# Machine otherwise idle by here: benchmarking is done, calling is done.
step "thread A/B"
C=chr20 REPS=2 THREAD_SET="1 2 5" OUT=/tmp/thread_ab.tsv \
    bash scripts/wgs/thread_ab.sh 2>&1 | tail -12

# --- 4. Scheduler experiment -------------------------------------------------
# Same contigs, serial times already measured in this run, re-run packed under a memory budget.
# Six contigs spanning the size range, so the comparison is not just about the small easy ones.
step "scheduler experiment"
SUBSET="chr3 chr8 chr14 chr18 chr20 chr21"
$PY - "$SUBSET" <<'PY'
import re, sys, pathlib, json
# Serial baseline: each contig's own time from its call log.
out = {}
for c in sys.argv[1].split():
    log = pathlib.Path(f"work/wgs/{c}/{c}.log")
    if not log.exists():
        continue
    m = re.search(r"^ *([0-9.]+) real", log.read_text(), re.M)
    if m:
        out[c] = float(m.group(1))
pathlib.Path("/tmp/serial_times.json").write_text(json.dumps(out, indent=2))
print("serial baseline (s):", {k: round(v) for k, v in out.items()},
      "total", round(sum(out.values())))
PY

for c in $SUBSET; do rm -f "work/wgs/$c/$c.done"; done
step "scheduler dry run"
$PY scripts/wgs/schedule_wgs.py --contigs $SUBSET --dry-run --threads 2 2>&1 | tail -14

step "scheduler timed run"
SCHED_START=$(date +%s)
$PY scripts/wgs/schedule_wgs.py --contigs $SUBSET --threads 2 --budget-gb 24 --max-jobs 6 2>&1 | tail -20
SCHED_END=$(date +%s)
echo "SCHEDULED_WALL_SECONDS=$((SCHED_END - SCHED_START))" | tee /tmp/sched_wall.txt

step "speed comparison"
$PY - <<'PY'
import json, pathlib
serial = json.load(open("/tmp/serial_times.json"))
wall = int(open("/tmp/sched_wall.txt").read().split("=")[1])
tot = sum(serial.values())
print(f"serial, summed:    {tot:8.0f} s  ({tot/60:.1f} min)")
print(f"scheduled, wall:   {wall:8d} s  ({wall/60:.1f} min)")
print(f"speedup:           {tot/wall:8.2f}x")
PY
echo "AFTER_RUN_DONE"
