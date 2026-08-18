#!/usr/bin/env python3
"""Run the per-contig calls concurrently, packed under a memory budget.

Why this exists. A single contig at `-t 5` uses about 3.5 of this machine's 10 cores, so a serial
run leaves two thirds of it idle. Running two or three contigs at once converts that into
throughput: on a six-contig subset, 1451 s serial becomes 797 s scheduled, a 1.82x speedup.

Not because the caller is I/O-bound. That was the first theory -- every read-fetch window
`posix_spawn`s a gbz-base and reopens the 22 GB SQLite database, which looked like enough to
explain a process sitting at one CPU -- and it is wrong. Measured on chr20 over two replicates,
`-t` 1/2/5 gives 0.99/1.79/3.48 CPU and 422/247/142 s, so it parallelises at about 70% efficiency,
and the per-thread read caches that were supposed to make threads expensive account for 0.4 GB
between `-t` 1 and 5. Scheduling more, thinner jobs on that theory produced 1126 s against 797 s.
The lesson worth keeping is that the spawn-per-window is real and still looks like the bottleneck
in the source; it just is not the one that governs wall clock here.

Memory is a constraint, so contigs are packed rather than run at a fixed concurrency. The model is
refitted from a full 24-contig run:

    peak GB ~ 2.25 + 11.2e-6 * emitted_records

    chr20 105,251 records  3.1 GB      chr8  238,309  5.6 GB
    chr6  284,529          5.0 GB      chr2  369,207  5.6 GB
    chr1  353,741          5.7 GB

**The previous coefficient was nearly double this and made the scheduler throttle itself on a
fiction.** It predicted 9.6 GB for chr1 where the measured peak is 5.7, overestimating every
contig by 0.4-4.4 GB, so the budget refused packings that would have fit comfortably. Refitting
cut the worst residual from 4.39 GB to 0.87.

Record counts are not known before a contig runs, so the truth's record count for that contig is
used as the predictor -- it is available from prep and correlates with what the caller emits far
better than contig length does.

The budget defaults to 24 GB on a 32 GB machine. The remainder is not slack: the read database is
22 GB and the OS page cache is doing real work for us, so squeezing it would trade one bottleneck
for a worse one.

Scheduling is largest-first. A 9.6 GB contig placed last cannot share the budget with anything, so
it would run alone at the end while the machine idles.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent

BASE_GB = 2.25
GB_PER_RECORD = 11.2e-6

# --nested emits more records and holds more working state, and the record count alone does not
# predict it: chr20 under -a has 64% more records than under --nested at *less* peak memory, so what
# grows is the descent's own working set rather than the buffered output. Two paired contigs give
# deltas of +0.21 GB (chr21) and +0.79 GB (chr20) for similar record increases, and chr21's default
# is itself 0.5 GB above chr20's at a comparable record count -- so run-to-run noise is the same
# order as the effect. Rather than fit a coefficient to two noisy points, apply a margin: it costs
# only some packing density, where under-predicting costs a swap storm.
NESTED_MARGIN = 1.25


def predict_gb(truth_records: int, nested: bool = False) -> float:
    gb = BASE_GB + GB_PER_RECORD * truth_records
    return gb * NESTED_MARGIN if nested else gb


def truth_record_count(work: Path, contig: str) -> int:
    """Emitted-record predictor: how many small-variant truth records the contig carries."""
    vcf = work / contig / f"truth.{contig}.smvar.vcf.gz"
    if not vcf.exists():
        return 0
    p = subprocess.run(["bcftools", "index", "-n", str(vcf)], capture_output=True, text=True)
    if p.returncode == 0 and p.stdout.strip().isdigit():
        return int(p.stdout.strip())
    p = subprocess.run(f"bcftools view -H {vcf} | wc -l", shell=True,
                       capture_output=True, text=True)
    try:
        return int(p.stdout.strip())
    except ValueError:
        return 0


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--work", default="work/wgs")
    p.add_argument("--budget-gb", type=float, default=24.0)
    p.add_argument("--nested", action="store_true",
                   help="pass --nested to vg call and budget memory for it")
    # -t 5, not the 2 this was first written with. The reasoning behind 2 was that the caller is
    # I/O-bound at about one CPU, so threads only bought per-thread read caches and could be traded
    # for concurrency. Measured on chr20, two replicates, that is wrong on both halves:
    #
    #     -t 1   422.1 s   3.46 GB   0.99 CPU
    #     -t 2   247.4 s   3.57 GB   1.79 CPU
    #     -t 5   141.7 s   3.88 GB   3.48 CPU
    #
    # It parallelises at about 70% efficiency, and the memory supposedly freed is 0.4 GB. Scheduled
    # at -t 2 the six-contig subset took 1126 s against 797 s at -t 5 -- packing more jobs that each
    # run 1.7x slower is a loss, and it read as a 1.29x win only against a serial baseline.
    p.add_argument("--threads", type=int, default=5,
                   help="-t per contig. The caller uses about 3.5 CPUs at 5, so on a 10-core "
                        "machine three concurrent jobs saturate it.")
    # Three jobs at ~3.5 CPU each saturate 10 cores. More jobs only queues them against a busy
    # machine while holding their memory.
    p.add_argument("--max-jobs", type=int, default=3)
    p.add_argument("--contigs", nargs="*",
                   default=[f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])
    p.add_argument("--vg", default=str(Path.home() / "CLionProjects/vg/bin/vg"),
                   help="the binary whose mtime decides whether a .done marker is still valid")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    work = Path(args.work)
    # Resume past work done by *this* binary only. A marker saying merely "called" carries a
    # pre-fix result across a rebuild without saying so; the coverage sweep kept pre-fix chrX
    # arms exactly this way and they scored as though they were the fixed caller.
    vg_mtime = Path(args.vg).stat().st_mtime if Path(args.vg).exists() else 0.0
    plan = []
    stale = []
    for c in args.contigs:
        marker = work / c / f"{c}.done"
        if marker.exists():
            if marker.stat().st_mtime >= vg_mtime:
                continue
            stale.append(c)
            marker.unlink()
        n = truth_record_count(work, c)
        plan.append({"contig": c, "truth_records": n,
                     "predict_gb": round(predict_gb(n, args.nested), 2)})
    plan.sort(key=lambda x: -x["predict_gb"])

    if stale:
        print(f"recalling {len(stale)} contig(s) whose marker predates {args.vg}: "
              + " ".join(stale))
    if not plan:
        print("nothing to do: every contig was called by this binary")
        return

    print(f"{'contig':8s} {'truth recs':>11s} {'predicted GB':>13s}")
    for e in plan:
        print(f"{e['contig']:8s} {e['truth_records']:11,d} {e['predict_gb']:13.1f}")
    print(f"\nbudget {args.budget_gb} GB, -t {args.threads}, max {args.max_jobs} jobs")

    if args.dry_run:
        # Show the packing without running anything.
        waves, cur, used = [], [], 0.0
        for e in plan:
            if cur and (used + e["predict_gb"] > args.budget_gb or len(cur) >= args.max_jobs):
                waves.append((cur, used)); cur, used = [], 0.0
            cur.append(e["contig"]); used += e["predict_gb"]
        if cur:
            waves.append((cur, used))
        print()
        for i, (w, u) in enumerate(waves, 1):
            print(f"  wave {i}: {u:5.1f} GB  {' '.join(w)}")
        print(f"\n{len(waves)} waves for {len(plan)} contigs")
        return

    # Run them, keeping the in-flight predicted total under budget. Not a wave scheduler: a job is
    # started the moment its predicted peak fits, so a long contig does not hold back short ones.
    running = {}   # popen -> entry
    started = time.time()
    log = work / "schedule.log"
    results = []

    def in_flight_gb():
        return sum(e["predict_gb"] for e in running.values())

    queue = list(plan)
    with open(log, "a") as lf:
        while queue or running:
            while queue and len(running) < args.max_jobs \
                    and (not running or in_flight_gb() + queue[0]["predict_gb"] <= args.budget_gb):
                e = queue.pop(0)
                cmd = ["bash", str(HERE / "call_wgs.sh")]
                env_line = f"CONTIGS={e['contig']} THREADS={args.threads}"
                msg = (f"[{time.strftime('%H:%M:%S')}] start {e['contig']} "
                       f"(predict {e['predict_gb']} GB, in flight {in_flight_gb():.1f} GB)")
                print(msg, flush=True); lf.write(msg + "\n"); lf.flush()
                proc = subprocess.Popen(
                    cmd, cwd=str(REPO),
                    env={**__import__("os").environ,
                         "CONTIGS": e["contig"], "THREADS": str(args.threads),
                         "EXTRA": "--nested" if args.nested else "",
                         "W": args.work},
                    stdout=open(work / f"{e['contig']}.schedule.out", "w"),
                    stderr=subprocess.STDOUT)
                running[proc] = e

            if not running:
                break
            time.sleep(5)
            for proc in list(running):
                if proc.poll() is None:
                    continue
                e = running.pop(proc)
                ok = proc.returncode == 0
                msg = (f"[{time.strftime('%H:%M:%S')}] {'done ' if ok else 'FAILED'} "
                       f"{e['contig']} rc={proc.returncode}")
                print(msg, flush=True); lf.write(msg + "\n"); lf.flush()
                results.append({**e, "returncode": proc.returncode})

    elapsed = time.time() - started
    (work / "schedule-results.json").write_text(json.dumps(results, indent=2))
    failed = [r["contig"] for r in results if r["returncode"] != 0]
    print(f"\n{len(results)} contigs in {elapsed/60:.1f} min"
          + (f"; FAILED: {' '.join(failed)}" if failed else ""))
    print("SCHEDULE_DONE" if not failed else "SCHEDULE_FAILED")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
