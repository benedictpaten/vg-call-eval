#!/usr/bin/env python3
"""Coverage sweep: chr20 and chr6, short reads and ONT, at nested subsampled depths.

Design, and why:
  * Nested subsamples (scripts/coverage/subsample_gaf.py): the 5x reads are a subset of 10x, so
    every step up the ladder is a paired comparison, not two independent draws.
  * One GAF-Base database per level, built against the WHOLE graph (as titrate.sh and the ONT
    genome recipe both do): a --gaf-reads source has no fetch window and silently disables the
    depth term, which is part of what coverage changes.
  * The call line is arm.py's, flag for flag, so the full-coverage point IS the tier-2 arm. Its VCF
    must reproduce the tier-2 VCF byte for byte (body, sample column aside); that is the control
    that says the sweep measures coverage and nothing else.
  * Short reads stop at the source depth (30.3x on chr20): subsampling cannot reach 40x.

Resumable: every artefact is skipped if already present and complete.

Usage: scripts/coverage/sweep.py [--vg=/path/to/vg]    (tabulate: coverage_table.py; plot: plot_coverage.py)
"""
import concurrent.futures as cf
import json
import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path.home() / "PycharmProjects/vg-call-eval"
os.chdir(REPO)
# The 2026-09-23 run used a pinned, re-signed copy of the binary and a frozen copy of scripts/ (another
# session was editing the repo in place). Defaults here are the repo's own; pass --vg to pin a binary.
VG = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--vg=")),
          str(Path.home() / "CLionProjects/vg/bin/vg"))
PIN = REPO
BIN = Path.home() / ".local/bin"
os.environ["PATH"] = f"{BIN}:{os.environ['PATH']}"
OUT = REPO / "work/cov"
LOG = REPO / "work/cov/run"

SERIES = [
    dict(tech="sr", contig="chr20", gbz="work/wgs-tt/chr20/chr20.gbz", gbzdb="work/graph.hap32.gbz.db",
         fulldb="work/reads.hap32.gaf.db", ref="data/hprc-v2.1-mc-chm13-eval.HG002.hap32.gbz",
         sample="HG002", extra=[], preset=[], levels=[5, 10, 15, 20, 25], ds="chr20-34hap",
         control="work/ont-preset/c20mm095.vcf.gz"),
    dict(tech="sr", contig="chr6", gbz="work/wgs-tt/chr6/chr6.gbz", gbzdb="work/graph.hap32.gbz.db",
         fulldb="work/reads.hap32.gaf.db", ref="data/hprc-v2.1-mc-chm13-eval.HG002.hap32.gbz",
         sample="HG002", extra=[], preset=[], levels=[5, 10, 15, 20, 25], ds="chr6-34hap",
         control="work/ont-preset/c6mm095.vcf.gz"),
    dict(tech="ont", contig="chr20", gbz="work/E821-chr20/chr20.gbz", gbzdb="work/graph.E821.gbz.db",
         fulldb="work/reads.E821.chr20.gaf.db", ref="data/E821-16-sampled.gbz",
         sample="SAMPLE", extra=["--preset", "ont"], preset=["--preset", "long"],
         levels=[5, 10, 15, 20, 25, 30, 40], ds="chr20-34hap", control="work/ont-preset/o20base.vcf.gz"),
    dict(tech="ont", contig="chr6", gbz="work/E821-chr6/chr6_0_chr6.gbz", gbzdb="work/graph.E821.gbz.db",
         fulldb="work/reads.E821.chr6.gaf.db", ref="data/E821-16-sampled.gbz",
         sample="SAMPLE", extra=["--preset", "ont"], preset=["--preset", "long"],
         levels=[5, 10, 15, 20, 25, 30, 40], ds="chr6-34hap", control="work/ont-preset/o6base.vcf.gz"),
]


def ts():
    return time.strftime("%H:%M:%S")


def log(msg):
    print(f"[{ts()}] {msg}", flush=True)


def sh(cmd, **kw):
    r = subprocess.run(cmd, **kw)
    if r.returncode:
        raise RuntimeError(f"failed ({r.returncode}): {' '.join(map(str, cmd))[:300]}")
    return r


def free_gb():
    st = os.statvfs(REPO)
    return st.f_bavail * st.f_frsize / 2**30


def wait_for_disk(need_gb):
    while free_gb() < need_gb:
        log(f"  disk guard: {free_gb():.0f} GB free < {need_gb} GB, waiting")
        time.sleep(120)


def source_gaf(s, d):
    gaf = d / f"{s['contig']}.gaf"
    if gaf.exists() and (d / "source.done").exists():
        return gaf
    wait_for_disk(25)
    log(f"{s['tech']} {s['contig']}: extract source GAF")
    if s["tech"] == "ont":
        sh(["gaf-base", "decompress", s["fulldb"], "-r", s["ref"], "-o", str(gaf)])
    else:
        # The node-ID range is NOT dense (chr20 spans 6.4M IDs for 2.78M nodes), and gbz-base errors on
        # an ID the graph lacks rather than skipping it. List the real nodes, as prep_chrx.sh does.
        nodes = d / "nodes.txt"
        sh(f"{VG} convert -f {s['gbz']} | awk '$1==\"S\"{{print $2}}' | sort -n > {nodes}", shell=True)
        sh([sys.executable, "scripts/tier2/extract_reads_from_db.py", "--nodes", str(nodes),
            "--gaf-base", s["fulldb"], "--gbz-base", s["gbzdb"], "--out", str(gaf),
            "--gbz-base-binary", str(BIN / "gbz-base"), "--tmp", str(d / "extract.tmp.gaf")],
           stdout=open(d / "extract.log", "w"), stderr=subprocess.STDOUT)
        nodes.unlink()
    (d / "source.done").touch()
    return gaf


def subsample_level(s, d, gaf, c):
    """One level per pass. The keep rule is hash(read name) < c/source_coverage, so levels made in
    separate passes are still nested -- and only one level's plain-text GAF exists at a time, which
    matters for ONT, whose GAF lines run to ~100 kb."""
    lgaf = d / f"{s['contig']}.{c}x.gaf"
    if (d / f"db.{c}x.done").exists() or (d / f"db.{c}x.skipped").exists() or lgaf.exists():
        return
    L = int(open(f"work/wgs/{s['contig']}/{s['contig']}.fa.fai").readline().split("\t")[1])
    wait_for_disk(30)
    sh([sys.executable, "scripts/coverage/subsample_gaf.py", "--gaf", str(gaf), "--contig-length", str(L),
        "--out-prefix", str(d / s["contig"]), "--levels", str(c)],
       stdout=open(d / f"subsample.{c}x.log", "w"), stderr=subprocess.STDOUT)


def build_db(s, d, c):
    db = d / f"reads.{c}x.gaf.db"
    if (d / f"db.{c}x.done").exists():
        return db
    lgaf = d / f"{s['contig']}.{c}x.gaf"
    if not lgaf.exists():
        # subsample_gaf.py skips a level at or above the source depth. Starting the pipeline anyway
        # would leave `construct` blocked forever on a FIFO nobody opens for writing.
        (d / f"db.{c}x.skipped").write_text((d / f"subsample.{c}x.log").read_text())
        log(f"{s['tech']} {s['contig']} {c}x: above source depth, skipped")
        return None
    wait_for_disk(10)
    fifo = d / f"sort.{c}x.fifo"
    if fifo.exists():
        fifo.unlink()
    os.mkfifo(fifo)
    srt = subprocess.Popen(["gaf-base", "sort", str(lgaf), *s["preset"], "-o", str(fifo), "-p"],
                           stdout=open(d / f"sort.{c}x.log", "w"), stderr=subprocess.STDOUT)
    con = subprocess.run(["gaf-base", "construct", str(fifo), "-r", s["ref"], "-o", str(db), "--overwrite"],
                         stdout=open(d / f"construct.{c}x.log", "w"), stderr=subprocess.STDOUT)
    srt.wait()
    fifo.unlink()
    if con.returncode or srt.returncode or not db.exists() or db.stat().st_size == 0:
        raise RuntimeError(f"db build failed {s['tech']} {s['contig']} {c}x")
    (d / f"db.{c}x.done").touch()
    lgaf.unlink()                           # the database is the artefact; the GAF is not
    return db


def call(s, d, tag, db):
    vcf = d / f"{tag}.vcf.gz"
    if (d / f"call.{tag}.done").exists():
        return vcf
    raw = d / f"{tag}.vcf"
    cmd = ["/usr/bin/time", "-l", VG, "call", s["gbz"], "-p", f"CHM13#0#{s['contig']}", "-s", s["sample"],
           "-d", "2", "-t", "6", "--progress", "--read-likelihood", "--phased",
           "--gaf-base", str(db), "--gbz-base", s["gbzdb"], *s["extra"]]
    with open(raw, "w") as o, open(d / f"{tag}.log", "w") as e:
        rc = subprocess.call(cmd, stdout=o, stderr=e)
    if rc:
        raise RuntimeError(f"vg call failed: {s['tech']} {s['contig']} {tag}")
    sh(["bgzip", "-f", str(raw)])
    sh(["tabix", "-f", "-p", "vcf", str(vcf)])
    (d / f"call.{tag}.done").touch()
    return vcf


def score(s, d, tag, vcf):
    label = f"cov-{s['tech']}-{s['contig']}-{tag}"
    res = d / f"score.{tag}.json"
    if res.exists():
        return json.loads(res.read_text())
    svcf = vcf
    if s["sample"] != "SAMPLE":             # score_vcf.py hardcodes --cSample SAMPLE
        svcf = d / f"{tag}.s.vcf.gz"
        names = d / "sample.txt"
        names.write_text("SAMPLE\n")
        sh(["bcftools", "reheader", "-s", str(names), "-o", str(svcf), str(vcf)])
        sh(["bcftools", "index", "-t", "-f", str(svcf)])
    sh([sys.executable, str(PIN / "scripts/tier2/score_vcf.py"), "--vcf", str(svcf), "--label", label,
        "--dataset", s["ds"], "--threads", "4"],
       stdout=open(d / f"score.{tag}.log", "w"), stderr=subprocess.STDOUT)
    phase = d / "phase"
    phase.mkdir(exist_ok=True)
    sh([sys.executable, str(PIN / "scripts/tier2/phasing_benchmark.py"), "--calls", str(vcf),
        "--truth", f"work/tier2-{s['contig']}-hap32/truth.{s['contig']}.smvar.vcf.gz",
        "--out", str(phase / tag)], stdout=open(phase / f"{tag}.txt", "w"), stderr=subprocess.STDOUT)
    out = {"label": label, "tag": tag}
    sj = REPO / "work/sv-atlas" / f"score-{label}.json"
    if sj.exists():
        j = json.loads(sj.read_text())
        out["sv"] = j.get("sv")
        out["smallvar"] = j.get("smallvar")
    for line in (phase / f"{tag}.txt").read_text().splitlines():
        if line.startswith("whatshap intersection"):
            f = line.split()
            out["phase"] = {"pairs": int(f[2]), "switches": int(f[3]), "switch_pct": f[4]}
    res.write_text(json.dumps(out, indent=1, default=str))
    return out


def run_series(s):
    d = OUT / f"{s['tech']}-{s['contig']}"
    d.mkdir(parents=True, exist_ok=True)
    gaf = source_gaf(s, d) if not all((d / f"db.{c}x.done").exists() or (d / f"db.{c}x.skipped").exists()
                                      for c in s["levels"]) else None
    for c in s["levels"]:                   # subsample -> build -> delete, one level at a time
        subsample_level(s, d, gaf, c)
        build_db(s, d, c)
    if gaf is not None and gaf.exists():
        gaf.unlink()
    log(f"{s['tech']} {s['contig']}: databases ready")
    jobs = [(f"{c}x", d / f"reads.{c}x.gaf.db") for c in s["levels"]
            if (d / f"db.{c}x.done").exists()] + [("full", Path(s["fulldb"]))]
    for tag, db in jobs:
        vcf = call(s, d, tag, db)
        score(s, d, tag, vcf)
        log(f"{s['tech']} {s['contig']} {tag}: called + scored")
    # control: the full-coverage call must reproduce the tier-2 arm's VCF body
    a = subprocess.run(f"bcftools view -H {d / 'full.vcf.gz'} | cut -f1-9 | md5",
                       shell=True, capture_output=True, text=True).stdout.strip()
    b = subprocess.run(f"bcftools view -H {s['control']} | cut -f1-9 | md5",
                       shell=True, capture_output=True, text=True).stdout.strip()
    ga = subprocess.run(f"bcftools query -f '[%GT]\\n' {d / 'full.vcf.gz'} | md5",
                        shell=True, capture_output=True, text=True).stdout.strip()
    gb = subprocess.run(f"bcftools query -f '[%GT]\\n' {s['control']} | md5",
                        shell=True, capture_output=True, text=True).stdout.strip()
    ok = (a == b and ga == gb)
    (d / "control.txt").write_text(f"sites+GT identical to {s['control']}: {ok}\n")
    log(f"{s['tech']} {s['contig']}: CONTROL full-coverage == tier-2 arm: {ok}")


if __name__ == "__main__":
    LOG.mkdir(parents=True, exist_ok=True)
    log(f"coverage sweep start, vg={VG}")
    # two series at a time: one short-read, one ONT, so their different bottlenecks overlap
    with cf.ThreadPoolExecutor(2) as ex:
        order = [SERIES[0], SERIES[2], SERIES[1], SERIES[3]]   # sr20 | ont20, then sr6 | ont6
        futs = {ex.submit(run_series, s): s for s in order}
        for f in cf.as_completed(futs):
            s = futs[f]
            try:
                f.result()
            except Exception as e:  # noqa: BLE001
                log(f"!! {s['tech']} {s['contig']} FAILED: {e}")
    (LOG / "ALL.done").touch()
    log("COVERAGE_SWEEP_DONE")
