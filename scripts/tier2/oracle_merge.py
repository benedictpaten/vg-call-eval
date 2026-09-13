#!/usr/bin/env python3
"""Site-wise oracle for a two-floor --mismap-min scheme on ONT.

The ONT floor sweep pulls SNVs and indels in opposite directions: HP>=5 false positives
fall monotonically as the floor rises, while non-HP FPs bottom out around 0.10 and then
climb. One scalar cannot sit at both optima. This asks what a two-floor scheme could win
BEFORE anyone implements one.

Method. VCF column 3 is the snarl id (`>node>node`), which is genotype-independent and
therefore identical across arms -- POS is not (a record's POS depends on its genotype).
So each snarl is assigned wholesale to one arm, all of its SB blocks together, and the
merged file is a coherent genotyping in which every site was decided by exactly one arm.

Two assignment rules:
  class  -- a snarl any arm spells with a length-changing ALT goes to the high floor.
            The ceiling of what a per-variant-class parameter could reach.
  hp     -- a snarl whose reference context is a homopolymer run of >= HPMIN goes to the
            high floor. This one is deployable: the run length is a property of the
            reference, knowable before genotyping, which the class rule is not.

What this CANNOT show: the linkage layer couples sites, so a real two-floor caller would
reach a different linkage solution than either arm did. This is an upper bound on the
per-site gain, not a prediction of the shipped number.
"""
import sys, gzip, subprocess
from pathlib import Path

# arm.py is the ONT harness and lives under work/, which is gitignored; this driver is
# tracked, so it resolves that path explicitly rather than by its own location.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "work/ont-preset"))
import arm  # noqa: E402  -- reuse its reference, scorer and reporting

HPMIN = int(__import__("os").environ.get("HPMIN", "5"))


def read_arm(path):
    """snarl id -> list of raw record lines, plus the header lines."""
    hdr, sites = [], {}
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                hdr.append(line)
                continue
            c = line.split("\t", 5)
            sites.setdefault(c[2], []).append(line)
    return hdr, sites


def length_changing(lines):
    for line in lines:
        c = line.split("\t", 5)
        ref = c[3]
        if any(len(a) != len(ref) for a in c[4].split(",") if a != "."):
            return True
    return False


def hp_context(lines):
    """Max reference homopolymer run length anchored just past each record's POS."""
    s = arm.ref_seq()
    best = 0
    for line in lines:
        c = line.split("\t", 5)
        i = int(c[1])  # 0-based index of POS+1, the first base after the anchor
        if 0 <= i < len(s):
            best = max(best, arm.hp_run(i + 1, s[i]))
    return best


def merge(lo_vcf, hi_vcf, rule, out):
    hdr, lo = read_arm(lo_vcf)
    _, hi = read_arm(hi_vcf)
    picked = {"lo": 0, "hi": 0}
    rows = []
    for sid in set(lo) | set(hi):
        both = lo.get(sid, []) + hi.get(sid, [])
        if rule == "class":
            want_hi = length_changing(both)
        elif rule == "hp":
            want_hi = hp_context(both) >= HPMIN
        else:
            raise SystemExit(f"unknown rule {rule}")
        src = (hi if want_hi else lo).get(sid)
        if src is None:          # that arm did not call this site at all: it is a no-call there
            picked["hi" if want_hi else "lo"] += 0
            continue
        picked["hi" if want_hi else "lo"] += 1
        rows.extend(src)
    rows.sort(key=lambda l: (int(l.split("\t", 2)[1]), l.split("\t", 3)[2]))
    with open(out, "w") as f:
        f.writelines(hdr)
        f.writelines(rows)
    print(f"[{out.stem}] rule={rule} sites from low floor {picked['lo']}, "
          f"high floor {picked['hi']}, {len(rows)} records", flush=True)
    subprocess.check_call(["bgzip", "-f", str(out)])
    subprocess.check_call(["tabix", "-f", "-p", "vcf", str(out) + ".gz"])
    return Path(str(out) + ".gz")


if __name__ == "__main__":
    tag, lo_vcf, hi_vcf, rule = sys.argv[1:5]
    arm.select_arm(__import__("os").environ.get("ORACLE_ARM", "ont-chr20"))
    vcf = merge(Path(lo_vcf), Path(hi_vcf), rule, arm.OUT / f"{tag}.vcf")
    adir = arm.aardvark(tag, vcf)
    row = arm.report(tag, adir, vcf)
    # arm.py writes this from main(); the oracle has its own entry point, so do it here too
    # or the numbers survive only in a log.
    import json
    (arm.OUT / f"{tag}.json").write_text(json.dumps(row, indent=1))
    print(f"[{tag}] wrote {arm.OUT / (tag + '.json')}", flush=True)
