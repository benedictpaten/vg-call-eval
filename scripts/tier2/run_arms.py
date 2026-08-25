#!/usr/bin/env python3
"""Tier 2: run the caller arms on HG002 chr20 and compare against the GIAB draft benchmark.

Implements plan §9.4 steps 5-7. Deliberately a standalone script rather than a
path through `vgcalleval.pipeline`, because that pipeline builds its dataset by
simulation; here the dataset is 28.8 GB of real reads prepared once by hand.
The aardvark engine module is reused as-is.

Six arms, chosen to separate two effects that would otherwise be confounded
(plan §9.2a). Panel enumeration can only propose alleles some haplotype walks,
so on a thin graph its ceiling is the panel's content, not the caller's:

    caller \\ enumeration   support (pack, FlowTraversalFinder)   panel (GBWT walks)
    Poisson                 poisson (its default)                poisson-z
    read-likelihood         readlik-support                      readlik (its default),
                                                                 readlik-nomismap,
                                                                 readlik-nolink

Comparing down a column isolates the caller; comparing across a row isolates
what the sampled graph costs. Without both, a graph limitation reads as a caller
limitation.

The two callers default to different columns, and the arm names track each
caller's own default rather than a single flag. `vg call` now enumerates from the
panel by default under `--read-likelihood` on a GBZ that carries one, because that
measured better on every small-variant class and on three of four SV sets; it does
not do so for the Poisson caller, where the same switch lost SV F1 on all four
datasets. So `readlik` is the plain invocation and `readlik-support` is the one
carrying a flag (`--enumerate-support`), while for Poisson it is the other way
round. The ablation arms sit in the panel column because they ablate against
`readlik`, and an ablation that also changed enumeration would measure two things.
"""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO / "src"))

from vgcalleval.engines import aardvark  # noqa: E402

# Set from --contig; the constants remain as defaults so existing chr20 invocations
# keep working unchanged.
REF_SAMPLE = "CHM13"
REF_PATH = "CHM13#0#chr20"
CONTIG = "chr20"


@dataclass
class Arm:
    name: str
    extra_args: list[str]
    needs_pack: bool
    needs_reads: bool
    description: str = ""
    # filled in by the run
    seconds: float = 0.0
    peak_rss_gb: float = 0.0
    variants: int = 0
    metrics: dict = field(default_factory=dict)


def arms(readlik_extra: list[str] | None = None,
         readlik_panel_extra: list[str] | None = None) -> list[Arm]:
    """The six fixed arms.

    `readlik_extra` is appended to the read-likelihood arms only, so a caller-side change
    under evaluation can be measured across the whole matrix without editing this list --
    and without touching the Poisson arms, which it must not affect. Passing something that
    *does* affect them would break the comparison silently, so keep it to read-likelihood
    flags.

    `readlik_panel_extra` goes to the panel-enumeration arms alone, for flags that need a
    panel to act on. `--linkage-weight` is the case that forced this: `vg call` refuses it
    where enumeration is not from haplotypes, so putting it in `readlik_extra` makes
    `readlik-support` exit immediately. It then scores as zero variants with empty metrics,
    the matrix completes "successfully", and the page build is what finally fails -- several
    arms and forty minutes downstream of the actual error. `readlik-nolink` is excluded too,
    since it sets `--linkage-weight` itself and a second one would silently override it.
    """
    extra = list(readlik_extra or [])
    panel_extra = list(readlik_panel_extra or [])
    return [
        Arm("poisson", [], True, False,
            "current default: Poisson genotyping, support enumeration"),
        # Needs the pack despite -z: the pack becomes optional only for --read-likelihood,
        # because Poisson *genotyping* consumes support even when enumeration comes from
        # haplotypes. -z is explicit here and has to stay that way: it is not the Poisson
        # caller's default, and on this benchmark it costs that caller SV F1.
        Arm("poisson-z", ["-z"], True, False,
            "Poisson genotyping, panel enumeration from the graph's haplotypes"),
        # No -z: panel enumeration is what --read-likelihood does by default on a GBZ that
        # carries a panel. Spelling it out would still work, but the arm is meant to be the
        # shipped default, so it is run the way a user would run it.
        Arm("readlik", ["--read-likelihood"] + extra + panel_extra, False, True,
            "read-level likelihoods as shipped: panel enumeration, no pack file"),
        Arm("readlik-nomismap",
            ["--read-likelihood", "--no-mismap-term"] + extra + panel_extra, False, True,
            "as readlik, MAPQ mismapping term disabled, to measure its contribution"),
        # The linkage ablation, and the reason it is a standing arm rather than a one-off
        # measurement: --linkage-weight defaults to 2, so readlik above carries the HMM and
        # nothing in the matrix would show what it contributes. Measured once, the transition
        # model is about 12% of the genotype-F1 gain over no linkage at all -- worth keeping
        # visible in the same table as the rest rather than buried in a planning note.
        Arm("readlik-nolink", ["--read-likelihood", "--linkage-weight", "0"] + extra,
            False, True,
            "as readlik, the whole linkage layer off -- transition model, phasing and the settled\n            genotypes that follow from it"),
        # The like-for-like caller comparison, and the only read-likelihood arm that needs a
        # flag to get its enumeration: it holds enumeration fixed against `poisson` so the
        # difference between them is the genotyper alone.
        Arm("readlik-support", ["--read-likelihood", "--enumerate-support"] + extra, True, True,
            "read-level likelihoods, support enumeration -- the like-for-like caller comparison"),
    ]


def sh(cmd: list[str], out_path: Path | None = None, log: Path | None = None,
       time_it: bool = False) -> tuple[int, float, float]:
    """Run a command, optionally under /usr/bin/time -l. Returns (rc, seconds, peak_rss_gb)."""
    full = ["/usr/bin/time", "-l"] + cmd if time_it else cmd
    started = time.time()
    with open(out_path, "wb") if out_path else open("/dev/null", "wb") as out:
        proc = subprocess.run(full, stdout=out, stderr=subprocess.PIPE)
    elapsed = time.time() - started
    stderr = proc.stderr.decode(errors="replace")
    if log:
        log.write_text(stderr)
    peak = 0.0
    for line in stderr.splitlines():
        if "maximum resident set size" in line:
            peak = int(line.split()[0]) / 1073741824
    return proc.returncode, elapsed, peak


def rename_contig(src: Path, dst: Path) -> int:
    """Rewrite CHM13#0#chr20 -> chr20 in CHROM and ##contig lines.

    The benchmark VCF uses bare `chr20`. Skipping this does not error -- it
    produces zero true positives, which looks like a catastrophic caller failure
    rather than a naming mismatch. Asserted on below.

    **Now usually a no-op, and deliberately kept anyway.** Upstream vg changed
    `vg call` to report the *base* path name in CHROM, so records already arrive as
    `chr20` rather than `CHM13#0#chr20`; the pass-through branch below handles that
    without comment. Removing this would make the harness silently dependent on
    which vg build produced the VCF, which is exactly the class of breakage it
    exists to prevent -- and the failure mode is zero true positives, not an error.
    """
    kept = 0
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if line.startswith("##contig="):
                fout.write(line.replace(f"ID={REF_PATH}", f"ID={CONTIG}"))
            elif line.startswith("#"):
                fout.write(line)
            else:
                fields = line.split("\t", 1)
                if fields[0] == REF_PATH:
                    fout.write(f"{CONTIG}\t{fields[1]}")
                else:
                    fout.write(line)
                kept += 1
    return kept


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--vg", default=str(Path.home() / "CLionProjects/vg/bin/vg"))
    p.add_argument("--graph", default=str(HERE / "chr20_0_chr20.gbz"))
    p.add_argument("--pack", default=str(HERE / "chr20.pack"))
    p.add_argument("--gaf-base", default=str(REPO / "work/reads.gaf.db"))
    p.add_argument("--gbz-base", default=str(REPO / "work/graph.gbz.db"))
    p.add_argument("--reference", default=str(HERE / "chr20.fa"))
    p.add_argument("--truth-vcf", default=str(HERE / "truth.chr20.smvar.vcf.gz"))
    p.add_argument("--truth-bed", default=str(HERE / "truth.chr20.smvar.bed"))
    p.add_argument("--aardvark", default="aardvark")
    p.add_argument("--threads", type=int, default=6)
    # Not passed unless set: vg picks a per-backend default (4096 for --gaf-base),
    # and the harness should track that rather than pin a value that silently
    # becomes wrong when the default moves.
    p.add_argument("--read-window", type=int, default=0)
    p.add_argument("--out", default=str(HERE / "results"))
    p.add_argument("--only", nargs="*", help="run only these arms")
    # One shell-quoted string, shlex-split here, rather than nargs="*": argparse stops
    # consuming a variadic list at the first token that looks like an option, so
    # `--readlik-extra --depth-term 0.1` silently became two stray positionals.
    p.add_argument("--readlik-extra", default="",
                   help="extra flags appended to the read-likelihood arms only, as "
                        "one quoted string, for evaluating a caller-side change across the "
                        "whole matrix")
    p.add_argument("--readlik-panel-extra", default="",
                   help="extra flags for the panel-enumeration read-likelihood arms alone, as "
                        "one quoted string, for flags that need a haplotype panel to act on "
                        "and would make readlik-support exit immediately")
    p.add_argument("--contig", default="chr20",
                   help="contig to call; sets both the reference path and the "
                        "name written into the output VCF")
    args = p.parse_args()

    global REF_PATH, CONTIG
    CONTIG = args.contig
    REF_PATH = f"{REF_SAMPLE}#0#{args.contig}"

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    selected = [a for a in arms(shlex.split(args.readlik_extra),
                                shlex.split(args.readlik_panel_extra))
                if not args.only or a.name in args.only]

    for arm in selected:
        print(f"\n=== {arm.name}: {arm.description}", flush=True)
        raw = out / f"{arm.name}.raw.vcf"
        cmd = [args.vg, "call", args.graph, "-p", REF_PATH,
               "-t", str(args.threads), "--progress"] + arm.extra_args
        if arm.needs_pack:
            cmd += ["-k", args.pack]
        if arm.needs_reads:
            cmd += ["--gaf-base", args.gaf_base, "--gbz-base", args.gbz_base]
            if args.read_window:
                cmd += ["--read-window", str(args.read_window)]

        rc, secs, peak = sh(cmd, out_path=raw, log=out / f"{arm.name}.call.log", time_it=True)
        arm.seconds, arm.peak_rss_gb = secs, peak
        if rc != 0:
            print(f"  FAILED rc={rc}; see {arm.name}.call.log", flush=True)
            tail = (out / f"{arm.name}.call.log").read_text().splitlines()[-6:]
            print("  " + "\n  ".join(tail), flush=True)
            continue

        renamed = out / f"{arm.name}.vcf"
        arm.variants = rename_contig(raw, renamed)
        print(f"  {arm.variants} variants, {secs:.0f}s, {peak:.1f} GB peak", flush=True)
        raw.unlink()

        gz = renamed.with_suffix(".vcf.gz")
        subprocess.run(["bgzip", "-f", "-c", str(renamed)], stdout=open(gz, "wb"), check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", str(gz)], check=True)

        sample = subprocess.run(["bcftools", "query", "-l", str(gz)],
                                capture_output=True, text=True).stdout.split()
        query_sample = sample[0] if sample else "SAMPLE"

        adir = out / f"aardvark-{arm.name}"
        try:
            aardvark.compare(
                aardvark=args.aardvark,
                reference=Path(args.reference),
                truth_vcf=Path(args.truth_vcf),
                query_vcf=gz,
                regions_bed=Path(args.truth_bed),
                out_dir=adir,
                truth_sample="HG002",
                query_sample=query_sample,
                label=arm.name,
                options=aardvark.AardvarkOptions(threads=args.threads),
            )
            arm.metrics = {"summary": aardvark.read_summary(adir)}
        except Exception as exc:  # noqa: BLE001
            print(f"  aardvark failed: {exc}", flush=True)

    # An arm that produced nothing is a failed run, not a result. Until this check existed a
    # dead arm was logged, skipped, and written out as zero variants with empty metrics -- the
    # matrix then "succeeded", and the first visible symptom was a KeyError in the page build,
    # forty minutes and several arms downstream of the flag that caused it. Fail here instead,
    # where the log that explains it is one line up.
    empty = [a.name for a in selected if a.variants == 0]
    if empty:
        print(f"\nERROR: no variants from {', '.join(empty)} -- see the .call.log for each",
              file=sys.stderr)
        sys.exit(1)

    # Which vg produced these numbers, recorded rather than assumed.
    #
    # "One build, one pass" is the rule this harness exists to enforce, and until now it rested on
    # procedure: nothing in the results said which binary ran, so a page could claim it and be
    # wrong. It has already happened once here. Recorded per arm, not per file, because a run with
    # --only merges into an existing arms.json and a single top-level field would then describe
    # whichever run wrote last.
    #
    # `vg version` is the git describe baked in at *build* time, so a binary built from a dirty
    # tree names the last commit rather than what was compiled. That happened on the run this was
    # added for. It is the right field to record anyway -- it is what the binary says about itself
    # -- but a page claiming a commit is claiming the build, not the working tree.
    try:
        vg_version = subprocess.run([args.vg, "version"], capture_output=True, text=True,
                                    timeout=60).stdout.strip().splitlines()[0]
    except Exception:
        vg_version = ""
    payload = [
        {"arm": a.name, "description": a.description, "variants": a.variants,
         "seconds": round(a.seconds, 1), "peak_rss_gb": round(a.peak_rss_gb, 2),
         "vg_version": vg_version, "metrics": a.metrics}
        for a in selected
    ]
    # Merge rather than overwrite. A run restricted to --only would otherwise delete
    # every arm it was not asked about, and the file would still look well-formed --
    # it would just silently describe a smaller experiment. The same bug bit
    # compare_sv.py, and it cost a re-run here before being fixed.
    out_path = out / "arms.json"
    merged = {}
    if out_path.exists():
        for entry in json.loads(out_path.read_text()):
            merged[entry["arm"]] = entry
    for entry in payload:
        merged[entry["arm"]] = entry
    out_path.write_text(json.dumps(list(merged.values()), indent=2))
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
