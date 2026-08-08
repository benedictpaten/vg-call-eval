#!/usr/bin/env python3
"""Call one contig at a time, so only one chromosome is ever in memory.

Prototype answering "how much extra code is needed in vg for per-chromosome
calling?" -- the answer is none. Each step already exists:

    gbz-base query --snarls   extract one contig's subgraph without loading the graph
    vg gbwt -G --gbz-format   rebuild a callable GBZ from that subgraph
    vg call                   call it, with reads streamed per window
    (this script)             loop the contigs and merge the VCFs

Verified to produce records byte-identical to calling the whole graph at once.

**Check `--snarls` works on your graph before relying on this.** It needs precomputed
top-level chains, and `gbz-base construct` finds them only for graphs matching its
assumptions. On the HPRC v2.1 MC CHM13 graph used for tier 2 it found chains for **2 of
46 components** and warned that `--snarls` queries "may not work correctly", which makes
this route unusable there without first building a distance index and supplying
`vg chains graph.gbz graph.dist > graph.chains`. That graph fits in memory anyway, so
`vg chunk --gbz --contig <name>` is both simpler and better there -- it emits GBZ, so it
keeps the GBWT haplotypes that this route flattens to `unknown#N#...`.

Peak memory is therefore one contig's subgraph rather than the whole pangenome.
That is traded against time: the per-contig GBZ rebuild is paid N times, where a
whole-graph run loads once. Not measured at pangenome scale -- see the caveats in
the design doc before trusting the cost side of that trade.

Contig lengths come from a FASTA .fai because gbz-base cannot supply them: an
oversized --interval fails outright ("No successor for GBWT position"), and the
ReferenceIndex table is a sparse position index, not a length table.
"""

from __future__ import annotations

import argparse
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path


def contigs_from_db(db: Path) -> list[str]:
    """Reference contigs, straight from GBZ-Base metadata.

    Reads the Paths table rather than running a query, because there is no CLI
    for listing paths and `vg paths` would load the whole graph -- defeating the
    point. This touches only metadata columns, not the column-compressed
    alignment format that upstream documents as unstable.
    """
    with sqlite3.connect(f"file:{db}?mode=ro", uri=True) as conn:
        rows = conn.execute(
            "SELECT DISTINCT contig FROM Paths WHERE is_indexed = 1 ORDER BY contig"
        ).fetchall()
    return [r[0] for r in rows]


def lengths_from_fai(fai: Path) -> dict[str, int]:
    lengths = {}
    for line in fai.read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        lengths[fields[0]] = int(fields[1])
    return lengths


def run(cmd: list[str], stdout=None) -> None:
    result = subprocess.run(cmd, stdout=stdout, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        sys.exit(f"failed ({result.returncode}): {' '.join(cmd)}\n{result.stderr}")


def call_one_contig(args, contig: str, length: int, work: Path) -> Path:
    gfa, gbz, vcf = work / f"{contig}.gfa", work / f"{contig}.gbz", work / f"{contig}.vcf"

    # --context 0 --snarls is the load-bearing combination. The default 100bp
    # context pulls in neighbouring nodes; context 0 alone *truncates variation*
    # (144 of 215 nodes in testing), which would silently lose alleles. --snarls
    # extends to whole top-level snarls, which is what makes the subgraph
    # callable, and is why the database needs precomputed top-level chains.
    with open(gfa, "w") as out:
        run([args.gbz_base, "query", str(args.gbz_base_db),
             "--contig", contig, "--interval", f"0..{length}",
             "--context", "0", "--snarls"], stdout=out)

    run([args.vg, "gbwt", "-G", str(gfa), "--gbz-format", "-g", str(gbz)])

    # Note the haplotypes come back as `unknown#N#contig#0`: gbz-base keeps
    # metadata only for the queried path. Harmless here -- allele enumeration
    # needs the haplotype *walks*, not their sample names -- but it means the
    # rebuilt graph is not a faithful copy, only a call-equivalent one.
    cmd = [args.vg, "call", str(gbz), "-z", "--read-likelihood", "-t", str(args.threads)]
    if args.gaf_base:
        cmd += ["--gaf-base", str(args.gaf_base), "--gbz-base", str(args.gbz_base_db)]
    else:
        cmd += ["--gam", str(args.gam)]
        if args.gam_index:
            cmd += ["--gam-index", str(args.gam_index)]
    with open(vcf, "w") as out:
        run(cmd, stdout=out)
    return vcf


def merge_vcfs(vcfs: list[Path], out_path: Path) -> int:
    """Union the ##contig lines, keep one copy of everything else, concatenate bodies.

    Needed because each per-contig run only knows its own contig, so every output
    carries a one-contig header.
    """
    preamble, contig_lines, chrom_line, body = [], [], None, []
    for vcf in vcfs:
        for line in vcf.read_text().splitlines():
            if line.startswith("##contig"):
                if line not in contig_lines:
                    contig_lines.append(line)
            elif line.startswith("#CHROM"):
                chrom_line = line
            elif line.startswith("##"):
                if line not in preamble:
                    preamble.append(line)
            elif line.strip():
                body.append(line)

    with open(out_path, "w") as out:
        out.write("\n".join(preamble + contig_lines + [chrom_line] + body) + "\n")
    return len(body)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("gbz_base_db", type=Path, help="GBZ-Base (gbz-base construct)")
    p.add_argument("--fai", type=Path, required=True, help="reference .fai, for contig lengths")
    p.add_argument("-o", "--output", type=Path, required=True)
    p.add_argument("--gam", type=Path)
    p.add_argument("--gam-index", type=Path)
    p.add_argument("--gaf-base", type=Path)
    p.add_argument("--contigs", nargs="*", help="restrict to these contigs")
    p.add_argument("--work", type=Path, default=Path("percontig-work"))
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--vg", default="vg")
    p.add_argument("--gbz-base", default="gbz-base")
    p.add_argument("--keep", action="store_true", help="keep per-contig intermediates")
    args = p.parse_args()

    if not args.gam and not args.gaf_base:
        p.error("need reads: --gam (optionally with --gam-index) or --gaf-base")
    for tool in (args.vg, args.gbz_base):
        if shutil.which(tool) is None:
            sys.exit(f"not on PATH: {tool}")

    args.work.mkdir(parents=True, exist_ok=True)
    lengths = lengths_from_fai(args.fai)
    contigs = args.contigs or contigs_from_db(args.gbz_base_db)

    missing = [c for c in contigs if c not in lengths]
    if missing:
        sys.exit(f"no length in {args.fai} for: {', '.join(missing)}")

    vcfs = []
    for i, contig in enumerate(contigs, 1):
        print(f"[{i}/{len(contigs)}] {contig} ({lengths[contig]:,} bp)", file=sys.stderr)
        vcfs.append(call_one_contig(args, contig, lengths[contig], args.work))

    n = merge_vcfs(vcfs, args.output)
    print(f"{n} variants across {len(contigs)} contigs -> {args.output}", file=sys.stderr)

    if not args.keep:
        shutil.rmtree(args.work, ignore_errors=True)


if __name__ == "__main__":
    main()
