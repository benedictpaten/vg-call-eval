#!/usr/bin/env python3
"""H4 of plan §9.22: can the scorer tell nested insertion alleles apart?

`score_read_against_allele` charges every read base but no allele sequence falling
before a read's first anchor or after its last. The comment says that is deliberate --
charging it "would penalise a read for being short" -- and on reflection that is
information-theoretically right: a read ending inside a long allele says nothing about
sequence beyond its reach. Reads that *span* the site still charge the whole allele as
an internal gap, so a wrong long allele should be penalised by those.

The question is therefore not whether the rule is defensible but whether enough
spanning signal survives when alleles are long and nested. This builds the adversarial
case for it: four alleles sharing prefixes, so telling them apart depends entirely on
reads that reach the 3' junction.

    REF   flank1 ------------------------------- flank2
    A1    flank1 -> p(100) ---------------------- flank2
    A2    flank1 -> p(100) -> q(200) ------------ flank2
    A3    flank1 -> p(100) -> q(200) -> r(700) -- flank2

A read that ends inside p cannot distinguish A1, A2 and A3 -- correctly, since it has
seen no evidence either way. Only reads crossing each allele's own 3' junction can.
With 150 bp reads, A3's junction is 1000 bp from A1's, so the alleles differ in how much
spanning evidence exists for them. If the model still recovers each truth haplotype,
H4 is refuted as a source of bias; if it drifts toward the longest allele, it is not.
"""

from __future__ import annotations

import argparse
import random
import subprocess
import sys
from pathlib import Path

BASES = "ACGT"


def seq(n: int, rng: random.Random) -> str:
    return "".join(rng.choice(BASES) for _ in range(n))


def build_gfa(path: Path, rng: random.Random, flank: int) -> dict:
    f1, f2 = seq(flank, rng), seq(flank, rng)
    p, q, r = seq(100, rng), seq(200, rng), seq(700, rng)
    nodes = {1: f1, 2: p, 3: q, 4: r, 5: f2}
    links = [(1, 2), (2, 3), (3, 4), (4, 5), (1, 5), (2, 5), (3, 5)]
    walks = {
        "REF": ">1>5",
        "A1": ">1>2>5",
        "A2": ">1>2>3>5",
        "A3": ">1>2>3>4>5",
    }
    with open(path, "w") as fh:
        fh.write("H\tVN:Z:1.1\n")
        for nid, s in nodes.items():
            fh.write(f"S\t{nid}\t{s}\n")
        for a, b in links:
            fh.write(f"L\t{a}\t+\t{b}\t+\t0M\n")
        fh.write(f"P\tref\t1+,5+\t*\n")
        for i, (name, walk) in enumerate(walks.items()):
            ln = sum(len(nodes[int(x)]) for x in walk.replace(">", " ").split())
            fh.write(f"W\tsim\t{i}\tref\t0\t{ln}\t{walk}\n")
    return {"lengths": {"REF": 0, "A1": 100, "A2": 300, "A3": 1000}}


def sh(cmd: list[str], out=None) -> None:
    r = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        sys.exit(f"failed: {' '.join(cmd)}\n{r.stderr[-1500:]}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vg", default="/Users/benedictpaten/CLionProjects/vg/bin/vg")
    ap.add_argument("--work", required=True)
    ap.add_argument("--flank", type=int, default=600)
    ap.add_argument("--depth", type=int, default=30)
    ap.add_argument("--readlen", type=int, default=150)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    W = Path(args.work)
    W.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    info = build_gfa(W / "site.gfa", rng, args.flank)
    sh([args.vg, "gbwt", "-G", str(W / "site.gfa"), "--gbz-format", "-g", str(W / "site.gbz")])

    print(f"{'truth':<7}{'true len':>9}   called GT / ALT length     verdict")
    for truth, tlen in info["lengths"].items():
        # Simulate from this haplotype's walk only, so the sample is homozygous for it.
        # The GBZ names walks <sample>#<hap>#<contig>#<frag>.
        hap_path = f"sim#{list(info['lengths']).index(truth)}#ref#0"
        span = 2 * args.flank + tlen
        n = max(1, span * args.depth // args.readlen)
        gam = W / f"{truth}.gam"
        with open(gam, "wb") as fh:
            sh([args.vg, "sim", "-x", str(W / "site.gbz"), "-n", str(n), "-l", str(args.readlen),
                "-a", "-s", str(args.seed), "-P", hap_path], out=fh)
        vcf = W / f"{truth}.vcf"
        with open(vcf, "w") as fh:
            sh([args.vg, "call", str(W / "site.gbz"), "-z", "--read-likelihood",
                "--gam", str(gam), "-t", "2"], out=fh)
        got = []
        for line in open(vcf):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            gt = f[9].split(":")[0]
            idx = {int(x) for x in gt.replace("|", "/").split("/") if x.isdigit() and int(x) > 0}
            alts = f[4].split(",")
            for i in sorted(idx):
                if i <= len(alts):
                    got.append((gt, len(alts[i - 1]) - len(f[3])))
        if not got:
            called = "0/0 (no ALT)"
            ok = (tlen == 0)
        else:
            called = ", ".join(f"{g} len {L}" for g, L in got)
            ok = any(L == tlen for _, L in got)
        print(f"{truth:<7}{tlen:>9}   {called:<28}{'OK' if ok else 'WRONG'}")


if __name__ == "__main__":
    main()
