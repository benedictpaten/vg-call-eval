#!/usr/bin/env python3
"""Follow-on to H4: can the scorer tell *similar* alleles apart, and how similar is too similar?

The first H4 test used random, mutually distinguishable insertion sequences and the
caller recovered every truth haplotype exactly -- so the free-unanchored-sequence rule
is not a length bias. But that is not the case the 34-haplotype graph presents. There the
competing alleles come from 32 *recombinant* haplotypes, so they are variations on each
other, and §9.23 found the false-call rate reaching 0.768 at sites offering ten or more.

This builds that case directly: N insertion alleles of the same length differing by k
substitutions each, one of them true, reads simulated from it with a realistic error rate.
Sweeping N and k separates two explanations that the real data cannot:

  - if accuracy falls with N at fixed k, the failure is *multiplicity* -- more candidates,
    more chances for noise to win, which is winner's curse in the evidence rather than
    the decision;
  - if accuracy falls with decreasing k at fixed N, the failure is *resolution* -- the
    reads cannot separate alleles this similar, and no amount of reweighting will help.

The two have different remedies, which is why it is worth telling them apart.
"""

from __future__ import annotations

import argparse
import random
import subprocess
import sys
from pathlib import Path

BASES = "ACGT"


def mutate(s: str, k: int, rng: random.Random) -> str:
    """k substitutions at distinct positions, each to a different base."""
    out = list(s)
    for i in rng.sample(range(len(s)), k):
        out[i] = rng.choice([b for b in BASES if b != out[i]])
    return "".join(out)


def sh(cmd: list[str], out=None) -> None:
    r = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        sys.exit(f"failed: {' '.join(cmd)}\n{r.stderr[-1500:]}")


def build(work: Path, n_alleles: int, k: int, ins_len: int, flank: int,
          rng: random.Random) -> None:
    f1, f2 = ("".join(rng.choice(BASES) for _ in range(flank)) for _ in range(2))
    base = "".join(rng.choice(BASES) for _ in range(ins_len))
    alleles = [base] + [mutate(base, k, rng) for _ in range(n_alleles - 1)]

    # Each allele is its own node between the flanks: a bubble of parallel alternates,
    # which is what a snarl over related haplotypes looks like once decomposed.
    nodes = {1: f1, 2: f2}
    nid = 3
    ids = []
    for a in alleles:
        nodes[nid] = a
        ids.append(nid)
        nid += 1
    with open(work / "site.gfa", "w") as fh:
        fh.write("H\tVN:Z:1.1\n")
        for i, s in nodes.items():
            fh.write(f"S\t{i}\t{s}\n")
        fh.write("L\t1\t+\t2\t+\t0M\n")          # reference: no insertion
        for i in ids:
            fh.write(f"L\t1\t+\t{i}\t+\t0M\n")
            fh.write(f"L\t{i}\t+\t2\t+\t0M\n")
        fh.write("P\tref\t1+,2+\t*\n")
        for h, i in enumerate(ids):
            fh.write(f"W\tsim\t{h}\tref\t0\t{2*flank+ins_len}\t>1>{i}>2\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vg", default="/Users/benedictpaten/CLionProjects/vg/bin/vg")
    ap.add_argument("--work", required=True)
    ap.add_argument("--ins-len", type=int, default=300)
    ap.add_argument("--flank", type=int, default=600)
    ap.add_argument("--depth", type=int, default=30)
    ap.add_argument("--readlen", type=int, default=150)
    ap.add_argument("--error", type=float, default=0.01)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--truth-ref", action="store_true",
                    help="simulate from the reference (no insertion) while the graph "
                         "offers N insertion alleles. This is the false-positive mode: "
                         "the sample carries nothing the alternates describe, so every "
                         "non-reference call is an error. On real data most sites are "
                         "like this -- the graph offers haplotypes the sample lacks.")
    ap.add_argument("--map", action="store_true",
                    help="map the simulated reads with giraffe instead of using vg sim's "
                         "true paths -- which is the whole point, see below")
    args = ap.parse_args()

    W = Path(args.work)
    W.mkdir(parents=True, exist_ok=True)
    print(f"insertion {args.ins_len} bp, flank {args.flank}, depth {args.depth}x, "
          f"{args.readlen} bp reads, error {args.error}, {args.reps} reps, "
          f"reads {'MAPPED with giraffe' if args.map else 'from vg sim true paths'}")
    if args.truth_ref:
        print("TRUTH = REFERENCE: every non-reference call is a false positive")
    print(f"{'alleles':>8}{'k subs':>8}{'correct':>9}{'as ref':>8}{'wrong alt':>11}{'GQ med':>8}")
    for n_alleles in (2, 4, 8, 16):
        for k in (20, 5, 2, 1):
            correct = as_ref = wrong = 0
            gqs = []
            for rep in range(args.reps):
                rng = random.Random(1000 * rep + 7)
                build(W, n_alleles, k, args.ins_len, args.flank, rng)
                sh([args.vg, "gbwt", "-G", str(W / "site.gfa"), "--gbz-format",
                    "-g", str(W / "site.gbz")])
                n = max(1, (2 * args.flank + args.ins_len) * args.depth // args.readlen)
                with open(W / "sim.gam", "wb") as fh:
                    sh([args.vg, "sim", "-x", str(W / "site.gbz"), "-n", str(n),
                        "-l", str(args.readlen), "-a", "-s", str(rep + 1),
                        "-e", str(args.error),
                        "-P", "ref" if args.truth_ref else "sim#0#ref#0"], out=fh)
                reads = W / "sim.gam"
                if args.map:
                    # vg sim -a emits each read already on its true graph path, so the
                    # scorer is handed the answer and discrimination is trivial. Real
                    # reads arrive from a mapper that had to choose among the alleles,
                    # and it is that choice -- not the scoring -- that similar
                    # haplotypes make hard. Round-trip through giraffe to restore it.
                    # Clear stale index artefacts: the graph is rebuilt every iteration,
                    # and giraffe refuses to run when a derived index is older than the
                    # distance index it depends on.
                    for stale in W.glob("site.*"):
                        if stale.suffix not in (".gfa", ".gbz"):
                            stale.unlink(missing_ok=True)
                    sh([args.vg, "index", "-j", str(W / "site.dist"), str(W / "site.gbz")])
                    sh([args.vg, "minimizer", "-d", str(W / "site.dist"),
                        "-o", str(W / "site.min"), str(W / "site.gbz")])
                    with open(W / "reads.fq", "w") as fh:
                        sh([args.vg, "view", "-X", str(W / "sim.gam")], out=fh)
                    with open(W / "mapped.gam", "wb") as fh:
                        sh([args.vg, "giraffe", "-Z", str(W / "site.gbz"),
                            "-d", str(W / "site.dist"), "-m", str(W / "site.min"),
                            "-f", str(W / "reads.fq"), "-t", "2"], out=fh)
                    reads = W / "mapped.gam"
                with open(W / "out.vcf", "w") as fh:
                    sh([args.vg, "call", str(W / "site.gbz"), "-z", "--read-likelihood",
                        "--gam", str(reads), "-t", "2"], out=fh)
                hit = False
                any_call = False
                for line in open(W / "out.vcf"):
                    if line.startswith("#"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    d = dict(zip(f[8].split(":"), f[9].split(":")))
                    gt = d.get("GT", "./.")
                    idx = {int(x) for x in gt.replace("|", "/").split("/")
                           if x.isdigit() and int(x) > 0}
                    if d.get("GQ", "").isdigit():
                        gqs.append(int(d["GQ"]))
                    if not idx:
                        continue
                    any_call = True
                    alts = f[4].split(",")
                    # truth is allele index 0 of the built set, i.e. the unmutated base
                    for i in idx:
                        if i <= len(alts) and len(alts[i - 1]) - len(f[3]) == args.ins_len:
                            hit = True
                if args.truth_ref:
                    # Truth is reference: any non-reference call is a false positive.
                    if not any_call:
                        correct += 1
                    else:
                        wrong += 1
                elif not any_call:
                    as_ref += 1
                elif hit:
                    correct += 1
                else:
                    wrong += 1
            gm = sorted(gqs)
            print(f"{n_alleles:>8}{k:>8}{correct:>9}{as_ref:>8}{wrong:>11}"
                  f"{gm[len(gm)//2] if gm else 0:>8}")


if __name__ == "__main__":
    main()
