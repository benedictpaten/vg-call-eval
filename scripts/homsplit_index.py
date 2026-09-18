#!/usr/bin/env python3
"""One pass over an anchors file -> a compact index for the hom-split switch analysis.

Emits a pickle holding, per SITE (a snarl):
    kind     'het' (2 slots, different alleles), 'homsplit' (2 slots, SAME allele),
             'one' (a single slot: an unsplit homozygote or a haploid -- the file cannot
             tell them apart, and the sweep never needs to, because a site unsplit at the
             defaults stays unsplit at every stricter setting)
    node     the smallest boundary node id, as a positional surrogate
    reads    the read ids pinned there, deduplicated across both pins and both slots
    slot_of  read id -> slot, for the split sites (both pins agree, so last writer is fine)

and per READ the sites it touches. Read ids are FILE-LOCAL indices into this file's own
`#read` table, so the names are kept too: nothing may be joined across anchor files by id.
"""
import argparse, collections, pickle, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    names = []
    # snarl -> {slot: allele}, snarl -> node, snarl -> {read: slot}
    slots = collections.defaultdict(dict)
    node_of = {}
    reads_of = collections.defaultdict(dict)
    read_sites = collections.defaultdict(set)

    cur_snarl, cur_slot = None, None
    n_a = n_r = 0
    for line in open(a.anchors):
        if line[0] == "#":
            if line.startswith("#read\t"):
                c = line.rstrip("\n").split("\t")
                rid = int(c[1])
                while len(names) <= rid:
                    names.append(None)
                names[rid] = c[2]
            continue
        c = line.rstrip("\n").split("\t")
        if c[0] == "A":
            node, snarl, slot, allele = int(c[1]), c[2], int(c[3]), int(c[4])
            slots[snarl][slot] = allele
            if snarl not in node_of or node < node_of[snarl]:
                node_of[snarl] = node
            cur_snarl, cur_slot = snarl, slot
            n_a += 1
        elif c[0] == "R" and cur_snarl is not None:
            rid = int(c[1])
            reads_of[cur_snarl][rid] = cur_slot
            read_sites[rid].add(cur_snarl)
            n_r += 1

    kind = {}
    for snarl, sl in slots.items():
        if len(sl) == 2:
            kind[snarl] = "het" if sl[0] != sl[1] else "homsplit"
        else:
            kind[snarl] = "one"

    out = {
        "names": names,
        "kind": kind,
        "node_of": node_of,
        "reads_of": {s: dict(d) for s, d in reads_of.items()},
        "read_sites": {r: sorted(s) for r, s in read_sites.items()},
    }
    with open(a.out, "wb") as f:
        pickle.dump(out, f, protocol=4)
    counts = collections.Counter(kind.values())
    print(f"{n_a} A rows, {n_r} R rows, {len(names)} reads", file=sys.stderr)
    print(f"sites: {dict(counts)}", file=sys.stderr)


if __name__ == "__main__":
    main()
