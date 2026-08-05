#!/usr/bin/env python3
"""Extract one chromosome's reads from the GAF-Base, instead of streaming the whole GAF.

Measured ~25x faster than the streaming filter: the database is sorted by node ID
and randomly accessible, so chr20's reads are fetched directly rather than found by
reading 28.8 GB of gzip. ~0.25 s per 4096-node query returning ~29,000 reads, so
chr20's 2,382,533 nodes cost roughly 2.5 minutes against well over an hour.

De-duplication is required, and the obvious key is wrong. A read overlapping nodes
in two adjacent chunks is returned by both queries, and feeding duplicates to
`vg pack` would double-count coverage at chunk boundaries -- inflating exactly the
support the Poisson caller genotypes from.

But **paired mates share a read name**: 20,000 records of this GAF carry 10,000
distinct names. De-duplicating on the name alone silently drops one mate of every
pair, which showed up here as an implausible ~50% duplicate rate. The key is
therefore name plus start position; two records agreeing on both are the same
alignment returned twice, while mates differ in position and are both kept.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--nodes", default=str(HERE / "chr20_all_nodes.txt"))
    p.add_argument("--gaf-base", default=str(HERE.parent / "reads.gaf.db"))
    p.add_argument("--gbz-base", default=str(HERE.parent / "graph.gbz.db"))
    p.add_argument("--out", default=str(HERE / "chr20.reads.fromdb.gaf"))
    p.add_argument("--chunk", type=int, default=4096)
    p.add_argument("--gbz-base-binary", default="gbz-base")
    p.add_argument("--tmp", default="/tmp/gafbase_extract.gaf")
    args = p.parse_args()

    nodes = [l.strip() for l in open(args.nodes) if l.strip()]
    total_chunks = (len(nodes) + args.chunk - 1) // args.chunk
    print(f"{len(nodes):,} nodes in {total_chunks} chunks of {args.chunk}", flush=True)

    seen: set[tuple] = set()
    written = dupes = 0
    started = time.time()

    with open(args.out, "w") as out:
        for i in range(total_chunks):
            sl = nodes[i * args.chunk:(i + 1) * args.chunk]
            cmd = [args.gbz_base_binary, "query", args.gbz_base]
            for n in sl:
                cmd += ["-n", n]
            cmd += ["--context", "0", "--gaf-base", args.gaf_base,
                    "--gaf-output", args.tmp, "--alignments", "overlapping"]
            r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
            if r.returncode != 0:
                sys.exit(f"chunk {i} failed: {r.stderr.strip()[:300]}")
            with open(args.tmp) as fh:
                for line in fh:
                    if line.startswith("@"):
                        continue
                    f = line.split("\t", 6)
                    # name + strand + path + query start: enough to tell mates apart
                    key = (f[0], f[4], f[5], f[2]) if len(f) >= 6 else (line,)
                    if key in seen:
                        dupes += 1
                        continue
                    seen.add(key)
                    out.write(line)
                    written += 1
            if (i + 1) % 50 == 0 or i + 1 == total_chunks:
                el = time.time() - started
                rate = (i + 1) / el
                print(f"  {i+1}/{total_chunks} chunks, {written:,} reads, {dupes:,} dupes, "
                      f"{el:.0f}s ({(total_chunks-i-1)/rate:.0f}s left)", flush=True)

    el = time.time() - started
    print(f"done: {written:,} unique reads, {dupes:,} duplicates dropped "
          f"({100*dupes/(written+dupes):.1f}%), {el:.0f}s")


if __name__ == "__main__":
    main()
