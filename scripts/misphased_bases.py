#!/usr/bin/env python3
"""Mis-phased bases from a whatshap switch-error BED.

A run of k contiguous BED intervals changes parity iff k is ODD -- an even run goes out and comes
straight back. Counting any run of length != 2 as a switch (the obvious reading) miscounts every
run of 4, which is how chr20 came out at 28 events instead of 27.
"""
import sys
raw = sorted((int(l.split()[1]), int(l.split()[2])) for l in open(sys.argv[1])
             if not l.startswith(("track", "#")))
span = int(sys.argv[2])
runs, cur = [], []
for r in raw:
    if cur and r[0] == cur[-1][1]:
        cur.append(r)
    else:
        if cur:
            runs.append(cur)
        cur = [r]
if cur:
    runs.append(cur)
ch = [r[0][0] for r in runs if len(r) % 2 == 1]
edges = [0] + ch + [span]
wrong = sum(edges[i + 1] - edges[i] for i in range(len(edges) - 1) if i % 2 == 1)
wrong = min(wrong, span - wrong)
print(f"{len(ch)}\t{wrong}\t{100*wrong/span:.2f}\t{span/1e6/(len(ch)+1):.2f}")
