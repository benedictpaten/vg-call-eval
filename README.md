# vg-call-eval

Concordance and performance testing for `vg call`, kept **outside the vg tree** so that none of its
dependencies (aardvark, truvari, a Python analysis stack) land in vg's build.

It exists to answer one question that vg's own test suite cannot: **is a change to `vg call` actually
more accurate?** vg's in-tree harness has exactly one truth-based concordance assertion, and it covers
the `-v` re-genotyping path — the default de novo path has never been measured against truth.

Implements stages 3b, 4 and 4b of the read-likelihood design.

## The one thing to read before quoting a number

**Tier-0 numbers are optimistic and are not absolute performance.** Reads are simulated from the graph
and mapped back to that same graph, so mapping is unrealistically easy. Tier 0 exists to compare
callers *to each other* and to calibrate `GQ`. Absolute numbers need tier 2 (real data), which is
deferred. See [docs/simulation.md](docs/simulation.md), which also documents exactly how the
simulation works and what it cannot tell you.

A worked example of why this matters: at 20 kb and 20x, **every caller scores F1 = 1.0000**. The task
is simply too easy to discriminate. At 4x with 100 bp reads the same harness separates them. If a
configuration gives everything a perfect score, that is a statement about the configuration, not the
callers.

## Sanity controls are not optional

A comparison harness that is subtly wrong is worse than none, because it produces confident numbers.
Two controls gate everything, and they run as tests:

```bash
PYTHONPATH=src python3 -m pytest tests/ -q
```

* **Positive** — identical inputs must score 1.0.
* **Negative** — deliberately dropping and corrupting calls must be detected. If this passes at 1.0,
  the harness cannot see errors and every number it has produced is meaningless.

## Install

Needs `vg`, `bcftools`, `samtools`/`tabix`/`bgzip`, and `aardvark`. See [docs/install.md](docs/install.md)
— note aardvark ships **only an x86_64 Linux binary**, so on macOS/ARM it must be built from source
(it is pure Rust and builds cleanly).

## Run

```bash
PYTHONPATH=src python3 -m vgcalleval.cli run \
    --out /tmp/eval --vg /path/to/vg --ref-length 60000 --depth 4 --read-length 100
```

Add `--vg-depthfix /path/to/patched/vg` to include the `poisson-depthfix` arm. The harness takes a
**binary path per arm**, so two vg builds can be compared in a single matrix — which is how the
`depth_err` bug's effect on the baseline gets quantified.

### Arms

| Arm | Purpose |
|---|---|
| `poisson` | the current default, as shipped |
| `poisson-depthfix` | with the `depth_err` one-liner patched. **Verified byte-identical output** - that bug is inert, since its only consumer in the likelihood is commented out. Retained as a control: if this arm ever diverges from `poisson`, someone has re-enabled the depth term. |
| `readlik` | the read-level likelihood caller |
| `readlik-nomismap` | `--no-mismap-term`, to measure what the mismapping term contributes |
| `readlik-gbwt-nopack` | `-z` haplotype enumeration with no pack file |

## Status

Working: tier-0 simulation, the caller matrix, aardvark comparison, sanity controls, per-arm timing.

Not yet built: `GQ`-sweep PR curves, the `read_weight` calibration fit (stage 4b), the truvari SV
cross-check, tier 1 (vg's HGSVC fixture) and tier 2 (real data). Tier 2 is additionally blocked on
`vg call`'s read source becoming practical at scale — the current in-memory backend loads the whole
GAM.
