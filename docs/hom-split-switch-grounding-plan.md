# Grounding hom splitting in switch events

`--anchors-hom-split` splits 87,456 of chr20's 92,935 homozygous sites (94.1%). The gate
(`|lo| >= --split-min-q 0.5` on each read, then `side0 >= 2 && side1 >= 2`) is a crash detector, not
a quality filter: at 40x with reads from both haplotypes, both sides clear 2 with enormous margin and
98.4% of reads clear 0.5 nats. The question is not "is 94.1% too many" in the abstract -- it is
**whether the splits bridge the phasing's switch events**, and that is what this measures.

## The harm, stated precisely

A hom split labels its two slots by `sign(lo)`, the read's cross-site strand against the settled
chain. The chain is internally consistent everywhere, so the split always *faithfully follows the
chain*. It cannot invent a switch. What it can do is **carry the chain's existing switch across a
junction the anchor graph would otherwise have refused to link.**

Where the phasing switches, het density is usually low -- that is often why it switched. Without hom
splits the graph has no haplotype-carrying node through such a region: the two haplotype paths
collapse into a shared node and separate again, and the graph honestly declines to say which side
joins which. With hom splits it states a confident, wrong pairing. **That is the whole cost, and it
is binary per switch: bridged or not.**

So `--split-min-q` and a coverage-scaled `--split-min-side` are worth testing not because they raise
average placement accuracy -- they cannot, see below -- but because collapsing a site *restores the
break*.

## Why the previous sweep found nothing

`split-min-q-cannot-reduce-switches.md` swept `--split-min-q` 0.5 -> 8 and reported agreement
identical to four decimals at every threshold (2,939,707/3,117,403 on every row). That number is the
**het-site self-check**, which does not depend on `--split-min-q` at all -- a constant, not a null
result. The parameter did move the output: 87,456 -> 74,176 sites split. It was simply never
measured. This plan measures the thing that moves.

The doc's other finding stands and bounds the scope: **`--split-min-q` never reaches placement.** It
gates which sites may split; every read at a site that does split is placed by `sign(lo)` with no
threshold. So both levers are *site* gates, and their only action is collapse.

## Run budget: 3 whole-chromosome runs

The gate is anchor-emission-only -- it cannot touch the VCF, so F1 and whatshap switch counts are
invariant across the entire sweep, and every arm shares one phasing. Combined with `lo` being one
number per read, **every (`--split-min-q`, `--split-min-side`) combination is reconstructable offline
from a single run**:

- an unsplit hom site emits one slot holding *every* pinned read, and a split one emits the same
  reads across two slots, so the pinned read set is in the file either way, whatever the gate did;
- `lo` is constant per read across all hom sites;
- `reliability` is computed over the site's one distinct allele, not over the two slots, so it is
  invariant to splitting and the reconstruction is faithful, not approximate.

| run | what | why it is needed |
|---|---|---|
| **A** chr20, shipping defaults, instrumented | `--anchors-out --anchors-hom-split` + a per-read `lo` dump | baseline *and* the substrate for the entire offline sweep |
| **B** chr20, shipping binary, chosen setting | the implemented coverage-scaled min-side | proves the offline reconstruction is exact, and that the change is anchors-only (VCF byte-identical to A) |
| **C** chr6, instrumented, defaults | hold-out | both arms reconstructed offline, through the pipeline B validated |

Run C does not wait on the chr20 result, because it too is run at defaults and swept offline. The
instrumentation is an env-gated dump of `name<TAB>lo`, ~79k lines, deleted before the PR is touched;
B is what proves it changed no behaviour.

## Measurements

**Primary -- bridged switches.** whatshap `--switch-error-bed` gives chr20's 27 true switches (the
27/76 in `all_switchflips`; flips are local and the standing instruction is that they are minor).
For each, ask whether the anchor graph carries a haplotype-resolved path continuously across the
junction. Hom sites have no coordinate for 60.6% of splits (no VCF line), so position comes from
node-id interpolation between the positioned sites, which are dense enough -- one VCF record per
~570 bp against a ~15 kb read -- for the bracketing to be unambiguous.

**Cost side -- run N50.** Haplotype-resolved run length over the compacted graph, via
`anchor_graph_gfa.py`'s unitig machinery. A setting that unbridges switches by collapsing the graph
everywhere has bought nothing.

**Secondary, higher-powered -- enrichment.** ~1.2% of splits lie within a read length of a switch
(27 switches x ~40 splits per +/-15 kb window, against 87,456). Does the set a stricter gate
collapses over-represent switch neighbourhoods? At 13,280 collapsed, chance gives ~160; a 2x
enrichment gives ~320 against a Poisson sd of 13. Good power, and it is the test of *mechanism*: if
the collapsed set is not enriched, the lever is removing low-coverage sites, not switch-bridging
ones, and should be described that way.

**A stated prediction, so the result can falsify it.** A coverage-scaled `--split-min-side` may well
be orthogonal to switches: at a switch-spanning site the reads still divide ~50/50 by `sign(lo)` --
both sides large -- so the gate sees exactly what it sees at a clean site. The mechanism by which
either lever could bite is that reads *spanning* a switch accumulate contradictory evidence and so
carry depressed `|lo|`. **Whether `|lo|` is depressed near truth switches is therefore the single
diagnostic that decides whether either parameter is the right lever at all**, and it is one cheap
plot off run A. If it is flat, both levers are noise filters and tuning them against 27 events would
be fitting noise -- which is the answer, not a failure.

## The two levers

**A. coverage-scaled `--split-min-side`** -- `max(2, ceil(f * n_confident))`, floor of 2 as
specified, `f` in {0.10, 0.15, 0.20, 0.25, 0.33}. Bites on lopsided sites -- 30 reads one way and 2
the other -- which is the "relabelled, not partitioned" case the existing code comment names but the
fixed 2 does not catch.

**B. `--split-min-q`** in {0.5, 1, 2, 4, 8}, swept jointly with A over the same offline substrate.

Not pursued: a **placement**-side confidence bar. The banded table shows sub-threshold reads are
71.7% correct against a 50% coin, so coining them is worse, and dropping them produced the
deletion-shaped holes the `lo == 0` coin exists to avoid. The defensible version is a third "weak
claim" state the format has no room for.

## Gates

- Run B's VCF byte-identical to run A's: the change is anchors-only, as claimed.
- Run B's anchors file identical to the offline reconstruction at the same setting: the sweep is
  exact, not approximate.
- chr6 never fits a parameter. The setting is chosen on chr20 and reported on chr6.
- If the `|lo|`-near-switch diagnostic is flat, stop and report that, rather than tuning to 27 events.
