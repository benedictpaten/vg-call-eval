# `vg call` refactor — document set

Read-level genotype likelihoods for `vg call`: an explicit `P(reads | genotype)` model offered
alongside the existing aggregate-depth Poisson caller. Built (PR vgteam/vg#4990) and evaluated on
real HG002 data over two chromosomes and two graphs.

## Start here

| Document | What it is | Read it when |
|---|---|---|
| [vg-read-likelihood-design.md](vg-read-likelihood-design.md) | **The authoritative description of what exists** — model, architecture, scoring, read retrieval, output fields, correctness requirements | you need to understand or change the caller |
| [vg-call-eval-plan.md](vg-call-eval-plan.md) | The evaluation harness, the full investigation log, and the forward plan | you need a number or the derivation behind one (**start at "Where this stands"**; Appendix A is the log, with an index), or you want to know what is left to do (**Appendix C**) |
| [vg-call-characterization.md](vg-call-characterization.md) | How `vg call` works today, independent of this change | you are new to `vg call` |
| [read-likelihood-genotyping-plan.md](read-likelihood-genotyping-plan.md) | Source reading, prior art in the tree, and the reasoning behind each settled decision | you want to know *why* a decision went the way it did |

Outbound drafts, standalone:

- [gbz-base-c-api-request.md](gbz-base-c-api-request.md) — the reads-only-query ask for gbz-base upstream
- [vg-issue-draft-gaf-trailing-node.md](vg-issue-draft-gaf-trailing-node.md) — `vg convert -G` drops a trailing zero-length path node

And one prototype:

- [call-per-contig.py](call-per-contig.py) — calls one contig at a time so only one chromosome is
  ever resident. Verified byte-identical to a whole-graph run, and **not used for tier 2**: on the
  HPRC v2.1 MC CHM13 graph `gbz-base --snarls` found chains for 2 of 46 components, and that graph
  fits in memory anyway. Kept because the question it answers — how much new vg code per-chromosome
  calling needs — has a useful answer: none.

These documents live **in the harness repo** so they sit alongside the evaluations they explain.
Generated results are one directory over, in [`../docs/`](../docs/) — regenerated from run
artefacts rather than transcribed, so when a number here and a number there disagree, that one is
right.

| Generated | |
|---|---|
| [../docs/tier2-chr20-results.md](../docs/tier2-chr20-results.md), [../docs/tier2-chr6-results.md](../docs/tier2-chr6-results.md) | full five-arm accuracy tables per chromosome |
| [../docs/tier2-chr20-hap32.md](../docs/tier2-chr20-hap32.md), [../docs/tier2-chr6-hap32.md](../docs/tier2-chr6-hap32.md) | 4- vs 34-haplotype graph |
| [../docs/tier2-quality-signals.md](../docs/tier2-quality-signals.md) | how calls are ranked, and the filters that did not help |

**Source links of the form `vg/src/...` resolve against a `vg` checkout, not this repo.** They are
line-anchored citations into the caller, kept because the reasoning depends on them; the line
numbers are accurate as of PR vgteam/vg#4990 and will drift.

## The one-paragraph version

The model builds a reads × alleles matrix of `ln P(read | allele)` from the graph-implied alignment
(no fresh dynamic programming), normalises each row by the read's best allele, and scores every
genotype as a mixture over the haplotypes in it plus a MAPQ-derived "this read came from somewhere
else" term. On real data it beats the Poisson caller on genotype F1 on every variant class of every
dataset tested, and its margin *widens* as the graph gets richer — 4 haplotypes to 34 — because more
haplotypes offer more true alleles and more wrong ones, and only a caller that scores alleles
individually can tell them apart. Structural variants are the weak class for both callers. Cost is
now near parity after the read path was made 5.8× faster with byte-identical output.

## Conventions in these documents

- **Claims that were later withdrawn are kept, not deleted**, and gathered in
  [vg-call-eval-plan.md](vg-call-eval-plan.md) Appendix B. Several of the useful results here came
  from noticing that an earlier confident claim was wrong, and deleting the wrong claim would delete
  the reason the right one was looked for.
- **Section numbers are stable.** Code comments and commit messages cite them (`plan §9.20`,
  `design §5.3.3`), so sections are relabelled or marked superseded rather than renumbered.
- **Numbers in the harness repo's `docs/` are generated; numbers here are transcribed.** When they
  disagree, the generated ones are right.
