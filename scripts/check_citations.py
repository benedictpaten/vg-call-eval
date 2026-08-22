#!/usr/bin/env python3
"""Check that the documentation still matches the vg tree it cites.

Stage 12 of planning/decide-then-render.md asks for this to be a re-runnable script rather than a
one-time read-through, because the failure it guards against is silent: a section survives a deletion
it described, and nothing about the prose looks wrong.

Two checks, deliberately scoped differently.

1. `src/foo.cpp:123` citations resolve -- the file exists and is at least that long. Run over every
   markdown tree given. This catches a citation that has fallen off the end of a file or into a file
   that no longer exists. It does NOT catch a citation that still resolves but now points somewhere
   else: line numbers drift under any edit above them, and no cheap check can tell.

2. Backticked code symbols exist somewhere in vg's sources. Run over VG'S OWN doc/ ONLY, because that
   is the documentation that ships with the code and has to describe the code as it is. The eval
   repo's docs and plans are a record of what was done and when, so they cite deleted things on
   purpose -- `apply_linkage_change` and the `nested_*` FILTERs appear there precisely to say they are
   gone. Gating on those would either fail forever or force the history to be rewritten.

Both checks are coarse nets, and it is worth knowing where the holes are rather than trusting a green
result too far:

  - Check 1 cannot tell a citation that drifted from one that did not. Only that it still lands inside
    the file.
  - Check 2 greps the raw source text, comments included, so a symbol that survives only in a comment
    about its own deletion counts as present. Verified: run against the doc as it stood before stage 12
    it flags `nested_haploid` and `nested_unreachable` but not `nested_diploid`, because that one name
    is still written in a unit-test comment saying what the FILTER used to label.

Exit status 1 if any check fails, so it can gate.
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

CITE = re.compile(r"(src/[A-Za-z0-9_./-]+\.(?:cpp|hpp)):(\d+)(?:-(\d+))?")
# A backticked identifier that looks like code rather than prose: it has an underscore or a scope
# operator. `GQN` and `PASS` are neither, and are not symbols to look up.
SYMBOL = re.compile(r"`([A-Za-z_][A-Za-z0-9_]*(?:::[A-Za-z_][A-Za-z0-9_]*)*)`")

# Not vg symbols, and cited in vg's docs for other reasons: VCF/FORMAT keys, command-line and
# environment names, and field names of other tools' formats. Listed rather than pattern-matched so
# that adding one is a decision someone made.
NOT_SYMBOLS = {
    "read_likelihood", "min_confidence", "linkage_weight", "genotype_snarls", "ploidy_bed",
    "no_nested", "enumerate_support", "mosaic_out", "top_down", "gaf_base", "gbz_base",
    "freq_prior", "rho_min", "phase_set", "sample_name",
    # Quantities named in a derivation, not identifiers: the prose backticks them as terms.
    "mean_pair_mass",
}

# vg's wiki pages are upstream and are not this work's to keep current. They are full of environment
# variables and CMake commands, which are backticked, contain underscores, and are not vg symbols.
SKIP_DIRS = {"wiki"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vg", default=str(Path.home() / "CLionProjects/vg"))
    ap.add_argument("--paths", nargs="*", default=["docs", "planning"],
                    help="markdown trees to check citations in (check 1)")
    ap.add_argument("--vg-docs", default="doc",
                    help="vg's own doc dir, relative to --vg, for the symbol check (check 2)")
    ap.add_argument("--skip-symbols", action="store_true")
    args = ap.parse_args()

    vg = Path(args.vg)
    if not (vg / "src").is_dir():
        print(f"no vg checkout at {vg}", file=sys.stderr)
        return 2

    failures: list[str] = []

    # ---- check 1: citations resolve
    counts: dict[str, "int | None"] = {}

    def lines_in(rel: str):
        if rel not in counts:
            p = vg / rel
            counts[rel] = sum(1 for _ in p.open("rb")) if p.is_file() else None
        return counts[rel]

    cited = 0
    trees = [Path(p) for p in args.paths]
    vg_doc_dir = vg / args.vg_docs
    if vg_doc_dir.is_dir():
        trees.append(vg_doc_dir)
    for root in trees:
        for md in sorted(root.rglob("*.md")):
            for lineno, line in enumerate(md.read_text(errors="replace").splitlines(), 1):
                for m in CITE.finditer(line):
                    rel, start, end = m.group(1), int(m.group(2)), m.group(3)
                    cited += 1
                    n = lines_in(rel)
                    if n is None:
                        failures.append(f"{md}:{lineno}: {rel} does not exist")
                    elif start < 1 or (int(end) if end else start) > n:
                        failures.append(f"{md}:{lineno}: {rel}:{m.group(0).split(':',1)[1]}"
                                        f" is past end of file ({n} lines)")
    print(f"check 1: {cited} citations, {len(failures)} unresolved")

    # ---- check 2: symbols in vg's own docs still exist in vg's sources
    sym_failures: list[str] = []
    if not args.skip_symbols and vg_doc_dir.is_dir():
        blob = subprocess.run(
            ["bash", "-c",
             f"cat {vg}/src/*.cpp {vg}/src/*.hpp {vg}/src/subcommand/*.cpp "
             f"{vg}/src/unittest/*.cpp 2>/dev/null"],
            capture_output=True, text=True).stdout
        seen = 0
        for md in sorted(vg_doc_dir.rglob("*.md")):
            if SKIP_DIRS & set(md.relative_to(vg_doc_dir).parts):
                continue
            for lineno, line in enumerate(md.read_text(errors="replace").splitlines(), 1):
                for m in SYMBOL.finditer(line):
                    sym = m.group(1)
                    if ("_" not in sym and "::" not in sym) or sym in NOT_SYMBOLS:
                        continue
                    seen += 1
                    if sym.split("::")[-1] not in blob:
                        sym_failures.append(f"{md}:{lineno}: `{sym}` is not in vg's sources")
        print(f"check 2: {seen} symbols in {vg_doc_dir}, {len(sym_failures)} absent from src")
    failures += sym_failures

    for f in failures:
        print(f"  {f}")
    print(f"\n{len(failures)} failures")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
