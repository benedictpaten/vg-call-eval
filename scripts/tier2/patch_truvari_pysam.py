#!/usr/bin/env python3
"""Work around a pysam 0.24.0 defect that makes `truvari refine` fail on this machine.

Symptom. `truvari refine` dies in `phab.set_regions` with

    SamtoolsError: 'samtools returned with error 1:
                    stdout=, stderr=[faidx] Could not build fai index ....fa.fai'

Cause. `phab.set_regions` extracts the refine regions with

    fout.write(samtools.faidx(self.reference_fn, "-r", regions_file_name))

relying on pysam's dispatcher to *return* samtools' stdout. Under pysam 0.24.0 here
that capture is broken for the samtools dispatcher specifically: `faidx` returns an
empty string and writes its output to the process's real stdout instead. phab therefore
writes a zero-byte reference fasta and the very next line, which indexes it, fails.

The failure is doubly misleading. It names the *index* as the problem when the fasta is
what is empty, and the sequence that should have gone into the file appears in the
terminal, so a run looks like it produced enormous output and then failed for an
unrelated reason.

Scope. `pysam.bcftools` capture is unaffected (verified: `bcftools view -h` returns
1,608 characters), and only these two call sites need it, so the patch is two lines.

This is a patch to a third-party package inside `work/truvari-venv`, which is a local
analysis artefact rather than a checked-in dependency. It lives here as a script rather
than as a hand-edit so that rebuilding the venv and re-running this reproduces the
environment exactly. Idempotent, and `--revert` restores the original.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
VENV = REPO / "work/truvari-venv"

ORIGINAL = '''        # Pull sequences
        out_fn = truvari.make_temp_filename(suffix='.fa')
        with open(out_fn, 'w') as fout:
            fout.write(samtools.faidx(
                self.reference_fn, "-r", regions_file_name))
        # Facilitate fetching
        samtools.faidx(out_fn)'''

PATCHED = '''        # Pull sequences
        # PATCHED by vg-call-eval scripts/tier2/patch_truvari_pysam.py:
        # pysam 0.24.0's samtools dispatcher does not capture faidx stdout here -- it
        # returns "" and leaks the sequence to the real stdout -- so this wrote a
        # zero-byte fasta and the index step below failed. Shell out instead.
        import subprocess as _sp
        out_fn = truvari.make_temp_filename(suffix='.fa')
        with open(out_fn, 'w') as fout:
            _sp.run(["samtools", "faidx", str(self.reference_fn),
                     "-r", str(regions_file_name)], stdout=fout, check=True)
        # Facilitate fetching
        _sp.run(["samtools", "faidx", str(out_fn)], check=True)'''


def find_phab() -> Path:
    hits = sorted(VENV.glob("lib/python*/site-packages/truvari/phab.py"))
    if not hits:
        sys.exit(f"no truvari/phab.py under {VENV}")
    return hits[0]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--revert", action="store_true")
    args = ap.parse_args()

    phab = find_phab()
    backup = phab.with_suffix(".py.orig")
    text = phab.read_text()

    if args.revert:
        if backup.exists():
            shutil.copy(backup, phab)
            print(f"reverted {phab}")
        else:
            print("no backup to revert to")
        return

    if PATCHED in text:
        print(f"already patched: {phab}")
        return
    if ORIGINAL not in text:
        sys.exit(f"{phab} does not contain the expected block; truvari version changed?"
                 " Inspect set_regions by hand before patching.")

    if not backup.exists():
        shutil.copy(phab, backup)
    phab.write_text(text.replace(ORIGINAL, PATCHED))
    print(f"patched {phab} (backup at {backup.name})")


if __name__ == "__main__":
    main()
