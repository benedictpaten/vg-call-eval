#!/usr/bin/env bash
# Re-run the read-likelihood arms at the new default floor (0.01). poisson and
# poisson-z are unaffected -- they do not use the read-likelihood model at all --
# so their existing numbers stay valid and are not recomputed.
set -o pipefail
cd "$(dirname "$0")/../.."
export PATH="/private/tmp/claude-501/-Users-benedictpaten-My-Drive-papers-and-projects-2026-vg-planning/65260683-9f63-4117-b57b-52dccfd3b871/scratchpad/gbztools/bin:$PATH"
python3 work/tier2-chr20/run_arms.py --only readlik readlik-nomismap readlik-z --threads 6 --read-window 1024
mv work/tier2-chr20/results/arms.json work/tier2-chr20/results/arms.floor-0.01.json
echo ALLDONE
