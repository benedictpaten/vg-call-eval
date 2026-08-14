# Benchmarking PanGenie on chr6: what it needs, and what blocks it here

PanGenie (Ebler et al. 2022) is the natural comparator for this caller. Its structure is the one the
linkage layer follows — hidden states are ordered pairs of panel haplotypes, transitions are
Li–Stephens, genotypes come from posteriors over the states implying them — with one difference that
is the whole point of comparing: **its emission is a k-mer model over raw reads, ours is
`ln P(reads | genotype)` from alignments.** A head-to-head on the same graph, sample and truth set
would isolate exactly that substitution.

This note records what running it requires and why it did not run here, so the next attempt starts
from the blockers rather than rediscovering them.

## What PanGenie needs

Two steps, or one combined ([user guide](https://pangenie.readthedocs.io/en/latest/usage.html)):

```
PanGenie-index -v <bubbles.vcf> -r <reference.fa> -t <threads> -o <prefix>
PanGenie -f <prefix> -i <reads.fa/fq> -s <sample> -j <threads> -t <threads> -o <prefix>
```

Constraints that matter, and each of which bit here:

- **Reads must be raw FASTA/FASTQ** (or a Jellyfish `.jf`), **uncompressed**. PanGenie is
  alignment-free; it counts k-mers in the read set.
- **The reference FASTA and panel VCF must be uncompressed.**
- **The panel VCF must be fully phased**, with each sample phased in a *single block from start to
  end*, and must contain **no overlapping variants** — overlapping variation has to be collapsed
  into one multi-allelic record. Minigraph-Cactus VCFs produced with `--vcf` are already
  `vcfbub`-filtered and satisfy this.
- **Diploid only.**
- **Designed for whole-genome genotyping**, explicitly not a restricted region.

## Why it did not run here

### 1. It does not build on this machine (arm64 macOS)

Not a flag problem. The layers disagree about C++ ABI:

- PanGenie sets `CMAKE_CXX_STANDARD 20`. Homebrew's jellyfish 2.3.1 headers use
  `std::get_temporary_buffer`, which modern libc++ has removed outright — the
  `_LIBCPP_ENABLE_CXX20_REMOVED_TEMPORARY_BUFFER` escape hatch no longer works — and contain a CRTP
  `static_cast` clang rejects. Dropping to `-DCMAKE_CXX_STANDARD=17` changes nothing.
- Compiling instead with Homebrew **GCC 15 and `-fpermissive` succeeds**, and then the *link* fails:
  Homebrew's `libjellyfish-2.0` was built with clang/libc++ and exports libc++-mangled symbols,
  while GCC objects want `std::__cxx11::` libstdc++ ones. Unresolved symbols include
  `jellyfish::RectangularBinaryMatrix::pseudo_inverse`, `jellyfish::file_header::matrix`,
  `Json::Value::operator[]` and the zlib `gz*` family.

The way through is to rebuild **jellyfish itself from source with GCC 15**, repoint `pkg-config`, and
build PanGenie against that. Tractable, untried, and worth an hour only once the data blocker below
is solved. The supported path is Linux x86_64, where Singularity and conda images exist.

Two smaller build fixes are already found and are needed regardless:
`find_package(cereal ...)` is commented out in `CMakeLists.txt`, so `/opt/homebrew/include` must be
added to `include_directories`; and CMake needs
`-DCMAKE_POLICY_VERSION_MINIMUM=3.5` with current CMake.

### 2. There are no raw reads — and this is the decisive one

The harness holds **aligned GAF and a GBZ, not FASTQ**. The GAF has no sequence column; read
sequence is recoverable in principle from the path plus the `cs:Z:` string, and base qualities sit in
`bq:Z:`, so a GAF → GAM → FASTQ reconstruction is possible with `vg`.

**It should not be used for this benchmark.** The chr6 GAF contains only reads that `vg giraffe`
placed on chr6. Feeding PanGenie a read set pre-filtered by another aligner removes precisely the
problem its k-mer model exists to handle, and does so in a direction that flatters us: reads that
failed to map, or mapped elsewhere, carry k-mers PanGenie would otherwise weigh, and their absence
can only make its genotypes worse. A comparison built that way would report our caller ahead partly
by construction, which is worse than no comparison.

Doing it honestly needs the original HG002 short reads (GIAB 30x Illumina, public). They are not
chromosome-partitioned, so the whole set is required even for a chr6 panel — which is also what
PanGenie's whole-genome design wants.

### 3. Panel VCF preparation is unfinished

`vg deconstruct -a` (used for the offline linkage work) emits **nested, overlapping** records, which
violates the no-overlapping-variants requirement. The panel would need deconstructing without `-a`,
then `vcfbub` — a Rust tool, also not installed here. `vg deconstruct` also emits `.` genotypes where
a haplotype does not traverse a site, which has to be reconciled with "fully phased, one block per
sample".

## What a fair benchmark would look like

1. Linux x86_64; PanGenie via Singularity or conda rather than a source build.
2. Original HG002 FASTQs, whole read set, uncompressed.
3. Panel: `vg deconstruct` (no `-a`) on the chr6 graph → `vcfbub` → resolve missing genotypes →
   uncompressed VCF; reference FASTA extracted from the graph, uncompressed.
4. `PanGenie-index` then `PanGenie -s HG002`.
5. Score with the **same** aardvark and truvari pipeline against the **same** T2T-Q100 truth and
   confident regions, so the numbers sit directly beside the existing arms.
6. State the caveat that a chr6 panel is off-label for a tool designed for whole-genome input, and
   that PanGenie is a *re-genotyper*: it genotypes the variant set it is given, so it should be
   compared against our `-z`/panel-enumeration arms, which are also restricted to panel alleles,
   rather than against support enumeration.

That last point is worth stating plainly: PanGenie and `readlik` (panel enumeration) answer the same
question — genotype this panel's alleles in this sample — from different evidence. That is the
comparison worth having.
