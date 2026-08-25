# Request to gbz-base upstream: a reads-only query

**Status: REWRITTEN AND READY, NOT SENT.** Intended as a GitHub issue on
<https://github.com/jltsiren/gbz-base>, or as the basis for a message to Jouni. Nothing has been
posted. The ask below supersedes the C-interface draft kept at the bottom of this file, which was
written before we built the subprocess backend and measured it.

---

**Title:** A reads-only query mode for `gbz-base query`

Hi Jouni,

`vg call` now has a read-level genotyping model -- an explicit reads x alleles likelihood matrix per
snarl, scored as `P(reads | genotype)` rather than the aggregate-support Poisson model it used to
use. That needs random access to the reads overlapping a snarl, and GAF-Base is what we use for it,
through `gbz-base query`. It works well and we are not blocked on anything. This is a request for one
flag, plus a small robustness note.

**First, three things worth saying because they are good.** Run at scale on HG002, 596,017,764
alignments against a 100 M-node HPRC v2.1 MC CHM13 graph:

- **The build performs very well.** `gaf-base sort` did 596 M records in 47 min at 6.5 GB peak,
  single-threaded, and `gaf-base construct` in 51 min at 15 GB, running concurrently through a FIFO,
  so 51 min wall in total. The database is 20.74 GB -- **36.5 bytes per alignment**, about 3x smaller
  than the equivalent sorted GAM. For comparison `vg gamsort` needed 249 s on 8 M records where
  `gaf-base sort` needed 17 s.
- **Queries are fast.** Against that 22 GB database, 0.04 s for a 380-node query and 0.22 s for a
  4,000-node one, near-linear in reads returned. Nothing here needs apologising for.
- **The graph-consistency check is enforced, and that is the right call.** Querying the whole-genome
  GAF-Base with a GBZ-Base built from one extracted chromosome fails with `The graph is not a valid
  reference for the alignments`. That is exactly the mismatch we expected to have to police
  ourselves, and we would rather you kept policing it.

**The ask: a query that returns reads without building a subgraph.**

We only ever want the reads. `gbz-base query` is a subgraph tool that can also return the overlapping
alignments, so every query we make builds a GFA subgraph and we send it to `/dev/null`. Timing one
realistic query -- 1,024 nodes, 5,512 reads returned, on a small graph so the constant costs are
visible:

| component | cost |
|---|---|
| `fork`/`exec` of the subprocess | 3 ms |
| subgraph extraction, output to `/dev/null` | **20 ms** |
| read decode | 41 ms |

So about a third of every query is work we discard. Over a chromosome that is 1,446 queries for
chr20 alone and roughly 24,000 for a genome.

`--reference-only` looked like an existing route to this and is not: on a node-based query it exits
101 with `Reference-only output is not supported for node-based queries` and writes nothing. Either
supporting it for `-n` queries, or a dedicated reads-only mode, would give us what we want.

Whether it arrives as a CLI flag, a Rust API or an `extern "C"` symbol matters much less to us than
that it exists -- and a **CLI flag is the cheapest thing for you to provide and for us to adopt**,
since we already parse the GAF you write. We would then be staying on the binary by choice rather
than as a compromise, which is the position we would prefer to be in anyway.

**A second, smaller request: skip unknown node IDs instead of failing the query.**

`gbz-base query` errors on a node ID that is not in the graph rather than ignoring it. With a sparse
ID space -- 100 M nodes across a 224 M ID range here -- a caller has to pre-filter every query
against the graph, and one stray ID aborts the whole thing. A flag to skip unknown nodes would be a
small and genuinely useful addition.

**Why we are not reading the SQLite ourselves.** Not a complaint, just so the ask makes sense. The
integer codecs are already byte-identical to what vg links today -- `ByteCode`
(`MASK 0x7F / FLAG 0x80 / SHIFT 7`) and `RLE` (`value + sigma*(len-1)`) match `gbwt`'s `internal.h`
and `Run::encodeBasic` exactly, and `Nodes.edges || Nodes.bwt` parses directly as a
`gbwt::CompressedRecord`. But decoding an `AlignmentBlock` means reproducing the conditional column
layout, the LF-walk through `Nodes` to recover target paths, and rANS 4x16 for quality strings. The
formats are documented as able to change without warning and are version-checked exactly, so a
hand-maintained reader is something we would rather not own. Going through a supported interface is
what we would choose even if the format were frozen.

**One question that stands on its own:** is the on-disk schema intended to be third-party-readable at
any point, or should we assume it stays an internal detail? A clear "internal only" is a perfectly
good answer -- it just tells us to go through a supported interface, which is what we want regardless.

**Two notes, not requests.**

- **MAPQ and base qualities both survive**, so this question does not need asking: the GAF that comes
  back carries MAPQ in column 12 and per-base qualities as `bq:Z:`, plus a `cs:Z:` difference string.
  The quality-aware scoring path is the normal case for us and the no-quality fallback stays an edge
  case.
- **The `--alignments` default is a trap for a consumer like us.** It defaults to `clipped`, which can
  return one read as several fragments; anything computing per-read statistics needs `overlapping`.
  On a 9-node query: `overlapping` -> 290 records and 290 distinct names, `clipped` -> 289,
  `contained` -> 0. We pass `--alignments overlapping` everywhere for this reason.

Thanks -- and thanks for gaf-base generally. The snarl-boundary query is a really good fit for what
variant callers actually need, and it is the reason this model is practical at all.

---

What follows is the original draft, kept for the parts still worth reusing.

---

**Title:** C interface for querying GAF-Base from C++ (for `vg call`)

Hi Jouni,

We're working on a read-level genotyping model for `vg call` — building an explicit
reads × alleles likelihood matrix per snarl and scoring genotypes as `P(reads | genotype)`,
rather than the current aggregate-support Poisson model. That means `vg call` needs random
access to the read alignments overlapping a given snarl, which it has never needed before
(today it consumes only a `vg pack` coverage index).

GAF-Base is by some distance the best fit for this of the options we looked at. The reason is
your query surface rather than the storage: you already support extracting a subgraph between
two snarl boundary nodes and returning the reads overlapping it, with `Overlapping` /
`Contained` / `Clipped` modes. That is *exactly* the access pattern — a snarl is defined by its
two boundary nodes. By comparison vg's own GAM index (`StreamIndex`/`.gai`) is a node-ID range
scan that over-fetches and needs filtering, and tabix'd GAF is more awkward still.

The obstacle is that vg is C++ and there's currently no C-callable interface. We looked at
reading the SQLite directly and concluded we shouldn't:

- Encouragingly, the integer codecs are already byte-identical to what vg links today —
  `ByteCode` (`MASK 0x7F / FLAG 0x80 / SHIFT 7`) and `RLE` (`value + sigma·(len−1)`) match
  `gbwt`'s `internal.h` / `Run::encodeBasic` exactly, and `Nodes.edges || Nodes.bwt` is
  directly parseable as a `gbwt::CompressedRecord`. zstd we already have.
- But decoding an `AlignmentBlock` means reproducing the conditional column layout (which
  varints appear in `numbers` depending on per-alignment flag bits, block `read_length`, and
  whether a difference string was stored), plus the LF-walk through `Nodes` to recover target
  paths, plus rANS 4x16 for quality strings — and then tracking all of it across releases.
  Given that the formats are explicitly documented as able to change without warning, with an
  exact-match version check and (by our count) three format-breaking changes in about four
  months, a hand-maintained C++ reader would break more often than it works.

So rather than duplicate your decoder, we'd like to ask whether you'd consider exposing a small
`extern "C"` surface over the machinery that already exists. Everything needed appears to be
there in `ReadSet` — what's missing is only a flat wrapper. Something on this order would be
enough for us:

```c
// open / close (gbz path needed for node sequences in the default reference-based mode)
gafbase_t*  gafbase_open(const char* gaf_db_path, const char* gbz_db_path, char** error_out);
void        gafbase_close(gafbase_t*);

// the one we actually care about: reads overlapping the subgraph between a snarl's
// two boundary nodes; overlap_mode = Overlapping | Contained | Clipped
gafbase_reads_t* gafbase_query_snarl(gafbase_t*, int64_t start_node, int64_t end_node,
                                     int overlap_mode, char** error_out);

// useful to have as well, if it's cheap: reads overlapping an explicit set of node ids
gafbase_reads_t* gafbase_query_nodes(gafbase_t*, const int64_t* node_ids, size_t n,
                                     int overlap_mode, char** error_out);
void             gafbase_reads_free(gafbase_reads_t*);

// results as GAF text
size_t      gafbase_reads_count(const gafbase_reads_t*);
const char* gafbase_reads_line(const gafbase_reads_t*, size_t i);

const char* gafbase_schema_version(gafbase_t*);
```

Three deliberate choices there, all aimed at keeping this cheap for you to maintain:

0. **The snarl query is the one that matters**, not the node-id list — the boundary-node form is the
   whole reason gaf-base is the right fit for us, and a node-id list would just be the GAM-index
   access pattern with extra steps. If only one of the two is worth wrapping, it's the first.
1. **Return GAF text, not a parsed struct.** No struct layout to version, so the ABI stays
   stable across any internal format change — which is the whole point for us. And it costs us
   nothing, because vg already has GAF text → `Alignment` (`gafkluge::parse_gaf_record` →
   `vg::io::gaf_to_alignment`). It should also be close to what `gaf-base decompress` already
   produces.
2. **Nothing about the graph.** vg already holds the graph in memory, so we only need the
   alignment side.

Questions, in case any of this changes your view:

- Would you be open to this in principle? If the shape is wrong, we're happy to adapt to
  whatever you'd prefer to support.
- **Thread safety:** vg queries from multiple OpenMP threads. Would a single handle be safe for
  concurrent queries, or should we open one per thread?
- **MAPQ and base qualities:** does the GAF that comes back out carry both? Our model needs MAPQ to
  weight how much a read is trusted at a site, and per-base qualities to charge each mismatch its own
  error probability. We see that `--no-quality` databases exist — is that the common case in practice?
  It's not a blocker either way; it just changes which scoring path we treat as the default.
- **Build/distribution:** this is the part we're least sure about on our end. A Rust
  `staticlib` would put the Rust toolchain into every vg build environment including release
  CI, and a `cdylib` would break vg's fully-static build target. Do you have a view on whether
  prebuilt static archives are plausible, or should vg treat gaf-base support as an optional
  compile-time feature?
- Separately, and useful regardless of the above: **is the on-disk schema intended to be
  third-party-readable at any point**, or should we assume it stays an internal detail? A clear
  "internal only" is a perfectly good answer — it just tells us to go through a supported
  interface, which is what we'd prefer anyway.

For what it's worth, none of our work is blocked on this. The read source sits behind an
interface, with a backend that shells out to `gbz-base query` — not `gaf-base`, which has no
random-access query (`construct`/`decompress`/`sort` only) — so your binary does the decoding
and format changes cost us nothing. Swapping that for a library call later touches one function.

**Two of our own questions below turned out to be answerable by building it:**

- **MAPQ and base qualities both survive.** The GAF that comes back carries MAPQ in column 12
  and per-base qualities as `bq:Z:`, plus a `cs:Z:` difference string. So the quality-aware
  scoring path is the normal case and the no-quality fallback stays an edge case. No change
  needed; noting it so the question is not asked twice.
- **The `--alignments` default is a trap for a consumer like us.** It defaults to `clipped`,
  which can return one read as several fragments; anything computing per-read statistics needs
  `overlapping`. On a 9-node query: `overlapping` → 290 records / 290 distinct names,
  `clipped` → 289, `contained` → 0. Not a request, just a note.

Thanks — and thanks for gaf-base generally; the snarl-boundary query is a really good fit for
what variant callers actually need.
