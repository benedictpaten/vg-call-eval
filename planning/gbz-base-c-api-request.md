# Draft: request to gbz-base upstream for a C interface

**Status: DRAFT, still not sent.** Intended as a GitHub issue on
<https://github.com/jltsiren/gbz-base>, or as the basis for an email/Slack message to Jouni.
Nothing has been posted.

**IMPORTANT: the ask below is now wrong, and needs rewriting before it is sent.** The subprocess backend
described at the end as an interim measure has since been *built and measured*
([vgteam/vg#4990](https://github.com/vgteam/vg/pull/4990)), and the measurement contradicts the premise of
the letter.

Timing one realistic query (1024 nodes, 5512 reads returned) on a 400 kb graph:

| component | cost |
|---|---|
| `fork`/`execvp` | **3 ms** |
| subgraph extraction, output to `/dev/null` | **20 ms** |
| read decode | 41 ms |
| total | 62 ms |

So **the subprocess is 5% of the cost and is not worth a letter.** The thing that is worth asking about is
the 20 ms: `gbz-base query` is a subgraph tool that can also return overlapping reads, and we only ever
want the reads, so a third of every query builds a GFA subgraph we send straight to `/dev/null`.

**Update after optimising our side (design §6.5).** With vg's own read path fixed — the whole-chromosome
run went 570 s to 99 s — the subprocess is even less of the remaining cost than the table above
suggests. Two further data points for the letter:

- `gbz-base query` is **fast**: 0.04 s for a 380-node query and 0.22 s for 4,000 nodes against the
  22 GB whole-genome GAF-Base, near-linear in reads returned. Nothing here needs apologising for.
- **`--reference-only` panics on node-based queries**: "Reference-only output is not supported for
  node-based queries", exit 101, no output written. It looked at first like an existing route to the
  reads-only query below, and it is not. Either supporting it for `-n` queries, or a dedicated
  reads-only mode, would give us what we want.

**The ask should therefore become: a reads-only query that does no subgraph construction.** Whether it
arrives as an `extern "C"` symbol, a Rust API, or simply a flag on the existing CLI matters far less than
that it exists — and a CLI flag would be much cheaper for upstream to provide and for us to adopt, since
we already parse the GAF it writes. Most of the ABI-stability reasoning below then becomes irrelevant,
because we would stay on the binary by choice rather than as a compromise.

Two things below are still worth keeping: the explanation of why we will not hand-maintain a decoder of
the on-disk format, and the questions about MAPQ/base qualities and about whether the schema is meant to
be third-party-readable. Both stand on their own.

**Three further data points from running this at scale** (HG002, 596,017,764 alignments on a 100 M-node
HPRC v2.1 MC CHM13 graph), all worth mentioning to Jouni because two are compliments and one is a
usability trap:

- **It performs very well.** `gaf-base sort` did 596 M records in 47 min at 6.5 GB peak, single-threaded,
  and `gaf-base construct` in 51 min at 15 GB — running *concurrently* through a FIFO, so 51 min wall in
  total. The database is 20.74 GB, i.e. **36.5 bytes per alignment**, about 3× smaller than the
  equivalent sorted GAM. For comparison `vg gamsort` needed 249 s on 8 M records where `gaf-base sort`
  needed 17 s.
- **The graph-consistency check is enforced, and that is the right call.** Querying the whole-genome
  GAF-Base with a GBZ-Base built from a single extracted chromosome fails with
  `The graph is not a valid reference for the alignments`. Worth saying we regard that as correct
  behaviour — it is exactly the mismatch we were expecting to have to police ourselves.
- **`gbz-base query` errors on a node ID that is not in the graph**, rather than ignoring it. With sparse
  ID spaces (100 M nodes across a 224 M ID range here) a caller must pre-filter every query against the
  graph, and a single stray ID aborts the whole query. A flag to skip unknown nodes would be a small,
  genuinely useful addition — this is probably the one concrete request worth making alongside the
  reads-only query.

Also, when rewriting: the count of format-breaking releases reads as a complaint when its purpose is only
to explain why we will not hand-maintain a decoder. Keep the fact, lose the tallying tone.

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
