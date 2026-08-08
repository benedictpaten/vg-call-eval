# Install

## Tools from your package manager

```bash
brew install bcftools samtools        # or apt-get
```

`vg` is not installed by this repo; pass its path with `--vg`.

## aardvark

**There is no macOS or ARM release binary** — as of v0.10.5 the only published asset is
`x86_64-unknown-linux-gnu`. On Linux x86_64, download it from
[releases](https://github.com/PacificBiosciences/aardvark/releases) and put it on `PATH`.

Everywhere else, build from source. It is pure Rust (`noodles`, not `rust-htslib`), so there are no
system C library dependencies and it builds without special handling:

```bash
brew install rust                     # or rustup
git clone --depth 1 https://github.com/PacificBiosciences/aardvark.git
cd aardvark && cargo build --release
cp target/release/aardvark ~/.local/bin/
```

Verified on macOS arm64 with aardvark 0.10.5.

## truvari

Needed for the tier-2 structural-variant comparison. aardvark's `Sv*` categories are scored against
the *small-variant* truth set, which holds no record over 50 bp, so they cannot answer an SV question;
truvari is scored against the structural benchmark and is the SV metric.

Kept in its own virtualenv so its dependency stack (numpy, pysam, pandas) stays out of the harness,
which otherwise uses only the standard library:

```bash
python3 -m venv work/truvari-venv && work/truvari-venv/bin/pip install truvari
```

`scripts/tier2/truvari_sv.py` looks for `work/truvari-venv/bin/truvari` by default; pass `--truvari`
to point elsewhere.

## Python

```bash
python3 -m pip install pytest         # tests only; the harness itself uses the stdlib
```

The harness deliberately drives `bcftools`, `aardvark` and `truvari` as subprocesses rather than
binding to them, so there is no compiled Python dependency to manage in the harness itself.
