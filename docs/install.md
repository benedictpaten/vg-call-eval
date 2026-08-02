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

## Python

```bash
python3 -m pip install pytest         # tests only; the harness itself uses the stdlib
```

The harness deliberately drives `bcftools` and `aardvark` as subprocesses rather than binding to
them, so there is no compiled Python dependency to manage.
