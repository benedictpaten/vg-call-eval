"""Command line entry point."""

from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path

from .engines import aardvark
from .pipeline import build_tier0_dataset, default_arms, run_arm
from .simulate import SimParams


def _tool(name: str, override: str | None = None) -> str:
    path = override or shutil.which(name)
    if not path:
        raise SystemExit(f"required tool not found on PATH: {name}")
    return path


def cmd_run(args: argparse.Namespace) -> None:
    vg = _tool("vg", args.vg)
    av = _tool("aardvark", args.aardvark)

    params = SimParams(
        ref_length=args.ref_length,
        seed=args.seed,
        depth=args.depth,
        read_length=args.read_length,
    )

    root = Path(args.out).resolve()
    dataset = build_tier0_dataset(root / "dataset", params, vg, threads=args.threads)
    print(f"dataset ready: {dataset.directory}")
    print(f"  truth: {json.loads((dataset.directory / 'truth_counts.json').read_text())}")

    results = []
    for arm in default_arms(vg, args.vg_depthfix):
        print(f"running arm: {arm.name}")
        info = run_arm(arm, dataset, root / "calls", threads=args.threads)
        out_dir = root / "compare" / arm.name
        aardvark.compare(
            aardvark=av,
            reference=dataset.reference,
            truth_vcf=dataset.truth_vcf,
            query_vcf=Path(info["vcf"]),
            regions_bed=dataset.confident_bed,
            out_dir=out_dir,
            truth_sample=params.sample,
            query_sample=params.sample,
            label=arm.name,
            options=aardvark.AardvarkOptions(threads=args.threads),
        )
        info["compare_dir"] = str(out_dir)
        results.append(info)
        print(f"    {info['seconds']:.1f}s, compared -> {out_dir}")

    (root / "results.json").write_text(json.dumps(results, indent=2))
    print(f"\nwrote {root / 'results.json'}")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(prog="vgcalleval")
    sub = parser.add_subparsers(dest="command", required=True)

    run_p = sub.add_parser("run", help="build a tier-0 dataset and run the caller matrix")
    run_p.add_argument("--out", required=True)
    run_p.add_argument("--vg", default=None, help="path to the vg binary")
    run_p.add_argument("--vg-depthfix", default=None,
                       help="path to a vg built with the depth_err fix, to add that arm")
    run_p.add_argument("--aardvark", default=None)
    run_p.add_argument("--ref-length", type=int, default=200_000)
    run_p.add_argument("--depth", type=float, default=30.0)
    run_p.add_argument("--read-length", type=int, default=150)
    run_p.add_argument("--seed", type=int, default=1)
    run_p.add_argument("--threads", type=int, default=4)
    run_p.set_defaults(func=cmd_run)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
