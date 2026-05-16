"""CLI bridge to keep DeepLB behavior stable while exposing a Python entrypoint."""

from __future__ import annotations

import argparse

from .pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="deeplb",
        description="Thin Python CLI wrapper for Scripts/DeepLB_pipeline.sh",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Append -d to the underlying DeepLB shell command.",
    )
    parser.add_argument(
        "script_args",
        nargs=argparse.REMAINDER,
        help="Arguments forwarded directly to Scripts/DeepLB_pipeline.sh",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    script_args = list(args.script_args)
    if script_args and script_args[0] == "--":
        script_args = script_args[1:]

    result = run_pipeline(script_args, dry_run=args.dry_run, check=False)
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
