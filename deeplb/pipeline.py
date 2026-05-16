"""Thin wrappers for invoking the existing DeepLB pipeline script."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Iterable, Sequence


def _repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def _pipeline_script() -> Path:
    return _repo_root() / "Scripts" / "DeepLB_pipeline.sh"


def run_pipeline(
    script_args: Sequence[str] | None = None,
    *,
    dry_run: bool = False,
    check: bool = True,
    capture_output: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Run the existing DeepLB shell pipeline with pass-through arguments.

    Parameters
    ----------
    script_args:
        Arguments passed to ``Scripts/DeepLB_pipeline.sh``.
    dry_run:
        If True, appends ``-d`` unless it is already present.
    check:
        If True, raises ``CalledProcessError`` for non-zero exit status.
    capture_output:
        If True, captures stdout/stderr and returns them in the result.
    """
    args = list(script_args or [])
    if dry_run and "-d" not in args:
        args.append("-d")

    command: Iterable[str] = ["bash", str(_pipeline_script()), *args]
    return subprocess.run(
        command,
        cwd=_repo_root(),
        text=True,
        check=check,
        capture_output=capture_output,
    )
