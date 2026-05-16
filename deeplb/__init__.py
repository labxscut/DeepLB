"""Compatibility-first Python interface for DeepLB.

This package provides thin wrappers around the existing shell-driven workflow
so current project behavior is preserved.
"""

from .pipeline import run_pipeline

__all__ = ["run_pipeline"]
