# DeepLB Manual (Compatibility-First)

## Scope
This manual documents the current DeepLB workflow as implemented in existing
scripts. It is intentionally compatibility-first and does not redefine the
pipeline logic.

## Primary entrypoint
- Shell: `Scripts/DeepLB_pipeline.sh`
- Python wrapper: `python -m deeplb -- [DeepLB_pipeline.sh args]`
- Optional installed CLI: `deeplb -- [DeepLB_pipeline.sh args]`

## Typical workflow
1. Prepare data under `Predata/`.
2. Validate paths in `Scripts/Part2.Pseudo-fragment_Generation_by_mMTS/env_module.py`.
3. Run one module (`part1`, `part2`, `part3`) or all modules via `-u`.

## Safety notes
- Existing script behavior is treated as source of truth.
- New Python interfaces are thin wrappers over the shell pipeline.
- Avoid changing model internals during this phase.
