---
name: release
description: >-
  Steps for releasing a new version of py-acidgenomes to the private index
  and to PyPI. Use when asked to cut a release, bump the version, update the
  changelog for a release, or run koopa app python publish for this
  repository.
---

# Releasing py-acidgenomes

The full checklist lives in `CONTRIBUTING.md` at the repository root. Read
and follow it in order. Do not restate it here; read the file itself.

## Three traps that caused past incidents

- **Choose the version number only after the collision check in
  `CONTRIBUTING.md` step 3.** Do not pick a number by guessing the next one
  in sequence. Two independent lines of work claimed `0.3.0` this way on
  2026-09-25, and one of them was already released.
- **Never edit `version` or `current_version` in `pyproject.toml` on a
  feature commit.** Only the release step (`bumpver update`) touches them.
- **Never add a commit between the `bumpver update` commit and
  `koopa app python publish`.** The publish step tags whatever commit is
  checked out, so a commit in between ties the release tag to code that was
  never built or published.
