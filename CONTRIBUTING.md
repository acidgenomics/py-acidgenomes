# Contributing

## Rules for feature work

- Never edit `version` or `current_version` in `pyproject.toml` on a feature
  commit. Only the release step (below) edits them.
- Add notes for unreleased work under a single `## Unreleased` heading in
  `CHANGELOG.md`. Never invent a version number on a feature commit.

These two rules exist because of an incident on 2026-09-25: two lines of work
each wrote `0.3.0` into `pyproject.toml` and `CHANGELOG.md` independently, and
one of them had already been released and published to PyPI under that number.

## Release checklist

1. Run `git switch develop && git pull`. Start the release from an integrated
   `develop`.
2. Read the `## Unreleased` notes and pick the version part to bump. A new
   public function is a minor bump. A fix alone is a patch bump.
3. Check that the number is free, before you write it anywhere:
   - `git fetch --tags --prune`
   - `git ls-remote --tags origin | grep -F "vX.Y.Z"` must print nothing.
   - `curl -s https://pypi.org/pypi/acidgenomics-acidgenomes/json | python3 -c
     "import json,sys; print(sorted(json.load(sys.stdin)['releases']))"` must
     not list `X.Y.Z`.
   - Stop if either command finds the number. Pick the next free number
     instead.
4. Rename `## Unreleased` to `## X.Y.Z (YYYY-MM-DD)` in `CHANGELOG.md`. Do not
   add a new empty `## Unreleased` heading. The next feature commit adds one.
5. Run the checks. All of these must pass:
   - `ruff format --check .`
   - `ruff check .`
   - `ty check`
   - `pyright`
   - `interrogate`
   - `pytest`
6. Run `bumpver update --minor` (or `--patch`). This commits
   `Bump to vX.Y.Z.`. It does not tag and does not push (`tag` and `push` are
   `false` in `[tool.bumpver]`; see below for why).
7. **Warning: publish from the bump commit. Do not add any commit between
   step 6 and step 8.** The publish step reads the version from
   `pyproject.toml` and tags the current `HEAD`. A commit in between would
   tie the release tag to code that was never built or published.
8. Run `koopa app python publish .`. This one command builds the package,
   refuses to overwrite an already-published artifact whose content differs,
   uploads to `python.acidgenomics.com` and then to PyPI, reindexes the
   private index, and creates and pushes the `vX.Y.Z` git tag.
9. Push `develop`, open the pull request to `main`, and merge it.
10. Run `koopa app python publish-docs .`.

### Recovery

- If the private-index upload succeeded but the PyPI upload failed (for
  example, a rate limit), run `koopa app python publish --pypi-only .`. It
  uploads the bytes already on the private index and never rebuilds.
- `--force` skips the artifact-collision check. Use it only for a deliberate,
  already-decided in-place correction of a version that is already live.
  Every other case should bump the version instead.

### Why `bumpver` does not tag

`koopa app python publish` already creates and pushes the `vX.Y.Z` tag. One
tool owning the tag means it can never point at a commit that was not
actually published. Leave `[tool.bumpver] tag = false`.

### Why `CHANGELOG.md` is not in `[tool.bumpver.file_patterns]`

`bumpver` replaces every match of the version pattern in a file. The
changelog holds every past release heading, so a `## {version} (date)`
pattern would also rewrite an earlier release's heading. Step 4 stays
manual.
