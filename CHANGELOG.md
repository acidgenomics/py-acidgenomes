# Changelog

## 0.2.0 (2026-08-21)

### Bug Fixes

- Work around a bug in `genomicranges` 0.8.4, where
  `GenomicRanges.to_pandas()` silently doubles the row count. It sets the range
  DataFrame's index to `gr.names` before it concatenates `gr.mcols.to_pandas()`,
  whose index is still a positional `RangeIndex`. `pandas.concat(..., axis=1)`
  then outer-joins two disjoint indices instead of binding columns
  positionally, so every row loses either its coordinate columns or its
  metadata columns.
- Add `_genomicranges_to_pandas_safe()`, and use it in place of
  `GenomicRanges.to_pandas()` at both call sites: the `_GRangesWrapper.data`
  property, and `make_ensembl_genes_from_gtf()`. Every provider (Ensembl,
  GENCODE, RefSeq, UCSC, FlyBase, WormBase) routes through the shared
  `dataframe_to_granges` / `make_granges_from_gff` pipeline, so this fixes all
  of them, not Ensembl alone.
- The `0.2.0` package index artifact was rebuilt in place on 2026-08-21 to
  carry this fix. The `v0.2.0` git tag still points at the original
  2026-07-13 commit; the fix landed in a later commit on `main`, so
  `git checkout v0.2.0` will not reproduce the artifact currently served
  from the index. See git history for the exact commit if this matters.

### Tests

- Add `tests/test_genomicranges_to_pandas_safe.py`. It covers the upstream bug
  and the workaround through the real `dataframe_to_granges` construction path.

### Changes

- Handle `None`/NaN `biotype` values in `_apply_broad_class` that caused
  `AttributeError: 'float' object has no attribute 'lower'` with Ensembl
  release 116 GTF files.

## 0.1.0 (2026-06-19)

### Changes

- Switch license to Apache-2.0.
- Publish to `python.acidgenomics.com` (private PEP 503 index).
- Update installation instructions in README.

---

## 0.0.1 (initial)

Initial release.
