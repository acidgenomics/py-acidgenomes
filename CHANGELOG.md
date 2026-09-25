# Changelog

## Unreleased

### New Features

- Add `classify_curated_gene_groups()`, tagging genes into curated HGNC
  gene groups (`ribo_cyto`, `ribo_mito`, `hemoglobin`) sourced from HGNC's
  own `gene_group`/`gene_group_id` assignments, never a symbol regex. For
  Mus musculus, human HGNC groups are propagated via a fully
  identifier-based chain (HGNC `hgnc_id` -> JAX ortholog `mouse_mgi_id` ->
  MGI `ensembl_gene_id`), with no gene-symbol matching at any step. This is
  deliberately independent of `broad_class`, which is single-valued and
  already assigns every ribosomal/hemoglobin gene a value (`coding`,
  `pseudo`, etc).

### Bug Fixes

- Fix `make_jax_human_to_mouse()` silently renaming its own documented
  `mouse_mgi_id` column to `mouse_mgi_id_y` (with a meaningless,
  always-`NaN` `mouse_mgi_id_x` alongside it). The raw JAX report stacks
  one column per field across both species' rows; `_merge_jax_species()`
  never dropped the human-side copy before merging in the mouse-side
  value, so pandas resolved the name collision with its default `_x`/`_y`
  suffixing instead of the plain column name every consumer expects.

## 0.3.0 (2026-09-20)

### Bug Fixes

- Restore automatic Ensembl release detection using the machine-readable
  `VERSION` endpoint after Ensembl removed `current_README`.
- Replace the unavailable `useast.ensembl.org` endpoint with the canonical
  REST and website hosts for current and archived Ensembl queries.
- Support current `genomicranges` releases, which no longer re-export
  `GenomicRanges` from the package root and have fixed the historical
  `to_pandas()` row-count defect.

### Tests

- Add regression coverage for the Ensembl release endpoint, release-pinned
  HTTPS downloads, and archive URL resolution.

## 0.2.1 (2026-09-01)

### Changes

- Rename the PyPI distribution to `acidgenomics-acidgenomes`. The import name is
  unchanged: `import acidgenomes as ag` still works.
- Publish to PyPI instead of `python.acidgenomics.com`.
- Pin `[tool.uv.build-backend] module-name = "acidgenomes"` explicitly, since the
  build backend's default module name now derives from the distribution name
  (`acidgenomics_acidgenomes`), not the import name.

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
  carry this fix. The `v0.2.0` git tag was moved to match, so
  `git checkout v0.2.0` reproduces the artifact currently served from the
  index.

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
