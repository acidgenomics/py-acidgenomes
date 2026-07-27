# Changelog

## 0.2.0 (2026-07-13)

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
