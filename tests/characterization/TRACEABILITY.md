# Characterization Test Traceability

This document maps each characterization test to the pipeline producer it pins. It is the canonical answer to "if I touch <file>, which tests should I expect to flicker?"

This file is updated as `hpc`-tier and `local_integration`-tier tests are added. Comparator unit tests (`local_unit`) are NOT listed here — they test the comparators themselves, not the pipeline.

## Tests

| Test (file::function) | Producer (file::function) | Output pinned | Strategy | Notes |
|---|---|---|---|---|

(rows added as tests are written)

## Conventions

- Test naming: `test_<output>_from_<producer>` where producer is the `.py` file or function name.
- Producer column: file path relative to repo root, optionally with `::function` suffix.
- Output pinned: relative path from `reference_root`.
- Strategy: one of `CSV-EXACT`, `CSV-STRUCT`, `JSON-DEEP`, `JSON-MODULO-PATHS`, `PNG-PERCEPTUAL`, `PNG-EXISTS`.
