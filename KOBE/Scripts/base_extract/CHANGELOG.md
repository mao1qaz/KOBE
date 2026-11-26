# Changelog

## Portable CLI package

- Added the `base_extract` host-side command launcher.
- Added host-side output-directory creation, symlink-aware Apptainer binds, repeatable `--bind`, and a host-side `flock` lock.
- Replaced hardcoded laboratory software paths with container-side software discovery.
- Exposed the original analysis, thread, and memory settings as command-line options while retaining the original defaults, except REDItoolKnown threads default to 2 per sample as requested.
- Added multi-species GTF/site auto-matching from the sample-table `species` column.
- Added `--resume`, `--force_reditool`, and `--reuse_bam` modes (hyphenated aliases are also accepted).
- Added per-sample BAM preflight, symlink staging, and reference-dictionary validation for `--reuse_bam`.
- Preserved the original fastp -> STAR -> Q10 -> Picard -> REDItoolKnown -> per-species A/C/G/T merge workflow.


## Review fixes

- Added the missing `json` import required to read the `--reuse_bam` manifest in the Snakemake workflow.
- Standardized the `REDITOOL_KNOWN` environment-variable lookup.


## Current v1.0.0 behavior

- The public `--reuse_bam` and `--bam-dir` options are removed; `--force_reditool` is the only mode that activates automatic in-workdir Picard BAM validation/reuse.
- BAM preflight, reference-dictionary checks, staging, manifest creation, and Snakemake validation-marker creation are consolidated into the single `validate_reuse_bam.py` helper. The separate `prepare_reuse_bams.py` helper and per-sample Snakemake validation script/rule are no longer used.
- If a final Picard BAM is not reusable, `--force_reditool` can remap from the original FASTQ pair or, when those are unavailable, from the retained paired clean FASTQ files in `00_QC/<sample_id>/`.
