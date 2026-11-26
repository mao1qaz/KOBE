# Changelog

## 1.0.0

- Consolidated `base_extract` BAM reuse preflight and per-sample validation into one `validate_reuse_bam.py` helper; the helper now creates the Snakemake validation markers directly, eliminating the duplicate validation script/rule.
- When a Picard BAM cannot be reused, `base_extract --force_reditool` can remap from retained clean FASTQ files if the raw FASTQ pair is unavailable.
- `base_extract --force_reditool` can fall back to a retained complete clean FASTQ pair when raw FASTQ files are unavailable and the existing Picard BAM is not reusable.

- Retains the portable KOBE dispatcher, installer, and module command aliases.
- Renamed the public `run_pairwise_tests` command to `pairwise_tests`; the pairwise help example now demonstrates three groups and multiple comparisons.
- Updated the main KOBE title and command ordering; removed the retired `check_genomic_reference_bases` utility and its command registration/documentation.
- Expanded `find_back --help` with `--mapping_jobs`, per-step Mapping resource controls, defaults, and the REDItoolKnown even-thread constraint.
- Updated `find_back` sample scheduling to launch the next pending sample as soon as any active Mapping/Picard workflow completes.
- Clarified conditional `--snpeff-dir` requirements in `denovo_detect_anno` help and moved threshold parameters with defaults to Optional.
- Added standalone help examples for the six add utilities; their `--help` output no longer requires SIF configuration or Apptainer.
- `find_back` accepts optional repeatable `--gtf species=FILE` arguments; species without a GTF retain CDS-level recovery and Known/Novel results while genomic-coordinate mapping is skipped.
- Added Picard MarkDuplicates to `find_back`; REDItoolKnown consumes the Picard duplicate-removed BAM. Clean FASTQ files are retained by both `base_extract` and `find_back`.
- Removed public `--reuse_bam` and `--bam-dir` options. `--force_reditool` validates and reuses only the expected in-workdir Picard BAM and index; without it, a default run starts fresh.
- BAM reuse checks file presence/non-empty status, index presence/readability, `samtools quickcheck`, `samtools idxstats`, sample read-group identity, and reference sequence names/lengths. No MD5 comparison is performed. Invalid BAMs fall back to remapping when FASTQ inputs are available.
- `base_extract` Snakemake uses the supplied `--out-dir` as its explicit working directory, keeping `.snakemake` state and locks in the project output directory. Host-side output locks are removed on normal exit and catchable termination; stale Snakemake locks are handled on the next run.
- `denovo_detect_anno` defaults are RNA minimum coverage 5 and RNA FDR 0.05.
- Preserves existing candidate-site selection, REDItoolKnown counting parameters, statistical methods, and result-table schemas except for the requested workflow additions.
