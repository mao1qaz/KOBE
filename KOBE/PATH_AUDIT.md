# KOBE v1.0.0 portability/path-isolation audit

## Policy

The KOBE installation directory is used only to locate KOBE-owned executables,
launchers, Snakefiles, helper scripts, and bundled software resources. It must
not be used to infer a user's project root, input-data directory, or output
location.

Project input and output paths come from command-line options or from paths
explicitly recorded in a user-supplied sample manifest. Relative FASTQ paths
inside sample manifests are resolved relative to the manifest file itself.
Internal workflow subdirectories may be created beneath an explicitly supplied
work/output root.

Container paths such as `/opt/rna_editing/...` are intentional paths inside the
Apptainer image for bundled software, helper scripts, fonts, or reference
resources. They are not host-side project-data paths.

## Module audit

| Module/command | User-supplied project inputs | User-supplied output/work root | Conditional explicit path requirements |
|---|---|---|---|
| `base_extract` | `--sample-table`, required `--gtf-dir`, `--site-dir`; FASTQ/genome paths in the sample table | `--out-dir` | GTF is required for STAR genome-index generation; `--force_reditool` validates/reuses only final in-workdir Picard BAMs |
| `find_back` | `--sample_table`, `--cds_pos_dir`, optional repeatable `--gtf species=FILE`, `--edit_txt_dir`, `--candi_dir`; FASTQ/reference paths in the sample table | `--out_dir` | GTF is optional and per-species; missing GTF skips only genomic-coordinate mapping |
| `denovo_detect_map` | `--samples`, `--ref-dir`; FASTQ paths in the sample table | `--out-dir` | All listed paths are required |
| `denovo_detect_anno` | `--samples`, Module 01 result root via `--work-dir` | `--work-dir` | `--snpeff-dir` is mandatory when SnpEff annotation is needed |
| `snpeff_build` | `--gtf`, `--genome` | `--outdir` | All listed paths are required |
| `ortho_edit` | `--cds_dir`, `--edit_txt_dir` | `--work_dir` | All listed paths are required |
| `conserve_plot` | `--annotated`, `--cds_dir`, optional repeatable `--gtf` | `--outdir` | Required input/output paths are validated by the launcher |
| `add_new_site` | `--annotated`, one or more `--sites` | `--outdir` | All listed paths are required |
| `find_back_plot` | `--annotated`, `--cds_dir` | `--outdir` | All listed paths are required |
| `pairwise_tests` | `--input`, group/comparison definitions | `--outdir` | Input and output paths are required |
| `add_*` utilities | Required `--input`, `--gtf`, `--genome`, `--sites`, or `--reference_dir` as appropriate | Required `--output`, `--outdir`, `--output_dir`, or `--work_dir` | Python argument parsers mark project paths as required |

## v1.0.0 workflow changes

1. `find_back` accepts repeated `--gtf species=FILE` arguments. The species
   identifier must match the sample-table `species` column. GTFs are used only
   in the final CDS-to-genome mapping; species without a GTF retain all CDS-
   level recovery and Known/Novel results.
2. `find_back` adds Picard `MarkDuplicates` after BWA coordinate sorting.
   REDItoolKnown reads from the Picard duplicate-removed BAM. Both `find_back`
   and `base_extract` retain their clean FASTQ outputs.
3. The public `--reuse_bam` and `--bam-dir` options are removed. With
   `--force_reditool`, both modules inspect the in-workdir BAM/index and verify
   BAM integrity, index readability, sample read-group identity, and reference dictionary compatibility.
   Both modules check BAM `@SQ` sequence names and lengths against the current
   reference dictionary; no MD5 or sequence-content fingerprint comparison is used.
4. `base_extract` runs Snakemake with `--directory <out-dir>`. Its `.snakemake`
   state and locks therefore reside in the user-supplied project work directory,
   not the caller's current directory or HOME. A host-side `flock` protects the
   output directory; catchable exits remove the lock-file pathname. A
   subsequent run clears a stale Snakemake lock left by an uncatchable kill.
5. `find_back` also removes its host lock-file pathname on normal exit and
   catchable termination signals. As with any process, `SIGKILL` cannot run a
   cleanup trap; the kernel releases `flock` when the owning process exits.
6. `denovo_detect_anno` defaults are RNA minimum coverage 5 and RNA FDR 0.05.

## Verification performed for this source revision

- Bash syntax validation of 29 shell entry points/helpers.
- Python syntax compilation of 36 Python source files.
- `find_back` GTF mapper smoke test: one species with a supplied GTF receives
  genomic coordinates, while another species without a GTF is skipped without
  aborting the run.
- `find_back` BAM validator smoke test: valid BAM/index passes; a same-length
  reference sequence names/length compatibility is checked using the BAM header.
- `base_extract` BAM preflight smoke test: a complete BAM/index with matching
  reference dictionary is staged; an incompatible dictionary is rejected and
  falls back to FASTQ mapping when FASTQ inputs are available.
- CLI help checks confirm `find_back --gtf species=FILE`, `--mapping_jobs`,
  the REDItoolKnown even-thread constraint, and `denovo_detect_anno` defaults
  (RNA coverage 5, RNA FDR 0.05).
- `find_back` sample scheduling uses `wait -n -p` so any completed sample
  releases a slot immediately for the next pending sample (Bash 5.1+).
- The six add utilities display argparse help/examples without requiring SIF or
  Apptainer; the retired genomic-reference-base utility is absent.
- Host-launcher smoke tests with a stub Apptainer executable cover `find_back`
  with and without GTF arguments and `base_extract --force_reditool` with empty
  FASTQ fields; both launchers removed their output-directory lock files after
  the stubbed process exited.

These checks validate syntax, helper behavior, and selected path/validation
logic. They do not replace end-to-end testing with the target Apptainer SIF and
real project datasets on the target server. Snakemake is not installed in the
current build environment, so the modified Snakefile has not been executed
against the target SIF here.
