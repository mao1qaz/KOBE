# base_extract

A portable, multi-species command-line wrapper for the original Base Extraction Snakemake workflow.

## Analysis workflow

The original analysis order and count-table logic are retained:

1. fastp paired-end QC; the resulting paired clean FASTQ files are retained under `00_QC/<sample_id>/`.
2. STAR 2-pass alignment and coordinate sorting.
3. samtools MAPQ filtering (`MAPQ >= 10` by default), coordinate sorting and BAM indexing.
4. Picard MarkDuplicates on the Q10 BAM (`REMOVE_DUPLICATES=true` by default), retaining metrics and a BAM index.
5. REDItoolKnown base-count extraction using the Q10-filtered, Picard-deduplicated BAM.
6. Merge A/C/G/T counts from all samples of each species into one wide table; sites without reads are zero-filled.

The final table remains:

`<out-dir>/<species>/06_ALL_SAMPLES/<species>_all_samples_base_counts.tsv`

The per-species output structure remains:

```text
<out-dir>/
├── <species>/
│   ├── 00_REF/GENOME/
│   ├── 00_REF/STAR/
│   ├── 00_INPUT/
│   ├── 00_QC/<sample_id>/
│   ├── 01_MAP/<sample_id>/
│   ├── 02_Q10/<sample_id>/
│   ├── 03_PICARD/<sample_id>/
│   ├── 04_REDItools/<sample_id>/
│   └── 06_ALL_SAMPLES/<species>_all_samples_base_counts.tsv
├── logs/<species>/
├── tmp/STAR/<species>/<sample_id>/
├── tmp/Picard/<species>/<sample_id>/
├── .reuse_bam_inputs/       # internal staging of accepted Picard BAMs
├── .reuse_bam_validated/    # validation markers created by validate_reuse_bam.py
├── .reuse_bam_manifest.json # validation/reuse manifest
├── .base_extract.lock
└── base_extract.log
```

The original `prepare_known_sites.py`, `merge_all_sample_base_counts.py`, and `merge_reditools_counts.py` are included. `merge_all_sample_base_counts.py` is the script used by the current multi-sample Snakefile; `merge_reditools_counts.py` is retained as the original helper but is not called by the current Snakefile.

## Required input layout

### Sample table

A tab-delimited file with this header:

```text
sample_id    species    fq1    fq2    genome_fasta
```

An optional `type` column is accepted only when its values are `RNA`. The genome column can also be named `genome`, `reference_fasta`, or `ref`.

- `sample_id` and `species` must use letters, numbers, `.`, `_`, or `-` and begin with a letter or number.
- Relative FASTQ and genome paths are resolved relative to the sample-table directory.
- All samples of one species must point to the same genome FASTA.
- Normal runs require `fq1` and `fq2` values in the table. With `--force_reditool`, FASTQ fields may be empty when the sample has either a valid in-workdir Picard BAM/index or a retained complete clean FASTQ pair under `00_QC/<sample_id>/`.

Example: `samples.tsv.example`.

### GTF directory

For each species in the sample table, provide:

`<gtf-dir>/<species>.gtf`

### Known-site directory

For each species in the sample table, provide:

`<site-dir>/<species>_known_sites.tsv`

The known-site input is headerless and tab-delimited:

- Column 1: chromosome.
- Column 2: 1-based genomic position.
- Optional column 3: strand (`+`, `-`, `1`, `0`, `2`, `.`, `NA`, or empty).
- Additional columns are retained in the merged output.

The input is checked against the genome FASTA index. The prepared metadata retains the original input rows, while the current multi-sample merge creates one final row per unique chromosome-position coordinate. Duplicate coordinates with the same strand are therefore collapsed in the final wide table; conflicting strand assignments for one coordinate cause an error. Unknown strand is represented as `unknown_forward_reference`; it is not assumed to be plus strand.

## Basic command

```bash
base_extract \
    --sample-table "/path/to/project/base_extract_data/samples.tsv" \
    --out-dir "/path/to/project/base_extract_result" \
    --gtf-dir "/path/to/project/base_extract_data/gtf" \
    --site-dir "/path/to/project/base_extract_data/known_sites" \
    --threads 30 \
    --sif "/path/to/rna_editing_fonts.sif"
```

The GTF and known-site files are matched automatically using the species column. The launcher creates `--out-dir` on the host, resolves symlink targets, automatically binds the required host directories into Apptainer, and writes a run log to `<out-dir>/base_extract.log`.

The launcher holds a host-side `flock` lock for the entire Apptainer/Snakemake run. A second run targeting the same canonical output directory exits with an error. The launcher removes `.base_extract.lock` on normal exit and catchable termination signals. `SIGKILL` cannot be trapped; in that case the kernel releases the advisory lock when the process exits, and the next run safely reuses the lock path. Snakemake is always launched with `--directory <out-dir>`, so `.snakemake/` and its lock metadata stay inside the user-supplied project work directory rather than the caller's current directory or HOME.

For a `nohup` command, place the shell redirection log outside a not-yet-created output directory, or create that directory first. The launcher itself creates the output directory after it starts.

## Resume modes

### Default: fresh run

When neither `--resume` nor `--force_reditool` is supplied, the launcher uses Snakemake `--forceall`. Existing outputs are not treated as reusable results; all workflow rules are scheduled again.

### `--resume`

Reuse completed, up-to-date outputs and rerun missing, incomplete, or changed steps. The known-site file SHA256 is included in the `prepare_known_sites` rule parameters, so changing a site file triggers site preparation and dependent REDItoolKnown/count-merge jobs even if the file timestamp is unchanged. Upstream QC/Mapping outputs are retained when still valid.

If a FASTQ has been removed, `--resume` can continue only when the existing downstream outputs make that FASTQ unnecessary. If the missing FASTQ is needed to rebuild an incomplete step, Snakemake reports the missing input.

### `--force_reditool`

This is the only option that activates the dedicated in-workdir BAM validation/reuse path. A single `validate_reuse_bam.py` helper scans samples, validates BAM/index/reference compatibility, stages accepted BAMs, writes the manifest, and creates the validation markers consumed by Snakemake; no second BAM-validation script is used. Ordinary `--resume` remains Snakemake signature/dependency-based. It implies resume for upstream steps and forces `prepare_known_sites`, every sample-level `reditool_known` job, and `merge_all_samples`.

For each sample, KOBE checks the final BAM at:

`<out-dir>/<species>/03_PICARD/<sample_id>/<sample_id>.sort.q10.nodup.bam`

The BAM is reused only when:
- the BAM exists and is non-empty;
- its BAI index exists and is non-empty;
- `samtools quickcheck` succeeds;
- `samtools idxstats` can read the BAM/index;
- the BAM read-group `SM` value matches the sample ID;
- the BAM `@SQ` sequence names and lengths match the current reference FASTA dictionary.

A valid BAM bypasses fastp, STAR, Q10, and Picard for that sample and can be reused even when the original FASTQ files are no longer present. If the BAM/index is missing, incomplete, or incompatible, the affected mapping products are invalidated and the sample is remapped. KOBE first uses an available complete raw FASTQ pair (running fastp as needed); if raw FASTQ is unavailable, it can use the retained complete clean FASTQ pair under `00_QC/<sample_id>/` and skip fastp. If neither a valid BAM nor either complete FASTQ pair is available, the run stops with an error. The clean FASTQ pair is kept as a declared fastp output, while fastp HTML/JSON reports remain sidecar files.

### Removed options

The explicit `--reuse_bam` and `--bam-dir` options have been removed. BAM reuse is automatic only inside the current `--out-dir` and only when `--force_reditool` is supplied. Without `--force_reditool` or `--resume`, the default remains a fresh run and existing BAMs are not reused.

## Multi-process REDItoolKnown scheduling

REDItoolKnown is a separate Snakemake job for each sample. The default is `--reditools-threads 2`; Snakemake schedules independent sample jobs concurrently, subject to both `--threads` and the total `mem_mb` budget. There is no shared per-sample REDItoolKnown output directory.

Default resource settings:

- Total CPU budget: `--threads 30`.
- Total memory scheduling budget: `--max-memory-mb 90000`.
- REDItoolKnown per-sample thread request: `--reditools-threads 2`.
- REDItoolKnown per-sample memory request: `--reditools-memory-mb 12000`.
- Picard memory request: `--picard-memory-mb 81200`.

At 90,000 MB, REDItoolKnown alone could fit at most seven 12,000-MB jobs at once when no other jobs consume memory. Actual concurrency can be lower because other workflow jobs also request resources. Picard requests 81,200 MB and therefore normally runs without REDItoolKnown jobs concurrently under the default total-memory budget. The memory values are Snakemake scheduling requests, not hard process memory limits.

## All configurable analysis parameters

All values below can be changed on the command line. They are optional; ordinary runs need only the main paths and `--threads`.

### Thread requests

| Option | Default |
|---|---:|
| `--reference-prepare-threads` | 2 |
| `--fastp-threads` | 10 |
| `--star-index-threads` | 15 |
| `--star-threads` | 15 |
| `--samtools-threads` | 10 |
| `--picard-index-threads` | 2 |
| `--reditools-threads` | 2 |

### Memory requests (MB)

| Option | Default |
|---|---:|
| `--reference-prepare-memory-mb` | 6000 |
| `--fastp-memory-mb` | 4000 |
| `--star-index-memory-mb` | 64000 |
| `--star-map-memory-mb` | 16000 |
| `--q10-memory-mb` | 8000 |
| `--picard-memory-mb` | 81200 |
| `--reditools-memory-mb` | 12000 |
| `--max-memory-mb` | 90000 |

### fastp

| Option | Default |
|---|---:|
| `--fastp-n-base-limit` | 5 |
| `--fastp-qualified-quality` | 15 |
| `--fastp-unqualified-percent` | 50 |
| `--fastp-length-required` | 150 |
| `--fastp-disable-trim-poly-g` | enabled (`-G`) |

Use `--fastp-length-required none` to disable the minimum-length filter. Use `--fastp-enable-trim-poly-g` to enable poly-G trimming instead of the original `-G` behavior.

### STAR

| Option | Default |
|---|---:|
| `--star-sjdb-overhang` | 149 |
| `--star-limit-genome-generate-ram` | 60000000000 |
| `--star-limit-bam-sort-ram` | 12884901888 |
| `--star-genome-sa-index-nbases` | auto |
| `--star-multimap-max` | 1 |
| `--star-mismatch-max` | 10 |
| `--star-twopass-mode` | Basic |
| `--star-out-reads-unmapped` | Fastx |

### Q10, Picard and REDItoolKnown

| Option | Default |
|---|---:|
| `--q10-mapq` | 10 |
| `--picard-validation-stringency` | SILENT |
| `--picard-remove-duplicates` | enabled |
| `--picard-create-index` | enabled |
| `--picard-assume-sorted` | enabled |
| `--reditools-min-coverage` | 1 |
| `--reditools-min-base-quality` | 25 |
| `--reditools-min-mapq` | 10 |
| `--reditools-min-variant-reads` | 0 |
| `--reditools-min-editing-frequency` | 0.0 |

Boolean Picard options can be toggled with `--picard-keep-duplicates`, `--picard-no-create-index`, and `--picard-do-not-assume-sorted`.

## Container software

Software paths are discovered inside the SIF; no laboratory-server software paths are embedded in this package. The SIF must provide:

- fastp
- STAR
- samtools
- Python 3
- Snakemake
- Python 2 with pysam
- REDItoolKnown.py
- Picard JAR
- Java

The launcher searches the container PATH and the packaged `/opt/rna_editing` software tree for REDItools, Python 2, and Picard. `PICARD_JAR`, `REDI_PYTHON`, `REDITOOL_KNOWN`, `FASTP`, `STAR`, `SAMTOOLS`, `JAVA`, and `SNAKEMAKE` environment variables can override discovery inside the container if needed.

## Dry-run and help

```bash
base_extract --help
base_extract --sample-table samples.tsv --out-dir result --gtf-dir gtf --site-dir known_sites --sif rna_editing.sif --threads 30 --dry-run
```

The dry-run still validates the inputs and, with `--force_reditool`, stages/validates in-workdir BAM inputs before asking Snakemake to print the plan.
