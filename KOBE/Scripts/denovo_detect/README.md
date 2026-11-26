# denovo_detect CLI module

## Script layout

The shared `container_runner.sh` is placed directly in the `denovo_detect/`
directory alongside `denovo_detect_map` and `denovo_detect_anno`. Both launchers
source it from the same directory; no `_common/` subdirectory is required.


This package provides two user-facing commands for the existing RNA-editing de novo discovery workflow:

- `denovo_detect_map`: Module 01 QC and DNA/RNA mapping.
- `denovo_detect_anno`: Module 02 mpileup/VCF normalization, SnpEff annotation, candidate extraction, statistical filtering, and final site export.

The YAML configuration files are no longer required. Defaults below reproduce the uploaded `01_config.yaml` and `02_config.yaml` values. All user-adjustable parameters are accepted directly on the command line.

## Installation

Install the complete toolkit once using the package-level `install.sh` (see the top-level README). This registers `denovo_detect_map` and `denovo_detect_anno` alongside the other module commands. Set the SIF path in each command using `--sif`, or export `RNA_EDITING_SIF` once in the shell. The launchers use the tools bundled in the Apptainer image; they do not use host-specific software paths.

## Module 01: denovo_detect_map

### Required sample table

The tab-separated sample table must contain:

```text
sample_id  sample_name  species  type  fq1  fq2
```

`type` must be `RNA` or `DNA`. `sample_id` must be unique; `sample_name` is the biological matching label used by downstream analysis. Relative FASTQ paths are resolved relative to the sample-table directory.

The reference directory must contain one FASTA and one GTF per species:

```text
ref/
├── invicta.fa
├── invicta.gtf
├── mellifera.fa
└── mellifera.gtf
```

### Example command

```bash
BASE="/path/to/05_ant_orth"
SIF="/path/to/rna_editing_fonts.sif"

 denovo_detect_map \
    --samples "$BASE/denovo_detect_data/samples.tsv" \
    --ref-dir "$BASE/denovo_detect_data/ref" \
    --out-dir "$BASE/denovo_detect_result" \
    --cores 30 \
    --fastp-length-required 150 \
    --sif "$SIF"
```

The command creates `--out-dir` on the host before Apptainer starts. The launcher discovers FASTQ paths from the sample table, resolves symlinks, and adds the corresponding host paths to Apptainer binds. Additional binds can be passed repeatedly with `--bind PATH[:DEST[:OPTS]]`.

### Defaults and command-line overrides

| Parameter | CLI option | Default |
|---|---|---:|
| Total Snakemake CPU limit | `--cores` | Required |
| fastp threads | `--fastp-threads` | 10 |
| BWA threads | `--bwa-threads` | 15 |
| STAR threads | `--star-threads` | 15 |
| samtools sort threads | `--samtools-sort-threads` | 10 |
| Sambamba threads | `--sambamba-threads` | 20 |
| fastp N-base limit | `--fastp-n-base-limit` | 5 |
| fastp qualified quality | `--fastp-qualified-quality-phred` | 15 |
| fastp unqualified percent | `--fastp-unqualified-percent-limit` | 50 |
| fastp minimum length | `--fastp-length-required` | 150 |
| fastp disable length filtering | `--fastp-disable-length-filtering` | false |
| fastp average quality | `--fastp-average-qual` | 0 |
| fastp disable poly-G trimming | `--fastp-disable-trim-poly-g` | true |
| fastp disable adapter trimming | `--fastp-disable-adapter-trimming` | false |
| fastp detect adapter for PE | `--fastp-detect-adapter-for-pe` | false |
| STAR junction overhang | `--star-sjdb-overhang` | 149 |
| STAR index RAM | `--star-limit-genome-generate-ram` | 60000000000 |
| STAR BAM-sort RAM | `--star-limit-bam-sort-ram` | 12884901888 |
| Sambamba sort memory | `--sambamba-sort-memory` | 6G |
| Sambamba overflow list size | `--sambamba-overflow-list-size` | 10000000 |

For `denovo_detect_map`, `--fastp-length-required none` or `off` omits the
explicit fastp `-l` argument. The separate `--fastp-disable-length-filtering`
option passes fastp `-L`; these controls are not interchangeable.

Memory requests are also configurable using `--mem-fastp-mb`, `--mem-reference-prepare-mb`, `--mem-bwa-index-mb`, `--mem-star-index-mb`, `--mem-bwa-map-mb`, `--mem-star-map-mb`, `--mem-picard-mb`, `--mem-gatk-split-mb`, and `--mem-sambamba-sort-mb`. Defaults are respectively 4000, 6000, 8000, 64000, 8000, 16000, 81200, 65000, and 8000 MB.

`--cores` is the total Snakemake CPU limit. Rule-level thread defaults remain as listed above; Snakemake may reduce individual rule thread allocations to respect `--cores`.

The output layout remains unchanged:

```text
denovo_detect_result/
├── 00_QC/
├── 00_REF/
├── 01_MAP/
└── tmp/
```

The host-side log is `denovo_detect_result/denovo_detect_map.log`.

## Module 02: denovo_detect_anno

The command uses the same Module 01 work directory and sample table. It reads BAMs from `01_MAP/` and reference FASTA/FAI files from `00_REF/GENOME/`.

`--snpeff-dir` is required whenever SnpEff annotation must run (for example,
when a selected species has no annotated VCF yet, or when `--reannotate` is
used). It is not required for `--filter-only`, or when all selected species
already have reusable annotated VCFs.

### Example: first complete run

All six statistical parameters are deliberately explicit in this example:

```bash
BASE="/path/to/05_ant_orth"
SIF="/path/to/rna_editing_fonts.sif"

denovo_detect_anno \
    --work-dir "$BASE/denovo_detect_result" \
    --samples "$BASE/denovo_detect_data/samples.tsv" \
    --snpeff-dir "$BASE/snpeff_result" \
    --species invicta \
    --dna-min-coverage 10 \
    --dna-fdr 0.05 \
    --rna-min-coverage 5 \
    --rna-fdr 0.05 \
    --min-support 2 \
    --total-error-rate 0.005 \
    --cores 36 \
    --sif "$SIF"
```

If `--species` is omitted, all species in `samples.tsv` are selected. It can be repeated or contain a comma-separated list.

### Re-run only statistical filtering with new thresholds

Existing annotated VCFs and candidate tables are reused automatically. For an explicitly downstream-only run, add `--filter-only`:

```bash
denovo_detect_anno \
    --work-dir "$BASE/denovo_detect_result" \
    --samples "$BASE/denovo_detect_data/samples.tsv" \
    --species invicta \
    --dna-min-coverage 10 \
    --dna-fdr 0.05 \
    --rna-min-coverage 5 \
    --rna-fdr 0.05 \
    --min-support 2 \
    --total-error-rate 0.005 \
    --filter-only \
    --cores 36 \
    --sif "$SIF"
```

`--filter-only` requires existing `02_anno/<species>/genic.candidates.tsv` and `intergenic.candidates.tsv`. It does not run BAM validation, mpileup, VCF normalization, SnpEff, or candidate extraction.

### Re-annotate existing normalized VCFs

```bash
denovo_detect_anno \
    --work-dir "$BASE/denovo_detect_result" \
    --samples "$BASE/denovo_detect_data/samples.tsv" \
    --snpeff-dir "$BASE/snpeff_result" \
    --species invicta \
    --dna-min-coverage 10 \
    --dna-fdr 0.05 \
    --rna-min-coverage 5 \
    --rna-fdr 0.05 \
    --min-support 2 \
    --total-error-rate 0.005 \
    --reannotate \
    --cores 36 \
    --sif "$SIF"
```

`--reannotate` forces `snpeff_annotate` and `extract_regions_prefilter`. Existing mpileup and normalized VCF outputs are reused unless they are missing or stale according to Snakemake dependencies. The current threshold's downstream outputs are then regenerated. `--filter-only` and `--reannotate` are mutually exclusive.

When annotation must run (a required annotated VCF is missing, or `--reannotate` is selected), `--snpeff-dir` must explicitly identify the prebuilt SnpEff database root. If all annotated VCFs already exist, or `--filter-only` is used, this database directory is not needed. KOBE does not infer it from its installation directory or silently assume `$WORK_DIR/snpeff_result`.

### Module 02 defaults

| Parameter | CLI option | Default |
|---|---|---:|
| DNA minimum coverage | `--dna-min-coverage` | 10 |
| DNA FDR | `--dna-fdr` | 0.05 |
| RNA minimum coverage | `--rna-min-coverage` | 5 |
| RNA FDR | `--rna-fdr` | 0.05 |
| Minimum independent RNA sample support | `--min-support` | 2 |
| Total sequencing error rate | `--total-error-rate` | 0.005 |
| Genic extraction mode | `--genic-mode` | c_numeric |
| Accept simple SNV during normalization | `--accept-simple-snv` / `--no-accept-simple-snv` | enabled |
| mpileup minimum base quality | `--min-base-quality` | 30 |
| mpileup options | `--mpileup-options` | `-u -v -t ADF,ADR` |
| SnpEff options | `--snpeff-options` | `-canon -ud 0 -v -nodownload -noLog` |
| SnpEff Java heap | `--snpeff-heap` | 8g |
| mpileup memory request | `--mpileup-memory-mb` | 12000 MB |
| SnpEff memory request | `--snpeff-memory-mb` | 12000 MB |

The SnpEff JAR and Java 25 launcher default to the container paths:

- `/opt/rna_editing/jars/snpEff.jar`
- `/opt/rna_editing/bin/java25`

They can be overridden with `--snpeff-jar` and `--snpeff-java`. SnpEff database/config files are matched by species under `--snpeff-dir/database/<species>/` and `--snpeff-dir/config/<species>.config`.

### Threshold result naming

The existing five-field threshold tag is retained without adding a new run ID:

```text
DNA<coverage>_RNA<coverage>_NRNA<support>_DNAFDR<cutoff>_RNAFDR<cutoff>
```

`total_error_rate` is intentionally not included in the tag. Its value is written to `denovo_detect_anno.log` and passed to the existing R filtering script.

### Output layout (unchanged)

```text
denovo_detect_result/
├── 00_QC/
├── 00_REF/
├── 01_MAP/
├── 02_anno/
│   └── <species>/
│       ├── DNA_RNA.vcf30.vcf
│       ├── DNA_RNA.var.vcf
│       ├── DNA_RNA.var.snpeff.vcf
│       ├── genic.candidates.tsv
│       ├── intergenic.candidates.tsv
│       └── threshold-specific statistics and logs
├── edit_site_result/
│   └── <species>/
│       └── threshold-specific final site tables and PDFs
└── denovo_detect_anno.log
```

The threshold tag distinguishes different combinations of DNA coverage, RNA coverage, RNA support, DNA FDR, and RNA FDR. The output files for a repeated identical tag are updated by Snakemake when their dependencies are forced or newer.

## Runtime and locking

The intended usage is sequential: finish `denovo_detect_map` before starting `denovo_detect_anno`. No additional cross-module working-directory lock is introduced. Snakemake's own normal lock behavior is retained; do not use `--nolock`.

Both launchers create their output/log directories on the host, bind input/output paths into Apptainer, resolve input symlinks, support repeated `--bind`, and write logs under the work directory.
