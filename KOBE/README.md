# KOBE: Kit for orthology-based editing analysis – 8 modules in 1

A portable command-line toolkit for RNA-editing analysis. The package preserves the existing module workflows and exposes them through one dispatcher plus direct command aliases.

A module-by-module installation/data path review is documented in [PATH_AUDIT.md](PATH_AUDIT.md).

## Installation

Requirements on the host:

- Linux
- Bash
- Apptainer (or a compatible Singularity installation)
- Python 3 for launchers that inspect sample tables and paths
- `flock` for the `base_extract` and `find_back` output-directory locks
- `findmnt` is used when available to detect filesystem mount points; launchers fall back to binding the relevant input directory or file parent when it is unavailable

Unpack the archive to a stable location, then install the command links:

```bash
unzip KOBE_v1.0.0.zip -d "$HOME/software"
bash "$HOME/software/KOBE/install.sh"
```

The installer places command links in `~/.local/bin`. If that directory is not already in `PATH`, add this single line to `~/.bashrc`:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

Then run `source ~/.bashrc` or open a new shell.

You can also use the dispatcher without installing command aliases:

```bash
/path/to/KOBE/bin/KOBE --help
/path/to/KOBE/bin/KOBE base_extract --help
```

## Container image

The Apptainer SIF image is not bundled. Supply it either through `--sif FILE` where supported or by setting:

```bash
export RNA_EDITING_SIF="/path/to/rna_editing.sif"
```

The image must contain the module-specific tools and dependencies. Container-internal dependency locations such as `/opt/rna_editing` are intentional image paths, not host-server defaults. `RNA_EDITING_HOME` can override the container-internal RNA Editing dependency root where the module supports it.

## Available commands

Run `COMMAND --help` for the full options of each command.

| Command | Function |
|---|---|
| `snpeff_build` | Build a species-specific SnpEff database |
| `denovo_detect_map` | QC and map DNA/RNA sequencing data |
| `denovo_detect_anno` | VCF processing, SnpEff annotation, candidate filtering and final export |
| `base_extract` | Extract genomic base counts for known sites |
| `pairwise_tests` | Pairwise differential RNA-editing analysis |
| `ortho_edit` | Orthologous editing-site identification and cross-species mapping |
| `conserve_plot` | Conserved editing-site sequence extraction, filtering, plotting and coordinate export |
| `find_back` | Recover orthologous RNA-editing sites from RNA-only samples |
| `add_new_site` | Add newly recovered sites to the annotated table |
| `find_back_plot` | Plot conserved sites after adding recovered sites |
| `add_cds_sites_genomic_coordinates` | Convert CDS coordinates to genomic coordinates |
| `extract_cds_candidate_sites` | Format de novo editing sites for orthology input |
| `extract_longest_cds` | Extract the longest CDS transcript per gene |
| `genome_site_cds_position` | Map genomic candidate sites to CDS coordinates |
| `prepare_snpeff_gtf` | Prepare representative GTF and genome link for SnpEff |
| `supported_sites_to_editing_counts` | Convert supported-site tables to sample-level A/G count matrices |

## Pairwise differential analysis

The public command is `pairwise_tests` (the internal R script retains its
historical filename). Define each group with a repeated `--group` option and
list the requested pairwise comparisons in `--comparisons`:

```bash
pairwise_tests \
    --input /path/to/project/editing_counts.tsv \
    --group N=N1,N2,N3 \
    --group F=F1,F2,F3 \
    --group C=C1,C2,C3 \
    --comparisons N:F,N:C,F:C \
    --outdir /path/to/project/diff_result \
    --sif /path/to/rna_editing.sif \
    --bind /path/to/project
```

## Installation/data path separation

The KOBE installation directory is used only to locate KOBE's own launchers,
Snakefiles, helper scripts, and bundled software resources. It is never used
to infer a project-data directory. Project input files/directories and output
roots must be supplied through the relevant command-line options (or through
an explicitly documented sample manifest). If a required path is omitted,
the command stops with an error instead of searching relative to the KOBE
installation directory. Internal workflow subdirectories are created only
inside an output/work directory explicitly supplied by the user.

### GTF requirements differ by module

- `base_extract` requires `--gtf-dir` containing `<species>.gtf` for every
  species in the sample table. STAR genome-index generation uses this GTF via
  `--sjdbGTFfile`; it is not optional in this workflow.
- `find_back` accepts optional, repeatable `--gtf species=FILE` arguments.
  These GTFs are used only for the final CDS-to-genome coordinate conversion.
  A species without a supplied GTF still completes CDS-level site recovery,
  editing-level calculation, and Known/Novel classification; only genomic
  coordinate conversion is skipped for that species.

### BAM reuse with `--force_reditool`

`--force_reditool` is the explicit mode for reusing an existing final Picard
BAM and rerunning REDItoolKnown. It checks only the expected final Picard BAM
and its index for each sample, not a pre-Picard sorted BAM. The checks include
file presence/non-empty status, `samtools quickcheck`, `samtools idxstats`,
sample identity from BAM read groups, and exact reference sequence-name and
length dictionary compatibility. No MD5 or per-base sequence-content
comparison is performed. If the final Picard BAM is absent or fails validation,
the module rebuilds it from available FASTQ inputs; it never substitutes a
pre-Picard BAM for the required duplicate-removed BAM.

Expected final Picard BAM locations:

- `base_extract`: `<out-dir>/<species>/03_PICARD/<sample_id>/<sample_id>.sort.q10.nodup.bam`
- `find_back`: `<out-dir>/Mapping_BWA_to_CDS/<species>/<sample_id>/<sample_id>.sort.nodup.bam`

`base_extract` and `find_back` support explicit in-workdir BAM validation/reuse
when `--force_reditool` is supplied; their ordinary `--resume` behavior remains
signature/dependency-based. There is no separate `--reuse_bam` or `--bam-dir`
option. `base_extract` explicitly sets the Snakemake working directory to its
user-supplied `--out-dir`, so `.snakemake` state and locks stay in the project
work directory. `find_back` accepts optional repeated
`--gtf species=FILE` arguments for species-specific CDS-to-genome coordinate
mapping. `denovo_detect_anno` requires `--snpeff-dir` whenever SnpEff
annotation is actually needed. Container paths such as `/opt/rna_editing/...`
are intentional software-resource paths inside the SIF and are not host
project-data paths.

## Input paths, symlinks and Apptainer binds

The launchers automatically bind the paths they know about, including module scripts, declared input files/directories, output directories, sample-table FASTQ/reference paths, and resolved symlink targets where the module reads them. When `findmnt` is available, the containing filesystem mount is used; otherwise the launcher binds the relevant directory or file parent. This avoids assumptions about a particular cluster's `/data`, `/lab`, or home-directory layout.

Use repeatable `--bind` options when a workflow reads additional files not discoverable from its declared arguments or sample table. For example:

```bash
base_extract ... --bind /additional/storage:/additional/storage
```

The syntax of `--bind` follows each command's help. Some commands accept a simple directory path, while others also accept `HOST:CONTAINER` mappings.

Output directories are created on the host by launchers that manage output creation. Ensure the user has write permission to those locations.

## fastp minimum read length

The default fastp minimum read length is **150 bp** in the `find_back`, `base_extract`, and `denovo_detect_map` workflows. `find_back` exposes this as `--min_length` and allows the user to override it. `base_extract` and `denovo_detect_map` expose their existing fastp length options. This default does not change the separate fastp unqualified-base percentage setting (50%).

## Workflow logic

KOBE 1.0.0 includes the targeted workflow changes listed in `CHANGELOG.md`:
Picard duplicate removal in `find_back`, validated in-workdir BAM reuse under
`--force_reditool`, optional species-specific GTF coordinate mapping, and
updated `denovo_detect_anno` RNA thresholds. Other candidate-site selection,
REDItoolKnown count extraction, statistical methods, and existing result-table
schemas are retained. Validate each module with the target SIF and known test
datasets before production use; this package has not been end-to-end tested
against your cluster data in this build environment.


## Command-line type labels

Command help uses consistent value-type labels: `INT` (integer), `FLOAT`
(decimal number), `FILE` (file path), `DIR` (directory path), and `STRING`
(text or a structured string value such as `species=path`). For options that
also accept a special keyword, the accepted keyword is shown alongside the
type, for example `INT|none` or `INT|auto`.


## Boolean option handling

Boolean values passed by Snakemake command-line configuration are parsed
explicitly (`true`/`false`, `1`/`0`, `yes`/`no`, and `on`/`off`). This prevents
non-empty strings such as `"false"` from being interpreted as true in the
`denovo_detect_map` and `denovo_detect_anno` workflows.


## Optional fastp minimum-length value

In `denovo_detect_map`, `--fastp-length-required` accepts an integer, `none`,
or `off`. An integer is passed to fastp as its explicit `-l` value; `none`
and `off` omit that explicit argument. The separate
`--fastp-disable-length-filtering` switch passes fastp's `-L` option. These
controls are intentionally distinct.
