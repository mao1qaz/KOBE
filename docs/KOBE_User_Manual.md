# KOBE User Manual

**Version 1.0.0**

## Contents

1. [Introduction](#1-introduction)
2. [Installation and Setup](#2-installation-and-setup)
3. [De Novo RNA Editing Identification](#3-de-novo-rna-editing-identification)
4. [Orthology-based RNA Editing Analysis](#4-orthology-based-rna-editing-analysis)
5. [Differential RNA Editing Analysis](#5-differential-rna-editing-analysis)
6. [Genomic Base Extraction](#6-genomic-base-extraction)
7. [Utility Scripts](#7-utility-scripts)

---

# 1. Introduction

## 1.1. About KOBE

KOBE (K̲it for o̲rthology-b̲ased e̲diting analysis – 8 modules in 1) is a command-line toolkit for identifying and analyzing adenosine-to-inosine (A-to-I) RNA editing. It integrates de novo RNA editing site identification, orthology-based RNA editing analysis, orthology-guided site recovery, differential RNA editing analysis, and downstream visualization and data processing. Its modular design supports flexible analyses across samples and species.

## 1.2. Main Features

KOBE provides the following major functionalities:


- **De novo RNA Editing Site Identification:** Identifies candidate RNA editing sites from RNA sequencing data, with genomic DNA sequencing data used to distinguish RNA editing events from genomic variants. Identified sites can be further annotated to determine their genomic and functional characteristics.

- **Orthology-based RNA Editing Analysis:** Analyzes RNA editing sites across species based on orthologous gene or transcript relationships, enabling the identification and comparison of conserved and species-specific editing sites.

- **Conservation Visualization:** Generates visual summaries of conserved RNA editing sites and their distribution across species, facilitating comparisons of conserved and species-specific editing patterns.

- **Orthology-guided Site Recovery:** Recovers potential RNA editing sites in a target species using homologous or orthologous sequence information from other species, helping identify editing sites that may not be detected through conventional de novo approaches.

- **Post-recovery Integration and Visualization:** Integrates newly recovered sites into the orthology-based results and visualizes the updated conserved-site groups.

- **Differential RNA Editing Analysis:** Identifies differences in RNA editing levels between experimental conditions or sample groups, supporting the investigation of condition-associated RNA editing patterns.

- **Genomic Base Extraction:** Extracts nucleotide counts (A, C, G, and T) at specified genomic positions across samples, providing base-level information for downstream analyses and validation.

- **Utility Scripts:** Provides additional tools for data preparation, coordinate conversion, site annotation, result integration, and other auxiliary operations required during RNA editing analyses.

## 1.3. Intended Use and Scope

KOBE is intended for researchers studying RNA editing at the transcriptomic and comparative levels. Its modular design supports several types of research applications:

1. **Genome-wide RNA editing identification:** Identification and functional characterization of candidate RNA editing sites using RNA sequencing data and, where applicable, matched DNA sequencing data.

2. **Comparative RNA editing analysis:** Investigation of conserved and species-specific RNA editing patterns based on orthologous relationships across multiple species.

3. **Recovery of potential RNA editing sites:** Identification of additional candidate editing sites through orthology-guided analysis, complementing conventional de novo identification approaches.

4. **Differential RNA editing analysis:** Comparison of RNA editing levels between different experimental conditions, treatments, or sample groups.

5. **Downstream nucleotide-level analysis:** Extraction and integration of genomic base counts and other supporting information for further investigation of candidate editing sites.

KOBE is particularly useful for studies involving multiple samples, experimental conditions, or species, where consistent analysis procedures and comparable output formats are important.

## 1.4. Modular Design

KOBE is organized into independent modules, each responsible for a specific analytical task. Users can execute individual commands according to their data and research objectives without necessarily running the entire workflow.

For analyses involving multiple modules, the output of one module can serve as the input to subsequent modules, enabling a flexible and reproducible workflow.

Detailed descriptions of installation, data preparation, workflow organization, individual commands, and their corresponding input and output formats are provided in the following sections of this guide.

---

# 2. Installation and Setup

## 2.1. Download KOBE

Download the KOBE software package and its Apptainer container image.

- **KOBE software:** `KOBE_v1.0.0.zip`
- **Apptainer container image:** Download link will be provided here.

> Note: The container image download link will be added when the image is publicly available.

## 2.2. Install KOBE

Extract the downloaded software package and run the installation script.

```bash
mkdir -p /path/to/software
cd /path/to/software

# 2. Download KOBE_v1.0.0.zip into this directory
# (Download the package from the provided link)

# 3. Extract the package
mkdir -p KOBE_v1.0.0
unzip -o KOBE_v1.0.0.zip -d KOBE_v1.0.0

# 4. Run the installation script
bash KOBE_v1.0.0/KOBE/install.sh
```

Replace `/path/to/software` with your preferred installation directory.

## 2.3. Install Apptainer

KOBE uses Apptainer to run its containerized analysis environment. Apptainer 1.5.4 is recommended.

Create and activate a dedicated Conda environment:

```bash
# Create the Conda environment
conda create -y \
    -n apptainer-1.5.4 \
    --override-channels \
    -c conda-forge \
    apptainer=1.5.4 \
    python=3.11

# Activate the environment
conda activate apptainer-1.5.4
```

The Apptainer environment should be activated before running KOBE commands that require the container.

Download the KOBE container image and keep it in an accessible directory. 

## 2.4. PATH Configuration

The installation script automatically creates symbolic links for all KOBE commands in `$HOME/.local/bin/` by default.

After creating the command links, the installer checks whether the installation directory is included in the system `PATH`.

- If the directory is already in `PATH`, no additional configuration is required.
- If the directory is not in `PATH`, the installer displays the required configuration command.

In the latter case, add the displayed command to `~/.bashrc` and reload the shell configuration:

```bash
source ~/.bashrc
```

After configuration, all KOBE commands can be executed directly from any working directory.

> Note: The default command installation directory can be customized using the `RNA_EDITING_BIN_DIR` environment variable.

## 2.5. Verify Installation

**1. Check Apptainer installation**

```bash
apptainer --version
```

Expected output:

```text
apptainer version 1.5.4
```

**2. Check KOBE command availability**

```bash
command -v KOBE
```

Expected output (default installation):

```text
/home/username/.local/bin/KOBE
```

The displayed path depends on the command installation directory selected during installation.

**3. Check KOBE version**

```bash
KOBE --version
```

Expected output:

```text
KOBE 1.0.0
```

**4. Check container availability and write permissions**

Specify the paths to the downloaded container image and the project output directory:

```bash
SIF="/path/to/rna_editing.sif"
TEST_DIR="/path/to/project/output"
```

Run the following test:

```bash
mkdir -p "$TEST_DIR" &&
apptainer exec \
    --bind "$TEST_DIR:/tmp" \
    "$SIF" \
    sh -c 'f=$(mktemp /tmp/.kobe_write_test.XXXXXX) && rm "$f"' &&
echo "PASS: Container is accessible and the output directory is writable."
```

Expected output:

```text
PASS: Container is accessible and the output directory is writable.
```

Replace the example paths with the actual container image and project output directory. The test creates and removes a temporary file without modifying existing project data.

If all four checks pass, KOBE and its container environment are ready to use.

## 2.6. Running Long-running Commands in the Background

For long-running analyses, KOBE commands can be executed in the background to prevent interruption when the terminal session is closed.

**Recommended command format:**

```bash
nohup setsid --fork --wait KOBE <command> \
    [options] \
    > "/path/to/command.log" 2>&1 < /dev/null &
```

**Command components:**

| Component | Description |
|---|---|
| `nohup` | Prevents the process from terminating when the terminal session is closed. |
| `setsid` | Starts the command in a new session, detached from the current terminal session. |
| `--fork` | Creates a child process to execute the command. |
| `--wait` | Waits for the child process to finish and returns its exit status. |
| `>` | Redirects standard output to the specified log file. |
| `2>&1` | Redirects standard error to the same log file. |
| `< /dev/null` | Redirects standard input from `/dev/null`. |
| `&` | Runs the command in the background. |

**Example:**

```bash
nohup setsid --fork --wait KOBE denovo_detect_map \
    --samples "/path/to/denovo_detect_data/samples.tsv" \
    --ref-dir "/path/to/denovo_detect_data/ref" \
    --out-dir "/path/to/denovo_detect_result" \
    --cores 30 \
    --sif "/path/to/rna_editing.sif" \
    > "/path/to/denovo_detect_map.log" 2>&1 < /dev/null &
```

---

# 3. De Novo RNA Editing Identification

## 3.1. SnpEff Database Construction

To build a species-specific SnpEff database, provide a reference genome FASTA file and its corresponding full genome GTF annotation. KOBE prepares the representative-transcript GTF required by SnpEff and then builds the database.

### 3.1.1. Input Data

| Input | Format | Description |
|---|---|---|
| Reference genome | FASTA (`.fa`) | Genome assembly sequence. |
| Gene annotation | GTF (`.gtf`) | Full genome annotation corresponding to the reference genome; used for RNA-seq mapping and preparation of the representative annotation for SnpEff. |


### 3.1.2. Prepare SnpEff Input Files

The `prepare_snpeff_gtf` command prepares the files required for SnpEff database construction.

It locates the species-specific genome FASTA and full GTF files in the reference directory, creates a symbolic link to the genome FASTA in the output directory, and generates a representative GTF retaining the transcript with the longest CDS for each coding gene and the transcript with the longest spliced exon length for each noncoding gene. If exon features are unavailable, the transcript feature span is used.

**Command:**

```bash
KOBE prepare_snpeff_gtf \
    --reference_dir "/path/to/denovo_detect_data/ref" \
    --output_dir "/path/to/snpeff_data" \
    --species invicta \
    --force \
    --sif "/path/to/rna_editing.sif" \
    --bind "/path/to/project"
```

**Parameters:**

| Parameter | Description |
|---|---|
| `--reference_dir` | Directory containing the reference genome FASTA and GTF files. |
| `--output_dir` | Output directory for prepared SnpEff input files. |
| `--species` | Species identifier used to locate the corresponding reference files. |
| `--force` | Forces regeneration of existing preparation results. |
| `--sif` | Path to the Apptainer container image; may be omitted if `RNA_EDITING_SIF` is set. |
| `--bind` | Additional host directory to mount inside the container; may be repeated. |

**Output:**

```text
snpeff_data/
├── invicta.fa
└── invicta_longestCDS.gtf
```

- `invicta.fa`: Symbolic link to the original reference genome FASTA.
- `invicta_longestCDS.gtf`: Representative GTF retaining the longest-CDS transcript for each coding gene and the transcript with the longest spliced exon length for each noncoding gene.

### 3.1.3. Build SnpEff Database

The `snpeff_build` command constructs a species-specific SnpEff database using the prepared longest-CDS GTF annotation and reference genome FASTA.

**Input Data:**

- **GTF file:** Longest-CDS annotation generated by `prepare_snpeff_gtf`.
- **Genome FASTA:** Reference genome FASTA prepared in the previous step.

**Command:**

```bash
KOBE snpeff_build \
    --gtf "/path/to/snpeff_data/invicta_longestCDS.gtf" \
    --genome "/path/to/snpeff_data/invicta.fa" \
    --outdir "/path/to/snpeff_result" \
    --sif "/path/to/rna_editing.sif" \
    --bind "/path/to/project"
```

**Parameters:**

| Parameter | Description |
|---|---|
| `--gtf` | Longest-CDS GTF annotation file. |
| `--genome` | Reference genome FASTA file. |
| `--outdir` | Output directory for the SnpEff database. |
| `--species` | Optional species identifier. If omitted, KOBE infers it from the GTF and FASTA filenames and checks that the prefixes match. |
| `--cores` | Number of CPU cores for database construction (default: 6). |
| `--force` | Rebuilds an existing species database. Existing databases are protected by default. |
| `--sif` | Path to the Apptainer container image; may be omitted if `RNA_EDITING_SIF` is set. |
| `--bind` | Additional host directory to mount inside the container; may be repeated. |

**Important Outputs:**

```text
snpeff_result/
├── config/
│   └── invicta.config
├── logs/
│   └── invicta.build.log
└── database/
    └── invicta/
        ├── genes.gtf
        ├── sequences.fa
        ├── sequence.bin
        ├── sequence.NC_052664.1.bin
        ├── ...
        └── snpEffectPredictor.bin
```

- `config/invicta.config`: SnpEff configuration file.
- `logs/invicta.build.log`: Database construction log.
- `database/invicta/`: Species-specific SnpEff database files used for subsequent functional annotation.

The completed SnpEff database is used in the RNA editing site identification and annotation step.

---

## 3.2. Reads Mapping

The `denovo_detect_map` command performs quality control and genome alignment of paired-end RNA-seq and DNA-seq reads.

### 3.2.1. Input Data

The following inputs are required:

| Input | Format | Description |
|---|---|---|
| Sample information | TSV | `samples.tsv` containing sample IDs, sample names, species, sequencing types, and FASTQ paths. |
| Reference genome | FASTA | Species-specific reference genome. |
| Gene annotation | GTF | Full genome annotation corresponding to the reference genome, used by STAR for RNA-seq mapping and as input to prepare_snpeff_gtf for SnpEff input preparation. |
| Sequencing reads | FASTQ | Paired-end RNA-seq and DNA-seq reads. |

**Sample Information Table**

Sample information must be provided in a tab-delimited file (`samples.tsv`) with the following six columns:

| Column | Description |
|---|---|
| `sample_id` | Unique sample identifier. |
| `sample_name` | Biological sample or individual name. |
| `species` | Species identifier. |
| `type` | Sequencing data type: `RNA` or `DNA`. |
| `fq1` | Path to the forward FASTQ file. |
| `fq2` | Path to the reverse FASTQ file. |

Example:

```text
sample_id	sample_name	species	type	fq1	fq2
H3_worker_1	H3_worker	invicta	DNA	/path/to/denovo_detect_data/DNA/H3_worker_body_R1.fq.gz	/path/to/denovo_detect_data/DNA/H3_worker_body_R2.fq.gz
N1_worker_2	N1_worker	invicta	DNA	/path/to/denovo_detect_data/DNA/N1_worker_body_R1.fq.gz	/path/to/denovo_detect_data/DNA/N1_worker_body_R2.fq.gz
H3_worker_3	H3_worker	invicta	RNA	/path/to/denovo_detect_data/RNA/H3_worker_head_R1.fq.gz	/path/to/denovo_detect_data/RNA/H3_worker_head_R2.fq.gz
N1_worker_4	N1_worker	invicta	RNA	/path/to/denovo_detect_data/RNA/N1_worker_head_R1.fq.gz	/path/to/denovo_detect_data/RNA/N1_worker_head_R2.fq.gz
```

The `type` column distinguishes RNA and DNA libraries. The `sample_id` uniquely identifies each sequencing library, with each RNA library evaluated independently, whereas `sample_name` identifies the biological sample or individual, grouping corresponding DNA libraries for background validation.

**Input Directory Structure**

```text
denovo_detect_data/
├── samples.tsv
├── DNA/
│   ├── H3_worker_body_R1.fq.gz
│   ├── H3_worker_body_R2.fq.gz
│   ├── N1_worker_body_R1.fq.gz
│   └── N1_worker_body_R2.fq.gz
├── RNA/
│   ├── H3_worker_head_R1.fq.gz
│   ├── H3_worker_head_R2.fq.gz
│   ├── N1_worker_head_R1.fq.gz
│   └── N1_worker_head_R2.fq.gz
└── ref/
    ├── invicta.fa
    └── invicta.gtf
```

Reference files are matched using the `species` column in `samples.tsv`. For example, `invicta` is matched to `invicta.fa` and `invicta.gtf`.

### 3.2.2. Command

```bash
KOBE denovo_detect_map \
    --samples "/path/to/denovo_detect_data/samples.tsv" \
    --ref-dir "/path/to/denovo_detect_data/ref" \
    --out-dir "/path/to/denovo_detect_result" \
    --cores 30 \
    --fastp-length-required 150 \
    --sif "/path/to/rna_editing.sif" \
    --bind "/path/to/project"
```

### 3.2.3. Parameters

| Parameter | Description |
|---|---|
| `--samples` | Sample information table (`samples.tsv`). |
| `--ref-dir` | Directory containing species-specific FASTA and GTF files. |
| `--out-dir` | Output working directory for mapping results. Created automatically if it does not exist. |
| `--cores` | Total number of CPU cores available to Snakemake. |
| `--fastp-length-required` | Minimum read length required after fastp processing; passed to fastp as `-l` (default: 150). Use `none` or `off` to omit the explicit length filter. |
| `--sif` | Path to the KOBE Apptainer container image; may be omitted if `RNA_EDITING_SIF` is set. |
| `--bind` | Additional host directory to mount inside the container; may be repeated. |

### 3.2.4. Processing and Important Outputs

The mapping workflow includes the following major steps:

1. **Quality Control:** Uses fastp to process both RNA-seq and DNA-seq reads, generating cleaned FASTQ files and QC reports.
2. **Reference Preparation:** Generates the required genome indexes for BWA and STAR.
3. **RNA-seq Mapping:** Uses BWA and STAR for alignment, followed by alignment refinement, duplicate removal, and GATK SplitNCigarReads processing.
4. **DNA-seq Mapping:** Uses BWA for alignment, followed by Sambamba sorting, indexing, and duplicate removal.

The main output directory is organized as follows:

```text
denovo_detect_result/
├── 00_QC/
│   └── invicta/
│       ├── DNA/
│       │   ├── H3_worker_1/
│       │   └── N1_worker_2/
│       └── RNA/
│           ├── H3_worker_3/
│           └── N1_worker_4/
│
├── 00_REF/
│   ├── GENOME/
│   ├── BWA/
│   └── STAR/
│
├── 01_MAP/
│   └── invicta/
│       ├── RNA/
│       │   ├── BWA/
│       │   ├── STAR/
│       │   └── final/
│       └── DNA/
│
├── tmp/
└── denovo_detect_map.log
```

Important outputs include:

| Output | Description |
|---|---|
| `00_QC/` | Cleaned FASTQ files (`qc.R1.fastq.gz` and `qc.R2.fastq.gz`) and fastp QC reports (`qc.html` and `qc.json`). |
| `00_REF/` | Reference genome files and BWA/STAR indexes. |
| `01_MAP/` | RNA-seq and DNA-seq alignment results. |
| `denovo_detect_map.log` | Mapping workflow execution log. |

**Final BAM files:**

- RNA-seq:

```text
01_MAP/invicta/RNA/final/<sample_id>/split.bam
01_MAP/invicta/RNA/final/<sample_id>/split.bam.bai
```

- DNA-seq:

```text
01_MAP/invicta/DNA/<sample_id>/<sample_id>.nodup.bam
01_MAP/invicta/DNA/<sample_id>/<sample_id>.nodup.bam.bai
```

These final BAM files are used as input for the subsequent RNA editing site identification and annotation step.

---

## 3.3. RNA Editing Site Identification and Annotation

The `denovo_detect_anno` command performs candidate RNA editing site identification and functional annotation using the mapping results and the species-specific SnpEff database.

The workflow includes mpileup-based variant calling, VCF normalization, SnpEff functional annotation, candidate site extraction, statistical filtering, and final site export.

### 3.3.1. Input Data

| Input | Format | Description |
|---|---|---|
| Mapping results | BAM | RNA-seq and DNA-seq BAM files generated by `denovo_detect_map`. |
| Reference files | FASTA/FAI | Reference genome and index files generated during Mapping. |
| Sample information | TSV | The same `samples.tsv` used in the Mapping step. |
| SnpEff database | Database directory | Species-specific SnpEff database generated by `snpeff_build`. |

The `--work-dir` parameter must point to the same working directory specified by `denovo_detect_map`.

The SnpEff database directory should be the output directory generated by `snpeff_build`.

### 3.3.2. Command

```bash
nohup setsid --fork --wait KOBE denovo_detect_anno \
    --work-dir "/path/to/denovo_detect_result" \
    --samples "/path/to/denovo_detect_data/samples.tsv" \
    --snpeff-dir "/path/to/snpeff_result" \
    --species invicta \
    --dna-min-coverage 10 \
    --dna-fdr 0.05 \
    --rna-min-coverage 3 \
    --rna-fdr 0.01 \
    --min-support 2 \
    --total-error-rate 0.005 \
    --cores 20 \
    --sif "/path/to/rna_editing.sif" \
    --bind "/path/to/project" \
    > "/path/to/denovo_detect_result/denovo_detect_anno.nohup.log" 2>&1 < /dev/null &
```

### 3.3.3. Parameters

| Parameter | Description |
|---|---|
| `--work-dir` | Mapping working directory generated by `denovo_detect_map`. |
| `--samples` | Sample information table used in the Mapping step. |
| `--snpeff-dir` | Root directory of the prebuilt SnpEff databases, containing database/ and config/ subdirectories. Required when SnpEff annotation is needed, such as during the first annotation run or when --reannotate is specified. Not required when existing annotated VCFs can be reused or when --filter-only is used. |
| `--species` | Species identifiers to be analyzed, corresponding to the species column in samples.tsv. Accepts comma-separated values and can be repeated. Only specified species are processed. If omitted, all unique species in samples.tsv are analyzed. |
| `--dna-min-coverage` | Minimum DNA sequencing coverage required for candidate site filtering (default: 10). |
| `--dna-fdr` | FDR threshold for DNA-based filtering (default: 0.05). |
| `--rna-min-coverage` | Minimum RNA sequencing coverage required for candidate site filtering (default: 5). |
| `--rna-fdr` | FDR threshold for RNA-based filtering (default: 0.05). |
| `--min-support` | Minimum number of independent RNA-seq libraries supporting an editing site (default: 2). |
| `--total-error-rate` | Total sequencing error rate used during statistical filtering (default: 0.005). |
| `--cores` | Total number of CPU cores available to Snakemake. |
| `--filter-only` | Reuse existing candidate tables and run only threshold filtering and final export; does not require `--snpeff-dir`. |
| `--reannotate` | Force SnpEff annotation and candidate extraction before downstream filtering; requires `--snpeff-dir`. |
| `--sif` | Path to the KOBE Apptainer container image; may be omitted if `RNA_EDITING_SIF` is set. |
| `--bind` | Additional host directory to mount inside the container; may be repeated. |

### 3.3.4. Important Outputs

The annotation workflow generates intermediate results under `02_anno/` and final editing site results under `edit_site_result/`.

Example directory structure:

```text
denovo_detect_result/
├── 00_QC/
├── 00_REF/
├── 01_MAP/
│
├── 02_anno/
│   └── invicta/
│       ├── DNA_RNA.vcf30.vcf
│       ├── DNA_RNA.var.vcf
│       ├── DNA_RNA.var.snpeff.vcf
│       ├── DNA_RNA.snpeff_summary.html
│       ├── genic.candidates.tsv
│       ├── intergenic.candidates.tsv
│       └── logs/
│
├── edit_site_result/
│   └── invicta/
│       ├── genic.final_sites.DNA10_RNA3_NRNA2_DNAFDR0.05_RNAFDR0.01.tsv
│       ├── intergenic.final_sites.DNA10_RNA3_NRNA2_DNAFDR0.05_RNAFDR0.01.tsv
│       ├── all.final_sites.DNA10_RNA3_NRNA2_DNAFDR0.05_RNAFDR0.01.tsv
│       ├── genic.DNA10_RNA3_NRNA2_DNAFDR0.05_RNAFDR0.01.mutation_types.pdf
│       ├── intergenic.DNA10_RNA3_NRNA2_DNAFDR0.05_RNAFDR0.01.mutation_types.pdf
│       └── variant_type_proportion.DNA10_RNA3_NRNA2_DNAFDR0.05_RNAFDR0.01.tsv
│
├── denovo_detect_map.log
├── denovo_detect_anno.log
└── denovo_detect_anno.nohup.log
```

**Important intermediate files:**

| Output | Description |
|---|---|
| `DNA_RNA.vcf30.vcf` | Variant calls generated from mpileup. |
| `DNA_RNA.var.vcf` | Normalized VCF file. |
| `DNA_RNA.var.snpeff.vcf` | SnpEff-annotated VCF file. |
| `DNA_RNA.snpeff_summary.html` | SnpEff annotation summary report. |
| `genic.candidates.tsv` | Extracted genic candidate sites before statistical filtering. |
| `intergenic.candidates.tsv` | Extracted intergenic candidate sites before statistical filtering. |

**Final editing site tables:**

| Output | Description |
|---|---|
| `genic.final_sites.*.tsv` | Final filtered genic RNA editing sites with functional annotation and sample-level read information. |
| `intergenic.final_sites.*.tsv` | Final filtered intergenic RNA editing sites. |
| `all.final_sites.*.tsv` | Combined final genic and intergenic editing site results. |
| `*.mutation_types.pdf` | Mutation-type distribution plots for genic and intergenic regions. |
| `variant_type_proportion.*.tsv` | Summary of mutation-type proportions. |

The output filename includes a threshold tag representing the DNA coverage, RNA coverage, minimum RNA sample support, DNA FDR, and RNA FDR settings used in the analysis.

For example:

```text
DNA10_RNA3_NRNA2_DNAFDR0.05_RNAFDR0.01
```

The final genic editing site table can be used as input for downstream orthology analysis, while the combined table provides an overview of the final genic and intergenic editing sites.

---

The Mapping and SnpEff database construction steps can be performed independently. Both must be completed before running `denovo_detect_anno`, which requires the mapping results and the corresponding SnpEff database.

---

# 4. Orthology-based RNA Editing Analysis

## 4.1. Step 1. Orthology-based RNA Editing Analysis

### 4.1.1. Input Data

Two types of input data are required.

**1. Longest-CDS FASTA files**

The `data/` directory contains the longest CDS sequences for each species in FASTA format.

Example directory:

```text
data/
├── Dip_Drosophila_melanogaster_longestcds.fasta
├── Hym_Acromyrmex_echinatior_longestcds.fasta
├── Hym_Apis_mellifera_longestcds.fasta
└── Hym_Solenopsis_invicta_longestcds.fasta
```

Each FASTA file contains the longest CDS sequences for the corresponding species.

**2. Candidate CDS RNA editing site files**

The `edit_txt/` directory contains candidate CDS RNA editing sites for each species in tab-delimited text format (`.txt`).

Example directory:

```text
edit_txt/
├── Dip_Drosophila_melanogaster_candi_cds_sites.txt
├── Hym_Acromyrmex_echinatior_candi_cds_sites.txt
├── Hym_Apis_mellifera_candi_cds_sites.txt
└── Hym_Solenopsis_invicta_candi_cds_sites.txt
```

Each file contains six columns:

| Column | Description |
|---|---|
| `transcript` | Transcript identifier |
| `cds_pos` | Editing site position within the CDS |
| `anno` | Functional annotation |
| `AA_pos` | Amino acid position |
| `AA_change` | Amino acid substitution |
| `gene_id` | Gene identifier |

Example:

```text
transcript	cds_pos	anno	AA_pos	AA_change	gene_id
AECH19648-RA	35	missense_variant	12	Lys>Arg	AECH19648
AECH19648-RA	27	missense_variant	9	Ile>Met	AECH19648
AECH19648-RA	25	missense_variant	9	Ile>Val	AECH19648
```

Each row represents one candidate CDS RNA editing site.

The FASTA and candidate site files must correspond to the same species.

### 4.1.2. Command

```bash
nohup setsid --fork --wait KOBE ortho_edit \
    --cds_dir "/path/to/ortho_edit_data/data" \
    --edit_txt_dir "/path/to/ortho_edit_data/edit_txt" \
    --work_dir "/path/to/ortho_result" \
    --threads 20 \
    --sif "/path/to/rna_editing.sif" \
    --bind "/path/to/project" \
    > "/path/to/ortho_edit.log" 2>&1 < /dev/null &
```

### 4.1.3. Parameters

| Parameter | Description |
|---|---|
| `--cds_dir` | Directory containing longest-CDS FASTA files for all species. |
| `--edit_txt_dir` | Directory containing candidate CDS RNA editing site files for all species. |
| `--work_dir` | Output directory for orthology-based analysis. |
| `--threads` | CPU thread budget for the analysis (default: 8). |
| `--sif` | Path to the KOBE Apptainer container image; may be omitted if `RNA_EDITING_SIF` is set. |
| `--bind` | Additional Apptainer bind mount; may be repeated. |

### 4.1.4. Important Outputs

The main results include the following files and directories:

```text
ortho_result/
├── pep/
├── spe_cds_pos/
└── ortho_edit_result/
    ├── merged_edit_sites_codon_table_annotated.csv
    └── candi_per_species/
```

| Output | Description |
|---|---|
| `merged_edit_sites_codon_table_annotated.csv` | Merged and annotated orthology-based RNA editing site results. Used as input for subsequent conservation visualization and novel site integration. |
| `pep/` | Species-specific protein FASTA files generated from the longest-CDS sequences; can be used as input for optional OrthoFinder species-tree construction. |
| `spe_cds_pos/` | Species-specific CDS position information used in the site recovery step. |
| `candi_per_species/` | Candidate editing sites organized by species, used in the site recovery step. |

---

## 4.2. Step 2. Conservation Visualization

### 4.2.1. Input Data

This step requires the annotated results generated by Step 1 and reuses the longest-CDS FASTA files provided as input to Step 1. Species-specific GTF annotation files can also be optionally provided to retrieve the genomic coordinates of orthologous RNA editing sites in the corresponding species.

**1. Annotated orthology results**

The `--annotated` parameter specifies the merged and annotated CSV file generated by Step 1.

```text
ortho_result/
└── ortho_edit_result/
    └── merged_edit_sites_codon_table_annotated.csv
```

**2. Longest-CDS FASTA files**

The `--cds_dir` parameter specifies the directory containing the same longest-CDS FASTA files used as input to Step 1.

**3. Species-specific GTF annotation files (optional)**

The `--gtf` parameter accepts species-specific GTF annotation files in the format `species=GTF_file_path`.

When provided, these files are used to obtain the corresponding genomic coordinates of orthologous RNA editing sites in each species. GTF files can be supplied for individual species as needed.

### 4.2.2. Command

```bash
nohup setsid --fork --wait KOBE conserve_plot \
    --annotated "/path/to/ortho_result/ortho_edit_result/merged_edit_sites_codon_table_annotated.csv" \
    --cds_dir "/path/to/ortho_edit_data/data" \
    --min_species 2 \
    --species_order "melanogaster,mellifera,echinatior,invicta" \
    --source_priority "melanogaster,mellifera,echinatior,invicta" \
    --outdir "/path/to/conserve_plot_result" \
    --sif "/path/to/rna_editing.sif" \
    --gtf "invicta=/path/to/ortho_edit_data/edit_txt/gtf/invicta.gtf" \
    --gtf "mellifera=/path/to/ortho_edit_data/edit_txt/gtf/mellifera.gtf" \
    --col 8 \
    --bind "/path/to/project" \
    > "/path/to/conserve_plot.log" 2>&1 < /dev/null &
```

### 4.2.3. Parameters

| Parameter | Description |
|---|---|
| `--annotated` | Path to the annotated CSV file generated by `ortho_edit`. |
| `--cds_dir` | Directory containing the longest-CDS FASTA files for all species. |
| `--min_species` | Minimum number of species with recorded RNA editing sites required to define a conserved editing site (must be at least 2). For example, `2` requires records in at least two species. |
| `--species_order` | Comma-separated species identifiers specifying the plotting order, from top to bottom. |
| `--source_priority` | Optional comma-separated species priority order used to select representative transcripts. If omitted, the annotation column order is used. |
| `--outdir` | Output directory for conservation visualization results. |
| `--sif` | Path to the KOBE Apptainer container image; may be omitted if `RNA_EDITING_SIF` is set. |
| `--gtf` | Optional species-specific GTF annotation file, specified as `species=path`. Can be supplied multiple times; supplied GTFs are used to add genomic coordinates to the exported conserved-site table. |
| `--col` | Number of columns per row in the generated PDF figure (default: 8). |
| `--bind` | Additional Apptainer bind mount; may be repeated. |

### 4.2.4. Important Outputs

The main visualization result is:

```text
conserve_plot_result/
└── conserved_sites_ge2species_all_transcripts.pdf
```

`conserved_sites_ge2species_all_transcripts.pdf` presents the conserved RNA editing sites identified across species, with the minimum species requirement specified by `--min_species`.

---

## 4.3. Step 3. Optional Species Tree Construction

If a species tree is already available, this step can be skipped. Otherwise, use the species-specific protein FASTA files generated by `ortho_edit` in `ortho_result/pep` as input for OrthoFinder 3.1.5.

### 4.3.1. Download the OrthoFinder 3.1.5 SIF

The OrthoFinder 3.1.5 container is available from BioContainers on Quay.io:

```text
quay.io/biocontainers/orthofinder:3.1.5--hdfd78af_0
```

Download the image and save it as a SIF file:

```bash
conda activate apptainer-1.5.4
apptainer pull \
    "/path/to/software/orthofinder_3.1.5.sif" \
    "docker://quay.io/biocontainers/orthofinder:3.1.5--hdfd78af_0"
```

### 4.3.2. Run OrthoFinder

The following command uses DIAMOND for protein sequence searches and the MSA-based workflow for species-tree inference:

```bash
mkdir -p "/path/to/project/orthofinder_tree_result"

nohup setsid --fork --wait apptainer exec \
    --bind "/path/to/project:/path/to/project" \
    "/path/to/software/orthofinder_3.1.5.sif" \
    orthofinder \
    -f "/path/to/project/ortho_result/pep" \
    -S diamond \
    -t 8 \
    -a 4 \
    -M msa \
    -o "/path/to/project/orthofinder_tree_result/orthofinder_results" \
    > "/path/to/project/orthofinder_tree_result/orthofinder.log" \
    2>&1 < /dev/null &
```

### 4.3.3. Output

After OrthoFinder finishes, the rooted species tree is generated in Newick format at:

```text
orthofinder_tree_result/
└── orthofinder_results/
    └── Results_*/
        └── Species_Tree/
            └── SpeciesTree_rooted.txt
```

`SpeciesTree_rooted.txt` is the rooted tree file generated by OrthoFinder. The command above does not specify an external outgroup.

---

## 4.4. Step 4. Orthology-guided Site Recovery

### 4.4.1. Input Data

This step requires a sample information table, RNA-seq FASTQ files, and several input files and directories generated or prepared during the previous orthology analysis.

**1. Sample Information Table**

Sample information must be provided in a tab-delimited file (`samples.tsv`) with the following columns:

| Column | Description |
|---|---|
| `sample_id` | Unique identifier for each sequencing sample |
| `individual` | Biological individual identifier |
| `species` | Species identifier |
| `type` | Sequencing data type; this module currently accepts `RNA` samples only. |
| `fq1` | Path to the forward FASTQ file |
| `fq2` | Path to the reverse FASTQ file |
| `ref` | Path to the corresponding species-specific longest-CDS FASTA file |

Example:

```text
sample_id	individual	species	type	fq1	fq2	ref
Amel_N1_worker_head2	N1	mellifera	RNA	/path/to/find_back_data/Amel_N1_worker_head2_1.fq.gz	/path/to/find_back_data/Amel_N1_worker_head2_2.fq.gz	/path/to/ortho_edit_data/data/Hym_Apis_mellifera_longestcds.fasta
Amel_N2_worker_head2	N2	mellifera	RNA	/path/to/find_back_data/Amel_N2_worker_head2_1.fq.gz	/path/to/find_back_data/Amel_N2_worker_head2_2.fq.gz	/path/to/ortho_edit_data/data/Hym_Apis_mellifera_longestcds.fasta
Aech_C1_worker_head2	C1	echinatior	RNA	/path/to/find_back_data/Aech_C1_worker_head2_1.fq.gz	/path/to/find_back_data/Aech_C1_worker_head2_2.fq.gz	/path/to/ortho_edit_data/data/Hym_Acromyrmex_echinatior_longestcds.fasta
Aech_C2_worker_head2	C2	echinatior	RNA	/path/to/find_back_data/Aech_C2_worker_head2_1.fq.gz	/path/to/find_back_data/Aech_C2_worker_head2_2.fq.gz	/path/to/ortho_edit_data/data/Hym_Acromyrmex_echinatior_longestcds.fasta
Sinv_H3_worker_3	H3_worker	invicta	RNA	/path/to/find_back_data/Sinv_H3_worker_3_1.fq.gz	/path/to/find_back_data/Sinv_H3_worker_3_2.fq.gz	/path/to/ortho_edit_data/data/Hym_Solenopsis_invicta_longestcds.fasta
Sinv_N1_worker_4	N1_worker	invicta	RNA	/path/to/find_back_data/Sinv_N1_worker_4_1.fq.gz	/path/to/find_back_data/Sinv_N1_worker_4_2.fq.gz	/path/to/ortho_edit_data/data/Hym_Solenopsis_invicta_longestcds.fasta
```

The directory containing each `ref` FASTA must be writable because `find_back` creates BWA and FASTA index files beside the reference.

**2. Orthology Analysis Results**

The following directories are generated by Step 1 (`ortho_edit`):

| Input parameter | Directory | Description |
|---|---|---|
| `--cds_pos_dir` | `ortho_result/spe_cds_pos/` | Species-specific CDS position information. |
| `--candi_dir` | `ortho_result/ortho_edit_result/candi_per_species/` | Candidate editing sites organized by species. |

**3. Species-specific Annotation and Editing Site Files**

The following directories contain the corresponding species-specific input files:

| Input parameter | Directory | Description |
|---|---|---|
| `--edit_txt_dir` | `ortho_edit_data/edit_txt/` | Candidate CDS RNA editing site files used in Step 1. |
| `--gtf` | `ortho_edit_data/edit_txt/gtf/<species>.gtf` | Optional species-specific GTF files, passed individually as `--gtf species=FILE`. |

Pass each optional GTF file explicitly with `--gtf species=FILE`; the species identifier must match the corresponding value in the `species` column of `samples.tsv`.

When a GTF file is provided for a species, genomic coordinates of the corresponding orthologous RNA editing sites can also be generated for that species.

### 4.4.2. Command

```bash
nohup setsid --fork --wait KOBE find_back \
    --sample_table "/path/to/find_back_data/samples.tsv" \
    --out_dir "/path/to/find_back_result" \
    --cds_pos_dir "/path/to/ortho_result/spe_cds_pos" \
    --gtf "species1=/path/to/gtf/species1.gtf" \
    --gtf "species2=/path/to/gtf/species2.gtf" \
    --edit_txt_dir "/path/to/ortho_edit_data/edit_txt" \
    --candi_dir "/path/to/ortho_result/ortho_edit_result/candi_per_species" \
    --min_samples 1 \
    --threads 20 \
    --min_length 150 \
    --mapping_jobs 4 \
    --sif "/path/to/rna_editing.sif" \
    --bind "/path/to/project" \
    > "/path/to/find_back.log" 2>&1 < /dev/null &
```

### 4.4.3. Parameters

| Parameter | Description |
|---|---|
| `--sample_table` | Sample information table containing sample IDs, individual IDs, species, RNA-only sequencing type, FASTQ paths, and reference CDS paths. |
| `--out_dir` | Output directory for site recovery results. |
| `--cds_pos_dir` | Directory containing species-specific CDS position information generated by `ortho_edit`. |
| `--gtf` | Species-specific GTF annotation file, supplied as species=path. This parameter can be repeated for multiple species. Species without a supplied GTF will skip genomic coordinate conversion. |
| `--edit_txt_dir` | Directory containing the original species-specific RNA editing site files used as known-site references for Known/Novel classification. |
| `--candi_dir` | Directory containing species-specific candidate editing sites generated by `ortho_edit`. |
| `--min_samples` | Minimum number of RNA samples required to support a candidate editing site (default: 1). |
| `--mapping_jobs` | Maximum number of sample-level Mapping workflows running concurrently, including fastp, BWA alignment, BAM sorting, indexing, and Picard duplicate removal (default: 2). A new sample starts when any active sample workflow finishes. |
| `--threads` | Total CPU-thread budget for REDItoolKnown processing. Each sample job uses 2 threads, so the maximum concurrent jobs are `threads / 2`. Must be an even integer of at least 2 (default: 20). |
| `--min_length` | Minimum read length passed to fastp as `-l` (default: 150). |
| `--resume` | Reuse completed Mapping/Picard and REDItoolKnown results when their signatures and completeness checks are valid. |
| `--force_reditool` | Reuse compatible in-workdir Picard BAMs when available, but rerun candidate-site preparation, REDItoolKnown, and downstream site-recovery steps; if a BAM is not reusable, the required Mapping steps are rerun. |
| `--sif` | Path to the KOBE Apptainer container image; may be omitted if `RNA_EDITING_SIF` is set. |
| `--bind` | Additional Apptainer bind specification; may be repeated. |

For example, `--min_samples 1` requires a candidate site to have support in at least one sample.

### 4.4.4. Important Outputs

The main results are organized by species under the `statistics/` directory. A summary table covering all analyzed species is also generated.

**Output directory structure**

```text
find_back_result/
└── statistics/
    ├── editing_site_statistics_ge1_all_species.tsv
    │
    ├── echinatior/
    │   ├── echinatior_editing_site_statistics_ge1.tsv
    │   ├── echinatior_supported_sites_ge1.tsv
    │   └── echinatior_novel_supported_sites_ge1.tsv
    │
    ├── invicta/
    │   ├── invicta_editing_site_statistics_ge1.tsv
    │   ├── invicta_supported_sites_ge1.tsv
    │   └── invicta_novel_supported_sites_ge1.tsv
    │
    └── mellifera/
        ├── mellifera_editing_site_statistics_ge1.tsv
        ├── mellifera_supported_sites_ge1.tsv
        └── mellifera_novel_supported_sites_ge1.tsv
```

**1. All-species summary**

The file `editing_site_statistics_ge1_all_species.tsv` summarizes RNA editing site recovery statistics across all analyzed species. Each row represents one species.

| Column | Description |
|---|---|
| `species` | Species identifier. |
| `n_RNA_samples` | Number of RNA-seq samples analyzed for the species. |
| `min_samples_required` | Minimum number of samples required to support an editing site. |
| `n_summary_unique_sites` | Total number of unique candidate editing sites in the all-sample summary. |
| `n_known_supported_sites_ge_threshold` | Number of known editing sites meeting the support threshold. |
| `n_novel_supported_sites_ge_threshold` | Number of novel editing sites meeting the support threshold. |
| `n_original_known_unique_sites` | Total number of unique original known editing sites. |
| `n_original_known_sites_below_threshold` | Number of original known editing sites that do not meet the support threshold. |
| `known_site_file` | Path to the original known editing site file. |
| `summary_file` | Path to the all-sample candidate site summary file. |

**2. Species-specific statistics**

Each species has an individual statistics file:

`*_editing_site_statistics_ge1.tsv`

This file provides the site recovery statistics for the corresponding species.

**3. Supported editing sites**

The file `*_supported_sites_ge1.tsv` contains editing sites meeting the specified minimum sample support requirement.

The result table includes the following columns:

| Column | Description |
|---|---|
| `transcript` | Transcript identifier. |
| `cds_pos` | Editing site position within the CDS. |
| `anno` | Functional annotation (`Syn` or `Nonsyn`). |
| `codon` | Original codon sequence. |
| `AA_before` | Amino acid before editing. |
| `codon_pos` | Position of the edited nucleotide within the codon. |
| `edited_codon` | Codon sequence after editing. |
| `AA_after` | Amino acid after editing. |
| `handle` | Site identifier used to track the corresponding editing event. |
| `<sample_id>_A_reads` | Number of A reads supporting the site in the sample. |
| `<sample_id>_G_reads` | Number of G reads supporting the site in the sample. |
| `<sample_id>_editing_level` | RNA editing level in the sample. |
| `n_samples_editing_level_gt0` | Number of samples with an editing level greater than zero. |
| `known_status` | Classification of the site as `Known` or `Novel`. |

The three sample-specific columns are generated for each sample according to the sample IDs in `samples.tsv`.

Example:

```text
transcript	cds_pos	anno	codon	AA_before	codon_pos	edited_codon	AA_after	handle	Amel_N1_worker_head2_A_reads	Amel_N1_worker_head2_G_reads	Amel_N1_worker_head2_editing_level	Amel_N2_worker_head2_A_reads	Amel_N2_worker_head2_G_reads	Amel_N2_worker_head2_editing_level	n_samples_editing_level_gt0	known_status
XM_026442151.1	1168	Nonsyn	AGA	Arg	1	GGA	Gly	XM_026442151.1:1168:Arg>Gly	193	10	0.04926108	333	4	0.01186944	2	Known
```

**4. Novel supported editing sites**

The file `*_novel_supported_sites_ge1.tsv` contains newly recovered editing sites that are not included in the original known editing site set.

Its table format is the same as `*_supported_sites_ge1.tsv`.

Example:

```text
transcript	cds_pos	anno	codon	AA_before	codon_pos	edited_codon	AA_after	handle	Amel_N1_worker_head2_A_reads	Amel_N1_worker_head2_G_reads	Amel_N1_worker_head2_editing_level	Amel_N2_worker_head2_A_reads	Amel_N2_worker_head2_G_reads	Amel_N2_worker_head2_editing_level	n_samples_editing_level_gt0	known_status
XM_006567404.3	495	Syn	AGA	Arg	3	AGG	Arg	Sinv000964.1:450:Leu>Leu	681	0	0.00000000	681	3	0.00438596	1	Novel

```

The `*_novel_supported_sites_ge1.tsv` file excludes sites already present in the original known editing site set.

The suffix `ge1` indicates that the results use a minimum sample support threshold of at least one sample when `--min_samples 1` is specified.

---

## 4.5. Step 5. Integration of Newly Recovered Sites

### 4.5.1. Input Data

This step requires the annotated orthology results generated by Step 1 and the novel supported editing site files generated by Step 4.

**1. Annotated orthology results**

The `--annotated` parameter specifies the annotated CSV file generated by the `ortho_edit` module.

```text
ortho_result/
└── ortho_edit_result/
    └── merged_edit_sites_codon_table_annotated.csv
```

**2. Novel supported editing site files**

The `--sites` parameter specifies the species-specific novel supported editing site files generated by the `find_back` module in Step 4.

Example:

```text
find_back_result/
└── statistics/
    ├── invicta/
    │   └── invicta_novel_supported_sites_ge1.tsv
    ├── mellifera/
    │   └── mellifera_novel_supported_sites_ge1.tsv
    └── echinatior/
        └── echinatior_novel_supported_sites_ge1.tsv
```

Multiple `--sites` parameters can be supplied to include novel editing site files from different species.

### 4.5.2. Command

```bash
nohup setsid --fork --wait KOBE add_new_site \
    --annotated "/path/to/ortho_result/ortho_edit_result/merged_edit_sites_codon_table_annotated.csv" \
    --sites "/path/to/find_back_result/statistics/invicta/invicta_novel_supported_sites_ge1.tsv" \
    --sites "/path/to/find_back_result/statistics/mellifera/mellifera_novel_supported_sites_ge1.tsv" \
    --sites "/path/to/find_back_result/statistics/echinatior/echinatior_novel_supported_sites_ge1.tsv" \
    --outdir "/path/to/find_back_plot_result" \
    --sif "/path/to/rna_editing.sif" \
    --bind "/path/to/project" \
    > "/path/to/add_new_site.log" 2>&1 < /dev/null &
```

### 4.5.3. Parameters

| Parameter | Description |
|---|---|
| `--annotated` | Path to the annotated CSV file generated by `ortho_edit`. |
| `--sites` | Path to a species-specific novel supported editing site file generated by `find_back`. This parameter can be supplied multiple times to include results from multiple species. |
| `--outdir` | Output directory for the integrated results. |
| `--sif` | Path to the KOBE Apptainer container image; may be omitted if `RNA_EDITING_SIF` is set. |
| `--bind` | Additional Apptainer bind mount; may be repeated. |

### 4.5.4. Important Outputs

The main integrated result is:

```text
find_back_plot_result/
└── merged_edit_sites_codon_table_annotated_with_new_sites.csv
```

`merged_edit_sites_codon_table_annotated_with_new_sites.csv` contains the annotated orthology results integrated with newly recovered editing sites. It serves as the input for the visualization step that follows.

---

## 4.6. Step 6. Visualization of Integrated Editing Sites

### 4.6.1. Input Data

This step uses the integrated annotated results generated by Step 5 and reuses the longest-CDS FASTA files provided as input to Step 1.

**1. Integrated annotated results**

The `--annotated` parameter specifies the updated CSV file generated by `add_new_site` in Step 5.

```text
find_back_plot_result/
└── merged_edit_sites_codon_table_annotated_with_new_sites.csv
```

**2. Longest-CDS FASTA files**

The `--cds_dir` parameter specifies the directory containing the longest-CDS FASTA files used in Step 1.

### 4.6.2. Command

```bash
KOBE find_back_plot \
    --annotated "/path/to/find_back_plot_result/merged_edit_sites_codon_table_annotated_with_new_sites.csv" \
    --cds_dir "/path/to/ortho_edit_data/data" \
    --group_col new_id2_group \
    --min_species 2 \
    --species_order "melanogaster,mellifera,echinatior,invicta" \
    --source_priority "melanogaster,mellifera,echinatior,invicta" \
    --outdir "/path/to/find_back_plot_result/plots" \
    --sif "/path/to/rna_editing.sif" \
    --col 8 \
    --bind "/path/to/project"
```

### 4.6.3. Parameters

| Parameter | Description |
|---|---|
| `--annotated` | Path to the integrated annotated CSV file generated by `add_new_site`. |
| `--cds_dir` | Directory containing the longest-CDS FASTA files for all species. |
| `--group_col` | Grouping column in the annotated CSV; accepted values are `ID2_group` and `new_id2_group`. |
| `--min_species` | Minimum number of species with editing records required for a group to be retained (minimum: 1). |
| `--species_order` | Comma-separated species identifiers specifying the plotting order, from top to bottom. |
| `--source_priority` | Optional comma-separated species priority order used to select representative transcripts. If omitted, annotation column order is used. |
| `--outdir` | Output directory for the integrated visualization results. |
| `--sif` | Path to the KOBE Apptainer container image; may be omitted if `RNA_EDITING_SIF` is set. |
| `--col` | Number of columns per row in the generated PDF figures (default: 8). |
| `--bind` | Additional Apptainer bind specification; may be repeated. |

For example, `--min_species 2` requires an editing record in at least two species for a site to be considered conserved.

### 4.6.4. Important Outputs

The main visualization result is:

```text
plots/
└── sup_conserved_sites_new_id2_group_ge2species_all_transcripts.pdf
```

This PDF is written directly under the specified `--outdir` and presents conserved RNA editing sites after integrating newly recovered sites, using the orthologous groups specified by `--group_col`.

---

# 5. Differential RNA Editing Analysis

## 5.1. Input Data

The differential editing analysis requires a tab-delimited editing count table (`editing_counts.tsv`).

Each row represents one RNA editing site. The table contains a site identifier (`ID`) and A/G read counts for each sample.

**Example:**

```text
ID	Sinv_H3_worker_3_A_reads	Sinv_H3_worker_3_G_reads	Sinv_N1_worker_4_A_reads	Sinv_N1_worker_4_G_reads
Sinv009841.1:301	16	50	12	57
```

| Column | Description |
|---|---|
| `ID` | Editing site identifier, consisting of the transcript ID and CDS position in the format `transcript:cds_pos`. |
| `<sample_id>_A_reads` | Number of A reads supporting the editing site in the sample. |
| `<sample_id>_G_reads` | Number of G reads supporting the editing site in the sample. |

Each sample must have a corresponding pair of A and G read count columns. Sample identifiers must match those used in the group definitions for differential analysis.

## 5.2. Preparation of Differential Editing Input

If the editing count table has not been prepared, KOBE provides the `supported_sites_to_editing_counts` command to convert the site recovery results generated by the `find_back` module into the required format.

Both `*_supported_sites_ge1.tsv` and `*_novel_supported_sites_ge1.tsv` can be used as input because they share the same table format.

- `*_supported_sites_ge1.tsv`: Includes all supported editing sites, including known and novel sites.
- `*_novel_supported_sites_ge1.tsv`: Includes only newly recovered editing sites.

**Command:**

```bash
KOBE supported_sites_to_editing_counts \
    --input "/path/to/find_back_result/statistics/invicta/invicta_supported_sites_ge1.tsv" \
    --species invicta \
    --work_dir "/path/to/diff_edit_site_data" \
    --sif "/path/to/rna_editing.sif" \
    --bind "/path/to/project"
```

**Parameters:**

| Parameter | Description |
|---|---|
| `--input` | Supported editing site table generated by `find_back`. |
| `--species` | Species identifier used for the output filename. |
| `--work_dir` | Output working directory. |
| `--sif` | Path to the KOBE Apptainer container image; may be omitted if `RNA_EDITING_SIF` is set. |
| `--bind` | Additional Apptainer bind mount; may be repeated. |

**Output:**

```text
diff_edit_site_data/
└── invicta_editing_counts.tsv
```

The generated `*_editing_counts.tsv` file can be directly used for differential editing analysis.

## 5.3. Pairwise Differential Editing Analysis

The `pairwise_tests` command performs pairwise differential RNA editing analysis between user-defined sample groups.

### 5.3.1. Comparison Between Two Groups

Example: comparing H3 and N1 groups in *Solenopsis invicta*.

```bash
KOBE pairwise_tests \
    --input "/path/to/diff_edit_site_data/invicta_editing_counts.tsv" \
    --group H3=Sinv_H3_worker_3 \
    --group N1=Sinv_N1_worker_4 \
    --comparisons H3:N1 \
    --min_coverage 3 \
    --fdr_cutoff 0.01 \
    --outdir "/path/to/diff_result" \
    --bind "/path/to/project" \
    --sif "/path/to/rna_editing.sif"
```

### 5.3.2. Multiple-group Comparisons

Multiple samples can be assigned to each group, and multiple pairwise comparisons can be specified in a single command.

Example: comparing three groups (N, F, and M) in *Apis mellifera*.

```bash
KOBE pairwise_tests \
    --input "/path/to/diff_edit_site_data/mellifera_editing_counts.tsv" \
    --group N=N1_worker_head2,N2_worker_head2,N3_worker_head2,N4_worker_head2 \
    --group F=F1_worker_head2,F2_worker_head2,F3_worker_head2,F4_worker_head2 \
    --group M=N5_worker_head2,N6_worker_head2,F5_worker_head2,F6_worker_head2 \
    --comparisons N:F,N:M,F:M \
    --min_coverage 3 \
    --fdr_cutoff 0.01 \
    --outdir "/path/to/diff_result" \
    --bind "/path/to/project" \
    --sif "/path/to/rna_editing.sif"
```

### 5.3.3. Parameters

| Parameter | Description |
|---|---|
| `--input` | Editing count table generated by `supported_sites_to_editing_counts`. |
| `--group` | Defines a sample group in the format `GROUP_NAME=SAMPLE1,SAMPLE2,...`. This parameter can be supplied multiple times. |
| `--comparisons` | Specifies pairwise comparisons between groups. Use `:` to separate two group names and commas to separate multiple comparisons. |
| `--min_coverage` | Minimum A + G read coverage required in an individual sample for a site to be included in the analysis (default: 3). |
| `--fdr_cutoff` | BH-adjusted false discovery rate (FDR) significance threshold (default: 0.05). |
| `--outdir` | Output directory for differential editing analysis results. |
| `--force` | Replace existing output directories for the requested comparisons. |
| `--bind` | Additional Apptainer bind mount; may be repeated. |
| `--sif` | Path to the KOBE Apptainer container image; may be omitted if `RNA_EDITING_SIF` is set. |

For example:

- `--group N=N1_worker_head2,N2_worker_head2` defines group N using two samples.
- `--comparisons N:F,N:M,F:M` specifies three pairwise comparisons.
- `--min_coverage 3` requires at least three A + G reads in an individual sample at a given site.
- `--fdr_cutoff 0.01` sets the FDR significance threshold to 0.01.

### 5.3.4. Important Outputs

Each pairwise comparison generates a separate result directory containing statistical results and visualization files.

For example, comparing N and F:

```text
diff_result/
└── N_vs_F/
    ├── N_vs_F_hypothesis_tests.tsv
    ├── N_vs_F_summary.tsv
    ├── N_vs_F_scatter.pdf
    └── N_vs_F_pooled_scatter.pdf
```

| Output file | Description |
|---|---|
| `*_hypothesis_tests.tsv` | Site-level statistical testing results, including editing levels, differences between groups, P-values, FDR values, and significance information. |
| `*_summary.tsv` | Summary statistics for the pairwise comparison. |
| `*_scatter.pdf` | Scatter plots comparing mean RNA editing levels between the two groups. |
| `*_pooled_scatter.pdf` | Scatter plots comparing pooled RNA editing levels between the two groups. |

Four statistical methods are included in the differential analysis:

- Student's t-test
- Fisher's exact test
- REDIT-LLR
- REDIT-Regression

For multiple-group comparisons, a separate result directory is generated for each specified pairwise comparison.

For example, three comparisons (`N:F`, `N:M`, and `F:M`) generate:

```text
diff_result/
├── N_vs_F/
├── N_vs_M/
└── F_vs_M/
```

Each comparison directory contains its own statistical testing table, summary table, and PDF visualizations.

---

# 6. Genomic Base Extraction

## 6.1. Input Data

The base extraction module retrieves nucleotide counts (A, C, G, and T) from sequencing data at specified genomic positions.

The following input data are required.

### 6.1.1. Sample Information Table

Sample information must be provided in a tab-delimited file (`samples.tsv`) containing the following columns:

| Column | Description |
|---|---|
| `sample_id` | Unique sample identifier. |
| `species` | Species identifier. |
| `fq1` | Path to the forward FASTQ file. |
| `fq2` | Path to the reverse FASTQ file. |
| `genome_fasta` | Path to the corresponding reference genome FASTA file. |

Example:

```text
sample_id	species	fq1	fq2	genome_fasta
H3_worker_3	invicta	/path/to/base_extract_data/H3_worker_head_R1.fq.gz	/path/to/base_extract_data/H3_worker_head_R2.fq.gz	/path/to/ref/invicta.fa
N1_worker_4	invicta	/path/to/base_extract_data/N1_worker_head_R1.fq.gz	/path/to/base_extract_data/N1_worker_head_R2.fq.gz	/path/to/ref/invicta.fa
```

### 6.1.2. Reference Annotation Files

The `--gtf-dir` parameter specifies the directory containing species-specific GTF annotation files.

GTF files are matched using species identifiers from the sample table and should be named as `<species>.gtf`.

Example:

```text
ref/
├── invicta.fa
└── invicta.gtf
```

The reference genome FASTA path is specified separately in the `genome_fasta` column of `samples.tsv`.

### 6.1.3. Known Editing Site Files

The `--site-dir` parameter specifies the directory containing species-specific known editing site files. Files are matched using species identifiers.

The input files must be tab-delimited and **without a header**. The first two columns specify the chromosome and genomic position (1-based). An optional third column specifies the strand.

**Example (without strand information):**

```text
NC_052664.1	1688625
NC_052664.1	1692321
NC_052664.1	1748552
```

**Example (with strand information):**

```text
NC_052664.1	1688625	+
NC_052664.1	1692321	+
NC_052664.1	1748552	+
```

| Column | Description |
|---|---|
| 1 | Chromosome identifier. |
| 2 | Genomic position (1-based). |
| 3 (optional) | Strand information (`+` or `-`). |

When strand information is not provided, sites are extracted using the forward reference orientation and labeled `unknown_forward_reference` in the output.

Example directory:

```text
base_extract_data/
├── samples.tsv
├── invicta_known_sites.tsv
├── H3_worker_head_R1.fq.gz
├── H3_worker_head_R2.fq.gz
├── N1_worker_head_R1.fq.gz
└── N1_worker_head_R2.fq.gz
```

## 6.2. Command

### 6.2.1. Execution

```bash
nohup setsid --fork --wait KOBE base_extract \
    --sample-table "/path/to/base_extract_data/samples.tsv" \
    --out-dir "/path/to/base_extract_result" \
    --gtf-dir "/path/to/ref" \
    --site-dir "/path/to/base_extract_data" \
    --threads 30 \
    --sif "/path/to/rna_editing.sif" \
    --bind "/path/to/project" \
    --force_reditool \
    > "/path/to/base_extract.log" 2>&1 < /dev/null &
```

## 6.3. Parameters

| Parameter | Description |
|---|---|
| `--sample-table` | Sample information table containing sample IDs, species, FASTQ paths, and reference genome paths. |
| `--out-dir` | Output directory for base extraction results. |
| `--gtf-dir` | Directory containing species-specific GTF annotation files. |
| `--site-dir` | Directory containing species-specific known editing site files. |
| `--threads` | Total Snakemake CPU budget (default: 30). |
| `--resume` | Reuse completed outputs when their signatures and completeness checks are valid; rerun missing, incomplete, or changed steps. |
| `--force_reditool` | Reuse compatible in-workdir Picard BAMs when available, but force known-site preparation, all REDItoolKnown jobs, and the merge step to rerun. If a BAM cannot be reused, the required upstream mapping steps are rerun from available FASTQ files. |
| `--sif` | Path to the KOBE Apptainer container image; may be omitted if `RNA_EDITING_SIF` is set. |
| `--bind` | Additional Apptainer bind mount; may be repeated. |

The module uses the `species` column in `samples.tsv` to match species-specific GTF annotation and known editing site files (e.g., `invicta` is matched to `invicta.gtf` and `invicta_known_sites.tsv`). Without `--resume` or `--force_reditool`, the workflow runs from the beginning rather than reusing previous results.



## 6.4. Important Outputs

The main output directory is organized by species and processing stage.

Example:

```text
base_extract_result/
├── invicta/
│   ├── 00_INPUT/
│   │   ├── known_sites.redi.tab
│   │   └── site_metadata.tsv
│   │
│   ├── 00_QC/
│   │   ├── H3_worker_3/
│   │   │   ├── qc.R1.fastq.gz
│   │   │   ├── qc.R2.fastq.gz
│   │   │   ├── qc.html
│   │   │   └── qc.json
│   │   └── N1_worker_4/
│   │
│   ├── 00_REF/
│   │   ├── GENOME/
│   │   └── STAR/
│   │
│   ├── 01_MAP/
│   ├── 02_Q10/
│   ├── 03_PICARD/
│   ├── 04_REDItools/
│   ├── 05_MERGE/
│   └── 06_ALL_SAMPLES/
│       └── invicta_all_samples_base_counts.tsv
│
├── logs/
└── tmp/
```

The main processing stages are:

| Directory | Description |
|---|---|
| `00_INPUT/` | Prepared known editing site files and associated metadata. |
| `00_QC/` | FASTQ quality control results, including retained cleaned FASTQ files. |
| `00_REF/` | Reference genome and STAR index files. |
| `01_MAP/` | STAR alignment results. |
| `02_Q10/` | Mapping quality-filtered BAM files. |
| `03_PICARD/` | Duplicate-removed BAM files and Picard metrics. |
| `04_REDItools/` | REDItoolKnown results for individual samples. |
| `05_MERGE/` | Per-sample nucleotide count tables. |
| `06_ALL_SAMPLES/` | Merged nucleotide count table across all samples. |
| `logs/` | Execution logs for individual processing steps. |
| `tmp/` | Temporary files generated during analysis. |

### 6.4.1. Final Output: All-sample Base Counts

The main result is:

```text
invicta/06_ALL_SAMPLES/invicta_all_samples_base_counts.tsv
```

This file contains genomic nucleotide counts for all specified sites across all analyzed samples.

The output table contains the following columns:

| Column | Description |
|---|---|
| `chromosome` | Chromosome identifier. |
| `position` | Genomic position (1-based). |
| `strand` | Strand information (`+`, `-`, or `unknown_forward_reference`). |
| `<sample_id>_A_count` | Number of A bases in the sample. |
| `<sample_id>_C_count` | Number of C bases in the sample. |
| `<sample_id>_G_count` | Number of G bases in the sample. |
| `<sample_id>_T_count` | Number of T bases in the sample. |

Example:

```text
chromosome	position	strand	H3_worker_3_A_count	H3_worker_3_C_count	H3_worker_3_G_count	H3_worker_3_T_count	N1_worker_4_A_count	N1_worker_4_C_count	N1_worker_4_G_count	N1_worker_4_T_count
NC_052664.1	503815	unknown_forward_reference	0	0	6	3	0	0	185	1
NC_052664.1	1683835	unknown_forward_reference	1	0	3	0	2	0	5	0
```

Each row represents one genomic position, with four nucleotide counts reported for each sample.

The `strand` column indicates the orientation used for nucleotide extraction. Sites labeled `unknown_forward_reference` are extracted relative to the forward reference genome because strand information was not provided.

After further processing, the extracted base counts can be used for downstream analyses, including differential RNA editing site analysis and other analyses.

---

# 7. Utility Scripts

This section describes additional KOBE utility scripts for preparing input files, connecting analysis modules, and processing downstream results.

## 7.1. Prepare Orthology Analysis Inputs

The `ortho_edit` module requires two species-specific input files:

1. **Candidate CDS editing site table** (`*_candi_cds_sites.txt`): Contains editing site information in CDS coordinates.
2. **Longest-CDS FASTA file** (`*_longestcds.fasta`): Contains the longest CDS sequence for each gene.

KOBE provides different utility scripts to prepare these files from genomic coordinates or de novo RNA editing results.

### 7.1.1. Prepare Orthology Inputs from Genomic Coordinates

The `genome_site_cds_position` command uses genomic editing site coordinates, GTF annotation, and the reference genome to generate both required input files for orthology analysis.

**Input Data**

- **Site file**: Headerless, tab-delimited genomic coordinates, containing chromosome, genomic position (1-based), and optional strand information.
- **GTF file**: Full genome annotation.
- **Genome FASTA**: Corresponding uncompressed reference genome sequence.

Example site file:

```text
NC_052664.1	1688625
NC_052664.1	1692321
NC_052664.1	1748552
```

The third column can optionally provide strand information (`+` or `-`).

When strand information is absent, the script infers the strand using the GTF annotation. Genes are identified based on GTF CDS features, and the longest CDS transcript is selected for each gene without additional gene biotype filtering. However, if genes on both strands overlap the same genomic region, the ambiguous site is excluded from the candidate site output.

**Command**

```bash
KOBE genome_site_cds_position \
    --species Hym_Solenopsis_invicta \
    --sites "/path/to/invicta_known_sites.tsv" \
    --gtf "/path/to/invicta.gtf" \
    --genome "/path/to/invicta.fa" \
    --outdir "/path/to/genome_site_cds_position" \
    --sif "/path/to/rna_editing.sif" \
    --bind "/path/to/project"
```

**Parameters**

| Parameter | Description |
|---|---|
| `--species` | Species identifier used for output filenames. |
| `--sites` | Input genomic editing site coordinate file. |
| `--gtf` | Full genome GTF annotation file. |
| `--genome` | Corresponding uncompressed reference genome FASTA file. |
| `--outdir` | Output directory. |
| `--sif` | Path to the Apptainer container image; may be omitted if `RNA_EDITING_SIF` is set. |
| `--bind` | Additional Apptainer bind mount; may be repeated. |

**Outputs**

This command directly generates both input files required by `ortho_edit`, together with an execution log.

```text
genome_site_cds_position/
├── Hym_Solenopsis_invicta_candi_cds_sites.txt
├── Hym_Solenopsis_invicta_longestcds.fasta
└── Hym_Solenopsis_invicta_log.txt
```

- `*_candi_cds_sites.txt`: Candidate CDS editing site table.
- `*_longestcds.fasta`: Longest-CDS sequences for all genes.
- `*_log.txt`: Execution log.

The candidate site table has the following format:

```text
transcript	cds_pos	anno	AA_pos	AA_change	gene_id
AECH19648-RA	35	missense_variant	12	Lys>Arg	AECH19648
AECH19648-RA	27	missense_variant	9	Ile>Met	AECH19648
AECH19648-RA	25	missense_variant	9	Ile>Val	AECH19648
```

CDS sequences are retained regardless of whether their lengths are divisible by three.

### 7.1.2. Prepare Candidate CDS Editing Sites from De Novo Results

Alternatively, the `extract_cds_candidate_sites` command extracts CDS editing site information directly from annotated de novo RNA editing results.

**Input Data**

A functional annotation result generated by the `denovo_detect_anno` module.

Example:

```text
genic.final_sites.DNA10_RNA3_NRNA2_DNAFDR0.05_RNAFDR0.01.tsv
```

**Command**

```bash
KOBE extract_cds_candidate_sites \
    --input "/path/to/denovo_detect_result/edit_site_result/invicta/genic.final_sites.DNA10_RNA3_NRNA2_DNAFDR0.05_RNAFDR0.01.tsv" \
    --output "/path/to/Hym_Solenopsis_invicta_candi_cds_sites.txt" \
    --sif "/path/to/rna_editing.sif" \
    --bind "/path/to/project"
```

**Parameters**

| Parameter | Description |
|---|---|
| `--input` | Annotated de novo RNA editing site result. |
| `--output` | Output candidate CDS editing site table. |
| `--sif` | Path to the Apptainer container image; may be omitted if `RNA_EDITING_SIF` is set. |
| `--bind` | Additional Apptainer bind mount; may be repeated. |

**Output**

```text
Hym_Solenopsis_invicta_candi_cds_sites.txt
```

The output uses the same six-column format described in Section 1.1.

### 7.1.3. Extract Longest-CDS Sequences

When candidate editing sites are prepared from de novo results using `extract_cds_candidate_sites`, the longest-CDS FASTA file must be generated separately using `extract_longest_cds`.

This step is not required when using `genome_site_cds_position`, which generates both required files directly.

**Input Data**

- **GTF file**: Full genome annotation.
- **Genome FASTA**: Corresponding reference genome sequence.

**Command**

```bash
KOBE extract_longest_cds \
    --gtf "/path/to/invicta.gtf" \
    --genome "/path/to/invicta.fa" \
    --output "/path/to/Hym_Solenopsis_invicta_longestcds.fasta" \
    --sif "/path/to/rna_editing.sif" \
    --bind "/path/to/project"
```

**Parameters**

| Parameter | Description |
|---|---|
| `--gtf` | Full genome GTF annotation file. |
| `--genome` | Corresponding uncompressed reference genome FASTA file. |
| `--output` | Output longest-CDS FASTA file. |
| `--sif` | Path to the Apptainer container image; may be omitted if `RNA_EDITING_SIF` is set. |
| `--bind` | Additional Apptainer bind mount; may be repeated. |

**Output**

```text
Hym_Solenopsis_invicta_longestcds.fasta
```

The output contains the longest CDS sequence for each gene. CDS sequences are retained regardless of whether their lengths are divisible by three.

### 7.1.4. Inputs for Orthology Analysis

After preparing the two files, organize them into the corresponding directories for the `ortho_edit` module.

```text
ortho_edit_data/
├── data/
│   └── Hym_Solenopsis_invicta_longestcds.fasta
│
└── edit_txt/
    └── Hym_Solenopsis_invicta_candi_cds_sites.txt
```

The two files are supplied to `ortho_edit` through `--cds_dir` and `--edit_txt_dir`, respectively.

---

## 7.2. Add Genomic Coordinates to Orthology Results

The `add_cds_sites_genomic_coordinates` command adds genomic coordinates to the annotated orthology results generated by `ortho_edit`.

GTF annotation files can be provided for individual species. Only species with a supplied GTF file receive an additional genomic coordinate column in the output table. Coordinates are represented in the `chr:pos:strand` format.

### 7.2.1. Input Data

- **Annotated CSV**: Merged and annotated orthology results generated by `ortho_edit`.
- **GTF files**: Species-specific GTF annotation files. Each file is associated with its species identifier using the format `species=GTF_file_path`.

### 7.2.2. Command

```bash
KOBE add_cds_sites_genomic_coordinates \
    --input "/path/to/merged_edit_sites_codon_table_annotated.csv" \
    --gtf "melanogaster=/path/to/gtf/dmel-all-r6.06.gtf" \
    --gtf "mellifera=/path/to/gtf/mellifera.gtf" \
    --gtf "invicta=/path/to/gtf/invicta.gtf" \
    --output "/path/to/merged_edit_sites_with_genomic_coordinates.tsv" \
    --sif "/path/to/rna_editing.sif" \
    --bind "/path/to/project"
```

### 7.2.3. Parameters

| Parameter | Description |
|---|---|
| `--input` | Annotated orthology results generated by `ortho_edit`. |
| `--gtf` | Species-specific GTF file, specified as `species=GTF_file_path`. Can be supplied multiple times. |
| `--output` | Output annotated table with genomic coordinates. |
| `--sif` | Path to the Apptainer container image; may be omitted if `RNA_EDITING_SIF` is set. |
| `--bind` | Additional Apptainer bind mount; may be repeated. |

**Output**

```text
merged_edit_sites_with_genomic_coordinates.tsv
```

The output retains the original annotated results and adds genomic coordinate information for species with provided GTF files.

---

## 7.3. Prepare Differential Editing Input

The `supported_sites_to_editing_counts` command converts supported editing site results generated by `find_back` into the editing count table required by `pairwise_tests`.

Both known and novel supported editing site tables can be used. The input can contain sites with editing detected in at least one sample.

### 7.3.1. Command

```bash
KOBE supported_sites_to_editing_counts \
    --input "/path/to/find_back_result/statistics/invicta/invicta_novel_supported_sites_ge1.tsv" \
    --species invicta \
    --work_dir "/path/to/diff_edit_site_data" \
    --sif "/path/to/rna_editing.sif" \
    --bind "/path/to/project"
```

**Note:** `--species` is used only to determine the output filename; it does not filter the input sites by species.

For more details about the input format and differential analysis, see [Differential RNA Editing Analysis](#5-differential-rna-editing-analysis).

---

## 7.4. Prepare SnpEff Input Files

The `prepare_snpeff_gtf` command prepares reference files for SnpEff database construction. It retains the transcript with the longest CDS for each coding gene and the transcript with the longest spliced exon length for each noncoding gene (using transcript span if exon features are unavailable), then creates a symbolic link to the corresponding reference genome FASTA file in the output directory.

### 7.4.1. Command

```bash
KOBE prepare_snpeff_gtf \
    --reference_dir "/path/to/denovo_detect_data/ref" \
    --output_dir "/path/to/snpeff_data" \
    --species invicta \
    --force \
    --sif "/path/to/rna_editing.sif" \
    --bind "/path/to/project"
```

**Note:** The input reference directory should contain species-specific genome FASTA and full GTF files. The `--species` parameter specifies the target species, and `--force` regenerates existing preparation results. The Apptainer image can be supplied with `--sif` or through `RNA_EDITING_SIF`.

For more details about SnpEff preparation and database construction, see [De Novo Identification](#3-de-novo-rna-editing-identification).
