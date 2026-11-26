# RNA editing add utilities

All commands run through the configured Apptainer image. Set `RNA_EDITING_SIF` once, or pass `--sif` to each command. All input/output paths are explicit command-line arguments. `--bind` / `-B` can be repeated for additional filesystem locations. Launchers automatically create output directories and bind input paths, output parents, and resolved symlink targets. Each utility supports `--help` without requiring an SIF image or Apptainer installation.

```bash
bash /path/to/KOBE/install.sh
export RNA_EDITING_SIF="/path/to/rna_editing.sif"
```

The package installer registers all add-tool commands together with the other modules. See the package-level README for details.

## 1. genome_site_cds_position

Maps genomic candidate sites to CDS coordinates and produces its result files and log under `--outdir`.

```bash
genome_site_cds_position \\
  --species Hym_Solenopsis_invicta \\
  --sites /project/data/invicta_known_sites.tsv \\
  --gtf /project/ref/invicta.gtf \\
  --genome /project/ref/invicta.fa \\
  --outdir /project/results/genome_site_cds_position
```

## 2. add_cds_sites_genomic_coordinates

Maps species-specific CDS positions in the annotated table to genomic coordinates. Repeat `--gtf SPECIES=FILE` for each species.

```bash
add_cds_sites_genomic_coordinates \\
  --input /project/data/merged_edit_sites_codon_table_annotated.csv \\
  --gtf melanogaster=/project/ref/dmel-all-r6.06.gtf \\
  --gtf mellifera=/project/ref/mellifera.gtf \\
  --gtf invicta=/project/ref/invicta.gtf \\
  --output /project/results/merged_edit_sites_with_genomic_coordinates.tsv
```

## 3. extract_cds_candidate_sites

Extracts synonymous and missense CDS candidate sites for downstream orthology input.

```bash
extract_cds_candidate_sites \\
  --input /project/results/invicta/genic.final_sites.tsv \\
  --output /project/results/orthology/invicta_candi_cds_sites.txt
```

## 4. supported_sites_to_editing_counts

Converts a `supported_sites_ge1` table into an A/G count matrix. The output filename remains `<species>_editing_counts.tsv` under `--work_dir`.

```bash
supported_sites_to_editing_counts \\
  --input /project/results/mellifera_supported_sites_ge1.tsv \\
  --species mellifera \\
  --work_dir /project/results/diff_edit_site_data
```

## 5. extract_longest_cds

Extracts one longest CDS transcript per gene. If the FASTA `.fai` index is absent, the script uses `samtools faidx` from the container.

```bash
extract_longest_cds \\
  --gtf /project/ref/invicta.gtf \\
  --genome /project/ref/invicta.fa \\
  --output /project/results/Hym_Solenopsis_invicta_longestcds.fasta
```

## 6. prepare_snpeff_gtf

Prepares a species-specific representative GTF and genome link for SnpEff.

```bash
prepare_snpeff_gtf \\
  --reference_dir /project/ref \\
  --output_dir /project/snpeff_data \\
  --species invicta \\
  --force
```

Expected reference files: `<species>.fa` and `<species>.gtf` under `--reference_dir`.

## Common options

- `--sif /path/to/image.sif`: override the `RNA_EDITING_SIF` environment variable.
- `--bind /path/to/extra/location` or `-B /path/to/extra/location`: add an extra Apptainer bind. Repeat as needed.
- For `--gtf SPECIES=PATH` and `--genome SPECIES=PATH`, the launcher recognizes the path after the first equals sign and mounts it, including resolved symlink targets. Legacy colon syntax is also accepted.

The Python analysis routines are retained from the supplied source files. Changes in this package are limited to command-launcher/container integration, explicit output-path handling where required, and standalone command help/examples.


### Bind performance note

Explicit `--bind` / `-B` directories are mounted as supplied and are not
recursively searched for symlinks. Input files passed through file arguments
have their own symlink targets resolved automatically. Only `--reference_dir`
is scanned for symlink targets, because it is a designated input directory.
Avoid using a broad bind such as `/data` unless it is actually required;
the individual input paths are mounted automatically.
