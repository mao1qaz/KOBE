# Conserved RNA editing plotting module

This package provides three independent commands:

1. `conserve_plot`: run the original homologous conserved-site extraction, filtering, plotting, summary-statistics, and coordinate-export workflow.
2. `add_new_site`: add recovered sites to the original homologous groups using the original `transcript:cds_pos` matching logic.
3. `find_back_plot`: run the supplementary extraction, filtering, and plotting workflow using the updated annotation CSV.

The analysis and site-matching logic is retained. Changes are limited to explicit input/output paths, Apptainer execution, and use of Liberation Serif fonts installed in the container.

## Requirements

- Apptainer
- An RNA Editing SIF containing the Python environment and the `fonts-liberation2` package. The four Liberation Serif font files are expected at:
  - `/usr/share/fonts/truetype/liberation2/LiberationSerif-Regular.ttf`
  - `/usr/share/fonts/truetype/liberation2/LiberationSerif-Bold.ttf`
  - `/usr/share/fonts/truetype/liberation2/LiberationSerif-Italic.ttf`
  - `/usr/share/fonts/truetype/liberation2/LiberationSerif-BoldItalic.ttf`
- The analysis dependencies used by the scripts, including pandas and matplotlib. The summary step no longer requires the `venn` package.

The launchers automatically bind the script directory and the parent directories of input files, plus the CDS and output directories, at their original absolute paths inside the container. No `--bind` argument is required for the listed inputs.

## 1. Original conserved-site plotting

```bash
/path/to/conserve_plot \
  --annotated /path/to/merged_edit_sites_codon_table_annotated.csv \
  --cds_dir /path/to/ortho_edit_data \
  --min_species 2 \
  --species_order melanogaster,mellifera,echinatior,invicta \
  --source_priority melanogaster,mellifera,echinatior,invicta \
  --outdir /path/to/conserve_plot_result \
  --sif /path/to/rna_editing.sif \
  --gtf invicta=/path/to/invicta.gtf \
  --gtf mellifera=/path/to/mellifera.gtf \\
  --col 8
```

`--gtf` is optional and may be repeated. `--col` controls the number of plot columns (default: 8). Only species with a supplied GTF are included in genome-coordinate mapping. All outputs and `conserve_plot.log` are written directly under `--outdir`.

## 2. Add recovered sites

```bash
/path/to/add_new_site \
  --annotated /path/to/merged_edit_sites_codon_table_annotated.csv \
  --sites /path/to/invicta_novel_supported_sites.tsv \
  --sites /path/to/mellifera_novel_supported_sites.tsv \
  --sites /path/to/echinatior_novel_supported_sites.tsv \
  --outdir /path/to/find_back_plot_result \
  --sif /path/to/rna_editing.sif
```

Each `--sites` TSV is read using its `transcript` and `cds_pos` columns. The original matching constructs `transcript:cds_pos` keys and compares them with the species-specific `*_cds_pos` columns in the annotated CSV. No species argument or filename-based species inference is used. Other columns in the site TSV are not used by this matching step.

The updated file is written to:

`<outdir>/merged_edit_sites_codon_table_annotated_with_new_sites.csv`

The log is written to `<outdir>/add_new_site.log`.

## 3. Plot recovered-site groups

Run this after command 2:

```bash
/path/to/find_back_plot \
  --annotated /path/to/find_back_plot_result/merged_edit_sites_codon_table_annotated_with_new_sites.csv \
  --cds_dir /path/to/ortho_edit_data \
  --group_col new_id2_group \
  --min_species 2 \
  --species_order melanogaster,mellifera,echinatior,invicta \
  --source_priority melanogaster,mellifera,echinatior,invicta \
  --outdir /path/to/find_back_plot_result \
  --sif /path/to/rna_editing.sif
```

`--group_col` accepts `ID2_group` or `new_id2_group`. `--ncols` is optional and defaults to 8.

The log is written directly under `--outdir` as `find_back_plot.log`. Results are written under:

`<outdir>/sup_conserve_plot_result/<group_col>_ge<N>species/`

## Notes

- All path inputs are explicit; the scripts do not derive project paths from a fixed server directory.
- The font directory is fixed inside the SIF at `/usr/share/fonts/truetype/liberation2`; no host font directory is mounted.
- Existing result files with the same names are overwritten as in the original scripts. Choose a new `--outdir` to retain previous results.
- `--help` is available for all three commands.


### Additional bind mounts

All three command-line launchers accept repeatable `--bind DIR` arguments.
Each supplied host directory is mounted at the same absolute path inside the
container. The shared container runner also detects file and directory
symlinks among the supplied inputs and binds both the supplied path's parent
and the resolved symlink target's parent/directory.

Example:

```bash
conserve_plot ...
```


## Conserved-site summary outputs

The original workflow keeps the summary statistics without drawing a Venn
diagram. It writes:

- `conserved_sites_ge<N>species_summary.tsv`: total sites and Nonsyn/Syn/Mixed/NA counts and percentages.
- `conserved_sites_ge<N>species_species_counts.tsv`: conserved-site counts for each species.
- `conserved_sites_ge<N>species_venn_regions.tsv`: exact mutually exclusive species-combination counts. The historical filename is retained for compatibility; no Venn graphic is produced.

The species-combination table enumerates all non-empty combinations, including
zero-count combinations. It is not restricted to six species.

After the full workflow succeeds, temporary intermediate tables are removed:
`conserved_sites_ge<N>species.tsv`,
`conserved_sites_ge<N>species_transcripts.tsv`, and
`conserved_site_all_species_unmatched.tsv`. The all-species sequence table
and its per-species QC table are retained.
