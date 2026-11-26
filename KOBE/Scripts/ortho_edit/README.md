# ortho_edit

## Workflow script names

| Previous name | New name |
|---|---|
| `01_final.py` | `build_besthit_matrix.py` |
| `02_final.sh` | `extract_besthit_protein_sequences.sh` |
| `02_2.sh` | `reorder_protein_sequences.sh` |
| `03_final.sh` | `mafft_align_protein_sequences.sh` |
| `04_final.py` | `back_translate_protein_alignment.py` |
| `05_final.py` | `merge_orthologous_editing_sites.py` |
| `info.R` | `group_shared_orthologous_sites.R` |

Other helper scripts retain their descriptive names.

## Run

The launcher uses the project-relative SIF path by default in the current
repository layout. It also accepts `--sif` or the `RNA_EDITING_SIF` environment
variable.

```bash
/path/to/Scripts/ortho_edit/ortho_edit \
  --cds_dir /path/to/project/ortho_edit_data \
  --edit_txt_dir /path/to/project/ortho_edit_txt \
  --work_dir /path/to/project/ortho_result \
  --threads 20
```

`--threads` is the maximum CPU resource budget for the parallel steps.
It defaults to 8 and must be a positive integer.

## Thread allocation

For `--threads INT`:

| Step | Allocation |
|---|---|
| 1. DIAMOND | One DIAMOND job using up to N threads |
| 2. Protein extraction | Up to N worker processes, one thread per process |
| 3. Protein sequence reordering | Up to N worker processes, one thread per process |
| 4. MAFFT | Up to floor(N/2) processes, two threads per process; if N=1, one process with one thread |
| 5. CDS back-translation | Up to N worker processes, one thread per process |
| 6-11 | Serial |

The number of active workers is capped by the number of input files or tasks.
For example, `--threads 20` gives up to 20 single-threaded workers for Steps 2,
3 and 5, and up to 10 MAFFT jobs with two threads each.

Steps execute sequentially. Parallel workers from one step finish before the
next step begins. The thread budget is therefore reused by each step rather
than added across steps.

## Full rebuild behavior

Every run removes only the following workflow-generated paths under `WORK_DIR`
before starting:

- `pep/`
- `besthit_results/`
- `besthit_matrix.tsv`
- `non_multiple_of_3.txt`
- `pep_every/`
- `re_pep/`
- `aln_pep/`
- `aln_cds/`
- `spe_cds_pos/`
- `ortho_edit_result/`

The workflow then rebuilds all these outputs from the supplied input files.
`ortho_edit.log` is not removed by the cleanup; the launcher overwrites it for
the current run. Do not store unrelated data in the listed workflow-owned paths.

The input directories are read-only from the workflow's perspective:
`*_longestcds.fasta` files are read from `--cds_dir`, and
`*_candi_cds_sites.txt` files are read from `--edit_txt_dir`.

## Output layout

```text
WORK_DIR/
├── ortho_edit.log
├── besthit_matrix.tsv
├── non_multiple_of_3.txt
├── pep/
├── besthit_results/
├── pep_every/
├── re_pep/
├── aln_pep/
│   └── logs/
├── aln_cds/
├── spe_cds_pos/
└── ortho_edit_result/
    ├── merged_edit_sites_codon_table.csv
    ├── merged_edit_sites_codon_table_uniq.csv
    ├── merged_edit_sites_codon_table_annotated.csv
    ├── missed_rows_debug.txt
    └── candi_per_species/
```

`missed_rows_debug.txt` is generated only when unmatched rows are detected.

## Logs

- Main log: `WORK_DIR/ortho_edit.log` (overwritten on each run).
- Per-alignment MAFFT logs: `WORK_DIR/aln_pep/logs/`.
- Protein extraction debug log: `WORK_DIR/pep_every/pep_extraction_debug.log`.

## Important note about previous runs

New runs no longer write result tables or `candi_per_species/` into
`--edit_txt_dir`. Any old result tables or `candi_per_species/` directory
already present in that input directory are left untouched for safety; move
or remove them manually after confirming they are previous workflow outputs.

## Install the command (optional)

```bash
chmod +x ortho_edit
ln -s /absolute/path/to/ortho_edit /usr/local/bin/ortho_edit
```
