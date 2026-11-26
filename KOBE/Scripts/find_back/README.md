# find_back (RNA-seq-guided homologous editing-site recovery)

## Invocation

Set the SIF once in the shell:

```bash
export RNA_EDITING_SIF=/path/to/rna_editing.sif
```

Run the module with GTF files supplied individually for the species that need
genomic-coordinate mapping:

```bash
find_back \
  --sample_table /path/to/samples.tsv \
  --out_dir /path/to/find_back_result \
  --cds_pos_dir /path/to/spe_cds_pos \
  --gtf invicta=/path/to/gtf/invicta.gtf \
  --gtf mellifera=/path/to/gtf/mellifera.gtf \
  --edit_txt_dir /path/to/ortho_edit_txt \
  --candi_dir /path/to/ortho_edit_txt/candi_per_species \
  --min_samples 1 \
  --threads 20 \
  --mapping_jobs 4 \
  --min_length 150 \
  --sif /path/to/rna_editing.sif \
  --bind /path/to/project
```

`--gtf` is repeatable and optional. The species key must exactly match the
`species` column in the sample table. A species without a supplied GTF is still
processed through read mapping, REDItoolKnown, count summarization, and
Known/Novel classification; only CDS-to-genome coordinate mapping is skipped
for that species.

The launcher creates `--out_dir` on the host, resolves all explicit paths
independently of the caller's current directory, and binds sample-table inputs,
the supplied GTF files, project inputs, and output directory into Apptainer.
`--edit_txt_dir` supplies the original known-site references for Known/Novel
classification, while `--candi_dir` supplies species-specific `*_candi.txt`
annotations for REDItoolKnown. Neither is inferred from the KOBE installation
directory.

The in-container sample scheduler uses Bash 5.1 or newer (`wait -n -p`) so it can detect any completed child process and immediately refill the Mapping slot.\n\nThe launcher holds a host-side `flock` on `<out_dir>/.find_back.lock` for the
full Apptainer run. It removes the lock-file pathname on normal exit and
catchable termination signals. `SIGKILL` cannot be trapped; the kernel
releases the advisory lock when the process exits, and a subsequent run can
safely acquire the same lock path. A second simultaneous run against the same
canonical output directory exits with an error.

`--mapping_jobs` controls the maximum number of complete sample-level Mapping
workflows running concurrently, including fastp, BWA, sorting, indexing, and
Picard. Its default is 2. When any running sample workflow finishes, the
scheduler immediately starts the next pending sample rather than waiting for
the oldest-started sample. Within each sample, the steps remain sequential.

Mapping defaults are fastp 2 threads, BWA 5 threads, Sambamba sort 5 threads,
and sort memory 12G. They can be overridden independently with
`--fastp_threads`, `--bwa_threads`, `--sort_threads`, and `--sort_mem`.

`--threads` is the REDItoolKnown-stage CPU budget. Each REDItoolKnown process
uses 2 threads; the maximum concurrent sample jobs are `threads / 2`, limited
by the number of pending RNA samples. It must be an even integer >= 2; default
20. `--min_length` is passed to fastp `-l`; default 150.

## Mapping and read processing

For each RNA-seq sample, the workflow performs:

1. fastp paired-end quality control.
2. BWA-MEM alignment against the species-specific longest-CDS reference.
3. SAM-to-BAM conversion and coordinate sorting with Sambamba.
4. Picard `MarkDuplicates` with `REMOVE_DUPLICATES=true`.
5. BAM indexing and validation.
6. REDItoolKnown using the Picard-deduplicated BAM to extract candidate-site A/G evidence.

The clean paired FASTQ files are retained under:

`<out_dir>/QC/<species>/<sample_id>/`

The final BAM used by REDItoolKnown is:

`<out_dir>/Mapping_BWA_to_CDS/<species>/<sample_id>/<sample_id>.sort.nodup.bam`

The coordinate-sorted pre-Picard BAM is retained as an intermediate. The
Picard metrics file and final BAM index are retained alongside the final BAM.

## Resume and BAM reuse controls

- No resume flag: rebuild Mapping and rerun all REDItoolKnown sample tasks.
- `--resume`: reuse completed Mapping/Picard outputs when the stored input and
  parameter signature matches; regenerate candidate-site inputs and reuse
  REDItoolKnown outputs only when their candidate-file hash and output
  completeness permit it.
- `--force_reditool`: inspect each sample's in-workdir Picard BAM and index.
  Reuse the BAM only if the BAM is non-empty, `samtools quickcheck` succeeds,
  `samtools idxstats` can read the index, BAM read-group `SM` matches the
  sample ID, and BAM `@SQ` sequence names and lengths match the current CDS
  reference dictionary exactly by sequence names and lengths; no MD5 or
  per-base sequence-content comparison is performed. If validation fails,
  the sample is remapped from FASTQ (or retained clean FASTQ when available).
  If no valid BAM and no
  complete raw/clean FASTQ pair exist, the run stops with an error.
  Regardless of whether BAMs are reused, candidate-site preparation is
  refreshed, every REDItoolKnown sample task is forced, and downstream
  summaries/classification/coordinate mapping are rebuilt.

The explicit `--reuse_bam` option has been removed. `--force_reditool` is the
only option that enables automatic in-workdir BAM reuse. Without it, the
default remains a fresh run; KOBE does not silently reuse an old BAM.

## GTF-based genomic-coordinate mapping

GTF files are provided as repeated `--gtf species=path` options. For example:

```bash
--gtf invicta=/project/gtf/invicta.gtf \
--gtf mellifera=/project/gtf/mellifera.gtf
```

The GTF species identifier must match the sample table exactly. Missing GTFs
are not errors: the corresponding species retains all CDS-level results, but
no `chrom`, `genome_pos`, or `genome_strand` columns are appended for that
species. Supplied GTFs are validated as non-empty files and used only in the
final CDS-to-genome coordinate mapping step.

## Container dependencies

Inside the SIF the workflow expects `fastp`, `bwa`, `samtools`, `sambamba`,
`java`, Picard, `python3`, `perl`, a Python 2 executable named `python2` with
`pysam`, and `REDItoolKnown.py` somewhere under `$RNA_EDITING_HOME` (default
`/opt/rna_editing`). The Perl helper is expected at:

`/opt/rna_editing/helpers/add.intersect.vcf.info.pl`

The main script discovers REDItoolKnown, the Perl helper, and Picard under the
container software root. `JAVA`, `PICARD_JAR`, `REDI_PYTHON`, `REDITool`,
`FASTP`, `BWA`, `SAMTOOLS`, and `SAMBAMBA` can override discovery inside the
container if needed.

## Output/logging

- Main log: `<out_dir>/find_back.log`
- Mapping/fastp/Picard sample logs: `<out_dir>/logs/`
- REDItoolKnown per-sample logs: `<out_dir>/logs/reditool_known/`
- Clean FASTQ: `<out_dir>/QC/<species>/<sample_id>/`
- Mapping and final Picard BAMs: `<out_dir>/Mapping_BWA_to_CDS/<species>/<sample_id>/`
- Counts and summaries: `<out_dir>/candidate_counts/`
- Known/Novel statistics: `<out_dir>/statistics/`
- Optional genomic-coordinate tables: `<out_dir>/genomic_coordinates/`

The original candidate preparation, physical-CDS deduplication, handle
merging, G/(A+G) calculation, multi-sample summary, and Known/Novel
classification algorithms are retained.
