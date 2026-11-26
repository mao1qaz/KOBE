#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SAMPLE_TABLE=""
OUT_DIR=""
GTF_DIR=""
SITE_DIR=""
THREADS=30
REFERENCE_PREPARE_THREADS=2
FASTP_THREADS=10
STAR_INDEX_THREADS=15
STAR_THREADS=15
SAMTOOLS_THREADS=10
PICARD_INDEX_THREADS=2
REDITools_THREADS=2
REFERENCE_PREPARE_MEMORY_MB=6000
FASTP_MEMORY_MB=4000
STAR_INDEX_MEMORY_MB=64000
STAR_MAP_MEMORY_MB=16000
Q10_MEMORY_MB=8000
PICARD_MEMORY_MB=81200
REDITools_MEMORY_MB=12000
MAX_MEMORY_MB=90000
FASTP_N_BASE_LIMIT=5
FASTP_QUALIFIED_QUALITY=15
FASTP_UNQUALIFIED_PERCENT=50
FASTP_LENGTH_REQUIRED=150
FASTP_DISABLE_TRIM_POLY_G=true
STAR_SJDB_OVERHANG=149
STAR_LIMIT_GENOME_GENERATE_RAM=60000000000
STAR_LIMIT_BAM_SORT_RAM=12884901888
STAR_GENOME_SA_INDEX_NBASES=auto
STAR_MULTIMAP_MAX=1
STAR_MISMATCH_MAX=10
STAR_TWOPASS_MODE=Basic
STAR_OUT_READS_UNMAPPED=Fastx
Q10_MAPQ=10
REDITools_MIN_COVERAGE=1
REDITools_MIN_BASE_QUALITY=25
REDITools_MIN_MAPQ=10
REDITools_MIN_VARIANT_READS=0
REDITools_MIN_EDITING_FREQUENCY=0.0
PICARD_REMOVE_DUPLICATES=true
PICARD_CREATE_INDEX=true
PICARD_ASSUME_SORTED=true
PICARD_VALIDATION_STRINGENCY=SILENT
RESUME=0
FORCE_REDITOOL=0
REUSE_BAM=0
DRY_RUN=0

usage() {
    cat <<'USAGE'
Usage:
  base_extract --sample-table FILE --out-dir DIR --gtf-dir DIR --site-dir DIR [--sif FILE] [--bind DIR ...] [options]

Required input:
  --sample-table FILE          TSV: sample_id, species, fq1, fq2, genome_fasta
  --out-dir DIR                Output root; created on the host by launcher
  --gtf-dir DIR                Contains <species>.gtf
  --site-dir DIR               Contains <species>_known_sites.tsv
  --sif FILE                   RNA editing Apptainer SIF (or RNA_EDITING_SIF)

Main resources:
  --threads INT                Total Snakemake CPU budget (default: 30)
  --max-memory-mb INT          Total Snakemake mem_mb budget (default: 90000)

Per-rule threads (defaults):
  --reference-prepare-threads INT (default: 2)
  --fastp-threads INT           (default: 10)
  --star-index-threads INT      (default: 15)
  --star-threads INT            (default: 15)
  --samtools-threads INT        (default: 10)
  --picard-index-threads INT    (default: 2; samtools fallback indexing)
  --reditools-threads INT       (default: 2 threads per REDItoolKnown sample job)

Per-rule memory requests in MB (defaults):
  --reference-prepare-memory-mb INT (default: 6000)
  --fastp-memory-mb INT             (default: 4000)
  --star-index-memory-mb INT        (default: 64000)
  --star-map-memory-mb INT          (default: 16000)
  --q10-memory-mb INT               (default: 8000)
  --picard-memory-mb INT            (default: 81200)
  --reditools-memory-mb INT         (default: 12000)

fastp options (defaults preserve the original workflow):
  --fastp-n-base-limit INT        (default: 5)
  --fastp-qualified-quality INT   (default: 15)
  --fastp-unqualified-percent INT (default: 50)
  --fastp-length-required INT|none (default: 150)
  --fastp-enable-trim-poly-g      Enable poly-G trimming (default is disabled via -G)
  --fastp-disable-trim-poly-g     Disable poly-G trimming (default)

STAR options:
  --star-sjdb-overhang INT                (default: 149)
  --star-limit-genome-generate-ram INT    (default: 60000000000)
  --star-limit-bam-sort-ram INT           (default: 12884901888)
  --star-genome-sa-index-nbases INT|auto  (default: auto)
  --star-multimap-max INT                 (default: 1)
  --star-mismatch-max INT                 (default: 10)
  --star-twopass-mode STRING              (default: Basic)
  --star-out-reads-unmapped STRING        (allowed: Fastx|none; default: Fastx)

Other analysis options:
  --q10-mapq INT                          (default: 10)
  --reditools-min-coverage INT            (default: 1)
  --reditools-min-base-quality INT        (default: 25)
  --reditools-min-mapq INT                (default: 10)
  --reditools-min-variant-reads INT       (default: 0)
  --reditools-min-editing-frequency FLOAT (default: 0.0)
  --picard-keep-duplicates                Keep duplicates (default removes them)
  --picard-remove-duplicates              Remove duplicates (default)
  --picard-no-create-index                Use samtools to create index instead
  --picard-create-index                   Let Picard create index (default)
  --picard-do-not-assume-sorted           Disable ASSUME_SORTED (default enabled)
  --picard-assume-sorted                  Enable ASSUME_SORTED (default)
  --picard-validation-stringency STRING   (default: SILENT)

Reuse modes:
  --resume                               Reuse valid completed outputs; rerun changed/incomplete steps
  --force_reditool (alias: --force-reditool) Validate/reuse in-workdir Picard BAMs; force site preparation and all REDItoolKnown jobs

Other:
  --bind DIR              Additional Apptainer bind; may be repeated
  --dry-run                              Show Snakemake plan without executing jobs
  --help                                 Show this help

Software paths are discovered inside the SIF; no laboratory-server paths are embedded.
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample-table) SAMPLE_TABLE="$2"; shift 2 ;;
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        --gtf-dir) GTF_DIR="$2"; shift 2 ;;
        --site-dir) SITE_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --reference-prepare-threads) REFERENCE_PREPARE_THREADS="$2"; shift 2 ;;
        --fastp-threads) FASTP_THREADS="$2"; shift 2 ;;
        --star-index-threads) STAR_INDEX_THREADS="$2"; shift 2 ;;
        --star-threads) STAR_THREADS="$2"; shift 2 ;;
        --samtools-threads) SAMTOOLS_THREADS="$2"; shift 2 ;;
        --picard-index-threads) PICARD_INDEX_THREADS="$2"; shift 2 ;;
        --reditools-threads) REDITools_THREADS="$2"; shift 2 ;;
        --reference-prepare-memory-mb) REFERENCE_PREPARE_MEMORY_MB="$2"; shift 2 ;;
        --fastp-memory-mb) FASTP_MEMORY_MB="$2"; shift 2 ;;
        --star-index-memory-mb) STAR_INDEX_MEMORY_MB="$2"; shift 2 ;;
        --star-map-memory-mb) STAR_MAP_MEMORY_MB="$2"; shift 2 ;;
        --q10-memory-mb) Q10_MEMORY_MB="$2"; shift 2 ;;
        --picard-memory-mb) PICARD_MEMORY_MB="$2"; shift 2 ;;
        --reditools-memory-mb) REDITools_MEMORY_MB="$2"; shift 2 ;;
        --max-memory-mb) MAX_MEMORY_MB="$2"; shift 2 ;;
        --fastp-n-base-limit) FASTP_N_BASE_LIMIT="$2"; shift 2 ;;
        --fastp-qualified-quality) FASTP_QUALIFIED_QUALITY="$2"; shift 2 ;;
        --fastp-unqualified-percent) FASTP_UNQUALIFIED_PERCENT="$2"; shift 2 ;;
        --fastp-length-required) FASTP_LENGTH_REQUIRED="$2"; shift 2 ;;
        --fastp-enable-trim-poly-g) FASTP_DISABLE_TRIM_POLY_G=false; shift ;;
        --fastp-disable-trim-poly-g) FASTP_DISABLE_TRIM_POLY_G=true; shift ;;
        --star-sjdb-overhang) STAR_SJDB_OVERHANG="$2"; shift 2 ;;
        --star-limit-genome-generate-ram) STAR_LIMIT_GENOME_GENERATE_RAM="$2"; shift 2 ;;
        --star-limit-bam-sort-ram) STAR_LIMIT_BAM_SORT_RAM="$2"; shift 2 ;;
        --star-genome-sa-index-nbases) STAR_GENOME_SA_INDEX_NBASES="$2"; shift 2 ;;
        --star-multimap-max) STAR_MULTIMAP_MAX="$2"; shift 2 ;;
        --star-mismatch-max) STAR_MISMATCH_MAX="$2"; shift 2 ;;
        --star-twopass-mode) STAR_TWOPASS_MODE="$2"; shift 2 ;;
        --star-out-reads-unmapped) STAR_OUT_READS_UNMAPPED="$2"; shift 2 ;;
        --q10-mapq) Q10_MAPQ="$2"; shift 2 ;;
        --reditools-min-coverage) REDITools_MIN_COVERAGE="$2"; shift 2 ;;
        --reditools-min-base-quality) REDITools_MIN_BASE_QUALITY="$2"; shift 2 ;;
        --reditools-min-mapq) REDITools_MIN_MAPQ="$2"; shift 2 ;;
        --reditools-min-variant-reads) REDITools_MIN_VARIANT_READS="$2"; shift 2 ;;
        --reditools-min-editing-frequency) REDITools_MIN_EDITING_FREQUENCY="$2"; shift 2 ;;
        --picard-keep-duplicates) PICARD_REMOVE_DUPLICATES=false; shift ;;
        --picard-remove-duplicates) PICARD_REMOVE_DUPLICATES=true; shift ;;
        --picard-no-create-index) PICARD_CREATE_INDEX=false; shift ;;
        --picard-create-index) PICARD_CREATE_INDEX=true; shift ;;
        --picard-do-not-assume-sorted) PICARD_ASSUME_SORTED=false; shift ;;
        --picard-assume-sorted) PICARD_ASSUME_SORTED=true; shift ;;
        --picard-validation-stringency) PICARD_VALIDATION_STRINGENCY="$2"; shift 2 ;;
        --resume) RESUME=1; shift ;;
        --force_reditool|--force-reditool) FORCE_REDITOOL=1; RESUME=1; shift ;;
        --reuse_bam|--reuse-bam|--bam-dir)
            echo "ERROR: $1 has been removed. Use --force_reditool to validate/reuse in-workdir BAMs." >&2
            exit 1 ;;
        --dry-run|-n) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: Unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

# --force_reditool is the only mode that accepts/reuses validated in-workdir BAMs.
if (( FORCE_REDITOOL == 1 )); then
    REUSE_BAM=1
fi

for item in "$THREADS" "$REFERENCE_PREPARE_THREADS" "$FASTP_THREADS" "$STAR_INDEX_THREADS" "$STAR_THREADS" "$SAMTOOLS_THREADS" "$PICARD_INDEX_THREADS" "$REDITools_THREADS" \
    "$REFERENCE_PREPARE_MEMORY_MB" "$FASTP_MEMORY_MB" "$STAR_INDEX_MEMORY_MB" "$STAR_MAP_MEMORY_MB" "$Q10_MEMORY_MB" "$PICARD_MEMORY_MB" "$REDITools_MEMORY_MB" "$MAX_MEMORY_MB"; do
    [[ "$item" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: thread and memory values must be positive integers: $item" >&2; exit 1; }
done

if [[ -z "$SAMPLE_TABLE" || -z "$OUT_DIR" || -z "$GTF_DIR" || -z "$SITE_DIR" ]]; then
    echo "ERROR: --sample-table, --out-dir, --gtf-dir and --site-dir are required." >&2
    usage >&2
    exit 1
fi
if (( PICARD_MEMORY_MB < 1024 )); then
    echo "ERROR: --picard-memory-mb must be at least 1024 MB because the Picard heap is derived in whole GB." >&2
    exit 1
fi
for requested_mem in "$REFERENCE_PREPARE_MEMORY_MB" "$FASTP_MEMORY_MB" "$STAR_INDEX_MEMORY_MB" "$STAR_MAP_MEMORY_MB" "$Q10_MEMORY_MB" "$PICARD_MEMORY_MB" "$REDITools_MEMORY_MB"; do
    if (( requested_mem > MAX_MEMORY_MB )); then
        echo "ERROR: A per-rule memory request ($requested_mem MB) exceeds --max-memory-mb ($MAX_MEMORY_MB MB)." >&2
        exit 1
    fi
done

# Discover all tools from the container environment. No host/lab paths are embedded.
find_exe() { command -v "$1" 2>/dev/null || true; }
FASTP="${FASTP:-$(find_exe fastp)}"
STAR="${STAR:-$(find_exe STAR)}"
SAMTOOLS="${SAMTOOLS:-$(find_exe samtools)}"
PYTHON3="${PYTHON3:-$(find_exe python3)}"
SNAKEMAKE_BIN="${SNAKEMAKE:-$(find_exe snakemake)}"
JAVA="${JAVA:-$(find_exe java)}"
REDI_PYTHON="${REDI_PYTHON:-$(find_exe python2)}"
REDITOOL_KNOWN="${REDITOOL_KNOWN:-}"
PICARD_JAR="${PICARD_JAR:-}"

if [[ -z "$REDI_PYTHON" && -x /opt/rna_editing/bin/python2 ]]; then REDI_PYTHON=/opt/rna_editing/bin/python2; fi
if [[ -z "$FASTP" && -x /opt/rna_editing/bin/fastp ]]; then FASTP=/opt/rna_editing/bin/fastp; fi
if [[ -z "$STAR" && -x /opt/rna_editing/bin/STAR ]]; then STAR=/opt/rna_editing/bin/STAR; fi
if [[ -z "$SAMTOOLS" && -x /opt/rna_editing/bin/samtools ]]; then SAMTOOLS=/opt/rna_editing/bin/samtools; fi
if [[ -z "$SNAKEMAKE_BIN" && -x /opt/rna_editing/envs/analysis_py39/bin/snakemake ]]; then SNAKEMAKE_BIN=/opt/rna_editing/envs/analysis_py39/bin/snakemake; fi
if [[ -z "$REDITOOL_KNOWN" && -d /opt/rna_editing ]]; then
    REDITOOL_KNOWN="$(find /opt/rna_editing -type f -name REDItoolKnown.py -print -quit 2>/dev/null || true)"
fi
if [[ -z "$PICARD_JAR" && -d /opt/rna_editing ]]; then
    PICARD_JAR="$(find /opt/rna_editing -type f \( -iname 'picard.jar' -o -iname 'picard-*.jar' \) -print 2>/dev/null | sort -V | tail -n 1)"
fi

for item in FASTP STAR SAMTOOLS PYTHON3 SNAKEMAKE_BIN JAVA REDI_PYTHON REDITOOL_KNOWN PICARD_JAR; do
    value="${!item:-}"
    [[ -n "$value" && -f "$value" ]] || { echo "ERROR: Required container software not found: $item=${value:-[not found]}" >&2; exit 1; }
done
[[ -x "$FASTP" && -x "$STAR" && -x "$SAMTOOLS" && -x "$PYTHON3" && -x "$SNAKEMAKE_BIN" && -x "$JAVA" && -x "$REDI_PYTHON" ]] || {
    echo "ERROR: One or more discovered software paths are not executable." >&2; exit 1;
}
[[ -f "$REDITOOL_KNOWN" && -f "$PICARD_JAR" ]] || { echo "ERROR: REDItoolKnown or Picard JAR is missing." >&2; exit 1; }
"$REDI_PYTHON" -c 'import pysam' >/dev/null 2>&1 || { echo "ERROR: pysam is unavailable in REDItools Python: $REDI_PYTHON" >&2; exit 1; }

if (( REUSE_BAM == 1 )); then
    "$PYTHON3" "$SCRIPT_DIR/validate_reuse_bam.py" \
        --sample-table "$SAMPLE_TABLE" \
        --out-dir "$OUT_DIR" \
        --samtools "$SAMTOOLS"
fi

if (( FORCE_REDITOOL == 1 )); then
    MODE="force-reditool"
elif (( REUSE_BAM == 1 )); then
    MODE="reuse-bam"
elif (( RESUME == 1 )); then
    MODE="resume"
else
    MODE="fresh"
fi

CPU_REDI_JOBS=$((THREADS / REDITools_THREADS))
MEM_REDI_JOBS=$((MAX_MEMORY_MB / REDITools_MEMORY_MB))
if (( CPU_REDI_JOBS < MEM_REDI_JOBS )); then MAX_REDI_JOBS=$CPU_REDI_JOBS; else MAX_REDI_JOBS=$MEM_REDI_JOBS; fi

cat <<EOF
============================================================
Base Extraction workflow
============================================================
Mode: $MODE
Resume enabled: $RESUME
Force REDItoolKnown: $FORCE_REDITOOL
Reuse BAM: $REUSE_BAM
Sample table: $SAMPLE_TABLE
Output directory: $OUT_DIR
GTF directory: $GTF_DIR
Known-site directory: $SITE_DIR
Snakemake CPU budget: $THREADS
Snakemake memory budget: $MAX_MEMORY_MB MB
REDItoolKnown threads per sample: $REDITools_THREADS
REDItoolKnown memory request per sample: $REDITools_MEMORY_MB MB
REDItoolKnown concurrency upper bound by CPU/memory alone: $MAX_REDI_JOBS
Picard memory request: $PICARD_MEMORY_MB MB
fastp: $FASTP
STAR: $STAR
samtools: $SAMTOOLS
REDItoolKnown: $REDITOOL_KNOWN
REDItools Python: $REDI_PYTHON
Picard JAR: $PICARD_JAR
Java: $JAVA
============================================================
EOF

CONFIG_ARGS=(
    "work_dir=$OUT_DIR"
    "samples_tsv=$SAMPLE_TABLE"
    "gtf_dir=$GTF_DIR"
    "site_dir=$SITE_DIR"
    "reuse_bam=$REUSE_BAM"
    "resume_mode=$RESUME"
    "reuse_bam_manifest=$OUT_DIR/.reuse_bam_manifest.json"
    "tool_fastp=$FASTP"
    "tool_star=$STAR"
    "tool_samtools=$SAMTOOLS"
    "tool_redi_python=$REDI_PYTHON"
    "tool_reditool_known=$REDITOOL_KNOWN"
    "tool_picard=$PICARD_JAR"
    "tool_java=$JAVA"
    "reference_prepare_threads=$REFERENCE_PREPARE_THREADS"
    "fastp_threads=$FASTP_THREADS"
    "star_index_threads=$STAR_INDEX_THREADS"
    "star_threads=$STAR_THREADS"
    "samtools_threads=$SAMTOOLS_THREADS"
    "picard_index_threads=$PICARD_INDEX_THREADS"
    "reditools_threads=$REDITools_THREADS"
    "reference_prepare_memory_mb=$REFERENCE_PREPARE_MEMORY_MB"
    "fastp_memory_mb=$FASTP_MEMORY_MB"
    "star_index_memory_mb=$STAR_INDEX_MEMORY_MB"
    "star_map_memory_mb=$STAR_MAP_MEMORY_MB"
    "q10_memory_mb=$Q10_MEMORY_MB"
    "picard_memory_mb=$PICARD_MEMORY_MB"
    "reditools_memory_mb=$REDITools_MEMORY_MB"
    "fastp_n_base_limit=$FASTP_N_BASE_LIMIT"
    "fastp_qualified_quality=$FASTP_QUALIFIED_QUALITY"
    "fastp_unqualified_percent=$FASTP_UNQUALIFIED_PERCENT"
    "fastp_length_required=$FASTP_LENGTH_REQUIRED"
    "fastp_disable_trim_poly_g=$FASTP_DISABLE_TRIM_POLY_G"
    "star_sjdb_overhang=$STAR_SJDB_OVERHANG"
    "star_limit_genome_generate_ram=$STAR_LIMIT_GENOME_GENERATE_RAM"
    "star_limit_bam_sort_ram=$STAR_LIMIT_BAM_SORT_RAM"
    "star_genome_sa_index_nbases=$STAR_GENOME_SA_INDEX_NBASES"
    "star_multimap_max=$STAR_MULTIMAP_MAX"
    "star_mismatch_max=$STAR_MISMATCH_MAX"
    "star_twopass_mode=$STAR_TWOPASS_MODE"
    "star_out_reads_unmapped=$STAR_OUT_READS_UNMAPPED"
    "q10_mapq=$Q10_MAPQ"
    "reditools_min_coverage=$REDITools_MIN_COVERAGE"
    "reditools_min_base_quality=$REDITools_MIN_BASE_QUALITY"
    "reditools_min_mapq=$REDITools_MIN_MAPQ"
    "reditools_min_variant_reads=$REDITools_MIN_VARIANT_READS"
    "reditools_min_editing_frequency=$REDITools_MIN_EDITING_FREQUENCY"
    "picard_remove_duplicates=$PICARD_REMOVE_DUPLICATES"
    "picard_create_index=$PICARD_CREATE_INDEX"
    "picard_assume_sorted=$PICARD_ASSUME_SORTED"
    "picard_validation_stringency=$PICARD_VALIDATION_STRINGENCY"
)

CMD=("$SNAKEMAKE_BIN" -s "$SCRIPT_DIR/01_Snakefile" --directory "$OUT_DIR" --cores "$THREADS" --resources "mem_mb=$MAX_MEMORY_MB" --rerun-incomplete -p --config "${CONFIG_ARGS[@]}")
if (( RESUME == 0 )); then
    CMD+=(--forceall)
fi
if (( FORCE_REDITOOL == 1 )); then
    CMD+=(--forcerun prepare_known_sites reditool_known merge_all_samples)
fi
if (( DRY_RUN == 1 )); then
    CMD+=(-n)
fi

# A previous SIGKILL cannot run shell traps. Because the host-side flock
# guarantees exclusive access, clear only stale Snakemake locks in this workdir.
if (( DRY_RUN == 0 )); then
    "$SNAKEMAKE_BIN" -s "$SCRIPT_DIR/01_Snakefile" --directory "$OUT_DIR" --config "${CONFIG_ARGS[@]}" --unlock >/dev/null 2>&1 || true
fi

exec "${CMD[@]}"
