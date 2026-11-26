
#!/bin/bash
# ============================================================
# RNA-only Homologous Editing Recovery Mapping Pipeline
#
# Workflow:
#   1. Read sample table
#   2. Build BWA index and FASTA index for unique CDS references
#   3. fastp quality control
#   4. BWA-MEM mapping to species-specific CDS
#   5. BAM conversion
#   6. Coordinate sorting
#   7. BAM indexing
#   8. Picard MarkDuplicates on coordinate-sorted BAM; retain clean FASTQ
#   9. Prepare homologous candidate sites for REDItoolKnown
#  10. Run REDItoolKnown and extract per-sample A/G counts
#  11. Summarize candidate counts across samples for each species
#  12. Classify supported sites as known or novel for each species
#
# Usage:
#   bash main.sh \
#       --sample_table /path/to/samples.tsv \
#       --out_dir /path/to/find_back_result \
#       --cds_pos_dir /path/to/spe_cds_pos \
#       --gtf species=/path/to/species.gtf \
#       --edit_txt_dir /path/to/ortho_edit_txt \
#       --candi_dir /path/to/ortho_edit_txt/candi_per_species \
#       --min_samples 1
#
# Sample table columns:
#   sample_id individual species type fq1 fq2 ref
# ============================================================

set -Eeuo pipefail

# The any-completed-job scheduler uses Bash 5.1's `wait -n -p`.
if (( BASH_VERSINFO[0] < 5 || (BASH_VERSINFO[0] == 5 && BASH_VERSINFO[1] < 1) )); then
    echo "ERROR: find_back sample scheduling requires Bash 5.1 or newer (wait -n -p)." >&2
    exit 1
fi

# ============================================================
# 1. Software configuration
# ============================================================

FASTP="${FASTP:-$(command -v fastp || true)}"
BWA="${BWA:-$(command -v bwa || true)}"
SAMTOOLS="${SAMTOOLS:-$(command -v samtools || true)}"
SAMBAMBA="${SAMBAMBA:-$(command -v sambamba || true)}"
JAVA="${JAVA:-$(command -v java || true)}"

# ============================================================
# 2. Resource configuration
# ============================================================

FASTP_THREADS=2
BWA_THREADS=5
SORT_THREADS=5

MAX_JOBS=2
SORT_MEM="12G"
TOTAL_THREADS=20
REDITool_THREADS=2
REDITool_JOBS=$((TOTAL_THREADS / REDITool_THREADS))
MIN_LENGTH=150
RESUME=0
FORCE_REDITOOL=0
PICARD_HEAP=8g

# ============================================================
# 3. Parse arguments
# ============================================================

SAMPLE_TABLE=""
OUT_DIR=""
CDS_POS_DIR=""
EDIT_TXT_DIR=""
CANDI_DIR=""
declare -A GTF_FILES=()
MIN_SAMPLES=1

while [[ $# -gt 0 ]]; do

    key="$1"

    case "$key" in

        --sample_table)
            SAMPLE_TABLE="$2"
            shift 2
            ;;

        --out_dir)
            OUT_DIR="$2"
            shift 2
            ;;

        --cds_pos_dir)
            CDS_POS_DIR="$2"
            shift 2
            ;;

        --gtf)
            [[ $# -ge 2 && -n "$2" ]] || { echo "ERROR: --gtf requires species=FILE"; exit 1; }
            value="$2"
            [[ "$value" == *=* ]] || { echo "ERROR: --gtf requires species=FILE"; exit 1; }
            species="${value%%=*}"
            gtf_path="${value#*=}"
            [[ "$species" =~ ^[A-Za-z0-9_.-]+$ && -n "$gtf_path" ]] || { echo "ERROR: Invalid --gtf value: $value"; exit 1; }
            [[ -z "${GTF_FILES[$species]+x}" ]] || { echo "ERROR: --gtf supplied more than once for species: $species"; exit 1; }
            GTF_FILES["$species"]="$gtf_path"
            shift 2
            ;;

        --edit_txt_dir)
            EDIT_TXT_DIR="$2"
            shift 2
            ;;

        --candi_dir)
            CANDI_DIR="$2"
            shift 2
            ;;

        --min_samples)
            MIN_SAMPLES="$2"
            shift 2
            ;;

        --mapping_jobs)
            MAX_JOBS="$2"
            shift 2
            ;;

        --fastp_threads)
            FASTP_THREADS="$2"
            shift 2
            ;;

        --bwa_threads)
            BWA_THREADS="$2"
            shift 2
            ;;

        --sort_threads)
            SORT_THREADS="$2"
            shift 2
            ;;

        --sort_mem)
            SORT_MEM="$2"
            shift 2
            ;;

        --threads)
            TOTAL_THREADS="$2"
            shift 2
            ;;

        --min_length)
            MIN_LENGTH="$2"
            shift 2
            ;;

        --resume)
            RESUME=1
            shift
            ;;

        --force_reditool)
            FORCE_REDITOOL=1
            RESUME=1
            shift
            ;;

        --reuse_bam)
            echo "ERROR: --reuse_bam has been removed. Use --force_reditool for validated in-workdir BAM reuse."
            exit 1
            ;;

        -h|--help)
            cat <<'HELP'
Usage:
  find_back --sample_table FILE --out_dir DIR --cds_pos_dir DIR \
    --edit_txt_dir DIR --candi_dir DIR [--gtf species=FILE ...] \
    [--min_samples INT] [--mapping_jobs INT] [--fastp_threads INT] \
    [--bwa_threads INT] [--sort_threads INT] [--sort_mem STRING] \
    [--threads INT] [--min_length INT] [--resume | --force_reditool]

Resource defaults:
  --mapping_jobs 2; fastp 2 threads; BWA 5 threads; Sambamba sort 5 threads;
  sort memory 12G. Mapping/Picard workflows are sample-level and run in parallel.
  A new sample starts as soon as any running sample workflow completes.
  --threads 20 is the REDItoolKnown-stage budget; each REDItoolKnown job uses
  2 threads, so the maximum concurrent jobs are threads / 2. The value must
  be an even integer >= 2. --min_length defaults to 150.
HELP
            exit 0
            ;;

        *)
            echo "ERROR: Unknown option: $1"
            exit 1
            ;;

    esac

done

if [[ -z "$SAMPLE_TABLE" || -z "$OUT_DIR" || -z "$CDS_POS_DIR" || -z "$EDIT_TXT_DIR" || -z "$CANDI_DIR" ]]; then
    echo "ERROR: --sample_table, --out_dir, --cds_pos_dir, --edit_txt_dir and --candi_dir are required"
    exit 1
fi

if ! [[ "$MIN_SAMPLES" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --min_samples must be a positive integer"
    exit 1
fi

for resource_value in "$FASTP_THREADS" "$BWA_THREADS" "$SORT_THREADS" "$MAX_JOBS" "$TOTAL_THREADS" "$MIN_LENGTH"; do
    if ! [[ "$resource_value" =~ ^[1-9][0-9]*$ ]]; then
        echo "ERROR: thread/job resource values must be positive integers"
        exit 1
    fi
done
if (( TOTAL_THREADS < 2 || TOTAL_THREADS % 2 != 0 )); then
    echo "ERROR: --threads must be an even integer >= 2"; exit 1
fi
REDITool_JOBS=$((TOTAL_THREADS / REDITool_THREADS))

if [[ -z "$SORT_MEM" ]]; then
    echo "ERROR: --sort_mem cannot be empty"
    exit 1
fi

if [[ ! -f "$SAMPLE_TABLE" ]]; then
    echo "ERROR: Sample table not found: $SAMPLE_TABLE"
    exit 1
fi

if [[ ! -d "$CDS_POS_DIR" ]]; then
    echo "ERROR: CDS position directory not found: $CDS_POS_DIR"
    exit 1
fi

if [[ ! -d "$EDIT_TXT_DIR" ]]; then
    echo "ERROR: Original known-site directory not found: $EDIT_TXT_DIR"
    exit 1
fi

if [[ ! -d "$CANDI_DIR" ]]; then
    echo "ERROR: Candidate annotation directory not found: $CANDI_DIR"
    exit 1
fi

SAMPLE_TABLE="$(readlink -f "$SAMPLE_TABLE")"
CDS_POS_DIR="$(cd "$CDS_POS_DIR" && pwd)"
EDIT_TXT_DIR="$(cd "$EDIT_TXT_DIR" && pwd)"
CANDI_DIR="$(cd "$CANDI_DIR" && pwd)"

for species in "${!GTF_FILES[@]}"; do
    gtf_path="${GTF_FILES[$species]}"
    if [[ ! -f "$gtf_path" || ! -s "$gtf_path" ]]; then
        echo "ERROR: GTF missing, not a regular file, or empty for species $species: $gtf_path"
        exit 1
    fi
    GTF_FILES["$species"]="$(readlink -f "$gtf_path")"
done

mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"

# ============================================================
# 4. Output directories
# ============================================================

QC_DIR="$OUT_DIR/QC"
MAP_DIR="$OUT_DIR/Mapping_BWA_to_CDS"
LOG_DIR="$OUT_DIR/logs"

mkdir -p "$QC_DIR"
mkdir -p "$MAP_DIR"
mkdir -p "$LOG_DIR"

# Keep the complete module log inside the requested output directory.
exec > >(tee -a "$OUT_DIR/find_back.log") 2>&1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CANDIDATE_SCRIPT="$SCRIPT_DIR/prepare_candidate_sites.py"
REDI_RUN_SCRIPT="$SCRIPT_DIR/run_reditool_known.sh"
REDI_PARALLEL_SCRIPT="$SCRIPT_DIR/run_reditool_known_parallel.py"
SUMMARY_SCRIPT="$SCRIPT_DIR/summarize_candidate_counts.py"
CLASSIFY_SCRIPT="$SCRIPT_DIR/classify_known_novel_sites.py"
GTF_MAP_SCRIPT="$SCRIPT_DIR/map_cds_to_genome.py"
VALIDATE_BAM_SCRIPT="$SCRIPT_DIR/validate_bam_reference.py"

# Project data directories are supplied explicitly by the user. Never infer
# project-data locations from the KOBE installation directory.
STATISTICS_DIR="$OUT_DIR/statistics"
mkdir -p "$STATISTICS_DIR"

# REDItoolKnown configuration
RNA_HOME="${RNA_EDITING_HOME:-/opt/rna_editing}"
REDI_PYTHON="${REDI_PYTHON:-$(command -v python2 || true)}"
REDITool="${REDITool:-$(find "$RNA_HOME" -type f -name REDItoolKnown.py -print -quit 2>/dev/null)}"
INTERSECT_PL="${INTERSECT_PL:-$RNA_HOME/helpers/add.intersect.vcf.info.pl}"
PICARD_JAR="${PICARD_JAR:-$(find "$RNA_HOME" -type f \( -iname 'picard.jar' -o -iname 'picard-*.jar' \) -print 2>/dev/null | sort -V | tail -n 1)}"
REDI_THREADS=2

for required_file in \
    "$CANDIDATE_SCRIPT" \
    "$REDI_RUN_SCRIPT" \
    "$REDI_PARALLEL_SCRIPT" \
    "$SUMMARY_SCRIPT" \
    "$CLASSIFY_SCRIPT" \
    "$GTF_MAP_SCRIPT" \
    "$VALIDATE_BAM_SCRIPT" \
    "$REDITool" \
    "$INTERSECT_PL" \
    "$PICARD_JAR"; do

    if [[ ! -f "$required_file" ]]; then
        echo "ERROR: Required file not found: $required_file"
        exit 1
    fi

done

if [[ ! -d "$CANDI_DIR" ]]; then
    echo "ERROR: Candidate directory not found: $CANDI_DIR"
    exit 1
fi

if [[ ! -d "$EDIT_TXT_DIR" ]]; then
    echo "ERROR: Original known-site directory not found: $EDIT_TXT_DIR"
    exit 1
fi

if [[ ! -x "$REDI_PYTHON" ]]; then
    echo "ERROR: REDItools Python not executable: $REDI_PYTHON"
    exit 1
fi

# ============================================================
# 5. Check software
# ============================================================

echo "==================================================="
echo "Checking software"
echo "==================================================="

for exe in "$FASTP" "$BWA" "$SAMBAMBA" "$JAVA"; do

    if [[ ! -x "$exe" ]]; then
        echo "ERROR: Executable not found: $exe"
        exit 1
    fi

done

if ! command -v "$SAMTOOLS" >/dev/null 2>&1; then
    echo "ERROR: samtools not found"
    exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 not found"
    exit 1
fi

if ! command -v perl >/dev/null 2>&1; then
    echo "ERROR: perl not found"
    exit 1
fi

"$REDI_PYTHON" -c 'import pysam' || {
    echo "ERROR: pysam unavailable in REDItools Python: $REDI_PYTHON"
    exit 1
}

echo "fastp: $FASTP"
echo "BWA: $BWA"
echo "samtools: $(command -v "$SAMTOOLS")"
echo "sambamba: $SAMBAMBA"
echo "Java: $JAVA"
echo "Picard: $PICARD_JAR"
echo "REDItoolKnown: $REDITool"
echo "REDItools Python: $REDI_PYTHON"
echo "REDItools worker: $REDI_RUN_SCRIPT"
echo "REDItools scheduler: $REDI_PARALLEL_SCRIPT"
echo "Total thread budget: $TOTAL_THREADS"
echo "REDItools concurrent jobs: $REDITool_JOBS"
echo "Minimum read length: $MIN_LENGTH"
echo "Resume: $RESUME"
echo "Force REDItoolKnown: $FORCE_REDITOOL"
echo "Automatic BAM reuse enabled by --force_reditool: $FORCE_REDITOOL"
echo "REDItools threads per sample: $REDI_THREADS"
echo "Candidate annotations: $CANDI_DIR"
echo "Original known-site references: $EDIT_TXT_DIR"
echo "Known/novel classification script: $CLASSIFY_SCRIPT"
echo "GTF mappings supplied: ${#GTF_FILES[@]}"
for species in "${!GTF_FILES[@]}"; do echo "  GTF[$species]: ${GTF_FILES[$species]}"; done
echo "CDS-to-genome mapping script: $GTF_MAP_SCRIPT"
echo "Minimum supporting samples: $MIN_SAMPLES"
echo "Statistics output: $STATISTICS_DIR"
echo "Perl intersection: $INTERSECT_PL"
echo "REDI threads: $REDI_THREADS"

# ============================================================
# 6. Validate sample table
# ============================================================

MANIFEST="$LOG_DIR/samples.validated.tsv"

echo ""
echo "==================================================="
echo "Validating sample table"
echo "==================================================="

python3 - "$SAMPLE_TABLE" "$MANIFEST" "$FORCE_REDITOOL" <<'PY'

import csv
import os
import re
import sys
from pathlib import Path

table = Path(sys.argv[1])
output = Path(sys.argv[2])
allow_missing_fastq = sys.argv[3] == "1"

required = [
    "sample_id",
    "individual",
    "species",
    "type",
    "fq1",
    "fq2",
    "ref"
]

def resolve_file(value, base, label, line_number, allow_missing=False):

    value = (value or "").strip()
    if not value:
        if allow_missing:
            return ""
        raise FileNotFoundError(f"Line {line_number}: {label} is empty")

    p = Path(value).expanduser()
    if not p.is_absolute():
        p = base / p
    p = Path(os.path.abspath(str(p)))

    if not p.is_file():
        if allow_missing:
            return str(p)
        raise FileNotFoundError(f"Line {line_number}: {label} not found: {p}")
    return str(p)


seen = set()
rows = []

with open(table, newline="", encoding="utf-8-sig") as f:

    reader = csv.DictReader(f, delimiter="\t")

    if not reader.fieldnames:
        raise ValueError("Empty sample table")

    reader.fieldnames = [
        x.strip() for x in reader.fieldnames
    ]

    missing = [
        x for x in required
        if x not in reader.fieldnames
    ]

    if missing:
        raise ValueError(
            "Missing columns: " + ", ".join(missing)
        )

    for line_number, raw in enumerate(reader, start=2):

        row = {
            k.strip(): (v or "").strip()
            for k, v in raw.items()
            if k is not None
        }

        if None in raw:
            raise ValueError(
                f"Line {line_number}: malformed TSV row"
            )

        sample = row["sample_id"]
        individual = row["individual"]
        species = row["species"]
        sample_type = row["type"].upper()

        if not sample or not individual or not species:
            raise ValueError(
                f"Line {line_number}: missing sample information"
            )

        if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample):
            raise ValueError(
                f"Invalid sample_id: {sample}"
            )

        if not re.fullmatch(r"[A-Za-z0-9_.-]+", species):
            raise ValueError(
                f"Invalid species: {species}"
            )

        if sample_type != "RNA":
            raise ValueError(
                f"Sample {sample}: this module currently accepts RNA only"
            )

        key = (species, sample)

        if key in seen:
            raise ValueError(
                f"Duplicate sample: {species}/{sample}"
            )

        seen.add(key)

        fq1 = resolve_file(
            row["fq1"], table.parent, "fq1", line_number, allow_missing_fastq
        )

        fq2 = resolve_file(
            row["fq2"], table.parent, "fq2", line_number, allow_missing_fastq
        )

        ref = resolve_file(
            row["ref"], table.parent, "ref", line_number
        )

        rows.append([
            sample,
            individual,
            species,
            sample_type,
            fq1 or "__MISSING_FASTQ__",
            fq2 or "__MISSING_FASTQ__",
            ref
        ])

if not rows:
    raise ValueError("No samples found")

with open(output, "w", newline="") as f:

    writer = csv.writer(
        f,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow(required)
    writer.writerows(rows)

print(f"Samples validated: {len(rows)}")
print(f"Validated manifest: {output}")

PY

for gtf_species in "${!GTF_FILES[@]}"; do
    if ! awk -F '\t' -v species="$gtf_species" 'NR > 1 && $3 == species { found=1 } END { exit !found }' "$MANIFEST"; then
        echo "ERROR: --gtf species '$gtf_species' does not appear in the sample table."
        exit 1
    fi
done

# ============================================================
# 7. Build unique reference indexes
# ============================================================

echo ""
echo "==================================================="
echo "Building reference indexes"
echo "==================================================="

REFS=()

while IFS=$'\t' read -r sample individual species sample_type fq1 fq2 ref; do

    [[ -z "$sample" ]] && continue

    found=0

    for existing in "${REFS[@]}"; do

        if [[ "$existing" == "$ref" ]]; then
            found=1
            break
        fi

    done

    if (( found == 0 )); then
        REFS+=("$ref")
    fi

done < <(tail -n +2 "$MANIFEST")


ensure_reference_index() {

    local ref="$1"
    local missing=0
    local signature_dir="$OUT_DIR/state/reference"
    local ref_key
    local signature_file
    local current_signature

    mkdir -p "$signature_dir"
    ref_key="$(printf '%s' "$ref" | sha256sum | awk '{print $1}')"
    signature_file="$signature_dir/${ref_key}.sha256"
    current_signature="$(sha256sum "$ref" | awk '{print $1}')"

    local suffixes=(amb ann bwt pac sa)
    if [[ ! -s "$signature_file" || "$(cat "$signature_file")" != "$current_signature" ]]; then
        echo "Reference FASTA content changed or has no prior signature; rebuilding BWA index:"
        echo "$ref"
        missing=1
        rm -f "${ref}.amb" "${ref}.ann" "${ref}.bwt" "${ref}.pac" "${ref}.sa" "${ref}.fai"
    fi

    for suffix in "${suffixes[@]}"; do
        if [[ ! -s "${ref}.${suffix}" ]]; then
            missing=1
        fi
    done

    if (( missing == 1 )); then
        echo "Building BWA index:"
        echo "$ref"
        "$BWA" index "$ref"
    else
        echo "BWA index exists and matches reference signature:"
        echo "$ref"
    fi

    if [[ ! -s "${ref}.fai" ]]; then
        echo "Building FASTA index:"
        echo "$ref"
        "$SAMTOOLS" faidx "$ref"
    else
        echo "FASTA index exists:"
        echo "$ref"
    fi

    for suffix in "${suffixes[@]}"; do
        if [[ ! -s "${ref}.${suffix}" ]]; then
            echo "ERROR: BWA index incomplete: ${ref}.${suffix}"
            return 1
        fi
    done
    if [[ ! -s "${ref}.fai" ]]; then
        echo "ERROR: FASTA index missing: ${ref}.fai"
        return 1
    fi

    printf '%s\n' "$current_signature" > "$signature_file"
}


for ref in "${REFS[@]}"; do

    ensure_reference_index "$ref"

done

# ============================================================
# 8. Single sample mapping workflow
# ============================================================

run_sample() {

    local sample="$1"
    local individual="$2"
    local species="$3"
    local fq1="$4"
    local fq2="$5"
    local ref="$6"

    echo ""
    echo "==================================================="
    echo "START: $species / $sample"
    echo "Individual: $individual"
    echo "Reference: $ref"
    echo "==================================================="

    local QC_SAMPLE_DIR="$QC_DIR/$species/$sample"
    local MAP_SAMPLE_DIR="$MAP_DIR/$species/$sample"

    mkdir -p "$QC_SAMPLE_DIR" "$MAP_SAMPLE_DIR"

    local CLEAN_R1="$QC_SAMPLE_DIR/${sample}_1.clean.fq.gz"
    local CLEAN_R2="$QC_SAMPLE_DIR/${sample}_2.clean.fq.gz"
    local UNSORTED_BAM="$MAP_SAMPLE_DIR/${sample}.bam"
    local SORTED_BAM="$MAP_SAMPLE_DIR/${sample}.sort.bam"
    local SORTED_BAI="${SORTED_BAM}.bai"
    local SORTED_CSI="${SORTED_BAM}.csi"
    local PICARD_BAM="$MAP_SAMPLE_DIR/${sample}.sort.nodup.bam"
    local PICARD_BAI="${PICARD_BAM}.bai"
    local PICARD_METRICS="$MAP_SAMPLE_DIR/${sample}.picard.metrics.txt"
    local STATE_DIR="$OUT_DIR/state/mapping/$species"
    local QC_STATE="$OUT_DIR/state/qc/$species/${sample}.signature"
    local MAP_STATE="$STATE_DIR/${sample}.signature"
    local VALIDATE_BAM_SCRIPT="$SCRIPT_DIR/validate_bam_reference.py"
    local HAS_RAW_FASTQ=0
    local USE_EXISTING_CLEAN=0
    local INPUT_SIGNATURE=""
    local MAP_SIGNATURE=""

    mkdir -p "$STATE_DIR" "$(dirname "$QC_STATE")"

    if [[ -n "$fq1" && -s "$fq1" && -n "$fq2" && -s "$fq2" ]]; then
        HAS_RAW_FASTQ=1
        INPUT_SIGNATURE="$(printf '%s\n' \
            "$fq1:$(stat -c '%s:%Y' "$fq1")" \
            "$fq2:$(stat -c '%s:%Y' "$fq2")" \
            "$ref:$(stat -c '%s:%Y' "$ref")" \
            "fastp=$FASTP_THREADS;min_length=$MIN_LENGTH" \
            | sha256sum | awk '{print $1}')"
    fi

    if (( HAS_RAW_FASTQ == 1 )); then
        MAP_SIGNATURE="$(printf '%s\n' "$INPUT_SIGNATURE" \
            "bwa=$BWA_THREADS;sort=$SORT_THREADS;sort_mem=$SORT_MEM" \
            "picard=$PICARD_JAR;heap=$PICARD_HEAP;remove_duplicates=true" \
            | sha256sum | awk '{print $1}')"
    elif [[ -s "$CLEAN_R1" && -s "$CLEAN_R2" ]]; then
        USE_EXISTING_CLEAN=1
        MAP_SIGNATURE="$(printf '%s\n' \
            "$CLEAN_R1:$(stat -c '%s:%Y' "$CLEAN_R1")" \
            "$CLEAN_R2:$(stat -c '%s:%Y' "$CLEAN_R2")" \
            "$ref:$(stat -c '%s:%Y' "$ref")" \
            "bwa=$BWA_THREADS;sort=$SORT_THREADS;sort_mem=$SORT_MEM" \
            "picard=$PICARD_JAR;heap=$PICARD_HEAP;remove_duplicates=true" \
            | sha256sum | awk '{print $1}')"
    fi

    validate_final_bam() {
        [[ -s "$PICARD_BAM" ]] || return 1
        python3 "$VALIDATE_BAM_SCRIPT" \
            --bam "$PICARD_BAM" \
            --reference "$ref" \
            --samtools "$SAMTOOLS" \
            --sample-id "$sample"
    }

    # --force_reditool explicitly enables in-workdir BAM reuse. It accepts a
    # sample BAM only when the BAM, index, and current reference are compatible.
    if (( FORCE_REDITOOL == 1 )); then
        echo "Checking existing Picard BAM for --force_reditool: $species/$sample"
        if validate_final_bam; then
            echo "Reusing validated Picard BAM: $PICARD_BAM"
            return 0
        fi
        echo "No valid compatible Picard BAM found; mapping will be performed: $species/$sample"
        if (( HAS_RAW_FASTQ == 0 && USE_EXISTING_CLEAN == 0 )); then
            echo "ERROR: Existing BAM is missing/incomplete/incompatible and neither raw nor retained clean FASTQ is available for $species/$sample"
            return 1
        fi
    elif (( RESUME == 1 )) && [[ -s "$MAP_STATE" ]] && \
         [[ -n "$MAP_SIGNATURE" ]] && [[ "$(cat "$MAP_STATE")" == "$MAP_SIGNATURE" ]] && \
         validate_final_bam; then
        echo "Already completed; Mapping/Picard signature matches: $species/$sample"
        return 0
    fi

    # Remove incomplete or stale alignment products before remapping.
    rm -f "$UNSORTED_BAM" "$SORTED_BAM" "$SORTED_BAI" "$SORTED_CSI" \
          "$PICARD_BAM" "$PICARD_BAI" "${PICARD_BAM%.bam}.bai" \
          "$PICARD_METRICS"

    # --------------------------------------------------------
    # Step 1. fastp QC (or reuse retained clean FASTQ if raw FASTQ is absent)
    # --------------------------------------------------------

    echo ""
    echo "Step 1: fastp QC"

    if (( USE_EXISTING_CLEAN == 1 && HAS_RAW_FASTQ == 0 )); then
        echo "Raw FASTQ unavailable; using retained clean FASTQ files."
    elif (( RESUME == 1 )) && [[ -s "$CLEAN_R1" && -s "$CLEAN_R2" && -s "$QC_STATE" ]] && \
         [[ -n "$INPUT_SIGNATURE" ]] && [[ "$(cat "$QC_STATE")" == "$INPUT_SIGNATURE" ]]; then
        echo "Reusing existing QC FASTQ files; input and QC signature match."
    else
        [[ "$HAS_RAW_FASTQ" == 1 ]] || {
            echo "ERROR: Raw FASTQ is unavailable and retained clean FASTQ is incomplete for $species/$sample"
            return 1
        }
        rm -f "$CLEAN_R1" "$CLEAN_R2"
        "$FASTP" \
            -n 5 \
            -q 15 \
            -u 50 \
            -l "$MIN_LENGTH" \
            --thread "$FASTP_THREADS" \
            -d 4 \
            -i "$fq1" \
            -I "$fq2" \
            -o "$CLEAN_R1" \
            -O "$CLEAN_R2" \
            --disable_trim_poly_g \
            --html "$QC_SAMPLE_DIR/${sample}.fastp.html" \
            --json "$QC_SAMPLE_DIR/${sample}.fastp.json"
    fi

    if [[ ! -s "$CLEAN_R1" || ! -s "$CLEAN_R2" ]]; then
        echo "ERROR: clean FASTQ output missing for $sample"
        return 1
    fi
    if [[ -n "$INPUT_SIGNATURE" ]]; then
        printf '%s\n' "$INPUT_SIGNATURE" > "$QC_STATE"
    fi

    # --------------------------------------------------------
    # Step 2. BWA-MEM mapping
    # --------------------------------------------------------

    echo ""
    echo "Step 2: BWA-MEM mapping"

    local RG="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}"

    "$BWA" mem \
        -t "$BWA_THREADS" \
        -M \
        -R "$RG" \
        "$ref" \
        "$CLEAN_R1" \
        "$CLEAN_R2" \
    | "$SAMTOOLS" view \
        -b \
        -o "$UNSORTED_BAM" \
        -

    if [[ ! -s "$UNSORTED_BAM" ]]; then
        echo "ERROR: BWA BAM output missing for $sample"
        return 1
    fi

    # --------------------------------------------------------
    # Step 3. Coordinate sorting
    # --------------------------------------------------------

    echo ""
    echo "Step 3: Coordinate sorting"

    "$SAMBAMBA" sort \
        -t "$SORT_THREADS" \
        -m "$SORT_MEM" \
        -p \
        --tmpdir "$MAP_SAMPLE_DIR" \
        -o "$SORTED_BAM" \
        "$UNSORTED_BAM"

    if [[ ! -s "$SORTED_BAM" ]]; then
        echo "ERROR: Sorted BAM missing for $sample"
        return 1
    fi

    # --------------------------------------------------------
    # Step 4. Index sorted BAM (retained as an intermediate)
    # --------------------------------------------------------

    echo ""
    echo "Step 4: Index sorted BAM"
    "$SAMBAMBA" index "$SORTED_BAM"
    if [[ ! -s "$SORTED_BAI" && ! -s "$SORTED_CSI" ]]; then
        echo "ERROR: Sorted BAM index missing for $sample"
        return 1
    fi

    # --------------------------------------------------------
    # Step 5. Picard MarkDuplicates
    # --------------------------------------------------------

    echo ""
    echo "Step 5: Picard MarkDuplicates"

    local PICARD_TMP="$MAP_SAMPLE_DIR/picard_tmp"
    mkdir -p "$PICARD_TMP"
    "$JAVA" "-Xmx${PICARD_HEAP}" -jar "$PICARD_JAR" MarkDuplicates \
        TMP_DIR="$PICARD_TMP" \
        INPUT="$SORTED_BAM" \
        OUTPUT="$PICARD_BAM" \
        METRICS_FILE="$PICARD_METRICS" \
        VALIDATION_STRINGENCY=SILENT \
        ASSUME_SORTED=true \
        CREATE_INDEX=true \
        REMOVE_DUPLICATES=true

    if [[ ! -s "$PICARD_BAM" || ! -s "$PICARD_METRICS" ]]; then
        echo "ERROR: Picard MarkDuplicates output or metrics missing for $sample"
        return 1
    fi

    if [[ ! -s "$PICARD_BAI" && ! -s "${PICARD_BAM%.bam}.bai" ]] || \
       ! "$SAMTOOLS" idxstats "$PICARD_BAM" >/dev/null 2>&1; then
        rm -f "$PICARD_BAI" "${PICARD_BAM%.bam}.bai"
        "$SAMTOOLS" index "$PICARD_BAM"
    fi

    # Verify the final duplicate-removed BAM against the current reference
    # using BAM integrity, index readability, sample identity, and SN/LN.
    python3 "$VALIDATE_BAM_SCRIPT" \
        --bam "$PICARD_BAM" \
        --reference "$ref" \
        --samtools "$SAMTOOLS" \
        --sample-id "$sample"

    printf '%s\n' "$MAP_SIGNATURE" > "$MAP_STATE"
    rm -f "$UNSORTED_BAM"

    echo ""
    echo "==================================================="
    echo "DONE: $species / $sample"
    echo "Final Picard BAM: $PICARD_BAM"
    echo "Clean FASTQ retained: $CLEAN_R1 ; $CLEAN_R2"
    echo "==================================================="
}

# ============================================================
# 9. Parallel execution
# ============================================================

echo ""
echo "==================================================="
echo "Starting sample mapping"
echo "==================================================="

echo "MAX_JOBS: $MAX_JOBS"
echo "FASTP_THREADS: $FASTP_THREADS"
echo "BWA_THREADS: $BWA_THREADS"
echo "SORT_THREADS: $SORT_THREADS"
echo "SORT_MEM: $SORT_MEM"

PIDS=()
FAILED=0

cleanup_mapping_jobs() {
    local status=$?
    trap - EXIT INT TERM HUP
    for pid in "${PIDS[@]:-}"; do
        kill -TERM "$pid" 2>/dev/null || true
    done
    for pid in "${PIDS[@]:-}"; do
        wait "$pid" 2>/dev/null || true
    done
    exit "$status"
}
trap cleanup_mapping_jobs EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
trap 'exit 129' HUP

wait_any() {

    local completed_pid=""
    local wait_status=0
    local remaining_pids=()

    # Bash 5.1+ supports -p, which identifies the child reaped by wait -n.
    # This lets the scheduler release the exact completed slot immediately.
    if wait -n -p completed_pid "${PIDS[@]}"; then
        wait_status=0
    else
        wait_status=$?
    fi

    if [[ -z "$completed_pid" ]]; then
        echo "ERROR: wait -n did not identify a completed sample process."
        FAILED=1
        return 0
    fi

    if (( wait_status == 0 )); then
        echo "Job completed: $completed_pid"
    else
        echo "ERROR: Job failed: $completed_pid (exit status: $wait_status)"
        FAILED=1
    fi

    for pid in "${PIDS[@]}"; do
        [[ "$pid" == "$completed_pid" ]] || remaining_pids+=("$pid")
    done
    PIDS=("${remaining_pids[@]}")

}


while IFS=$'\t' read -r sample individual species sample_type fq1 fq2 ref; do

    [[ -z "$sample" ]] && continue
    [[ "$fq1" == "__MISSING_FASTQ__" ]] && fq1=""
    [[ "$fq2" == "__MISSING_FASTQ__" ]] && fq2=""

    run_sample \
        "$sample" \
        "$individual" \
        "$species" \
        "$fq1" \
        "$fq2" \
        "$ref" \
        > "$LOG_DIR/${species}__${sample}.log" \
        2> "$LOG_DIR/${species}__${sample}.error" &

    PIDS+=("$!")

    if (( ${#PIDS[@]} >= MAX_JOBS )); then

        wait_any

        if (( FAILED == 1 )); then
            echo "ERROR: Stopping new jobs because a sample failed."
            break
        fi

    fi

done < <(tail -n +2 "$MANIFEST")


# Wait for remaining jobs; each call reaps whichever sample finishes next.
while (( ${#PIDS[@]} > 0 )); do
    wait_any
done


if (( FAILED == 1 )); then

    echo ""
    echo "ERROR: Some samples failed."
    echo "Check logs in: $LOG_DIR"
    exit 1

fi

# ============================================================
# 10. Prepare species-specific candidate site files
# ============================================================

echo ""
echo "==================================================="
echo "Preparing candidate sites for REDItoolKnown"
echo "==================================================="

echo "CDS position directory: $CDS_POS_DIR"

echo "Candidate preparation script: $CANDIDATE_SCRIPT"

python3 "$CANDIDATE_SCRIPT" \
    --sample_table "$SAMPLE_TABLE" \
    --cds_pos_dir "$CDS_POS_DIR" \
    --out_dir "$OUT_DIR"

# ============================================================
# 11. Run REDItoolKnown sample jobs concurrently
# ============================================================

echo ""
echo "==================================================="
echo "Running sample-level REDItoolKnown parallel jobs"
echo "==================================================="

echo "Validated sample table: $MANIFEST"
echo "Candidate directory: $CANDI_DIR"
echo "REDItoolKnown: $REDITool"
echo "REDItools Python: $REDI_PYTHON"
echo "Concurrent sample processes: $REDITool_JOBS"
echo "Threads per REDItoolKnown process: $REDI_THREADS"
echo "Per-sample logs: $OUT_DIR/logs/reditool_known"

REDI_CMD=(python3 "$REDI_PARALLEL_SCRIPT" --sample_table "$MANIFEST" --out_dir "$OUT_DIR" --candi_dir "$CANDI_DIR" --reditool "$REDITool" --python "$REDI_PYTHON" --intersect_pl "$INTERSECT_PL" --worker_script "$REDI_RUN_SCRIPT" --jobs "$REDITool_JOBS" --threads "$REDITool_THREADS")
if (( RESUME == 1 )); then REDI_CMD+=(--resume); fi
if (( FORCE_REDITOOL == 1 )); then REDI_CMD+=(--force); fi
"${REDI_CMD[@]}"

# ============================================================
# 12. Summarize candidate counts across all samples
# ============================================================

echo ""
echo "==================================================="
echo "Summarizing candidate counts across all samples"
echo "==================================================="

echo "Summary script: $SUMMARY_SCRIPT"
echo "Sample manifest: $MANIFEST"
echo "Count directory: $OUT_DIR/candidate_counts"

python3 "$SUMMARY_SCRIPT" \
    --sample_table "$MANIFEST" \
    --count_dir "$OUT_DIR/candidate_counts"

# ============================================================
# 13. Classify known and novel sites for every species
# ============================================================

echo ""
echo "==================================================="
echo "Classifying known and novel RNA editing sites"
echo "==================================================="

echo "Classification script: $CLASSIFY_SCRIPT"
echo "Sample manifest: $MANIFEST"
echo "Species summary directory: $OUT_DIR/candidate_counts"
echo "Original known-site directory: $EDIT_TXT_DIR"
echo "Minimum supporting samples: $MIN_SAMPLES"
echo "Statistics output directory: $STATISTICS_DIR"

python3 "$CLASSIFY_SCRIPT" \
    --sample_table "$MANIFEST" \
    --summary_dir "$OUT_DIR/candidate_counts" \
    --edit_txt_dir "$EDIT_TXT_DIR" \
    --out_dir "$STATISTICS_DIR" \
    --min_samples "$MIN_SAMPLES"

# ============================================================
# 14. Append genomic coordinates to the three result tables
# ============================================================

echo ""
echo "==================================================="
echo "Appending genomic coordinates for species with supplied GTFs"
echo "==================================================="

if (( ${#GTF_FILES[@]} > 0 )); then
    GTF_MAP_ARGS=()
    for species in "${!GTF_FILES[@]}"; do
        GTF_MAP_ARGS+=(--gtf "$species=${GTF_FILES[$species]}")
    done
    python3 "$GTF_MAP_SCRIPT" \
        --sample_table "$MANIFEST" \
        --summary_dir "$OUT_DIR/candidate_counts" \
        --statistics_dir "$STATISTICS_DIR" \
        "${GTF_MAP_ARGS[@]}" \
        --min_samples "$MIN_SAMPLES" \
        --out_dir "$OUT_DIR/genomic_coordinates"
else
    echo "No --gtf arguments supplied; genomic-coordinate mapping skipped."
    echo "CDS-level candidate counts and Known/Novel results remain available."
fi
# ============================================================
# 15. Final summary
# ============================================================

echo ""
echo "==================================================="
echo "RNA-only homologous editing recovery completed!"
echo "==================================================="

echo "Sample table:"
echo "$SAMPLE_TABLE"

echo ""
echo "Output directory:"
echo "$OUT_DIR"

echo ""
echo "QC results:"
echo "$QC_DIR"

echo ""
echo "BWA mapping and Picard duplicate-removal results:"
echo "$MAP_DIR"

echo ""
echo "Candidate site files:"
echo "$OUT_DIR/known_sites"
echo "Sample-level candidate site copies:"
echo "$MAP_DIR/<species>/<sample>/knownEditingSites.tab"

echo ""
echo "Per-sample A/G count tables:"
find "$OUT_DIR/candidate_counts" \
    -type f \
    -name "*.candidate_counts.tsv" \
    | sort

echo ""
echo "Species-level multi-sample summary tables:"
find "$OUT_DIR/candidate_counts" \
    -type f \
    -name "*_all_samples_summary.tsv" \
    | sort

echo ""
echo "Raw Perl intersection results:"
find "$OUT_DIR" \
    -type f \
    -name "*_raw_*.txt" \
    | sort

echo ""
echo "Known/novel site statistics:"
find "$STATISTICS_DIR" \
    -type f \
    -name "*editing_site_statistics*.tsv" \
    | sort

echo ""
echo "Known/novel supported-site tables:"
find "$STATISTICS_DIR" \
    -type f \
    -name "*supported_sites*.tsv" \
    | sort

echo ""
echo "Genomic-coordinate result tables:"
find "$OUT_DIR" "$STATISTICS_DIR" \
    -type f \
    -name "*.genomic.tsv" \
    | sort

echo ""
echo "Logs:"
echo "$LOG_DIR"
echo "$OUT_DIR/find_back.log"
echo "$OUT_DIR/logs/reditool_known/"

echo ""
echo "Final BAM files:"

find "$MAP_DIR" \
    -type f \
    -name "*.sort.nodup.bam" \
    | sort

echo ""
echo "Total sorted BAM files:"
find "$MAP_DIR" -type f -name "*.sort.nodup.bam" | wc -l

echo ""
echo "ALL DONE"
