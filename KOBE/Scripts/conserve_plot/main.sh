#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="/opt/rna_editing/envs/analysis_py39/bin/python"

ANNOTATED=""
CDS_DIR=""
OUT_DIR=""
MIN_SPECIES=""
SPECIES_ORDER=""
SOURCE_PRIORITY=""
COL="8"
GTF_ARGS=()

usage() {
    cat <<'EOF'
Usage: bash main.sh \
  --annotated FILE --cds_dir DIR --outdir DIR \
  --min_species INT --species_order STRING \
  [--source_priority STRING] [--gtf STRING ... (SPECIES=FILE)] [--col INT]
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --annotated) ANNOTATED="$2"; shift 2 ;;
        --cds_dir) CDS_DIR="$2"; shift 2 ;;
        --outdir) OUT_DIR="$2"; shift 2 ;;
        --min_species) MIN_SPECIES="$2"; shift 2 ;;
        --species_order) SPECIES_ORDER="$2"; shift 2 ;;
        --source_priority) SOURCE_PRIORITY="$2"; shift 2 ;;
        --col) COL="$2"; shift 2 ;;
        --gtf) GTF_ARGS+=(--gtf "$2"); shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

if [[ -z "$ANNOTATED" || -z "$CDS_DIR" || -z "$OUT_DIR" || -z "$MIN_SPECIES" || -z "$SPECIES_ORDER" ]]; then
    echo "ERROR: --annotated, --cds_dir, --outdir, --min_species and --species_order are required." >&2
    usage >&2
    exit 1
fi

mkdir -p "$OUT_DIR"
ALL_SITE_FILE="$OUT_DIR/conserved_site_all_species_9nt_3aa.tsv"
FILTERED_SITE_FILE="$OUT_DIR/conserved_sites_ge${MIN_SPECIES}species.tsv"

printf '%s\n' "===================================================" "Conserved RNA editing analysis pipeline" "==================================================="
echo "ANNOTATED: $ANNOTATED"
echo "CDS_DIR: $CDS_DIR"
echo "OUT_DIR: $OUT_DIR"
echo "MIN_SPECIES: $MIN_SPECIES"
echo "SPECIES_ORDER: $SPECIES_ORDER"
echo "SOURCE_PRIORITY: ${SOURCE_PRIORITY:-annotation column order}"

echo "Step 1: Extract conserved 9nt and 3aa sequences"
EXTRACT_ARGS=(--annotated "$ANNOTATED" --cds_dir "$CDS_DIR" --out_dir "$OUT_DIR")
if [[ -n "$SOURCE_PRIORITY" ]]; then EXTRACT_ARGS+=(--source_priority "$SOURCE_PRIORITY"); fi
"$PYTHON" "$SCRIPT_DIR/01_extract_conserved_9nt_3aa.py" "${EXTRACT_ARGS[@]}"

echo "Step 2: Filter conserved sites"
"$PYTHON" "$SCRIPT_DIR/02_filter_conserved_sites.py" \
    --out_dir "$OUT_DIR" --min_species "$MIN_SPECIES"

echo "Step 3: Plot all conserved transcripts"
"$PYTHON" "$SCRIPT_DIR/03_plot_conserved_sites.py" \
    --all --site_file "$FILTERED_SITE_FILE" --annotated "$ANNOTATED" \
    --cds_dir "$CDS_DIR" --species_order "$SPECIES_ORDER" --out_dir "$OUT_DIR" --col "$COL"

echo "Step 4: Summarize conserved sites and exact species-combination regions"
"$PYTHON" "$SCRIPT_DIR/04_summary_conserved_sites_venn.py" \
    --annotated "$ANNOTATED" --cds_dir "$CDS_DIR" --out_dir "$OUT_DIR" \
    --min_species "$MIN_SPECIES" --species_order "$SPECIES_ORDER"

echo "Step 5: Export detailed conserved site coordinates"
"$PYTHON" "$SCRIPT_DIR/05_export_conserved_site_coordinates.py" \
    --out_dir "$OUT_DIR" --site_file "$FILTERED_SITE_FILE" \
    --annotated "$ANNOTATED" --min_species "$MIN_SPECIES" \
    --species_order "$SPECIES_ORDER" "${GTF_ARGS[@]}"

echo "Step 6: Remove intermediate files after successful completion"
for intermediate_file in \
    "$FILTERED_SITE_FILE" \
    "$OUT_DIR/conserved_sites_ge${MIN_SPECIES}species_transcripts.tsv" \
    "$OUT_DIR/conserved_site_all_species_unmatched.tsv" \
    "$OUT_DIR/conserved_sites_ge${MIN_SPECIES}species_venn.pdf" \
    "$OUT_DIR/conserved_sites_ge${MIN_SPECIES}species_venn.png"
do
    if [[ -f "$intermediate_file" ]]; then
        rm -f "$intermediate_file"
        echo "Removed intermediate/obsolete file: $intermediate_file"
    fi
done

echo "==================================================="
echo "Conserved RNA editing pipeline completed"
echo "All conserved sites: $ALL_SITE_FILE"
echo "Filtered sites: $FILTERED_SITE_FILE"
echo "Output directory: $OUT_DIR"
echo "==================================================="
