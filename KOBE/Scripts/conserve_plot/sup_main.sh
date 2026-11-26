#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="/opt/rna_editing/envs/analysis_py39/bin/python"

ANNOTATED=""
CDS_DIR=""
OUT_ROOT=""
MIN_SPECIES=""
GROUP_COL=""
SPECIES_ORDER=""
SOURCE_PRIORITY=""
NCOLS=8

usage() {
    cat <<'EOF'
Usage: bash sup_main.sh \
  --annotated FILE --cds_dir DIR --outdir DIR \
  --group_col STRING --min_species INT \
  --species_order STRING [--source_priority STRING] [--ncols INT]
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --annotated) ANNOTATED="$2"; shift 2 ;;
        --cds_dir) CDS_DIR="$2"; shift 2 ;;
        --outdir) OUT_ROOT="$2"; shift 2 ;;
        --group_col) GROUP_COL="$2"; shift 2 ;;
        --min_species) MIN_SPECIES="$2"; shift 2 ;;
        --species_order) SPECIES_ORDER="$2"; shift 2 ;;
        --source_priority) SOURCE_PRIORITY="$2"; shift 2 ;;
        --ncols) NCOLS="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: Unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

if [[ -z "$ANNOTATED" || -z "$CDS_DIR" || -z "$OUT_ROOT" || -z "$MIN_SPECIES" || -z "$GROUP_COL" || -z "$SPECIES_ORDER" ]]; then
    echo "ERROR: Required arguments are missing." >&2
    usage >&2
    exit 1
fi
if [[ "$GROUP_COL" != "ID2_group" && "$GROUP_COL" != "new_id2_group" ]]; then
    echo "ERROR: --group_col must be ID2_group or new_id2_group" >&2
    exit 1
fi

RUN_TAG="${GROUP_COL}_ge${MIN_SPECIES}species"
OUT_DIR="$OUT_ROOT"
mkdir -p "$OUT_DIR"
FILTERED_SITE_FILE="$OUT_DIR/sup_conserved_sites_${RUN_TAG}.tsv"

echo "=================================================="
echo "Recovered-site conserved RNA editing workflow"
echo "=================================================="
echo "ANNOTATED: $ANNOTATED"
echo "CDS_DIR: $CDS_DIR"
echo "GROUP_COL: $GROUP_COL"
echo "MIN_SPECIES: $MIN_SPECIES"
echo "SPECIES_ORDER: $SPECIES_ORDER"
echo "OUT_DIR: $OUT_DIR"

EXTRACT_ARGS=(--anno_file "$ANNOTATED" --fasta_dir "$CDS_DIR" --out_dir "$OUT_DIR")
if [[ -n "$SOURCE_PRIORITY" ]]; then EXTRACT_ARGS+=(--source_priority "$SOURCE_PRIORITY"); fi

echo "Step 1: Extract homologous CDS windows"
"$PYTHON" "$SCRIPT_DIR/sup_01_extract_conserved_9nt_3aa.py" "${EXTRACT_ARGS[@]}"

echo "Step 2: Filter sites using $GROUP_COL >= $MIN_SPECIES species"
"$PYTHON" "$SCRIPT_DIR/sup_02_filter_conserved_sites.py" \
    --out_dir "$OUT_DIR" --group_col "$GROUP_COL" --min_species "$MIN_SPECIES"

echo "Step 3: Plot filtered transcripts"
"$PYTHON" "$SCRIPT_DIR/sup_03_plot_conserved_sites.py" \
    --all --site_file "$FILTERED_SITE_FILE" --anno_file "$ANNOTATED" \
    --fasta_dir "$CDS_DIR" --species_order "$SPECIES_ORDER" \
    --out_dir "$OUT_DIR" --ncols "$NCOLS"

echo "=================================================="
echo "Recovered-site plotting workflow completed"
echo "Output directory: $OUT_DIR"
echo "=================================================="
