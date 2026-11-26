#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

CDS_DIR=""
EDIT_TXT_DIR=""
WORK_DIR=""
THREADS=8

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cds_dir)
            [[ $# -ge 2 ]] || { echo "ERROR: --cds_dir requires a value" >&2; exit 1; }
            CDS_DIR="$2"; shift 2 ;;
        --edit_txt_dir)
            [[ $# -ge 2 ]] || { echo "ERROR: --edit_txt_dir requires a value" >&2; exit 1; }
            EDIT_TXT_DIR="$2"; shift 2 ;;
        --work_dir)
            [[ $# -ge 2 ]] || { echo "ERROR: --work_dir requires a value" >&2; exit 1; }
            WORK_DIR="$2"; shift 2 ;;
        --threads)
            [[ $# -ge 2 ]] || { echo "ERROR: --threads requires a value" >&2; exit 1; }
            THREADS="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: main.sh --cds_dir DIR --edit_txt_dir DIR --work_dir DIR [--threads INT]"
            exit 0 ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            exit 1 ;;
    esac
done

if [[ -z "$CDS_DIR" || -z "$EDIT_TXT_DIR" || -z "$WORK_DIR" ]]; then
    echo "ERROR: --cds_dir, --edit_txt_dir and --work_dir are required" >&2
    exit 1
fi

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --threads must be a positive integer." >&2
    exit 1
fi

[[ -d "$CDS_DIR" ]] || { echo "ERROR: CDS directory not found: $CDS_DIR" >&2; exit 1; }
[[ -d "$EDIT_TXT_DIR" ]] || { echo "ERROR: Editing-site input directory not found: $EDIT_TXT_DIR" >&2; exit 1; }

CDS_DIR="$(cd "$CDS_DIR" && pwd -P)"
EDIT_TXT_DIR="$(cd "$EDIT_TXT_DIR" && pwd -P)"
mkdir -p "$WORK_DIR"
WORK_DIR="$(cd "$WORK_DIR" && pwd -P)"

path_is_within() {
    local parent="${1%/}"
    local child="${2%/}"
    [[ "$child" == "$parent" || "$child" == "$parent/"* ]]
}

paths_overlap() {
    path_is_within "$1" "$2" || path_is_within "$2" "$1"
}

if [[ "$WORK_DIR" == "/" || "$WORK_DIR" == "$HOME" ]]; then
    echo "ERROR: Unsafe work directory: $WORK_DIR" >&2
    exit 1
fi

for input_dir in "$CDS_DIR" "$EDIT_TXT_DIR"; do
    if paths_overlap "$WORK_DIR" "$input_dir"; then
        echo "ERROR: WORK_DIR must not be the same as, inside, or above an input directory." >&2
        echo "WORK_DIR: $WORK_DIR" >&2
        echo "Input:    $input_dir" >&2
        exit 1
    fi
done

shopt -s nullglob
CDS_FILES=("$CDS_DIR"/*_longestcds.fasta)
EDIT_FILES=("$EDIT_TXT_DIR"/*_candi_cds_sites.txt)
shopt -u nullglob

(( ${#CDS_FILES[@]} > 0 )) || { echo "ERROR: No *_longestcds.fasta files found in $CDS_DIR" >&2; exit 1; }
(( ${#EDIT_FILES[@]} > 0 )) || { echo "ERROR: No *_candi_cds_sites.txt files found in $EDIT_TXT_DIR" >&2; exit 1; }

CODON_TABLE="${RNA_EDITING_CODON_TABLE:-/opt/rna_editing/share/Codon_table.new.txt}"
[[ -s "$CODON_TABLE" ]] || { echo "ERROR: Codon table not found: $CODON_TABLE" >&2; exit 1; }

export PYTHONIOENCODING="${PYTHONIOENCODING:-utf-8}"
export MAFFT_BINARIES="${MAFFT_BINARIES:-/opt/rna_editing/mafft/libexec/mafft}"

for cmd in python3 Rscript diamond seqtk mafft xargs awk; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: Required command not found in container: $cmd" >&2
        exit 1
    }
done

for script in \
    build_besthit_matrix.py \
    extract_besthit_protein_sequences.sh \
    reorder_protein_sequences.sh \
    mafft_align_protein_sequences.sh \
    back_translate_protein_alignment.py \
    merge_orthologous_editing_sites.py \
    group_shared_orthologous_sites.R \
    check_cds_codon.py \
    extract_cds_pos.py \
    annotate_aa_change.py \
    generate_candi_per_species.py
do
    [[ -s "$SCRIPT_DIR/$script" ]] || {
        echo "ERROR: Required script is missing: $SCRIPT_DIR/$script" >&2
        exit 1
    }
done

EDIT_RESULT_DIR="$WORK_DIR/ortho_edit_result"

echo "==================================================="
echo "ORTHOLOGOUS RNA EDITING ANALYSIS"
echo "==================================================="
echo "CDS directory:       $CDS_DIR"
echo "Editing-site input:  $EDIT_TXT_DIR"
echo "Working directory:   $WORK_DIR"
echo "Editing results:     $EDIT_RESULT_DIR"
echo "CPU thread budget:   $THREADS"
echo "Codon table:         $CODON_TABLE"
echo "==================================================="

echo
echo "Cleaning previous workflow-generated results..."
for path in \
    "$WORK_DIR/pep" \
    "$WORK_DIR/besthit_results" \
    "$WORK_DIR/besthit_matrix.tsv" \
    "$WORK_DIR/non_multiple_of_3.txt" \
    "$WORK_DIR/pep_every" \
    "$WORK_DIR/re_pep" \
    "$WORK_DIR/aln_pep" \
    "$WORK_DIR/aln_cds" \
    "$WORK_DIR/spe_cds_pos" \
    "$EDIT_RESULT_DIR"
do
    rm -rf -- "$path"
done

mkdir -p "$EDIT_RESULT_DIR"

echo "Cleanup completed. All workflow steps will be rerun."

echo
echo "Step 1: Build best-hit matrix (DIAMOND threads=$THREADS)"
python3 "$SCRIPT_DIR/build_besthit_matrix.py" \
    --cds_dir "$CDS_DIR" \
    --edit_txt_dir "$EDIT_TXT_DIR" \
    --work_dir "$WORK_DIR" \
    --threads "$THREADS"

echo
echo "Step 2: Extract homologous protein sequences (up to $THREADS processes)"
bash "$SCRIPT_DIR/extract_besthit_protein_sequences.sh" "$WORK_DIR" "$THREADS"

echo
echo "Step 3: Reorder protein sequences (up to $THREADS processes)"
bash "$SCRIPT_DIR/reorder_protein_sequences.sh" "$WORK_DIR" "$THREADS"

echo
echo "Step 4: Align protein sequences with MAFFT"
bash "$SCRIPT_DIR/mafft_align_protein_sequences.sh" "$WORK_DIR" "$THREADS"

echo
echo "Step 5: Back-translate protein alignments to CDS (up to $THREADS processes)"
python3 "$SCRIPT_DIR/back_translate_protein_alignment.py" \
    --work_dir "$WORK_DIR" \
    --cds_dir "$CDS_DIR" \
    --threads "$THREADS"

echo
echo "Step 6: Project editing sites onto aligned CDS"
python3 "$SCRIPT_DIR/merge_orthologous_editing_sites.py" \
    --work_dir "$WORK_DIR" \
    --edit_txt_dir "$EDIT_TXT_DIR" \
    --result_dir "$EDIT_RESULT_DIR"

echo
echo "Step 7: Group and deduplicate shared orthologous sites"
Rscript "$SCRIPT_DIR/group_shared_orthologous_sites.R" \
    "$EDIT_RESULT_DIR/merged_edit_sites_codon_table.csv" \
    "$EDIT_RESULT_DIR/merged_edit_sites_codon_table_uniq.csv"

echo
echo "Step 8: Check CDS codon consistency"
python3 "$SCRIPT_DIR/check_cds_codon.py" \
    --work_dir "$WORK_DIR" \
    --cds_dir "$CDS_DIR" \
    --result_dir "$EDIT_RESULT_DIR"

echo
echo "Step 9: Extract CDS positions"
python3 "$SCRIPT_DIR/extract_cds_pos.py" \
    --work_dir "$WORK_DIR" \
    --result_dir "$EDIT_RESULT_DIR"

echo
echo "Step 10: Annotate amino-acid changes"
python3 "$SCRIPT_DIR/annotate_aa_change.py" \
    --result_dir "$EDIT_RESULT_DIR" \
    --codon_table "$CODON_TABLE"

echo
echo "Step 11: Generate candidate editing sites per species"
python3 "$SCRIPT_DIR/generate_candi_per_species.py" \
    --result_dir "$EDIT_RESULT_DIR" \
    --codon_table "$CODON_TABLE"

echo
echo "==================================================="
echo "ORTHOLOGOUS RNA EDITING ANALYSIS COMPLETED"
echo "==================================================="
echo "Best-hit matrix:      $WORK_DIR/besthit_matrix.tsv"
echo "Protein alignments:   $WORK_DIR/aln_pep"
echo "Aligned CDS:          $WORK_DIR/aln_cds"
echo "CDS positions:        $WORK_DIR/spe_cds_pos"
echo "Editing-site results: $EDIT_RESULT_DIR"
echo "Merged table:         $EDIT_RESULT_DIR/merged_edit_sites_codon_table.csv"
echo "Unique table:         $EDIT_RESULT_DIR/merged_edit_sites_codon_table_uniq.csv"
echo "Annotated table:      $EDIT_RESULT_DIR/merged_edit_sites_codon_table_annotated.csv"
echo "Candidate sites:      $EDIT_RESULT_DIR/candi_per_species"
echo "==================================================="
