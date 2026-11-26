#!/bin/bash
# ===================================================
# Extract sequences from besthit_matrix.tsv
# FINAL VERSION:
# - NO query row output
# - matrix order preserved
# - query species column replaced by query sequence
# ===================================================

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 DIR [INT]"
    exit 1
fi

work_dir="$1"
THREADS="${2:-8}"

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: threads must be a positive integer."
    exit 1
fi
pep_dir="${work_dir}/pep"
tsv_file="${work_dir}/besthit_matrix.tsv"
output_dir="${work_dir}/pep_every"

mkdir -p "$output_dir"

debug_log="${output_dir}/pep_extraction_debug.log"
> "$debug_log"

MAX_JOBS="$THREADS"


# ===================================================
# STEP 1: build fast index (protein_id -> species)
# ===================================================
species_map="${output_dir}/query_species.map"
> "$species_map"

echo "Building index..." >> "$debug_log"

for fasta in "$pep_dir"/*_longestpep.fasta; do

    species=$(basename "$fasta" | cut -d'_' -f1)

    awk -v sp="$species" '
    /^>/{
        split($0,a,"|");
        id=a[2];
        if(id!="") print id"\t"sp;
    }' "$fasta" >> "$species_map"

done


# ===================================================
# READ HEADER
# ===================================================
IFS=$'\t' read -r -a headers < <(head -n 1 "$tsv_file")


# ===================================================
# PROCESS FUNCTION
# ===================================================
process_line() {

    local query_id="$1"
    shift
    local fields=("$@")

    # find query species
    local query_species
    query_species=$(awk -v id="$query_id" '$1==id{print $2; exit}' "$species_map")

    local out_file="${output_dir}/${query_id}.pep.fa"
    > "$out_file"

    local query_seq=""

    # ===================================================
    # STEP 1: get query sequence (NO OUTPUT)
    # ===================================================
    if [[ -n "$query_species" ]]; then

        local query_pep="${pep_dir}/${query_species}_longestpep.fasta"

        full_id="${query_species}|${query_id}"

        query_seq=$(seqtk subseq "$query_pep" <(echo "$full_id") 2>/dev/null \
            | grep -v '^>' | tr -d '\n')

    fi


    # ===================================================
    # STEP 2: matrix traversal (ORDERED)
    # ===================================================
    for ((i=1;i<${#headers[@]};i++)); do

        local protein_entry="${fields[$i]}"

        # skip empty
        if [[ "$protein_entry" == "nohit" || -z "$protein_entry" ]]; then
            echo ">NO_HIT" >> "$out_file"
            echo "XXXXXXXXXXXXXXXX" >> "$out_file"
            continue
        fi

        local short_name="${protein_entry%%|*}"

        # ===================================================
        # REPLACE query species column
        # ===================================================
        if [[ "$short_name" == "$query_species" ]]; then

            echo ">${query_species}|${query_id}" >> "$out_file"
            echo "$query_seq" >> "$out_file"

            continue
        fi

        local pep_file="${pep_dir}/${short_name}_longestpep.fasta"

        if [[ ! -f "$pep_file" ]]; then
            echo ">${protein_entry}|PEP_FILE_MISSING" >> "$out_file"
            echo "XXXXXXXXXXXXXXXX" >> "$out_file"
            continue
        fi

        local seq
        seq=$(seqtk subseq "$pep_file" <(echo "$protein_entry") 2>/dev/null \
            | grep -v '^>' | tr -d '\n')

        if [[ -z "$seq" ]]; then
            echo ">${protein_entry}|NOT_FOUND" >> "$out_file"
            echo "XXXXXXXXXXXXXXXX" >> "$out_file"
        else
            echo ">${protein_entry}" >> "$out_file"
            echo "$seq" >> "$out_file"
        fi

    done

    echo "Processed $query_id" >> "$debug_log"
}


export -f process_line
export pep_dir output_dir debug_log species_map


# ===================================================
# PARALLEL EXECUTION
# ===================================================
{
    read -r header

    while IFS=$'\t' read -r -a fields; do

        [[ -z "${fields[0]}" ]] && continue

        query_id="${fields[0]}"

        process_line "$query_id" "${fields[@]}" &

        while (( $(jobs | wc -l) >= MAX_JOBS )); do
            sleep 0.1
        done

    done

} < "$tsv_file"

wait

echo "DONE: output in $output_dir"