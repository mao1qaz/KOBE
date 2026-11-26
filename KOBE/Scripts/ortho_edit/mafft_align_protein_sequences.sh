#!/usr/bin/env bash
set -Eeuo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 DIR [INT]" >&2
    exit 1
fi

work_dir="$1"
THREADS="${2:-8}"

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: threads must be a positive integer." >&2
    exit 1
fi

input_dir="${work_dir}/re_pep"
output_dir="${work_dir}/aln_pep"
log_dir="${output_dir}/logs"

if (( THREADS == 1 )); then
    threads_per_mafft=1
    parallel_jobs=1
else
    threads_per_mafft=2
    parallel_jobs=$((THREADS / 2))
fi

mkdir -p "$output_dir" "$log_dir"

shopt -s nullglob
pep_files=("$input_dir"/*.pep.fa)
shopt -u nullglob

total=${#pep_files[@]}
if (( total == 0 )); then
    echo "ERROR: No .pep.fa files found in $input_dir" >&2
    exit 1
fi

if (( parallel_jobs > total )); then
    parallel_jobs="$total"
fi

process_file() {
    set -Eeuo pipefail

    local pep_file="$1"
    local seq_id
    local safe_id
    local aln_file
    local log_file
    local tmp_in
    local tmp_out
    local tmp_aln

    seq_id=$(basename "$pep_file" .pep.fa)
    safe_id=$(printf '%s' "$seq_id" | sed 's/|/_/g; s/\//_/g')

    aln_file="${output_dir}/${safe_id}.aln.fasta"
    log_file="${log_dir}/${safe_id}.mafft.log"

    if [[ -s "$aln_file" ]]; then
        echo "Skip $seq_id (alignment exists)"
        return 0
    fi

    tmp_in=$(mktemp)
    tmp_out=$(mktemp)
    tmp_aln="${aln_file}.tmp.$$"
    trap 'rm -f -- "$tmp_in" "$tmp_out" "$tmp_aln"' EXIT

    awk '{
        if(/^>/) {
            print
            next
        }
        if(/^X+$/) {
            print "XXXXXXXXXX"
        } else {
            print
        }
    }' "$pep_file" > "$tmp_in"

    if [[ ! -s "$tmp_in" ]]; then
        echo "$seq_id: empty after preprocessing, skip"
        return 0
    fi

    echo "Aligning $seq_id (MAFFT threads=$threads_per_mafft)"

    if ! mafft --thread "$threads_per_mafft" --maxiterate 1000 --localpair "$tmp_in" \
        > "$tmp_out" 2> "$log_file"; then
        echo "ERROR: MAFFT failed for $seq_id; see $log_file" >&2
        return 1
    fi

    awk '{
        if(/^>.*NO_HIT/) {
            print
            getline
            gsub(/[^-]/, "X")
            print
        } else {
            print
        }
    }' "$tmp_out" > "$tmp_aln"

    mv -f -- "$tmp_aln" "$aln_file"
    rm -f -- "$tmp_in" "$tmp_out"
    trap - EXIT

    echo "Done $seq_id"
}

export -f process_file
export output_dir log_dir threads_per_mafft

echo "Starting $total MAFFT alignments: up to $parallel_jobs processes x $threads_per_mafft threads."

printf '%s\0' "${pep_files[@]}" |
    xargs -0 -r -n 1 -P "$parallel_jobs" bash -c 'process_file "$1"' _

echo "All MAFFT alignments completed."
echo "Results: $output_dir"
echo "Logs: $log_dir"
