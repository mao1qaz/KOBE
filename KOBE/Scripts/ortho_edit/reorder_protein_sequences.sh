#!/usr/bin/env bash
set -Eeuo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 DIR [INT]" >&2
    exit 1
fi

workspace="$1"
THREADS="${2:-8}"

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: threads must be a positive integer." >&2
    exit 1
fi

input_dir="${workspace}/pep_every"
output_dir="${workspace}/re_pep"

[[ -d "$input_dir" ]] || {
    echo "ERROR: Input directory not found: $input_dir" >&2
    exit 1
}

mkdir -p "$output_dir"

shopt -s nullglob
pep_files=("$input_dir"/*.pep.fa)
shopt -u nullglob

total=${#pep_files[@]}
if (( total == 0 )); then
    echo "ERROR: No .pep.fa files found in $input_dir" >&2
    exit 1
fi

MAX_JOBS="$THREADS"
if (( MAX_JOBS > total )); then
    MAX_JOBS="$total"
fi

process_file() {
    set -Eeuo pipefail

    local f="$1"
    local base
    local out
    local tmp_out

    base=$(basename "$f")
    out="$output_dir/$base"
    tmp_out="${out}.tmp.$$"

    trap 'rm -f -- "$tmp_out"' EXIT

    awk '
    BEGIN {RS=">"; ORS=""}
    NR>1 {
        header = $1
        seq = ""
        sub(/\n.*/,"",header)
        n = split($0, lines, "\n")
        for(i=2;i<=n;i++) seq = seq lines[i]
        gsub(/[\r\n]/,"",seq)
        if(seq ~ /[^X\-]/) real[++r] = ">"$0
        else xxx[++x] = ">"$0
    }
    END {
        for(i=1;i<=r;i++) print real[i] "\n"
        for(i=1;i<=x;i++) print xxx[i] "\n"
    }
    ' "$f" > "$tmp_out"

    mv -f -- "$tmp_out" "$out"
    trap - EXIT

    echo "Reordered $base -> $out"
}

export -f process_file
export output_dir

echo "Reordering $total protein files using up to $MAX_JOBS processes (1 thread per process)."

printf '%s\0' "${pep_files[@]}" |
    xargs -0 -r -n 1 -P "$MAX_JOBS" bash -c 'process_file "$1"' _

echo "All protein files reordered."
