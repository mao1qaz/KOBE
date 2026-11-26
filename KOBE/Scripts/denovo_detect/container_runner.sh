#!/usr/bin/env bash
# Shared Apptainer path and symlink handling for denovo_detect commands.
set -Eeuo pipefail

COMMON_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="$(cd "$COMMON_DIR/.." && pwd)"
declare -ag BIND_ARGS=()
declare -Ag _AUTO_BIND_SEEN=()
declare -ag USER_BIND_SPECS=()

require_host_tools() {
    command -v apptainer >/dev/null 2>&1 || {
        echo "ERROR: apptainer not found on the host." >&2
        exit 1
    }
    command -v python3 >/dev/null 2>&1 || {
        echo "ERROR: host python3 is required for path and sample-table handling." >&2
        exit 1
    }
    command -v readlink >/dev/null 2>&1 || {
        echo "ERROR: readlink is required on the host." >&2
        exit 1
    }
}

absolute_path() {
    python3 - "$1" <<'PY'
import os
import sys
print(os.path.abspath(os.path.expanduser(sys.argv[1])))
PY
}

_add_auto_bind() {
    local source="$1" destination="$2" key="$1:$2"
    [[ -d "$source" ]] || return 0
    if [[ -z "${_AUTO_BIND_SEEN[$key]+x}" ]]; then
        BIND_ARGS+=(--bind "$source:$destination")
        _AUTO_BIND_SEEN[$key]=1
    fi
}

add_bind_for_path() {
    local path="$1" lexical resolved candidate source
    [[ -n "$path" ]] || return 0
    lexical="$(absolute_path "$path")"
    resolved="$(readlink -f -- "$lexical" 2>/dev/null || true)"
    for candidate in "$lexical" "$resolved"; do
        [[ -n "$candidate" && "$candidate" == /* ]] || continue
        source=""
        if [[ "$candidate" == "$lexical" && -n "$resolved" && "$resolved" != "$lexical" ]]; then
            if [[ -d "$lexical" && ! -L "$lexical" ]]; then source="$lexical"; else source="$(dirname "$lexical")"; fi
        else
            if command -v findmnt >/dev/null 2>&1 && [[ -e "$candidate" ]]; then
                source="$(findmnt -T "$candidate" -n -o TARGET 2>/dev/null | head -n 1 || true)"
            fi
            if [[ -z "$source" || "$source" == "/" ]]; then
                if [[ -d "$candidate" ]]; then source="$candidate"; else source="$(dirname "$candidate")"; fi
            fi
        fi
        [[ "$source" != "/" && -d "$source" ]] || continue
        _add_auto_bind "$source" "$source"
    done
}

add_user_bind() {
    local specification="$1" source
    [[ -n "$specification" ]] || return 0
    source="${specification%%:*}"
    add_bind_for_path "$source"
    USER_BIND_SPECS+=("$specification")
}

finalize_user_binds() {
    local specification
    for specification in "${USER_BIND_SPECS[@]}"; do
        BIND_ARGS+=(--bind "$specification")
    done
}

add_symlink_targets_in_dir() {
    local directory="$1" link target
    [[ -d "$directory" ]] || return 0
    while IFS= read -r -d '' link; do
        target="$(readlink -f -- "$link" 2>/dev/null || true)"
        [[ -n "$target" ]] && add_bind_for_path "$target"
    done < <(find "$directory" -type l -print0 2>/dev/null || true)
}

collect_mapping_sample_paths() {
    local sample_table="$1" ref_dir="$2"
    python3 - "$sample_table" "$ref_dir" <<'PY'
import csv
import os
import sys

sample_table = os.path.abspath(sys.argv[1])
ref_dir = os.path.abspath(sys.argv[2])
base = os.path.dirname(sample_table)
with open(sample_table, "r", encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    required = {"sample_id", "sample_name", "species", "type", "fq1", "fq2"}
    if not reader.fieldnames or not required.issubset(reader.fieldnames):
        raise SystemExit("ERROR: samples.tsv requires sample_id, sample_name, species, type, fq1, fq2")
    seen = set()
    for line_no, row in enumerate(reader, 2):
        species = (row.get("species") or "").strip()
        if not species:
            raise SystemExit("ERROR: empty species at samples.tsv line %d" % line_no)
        for column in ("fq1", "fq2"):
            value = (row.get(column) or "").strip()
            if not value:
                raise SystemExit("ERROR: empty %s at samples.tsv line %d" % (column, line_no))
            path = value if os.path.isabs(value) else os.path.join(base, value)
            path = os.path.abspath(os.path.expanduser(path))
            if not os.path.isfile(path):
                raise SystemExit("ERROR: %s not found at samples.tsv line %d: %s" % (column, line_no, path))
            sys.stdout.buffer.write(os.fsencode(path) + b"\0")
        if species not in seen:
            seen.add(species)
            for suffix in (".fa", ".gtf"):
                path = os.path.join(ref_dir, species + suffix)
                if not os.path.isfile(path):
                    raise SystemExit("ERROR: reference file not found for %s: %s" % (species, path))
                sys.stdout.buffer.write(os.fsencode(path) + b"\0")
PY
}

collect_annotation_sample_paths() {
    local sample_table="$1" work_dir="$2"
    python3 - "$sample_table" "$work_dir" <<'PY'
import csv
import os
import sys

sample_table = os.path.abspath(sys.argv[1])
work_dir = os.path.abspath(sys.argv[2])
with open(sample_table, "r", encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    required = {"sample_id", "sample_name", "species", "type"}
    if not reader.fieldnames or not required.issubset(reader.fieldnames):
        raise SystemExit("ERROR: samples.tsv requires sample_id, sample_name, species, type")
    for row in reader:
        sid = (row.get("sample_id") or "").strip()
        species = (row.get("species") or "").strip()
        typ = (row.get("type") or "").strip().upper()
        if not sid or not species or typ not in ("DNA", "RNA"):
            raise SystemExit("ERROR: invalid sample row in samples.tsv")
        if typ == "DNA":
            bam = os.path.join(work_dir, "01_MAP", species, "DNA", sid, sid + ".nodup.bam")
        else:
            bam = os.path.join(work_dir, "01_MAP", species, "RNA", "final", sid, "split.bam")
        for path in (bam, bam + ".bai", os.path.join(work_dir, "00_REF", "GENOME", species + ".fa"),
                     os.path.join(work_dir, "00_REF", "GENOME", species + ".fa.fai")):
            sys.stdout.buffer.write(os.fsencode(path) + b"\0")
PY
}

run_in_container_logged() {
    local sif="$1" log_file="$2"
    shift 2
    finalize_user_binds
    echo "[denovo_detect] SIF: $sif"
    echo "[denovo_detect] Log: $log_file"
    echo "[denovo_detect] Automatic/user bind arguments: ${#BIND_ARGS[@]}"
    apptainer exec "${BIND_ARGS[@]}" "$sif" "$@" 2>&1 | tee -a "$log_file"
}
