#!/usr/bin/env bash

run_in_container_logged() {
    local sif="$1"
    local log_file="$2"
    shift 2

    local -a bind_paths=()
    while [[ $# -gt 0 && "$1" != "--" ]]; do
        bind_paths+=("$1")
        shift
    done
    if [[ $# -eq 0 ]]; then
        echo "Internal error: missing -- before container command" >&2
        return 2
    fi
    shift

    local -a bind_args=()
    local -A seen=()

    add_bind_dir() {
        local dir="$1"
        [[ -n "$dir" ]] || return 0
        dir="$(realpath -ms -- "$dir")"
        [[ -d "$dir" ]] || mkdir -p "$dir"
        if [[ -z "${seen[$dir]+x}" ]]; then
            seen["$dir"]=1
            bind_args+=(--bind "$dir:$dir")
        fi
    }

    local path original_path resolved_path bind_dir
    for path in "${bind_paths[@]}"; do
        [[ -n "$path" ]] || continue

        if [[ -d "$path" ]]; then
            # Bind the user-supplied directory path.
            add_bind_dir "$path"

            # If the directory itself is a symlink, also bind its real target.
            resolved_path="$(readlink -f -- "$path" 2>/dev/null || true)"
            if [[ -n "$resolved_path" && -d "$resolved_path" ]]; then
                add_bind_dir "$resolved_path"
            fi
        else
            # Bind the parent of the user-supplied file/path.
            original_path="$(realpath -ms -- "$path")"
            add_bind_dir "$(dirname "$original_path")"

            # Also bind the parent of a symlink's final target.
            if [[ -e "$path" || -L "$path" ]]; then
                resolved_path="$(readlink -f -- "$path" 2>/dev/null || true)"
                if [[ -n "$resolved_path" ]]; then
                    if [[ -d "$resolved_path" ]]; then
                        bind_dir="$resolved_path"
                    else
                        bind_dir="$(dirname "$resolved_path")"
                    fi
                    add_bind_dir "$bind_dir"
                fi
            fi
        fi
    done

    echo "Container bind mounts:"
    for ((i=0; i<${#bind_args[@]}; i+=2)); do
        echo "  ${bind_args[i+1]}"
    done

    apptainer exec "${bind_args[@]}" "$sif" "$@" 2>&1 | tee "$log_file"
}
