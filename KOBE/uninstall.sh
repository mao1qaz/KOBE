#!/usr/bin/env bash
set -Eeuo pipefail
SELF="$(readlink -f -- "${BASH_SOURCE[0]}")"
PKG_ROOT="$(cd -- "$(dirname -- "$SELF")" && pwd -P)"
BIN_DIR="${RNA_EDITING_BIN_DIR:-$HOME/.local/bin}"

# Remove every alias that resolves to this package's dispatcher, including
# aliases from older command sets that have since been renamed or retired.
for link_path in "$BIN_DIR"/*; do
  [[ -L "$link_path" ]] || continue
  link_target="$(readlink -f -- "$link_path" 2>/dev/null || true)"
  if [[ "$link_target" == "$PKG_ROOT/bin/KOBE" ]]; then
    rm -f -- "$link_path"
  fi
done

echo "KOBE command links removed from $BIN_DIR."
echo "The software directory itself was not deleted."
