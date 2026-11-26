#!/usr/bin/env bash
set -Eeuo pipefail
SELF="$(readlink -f -- "${BASH_SOURCE[0]}")"
PKG_ROOT="$(cd -- "$(dirname -- "$SELF")" && pwd -P)"
BIN_DIR="${RNA_EDITING_BIN_DIR:-$HOME/.local/bin}"
mkdir -p -- "$BIN_DIR"

COMMANDS=(
  KOBE
  snpeff_build denovo_detect_map denovo_detect_anno base_extract
  pairwise_tests ortho_edit conserve_plot find_back add_new_site find_back_plot
  add_cds_sites_genomic_coordinates extract_cds_candidate_sites
  extract_longest_cds genome_site_cds_position prepare_snpeff_gtf
  supported_sites_to_editing_counts
)

for command_name in "${COMMANDS[@]}"; do
  ln -sfn "$PKG_ROOT/bin/KOBE" "$BIN_DIR/$command_name"
done

# Remove stale aliases from earlier KOBE command sets (for example, renamed
# commands or retired utilities) when they still point to this dispatcher.
declare -A ACTIVE_COMMANDS=()
for command_name in "${COMMANDS[@]}"; do ACTIVE_COMMANDS["$command_name"]=1; done
for link_path in "$BIN_DIR"/*; do
  [[ -L "$link_path" ]] || continue
  link_target="$(readlink -f -- "$link_path" 2>/dev/null || true)"
  link_name="$(basename -- "$link_path")"
  if [[ "$link_target" == "$PKG_ROOT/bin/KOBE" && -z "${ACTIVE_COMMANDS[$link_name]+x}" ]]; then
    rm -f -- "$link_path"
  fi
done

# Remove the previous package-level command link when it points to the old
# RNA_editing_toolkit dispatcher. Module command aliases remain unchanged.
if [[ -L "$BIN_DIR/rna-edit" ]]; then
  old_target="$(readlink -f -- "$BIN_DIR/rna-edit" 2>/dev/null || true)"
  if [[ "$old_target" == */RNA_editing_toolkit/bin/rna-edit ]]; then
    rm -f -- "$BIN_DIR/rna-edit"
  fi
fi

chmod +x "$PKG_ROOT/bin/KOBE"
for f in "$PKG_ROOT"/Scripts/*/*; do
  [[ -f "$f" ]] || continue
  case "$f" in
    */*.sh|*/base_extract|*/find_back|*/denovo_detect_map|*/denovo_detect_anno|*/ortho_edit|*/conserve_plot|*/find_back_plot|*/add_new_site|*/pairwise_tests|*/snpeff_build|*/add_cds_sites_genomic_coordinates|*/extract_cds_candidate_sites|*/extract_longest_cds|*/genome_site_cds_position|*/prepare_snpeff_gtf|*/supported_sites_to_editing_counts)
      chmod +x "$f" ;;
  esac
done

echo "KOBE commands installed in: $BIN_DIR"
if [[ ":$PATH:" != *":$BIN_DIR:"* ]]; then
  echo
  echo "Add this single line to your shell startup file (e.g. ~/.bashrc):"
  printf 'export PATH="%s:$PATH"\n' "$BIN_DIR"
  echo
  echo "Then run: source ~/.bashrc"
else
  echo "The installation directory is already in PATH."
fi
