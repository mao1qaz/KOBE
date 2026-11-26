#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser(description="Generate candidate CDS editing sites for each species.")
parser.add_argument("--result_dir", required=True, help="Directory containing merged_edit_sites_codon_table_uniq.csv")
parser.add_argument("--codon_table", default=os.environ.get("RNA_EDITING_CODON_TABLE", "/opt/rna_editing/share/Codon_table.new.txt"), help="Codon table TSV (default: container resource)")
args = parser.parse_args()

result_dir = Path(args.result_dir)
csv_file = result_dir / "merged_edit_sites_codon_table_uniq.csv"
out_dir = result_dir / "candi_per_species"
out_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(csv_file)

cds_cols = [c for c in df.columns if c.endswith("_cds_pos")]

# Codon table
if not Path(args.codon_table).is_file():
    raise FileNotFoundError(f"Codon table not found: {args.codon_table}")
codon_df = pd.read_csv(args.codon_table, sep="\t")
codon_map = dict(zip(codon_df["Codon"].str.upper(), codon_df["AA"]))

def translate(codon):
    return codon_map.get(codon.upper(), "-")

def mutate_codon(codon, pos):
    codon = list(codon)
    codon[pos - 1] = "G"
    return "".join(codon)

for col in cds_cols:
    sp = col.replace("_cds_pos", "")
    codon_col = f"{sp}_codon"
    codon_pos_col = f"{sp}_codon_pos"
    out_file = out_dir / f"{sp}_candi.txt"
    lines = []

    for _, row in df.iterrows():
        cds_val = str(row[col])

        if "NOT_FOUND" in cds_val or ":---" in cds_val:
            continue

        try:
            seq_id, pos = cds_val.split(":")
            pos = int(pos)
        except:
            continue

        codon = str(row[codon_col]).upper()
        pos_on_codon = int(row[codon_pos_col])
        base = codon[pos_on_codon - 1]

        if base != "A":
            continue

        aa1 = translate(codon)
        mut_codon = mutate_codon(codon, pos_on_codon)
        aa2 = translate(mut_codon)
        anno = "Syn" if aa1 == aa2 else "Nonsyn"

        ID2 = str(row["ID2"])
        aa_raw = str(row["AA_change"])
        aa_part = aa_raw.split(":")[-1] if ":" in aa_raw else "-"
        final_id = f"{ID2}:{aa_part}"

        line = [
            seq_id,
            pos,
            anno,
            codon,
            aa1,
            pos_on_codon,
            mut_codon,
            aa2,
            final_id
        ]

        lines.append("\t".join(map(str, line)))

    with open(out_file, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"{sp} done: {len(lines)}")

print("ALL DONE")