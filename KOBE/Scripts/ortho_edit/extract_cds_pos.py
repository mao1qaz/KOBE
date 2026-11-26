#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser(description="Extract CDS positions for each species.")
parser.add_argument("--work_dir", required=True, help="Output working directory")
parser.add_argument("--result_dir", required=True, help="Directory containing merged_edit_sites_codon_table_uniq.csv")
args = parser.parse_args()

csv_file = os.path.join(args.result_dir, "merged_edit_sites_codon_table_uniq.csv")
out_dir = os.path.join(args.work_dir, "spe_cds_pos")
os.makedirs(out_dir, exist_ok=True)

df = pd.read_csv(csv_file)

cds_cols = [c for c in df.columns if c.endswith("_cds_pos")]

# =========================
# 1. Merged table
# =========================
new_df = df[["ID2"] + cds_cols].copy()
new_df.to_csv(os.path.join(out_dir, "merged_ID2_cds_pos.tsv"), sep="\t", index=False)

print("New table saved.")

# =========================
# 2. Per-species lists
# =========================
for col in cds_cols:

    sp = col.replace("_cds_pos", "")
    codon_pos_col = f"{sp}_codon_pos"
    codon_col = f"{sp}_codon"

    out_file = os.path.join(out_dir, f"{sp}_cds_pos.list")
    lines = []

    for _, row in df.iterrows():

        val = str(row[col])

        if "NOT_FOUND" in val:
            continue

        if ":---" in val:
            continue

        try:
            seq_id, pos = val.split(":")
        except:
            continue

        codon = str(row[codon_col])
        codon_pos = str(row[codon_pos_col])

        if codon in ["---", "nan", "NaN"]:
            continue

        try:
            codon_pos = int(codon_pos)
        except:
            continue

        codon = codon.upper()

        if codon_pos < 1 or codon_pos > 3:
            continue

        base = codon[codon_pos - 1]

        if base != "A":
            continue

        lines.append(f"{seq_id}\t{pos}")

    with open(out_file, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"{sp} list saved: {len(lines)}")