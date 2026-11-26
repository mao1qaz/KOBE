
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Annotate amino acid changes for each species.")
parser.add_argument("--result_dir", required=True, help="Directory containing merged_edit_sites_codon_table_uniq.csv")
parser.add_argument("--codon_table", default=os.environ.get("RNA_EDITING_CODON_TABLE", "/opt/rna_editing/share/Codon_table.new.txt"), help="Codon table TSV (default: container resource)")
args = parser.parse_args()
os.makedirs(args.result_dir, exist_ok=True)

csv_file = os.path.join(args.result_dir, "merged_edit_sites_codon_table_uniq.csv")
out_file = os.path.join(args.result_dir, "merged_edit_sites_codon_table_annotated.csv")

df = pd.read_csv(csv_file)

cds_cols = [c for c in df.columns if c.endswith("_cds_pos")]

# =========================
# Codon table
# =========================
codon_table_file = args.codon_table
if not os.path.isfile(codon_table_file):
    raise FileNotFoundError(f"Codon table not found: {codon_table_file}")
codon_df = pd.read_csv(codon_table_file, sep="\t")
codon_map = dict(zip(codon_df["Codon"].str.upper(), codon_df["AA"]))


def translate(codon):
    if codon in codon_map:
        return codon_map[codon]
    return "-"


def mutate_codon(codon, pos):
    codon = list(codon)
    codon[pos - 1] = "G"
    return "".join(codon)


# =========================
# Process each species
# =========================
for col in cds_cols:

    sp = col.replace("_cds_pos", "")
    codon_col = f"{sp}_codon"
    codon_pos_col = f"{sp}_codon_pos"

    aa_col_name = f"{sp}_AA_change"
    aa_changes = []

    for i, row in df.iterrows():

        val = str(row[col])

        # 1. Basic filtering
        if "NOT_FOUND" in val or ":---" in val:
            aa_changes.append("-")
            continue

        try:
            codon = str(row[codon_col]).upper()
            pos_on_codon = int(row[codon_pos_col])
        except:
            aa_changes.append("-")
            continue

        if codon in ["---", "XXX", "NAN"]:
            aa_changes.append("-")
            continue

        if pos_on_codon < 1 or pos_on_codon > 3:
            aa_changes.append("-")
            continue

        # 2. Check reference A
        base = codon[pos_on_codon - 1]

        if base != "A":
            aa_changes.append("-")
            continue

        # 3. AA1
        aa1 = translate(codon)

        # 4. Mutate codon (A->G)
        mut_codon = mutate_codon(codon, pos_on_codon)

        # 5. AA2
        aa2 = translate(mut_codon)

        if aa1 == "-" or aa2 == "-":
            aa_changes.append("-")
        else:
            aa_changes.append(f"{aa1}>{aa2}")

    df[aa_col_name] = aa_changes

# =========================
# Output
# =========================
df.to_csv(out_file, index=False)

print("Done:", out_file)
