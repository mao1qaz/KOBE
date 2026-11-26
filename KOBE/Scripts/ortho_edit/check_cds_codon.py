#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
from Bio import SeqIO

# =========================
# Path settings
# =========================
parser = argparse.ArgumentParser(description="Check CDS codon consistency.")
parser.add_argument("--work_dir", required=True, help="Output working directory")
parser.add_argument("--cds_dir", required=True, help="Directory containing *_longestcds.fasta files")
parser.add_argument("--result_dir", required=True, help="Directory containing merged_edit_sites_codon_table_uniq.csv")
args = parser.parse_args()

csv_file = os.path.join(args.result_dir, "merged_edit_sites_codon_table_uniq.csv")
fasta_dir = args.cds_dir
out_dir = os.path.join(args.work_dir, "spe_cds_pos")
os.makedirs(out_dir, exist_ok=True)

# =========================
# Species to FASTA mapping
# =========================
species_fasta = {}

for filename in sorted(os.listdir(fasta_dir)):
    if filename.endswith("_longestcds.fasta"):
        sp = filename.replace("_longestcds.fasta", "").split("_")[-1].lower()
        species_fasta[sp] = filename

# =========================
# Read CDS FASTA
# =========================
def load_fasta(fasta_path):
    seq_dict = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq_dict[rec.id] = str(rec.seq).upper()
    return seq_dict


# =========================
# Extract codon
# =========================
def get_codon(seq, pos):
    pos = int(pos)

    # 1-based -> 0-based codon grouping
    codon_index = (pos - 1) // 3
    start = codon_index * 3

    codon = seq[start:start + 3]

    if len(codon) < 3:
        return None

    return codon

# =========================
# Read CSV
# =========================
df = pd.read_csv(csv_file)

all_bad = []

# =========================
# Process each species
# =========================
for sp in species_fasta:

    cds_col = f"{sp}_cds_pos"
    codon_col = f"{sp}_codon"

    fasta_path = os.path.join(fasta_dir, species_fasta[sp])
    seqs = load_fasta(fasta_path)

    bad_rows = []

    for idx, row in df.iterrows():

        val = str(row[cds_col])
        expected_codon = str(row[codon_col])

        # -------------------------
        # Filter NOT_FOUND
        # -------------------------
        if "NOT_FOUND" in val:
            continue

        # -------------------------
        # Filter missing positions
        # -------------------------
        if ":---" in val:
            continue

        # -------------------------
        # Parse transcript and position
        # -------------------------
        try:
            seq_id, pos = val.split(":")
        except:
            continue

        if seq_id not in seqs:
            continue

        seq = seqs[seq_id]
        codon = get_codon(seq, pos)

        if codon is None:
            continue

        # -------------------------
        # Case-insensitive comparison
        # -------------------------
        codon = codon.upper()
        expected_codon = expected_codon.upper()

        # -------------------------
        # Compare codons
        # -------------------------
        if codon != expected_codon:
            bad_rows.append(val)

    # Output mismatched sites
    out_file = os.path.join(out_dir, f"{sp}_codon_mismatch.txt")
    pd.Series(bad_rows).to_csv(out_file, index=False, header=False)

    print(f"[{sp}] mismatch: {len(bad_rows)}")

print("Done.")