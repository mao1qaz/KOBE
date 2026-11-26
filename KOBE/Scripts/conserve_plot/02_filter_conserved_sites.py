#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

import pandas as pd

# ============================================================
# 1. Input / output
# ============================================================
parser = argparse.ArgumentParser(description="Select homologous RNA editing sites detected in at least N species.")
parser.add_argument("--out_dir", required=True, help="Conservation output directory")
parser.add_argument("--min_species", "-n", required=True, type=int, help="Minimum number of species with detected RNA editing")
args = parser.parse_args()

if args.min_species < 1:
    parser.error("--min_species must be at least 1")

out_dir = Path(args.out_dir)
in_file = out_dir / "conserved_site_all_species_9nt_3aa.tsv"
site_out = out_dir / f"conserved_sites_ge{args.min_species}species.tsv"
transcript_out = out_dir / f"conserved_sites_ge{args.min_species}species_transcripts.tsv"

if not in_file.is_file():
    raise FileNotFoundError(f"Input file not found: {in_file}")

# ============================================================
# 2. Read table
# ============================================================
df = pd.read_csv(in_file, sep="\t", dtype=str, keep_default_na=False)
required = ["ID2", "edited_species", "n_edited_species", "source_species", "transcript", "transcript_pos"]
missing = [col for col in required if col not in df.columns]
if missing:
    raise ValueError("Missing columns: " + ", ".join(missing))

# n_edited_species was derived from ID2_group in the preceding script.
df["n_edited_species"] = pd.to_numeric(df["n_edited_species"], errors="coerce")
invalid_count = int(df["n_edited_species"].isna().sum())
if invalid_count:
    print(f"WARNING: {invalid_count} rows have invalid n_edited_species and cannot pass the filter")

# ============================================================
# 3. Filter sites edited in at least N species
# ============================================================
selected = df.loc[df["n_edited_species"] >= args.min_species].copy()
selected["_sort_pos"] = pd.to_numeric(selected["transcript_pos"], errors="coerce")
selected = selected.sort_values(
    by=["n_edited_species", "transcript", "_sort_pos", "ID2"],
    ascending=[False, True, True, True], kind="stable"
).drop(columns="_sort_pos")

# Preserve every original site-level column, including sequences and states.
selected.to_csv(site_out, sep="\t", index=False)

# ============================================================
# 4. Transcript summary of selected sites
# ============================================================
invalid_transcripts = {"", "NA", "NAN", "NONE", "-"}
valid = selected.loc[~selected["transcript"].str.strip().str.upper().isin(invalid_transcripts)].copy()

summary_cols = [
    "transcript", "source_species", "n_sites", "max_n_edited_species",
    "mean_n_edited_species", "site_IDs", "transcript_positions", "edited_species_by_site"
]

if valid.empty:
    summary = pd.DataFrame(columns=summary_cols)
else:
    summary = valid.groupby("transcript", sort=False, as_index=False).agg(
        source_species=("source_species", lambda x: ";".join(dict.fromkeys(x))),
        n_sites=("ID2", "size"),
        max_n_edited_species=("n_edited_species", "max"),
        mean_n_edited_species=("n_edited_species", "mean"),
        site_IDs=("ID2", lambda x: ";".join(x)),
        transcript_positions=("transcript_pos", lambda x: ";".join(x)),
        edited_species_by_site=("edited_species", lambda x: " | ".join(x))
    )
    summary["mean_n_edited_species"] = summary["mean_n_edited_species"].round(3)
    summary = summary.sort_values(
        by=["n_sites", "max_n_edited_species", "transcript"],
        ascending=[False, False, True], kind="stable"
    )

summary.to_csv(transcript_out, sep="\t", index=False)

# ============================================================
# 5. Report
# ============================================================
print("\n========================================")
print(f"RNA editing sites shared by >= {args.min_species} species")
print("========================================")
print(f"Input sites: {len(df)}")
print(f"Selected sites: {len(selected)}")
print(f"Selected transcripts: {len(summary)}")
print(f"Selected sites without representative transcript: {len(selected) - len(valid)}")

if df["n_edited_species"].notna().any():
    print("\nDistribution of edited-species counts:")
    print(df["n_edited_species"].value_counts().sort_index().to_string())

if not summary.empty:
    print("\nTop transcripts by number of selected sites:")
    print(summary[["transcript", "source_species", "n_sites", "max_n_edited_species", "mean_n_edited_species"]].head(20).to_string(index=False))

print("\nOutputs:")
print(site_out)
print(transcript_out)
print("Done.")
