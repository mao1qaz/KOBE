#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Filter supplementary conserved sites by selected group definition and minimum species count."
    )
    p.add_argument("--out_dir", required=True, help="Supplementary output directory")
    p.add_argument("--min_species", required=True, type=int, help="Minimum number of edited species")
    p.add_argument(
        "--group_col",
        required=True,
        choices=["ID2_group", "new_id2_group"],
        help="Use the original or updated homologous editing group"
    )
    return p.parse_args()


args = parse_args()
if args.min_species < 1:
    raise SystemExit("--min_species must be >= 1")

out_dir = Path(args.out_dir)
in_file = out_dir / "sup_conserved_site_all_species_9nt_3aa.tsv"
if not in_file.is_file():
    raise FileNotFoundError(f"Input site table not found: {in_file}")

count_col = "n_ID2_group_species" if args.group_col == "ID2_group" else "n_new_id2_group_species"

df = pd.read_csv(in_file, sep="\t", dtype=str, keep_default_na=False)
required = [
    "ID2", "source_species", "transcript", "transcript_pos",
    "edited_species", "recovered_species", count_col
]
missing = [x for x in required if x not in df.columns]
if missing:
    raise ValueError("Input table is missing columns: " + ", ".join(missing))

df[count_col] = pd.to_numeric(df[count_col], errors="coerce")
invalid = int(df[count_col].isna().sum())
if invalid:
    print(f"WARNING: {invalid} rows have invalid values in {count_col}")

selected = df.loc[df[count_col] >= args.min_species].copy()
selected["n_edited_species"] = selected[count_col]
selected["filter_group_col"] = args.group_col
selected["_sort_pos"] = pd.to_numeric(selected["transcript_pos"], errors="coerce")
selected = selected.sort_values(
    by=["n_edited_species", "transcript", "_sort_pos", "ID2"],
    ascending=[False, True, True, True],
    kind="stable"
).drop(columns="_sort_pos")

tag = f"{args.group_col}_ge{args.min_species}species"
site_out = out_dir / f"sup_conserved_sites_{tag}.tsv"
transcript_out = out_dir / f"sup_conserved_sites_{tag}_transcripts.tsv"
selected.to_csv(site_out, sep="\t", index=False)

invalid_transcripts = {"", "NA", "NAN", "NONE", "-"}
valid = selected.loc[
    ~selected["transcript"].str.strip().str.upper().isin(invalid_transcripts)
].copy()

summary_cols = [
    "transcript", "source_species", "n_sites", "max_n_edited_species",
    "mean_n_edited_species", "site_IDs", "transcript_positions",
    "edited_species_by_site", "recovered_species_by_site"
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
        edited_species_by_site=("edited_species", lambda x: " | ".join(x)),
        recovered_species_by_site=("recovered_species", lambda x: " | ".join(x))
    )
    summary["mean_n_edited_species"] = summary["mean_n_edited_species"].round(3)
    summary = summary.sort_values(
        by=["n_sites", "max_n_edited_species", "transcript"],
        ascending=[False, False, True],
        kind="stable"
    )

summary.to_csv(transcript_out, sep="\t", index=False)

print("\n========== SUPPLEMENTARY FILTER SUMMARY ==========")
print(f"Group definition: {args.group_col}")
print(f"Minimum edited species: {args.min_species}")
print(f"Count column: {count_col}")
print(f"Input sites: {len(df)}")
print(f"Selected sites: {len(selected)}")
print(f"Selected transcripts: {len(summary)}")
print(f"Selected sites without representative transcript: {len(selected) - len(valid)}")
print("\nGroup-size distribution:")
print(df[count_col].value_counts(dropna=False).sort_index().to_string())
print("\nOutputs:")
print(site_out)
print(transcript_out)
print("Done.")
