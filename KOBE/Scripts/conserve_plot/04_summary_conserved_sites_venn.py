
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import itertools
from collections import Counter
from pathlib import Path

import pandas as pd


# ============================================================
# 1. Arguments and paths
# ============================================================

parser = argparse.ArgumentParser(
    description="Summarize conserved RNA editing sites and exact species-combination regions."
)

parser.add_argument("--annotated", required=True)
parser.add_argument("--cds_dir", required=True)
parser.add_argument("--out_dir", required=True)
parser.add_argument("--min_species", required=True, type=int)
parser.add_argument("--species_order", required=True)

args = parser.parse_args()

if args.min_species < 2:
    parser.error("--min_species must be >= 2")

out_dir = Path(args.out_dir)
fasta_dir = Path(args.cds_dir)
anno_file = Path(args.annotated)

site_file = out_dir / f"conserved_sites_ge{args.min_species}species.tsv"

species_order = [
    x.strip() for x in args.species_order.split(",") if x.strip()
]

if len(species_order) < 2:
    parser.error("--species_order must contain at least two species")

if len(species_order) != len(set(species_order)):
    parser.error("Duplicate species in --species_order")

out_dir.mkdir(parents=True, exist_ok=True)


# ============================================================
# 3. Helper functions
# ============================================================

def parse_species_list(value):
    return [
        x.strip()
        for x in str(value).split(";")
        if x.strip() and
        x.strip().upper() not in {"NA", "NAN", "-"}
    ]


def format_aa_change(value):
    s = str(value).strip()

    if s.upper() in {"", "NA", "NAN", "-", "NONE"}:
        return None

    if s.count(">") != 1:
        return None

    aa1, aa2 = [x.strip() for x in s.split(">", 1)]

    if not aa1 or not aa2:
        return None

    return f"{aa1}>{aa2}"


def classify_site(row):
    """
    Classification is based exclusively on species
    in which editing was actually detected.
    """

    edited_species = list(dict.fromkeys(
        parse_species_list(row.get("edited_species", ""))
    ))

    if not edited_species:
        return "NA"

    changes = []

    for sp in edited_species:

        change = format_aa_change(
            row.get(f"{sp}_AA_change", "")
        )

        if change is None:
            return "NA"

        changes.append(change)

    unique_changes = list(dict.fromkeys(changes))

    if len(unique_changes) > 1:
        return "Mixed"

    aa1, aa2 = unique_changes[0].split(">", 1)

    if aa1 == aa2:
        return "Syn"

    return "Nonsyn"


# ============================================================
# 4. Read filtered conserved sites
# ============================================================

if not site_file.is_file():
    raise FileNotFoundError(f"Site file not found: {site_file}")

if not anno_file.is_file():
    raise FileNotFoundError(f"Annotation file not found: {anno_file}")

df = pd.read_csv(
    site_file,
    sep="\t",
    dtype=str,
    keep_default_na=False
)

required = ["ID2", "edited_species", "n_edited_species"]

missing = [x for x in required if x not in df.columns]

if missing:
    raise ValueError("Missing columns: " + ", ".join(missing))

df["ID2"] = df["ID2"].str.strip()

df = df[df["ID2"] != ""].copy()

if df["ID2"].duplicated().any():
    raise ValueError(
        "Duplicate ID2 detected in filtered site file. "
        "Expected one row per homologous site."
    )

df["n_edited_species"] = pd.to_numeric(
    df["n_edited_species"],
    errors="coerce"
)

if (df["n_edited_species"] < args.min_species).any():
    raise ValueError(
        "Input contains sites below the requested minimum species number."
    )

print("Filtered sites:", len(df))
print("Species order:", ", ".join(species_order))


# ============================================================
# 5. Merge amino acid changes from annotated CSV
# ============================================================

try:
    anno = pd.read_csv(
        anno_file,
        dtype=str,
        keep_default_na=False,
        encoding="utf-8-sig"
    )
except UnicodeDecodeError:
    anno = pd.read_csv(
        anno_file,
        dtype=str,
        keep_default_na=False,
        encoding="gbk"
    )

if "ID2" not in anno.columns:
    raise ValueError("Annotated CSV missing ID2")

anno["ID2"] = anno["ID2"].str.strip()

aa_cols = [
    f"{sp}_AA_change"
    for sp in species_order
]

missing_aa = [
    col for col in aa_cols
    if col not in anno.columns
]

if missing_aa:
    raise ValueError(
        "Missing AA change columns in annotated CSV: " +
        ", ".join(missing_aa)
    )

missing_in_site = [
    col for col in aa_cols
    if col not in df.columns
]

if missing_in_site:

    anno_aa = anno[["ID2"] + missing_in_site].drop_duplicates(
        "ID2",
        keep="first"
    )

    df = df.merge(
        anno_aa,
        on="ID2",
        how="left",
        validate="one_to_one"
    )


# ============================================================
# 6. Validate edited species
# ============================================================

unknown_species = set()

for value in df["edited_species"]:
    unknown_species.update(
        set(parse_species_list(value)) - set(species_order)
    )

if unknown_species:
    raise ValueError(
        "Species in edited_species missing from --species_order: " +
        ", ".join(sorted(unknown_species))
    )

df["_edited_set"] = df["edited_species"].apply(
    lambda x: set(parse_species_list(x))
)

df["_actual_n"] = df["_edited_set"].apply(len)

mismatch = df[
    df["_actual_n"] != df["n_edited_species"]
]

if not mismatch.empty:
    raise ValueError(
        f"{len(mismatch)} sites have inconsistent n_edited_species "
        "and edited_species. Check ID2_group mapping."
    )


# ============================================================
# 7. Classify all conserved sites
# ============================================================

df["effect_class"] = df.apply(
    classify_site,
    axis=1
)

total_sites = len(df)

class_order = [
    "Nonsyn",
    "Syn",
    "Mixed",
    "NA"
]

class_counts = df["effect_class"].value_counts()

summary_rows = [
    {
        "category": "Total",
        "n_sites": total_sites,
        "percentage": 100.0 if total_sites else 0.0
    }
]

for category in class_order:

    count = int(class_counts.get(category, 0))

    percentage = (
        count / total_sites * 100
        if total_sites else 0
    )

    summary_rows.append({
        "category": category,
        "n_sites": count,
        "percentage": round(percentage, 3)
    })

summary = pd.DataFrame(summary_rows)

summary_file = out_dir / (
    f"conserved_sites_ge{args.min_species}species_summary.tsv"
)

summary.to_csv(
    summary_file,
    sep="\t",
    index=False
)


# ============================================================
# 8. Build species-specific editing sets
# ============================================================

site_ids = df["ID2"].tolist()

species_sets = {}

for sp in species_order:

    species_sets[sp] = {
        row["ID2"]
        for _, row in df.iterrows()
        if sp in row["_edited_set"]
    }

species_count_rows = []

for sp in species_order:

    species_count_rows.append({
        "species": sp,
        "n_conserved_sites": len(species_sets[sp])
    })

species_counts = pd.DataFrame(species_count_rows)

species_count_file = out_dir / (
    f"conserved_sites_ge{args.min_species}species_species_counts.tsv"
)

species_counts.to_csv(
    species_count_file,
    sep="\t",
    index=False
)


# ============================================================
# 9. Calculate exact Venn regions
# ============================================================

region_counter = Counter()

for _, row in df.iterrows():

    edited = row["_edited_set"]

    pattern = "".join(
        "1" if sp in edited else "0"
        for sp in species_order
    )

    region_counter[pattern] += 1


region_rows = []

n_species = len(species_order)

for bits in itertools.product("01", repeat=n_species):

    pattern = "".join(bits)

    if "1" not in pattern:
        continue

    combination = [
        sp for i, sp in enumerate(species_order)
        if pattern[i] == "1"
    ]

    region_rows.append({
        "pattern": pattern,
        "species_combination": ";".join(combination),
        "n_sites": region_counter.get(pattern, 0)
    })

region_df = pd.DataFrame(region_rows)

region_file = out_dir / (
    f"conserved_sites_ge{args.min_species}species_venn_regions.tsv"
)

region_df.to_csv(
    region_file,
    sep="\t",
    index=False
)


# ============================================================
# 10. Remove obsolete Venn graphics from previous runs
# ============================================================

# Exact species-combination statistics are retained in the TSV file.
# Venn diagram graphics are no longer generated.
for obsolete_file in (
    out_dir / f"conserved_sites_ge{args.min_species}species_venn.pdf",
    out_dir / f"conserved_sites_ge{args.min_species}species_venn.png",
):
    if obsolete_file.exists():
        obsolete_file.unlink()


# ============================================================
# 11. Print summary
# ============================================================

print("\n==============================================")
print("Conserved RNA editing site summary")
print("==============================================")

print(f"Minimum edited species: {args.min_species}")
print(f"Total conserved sites: {total_sites}")

for category in class_order:

    print(
        f"{category}: {class_counts.get(category, 0)}"
    )

print("\nSpecies-specific conserved editing sites:")

print(species_counts.to_string(index=False))

print("\nExact species-combination regions:")

print(region_df.to_string(index=False))

print("\nOutputs:")
print(summary_file)
print(species_count_file)
print(region_file)

print("\nVenn graphics: not generated")

print("\nALL DONE")
