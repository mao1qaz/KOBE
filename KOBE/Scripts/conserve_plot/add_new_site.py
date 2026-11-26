
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import re
import sys
from pathlib import Path


def clean(value):
    if value is None:
        return ""
    value = str(value).strip()
    if value.lower() in {"", "na", "nan", "none", "null"}:
        return ""
    return value


def normalize_position(value):
    value = clean(value)
    if re.fullmatch(r"\d+\.0+", value):
        value = value.split(".")[0]
    return value


def make_site_key(transcript, cds_pos):
    transcript = clean(transcript)
    cds_pos = normalize_position(cds_pos)

    if not transcript or not cds_pos:
        return ""

    return f"{transcript}:{cds_pos}"


def parse_group(value):
    value = clean(value)

    if not value:
        return []

    items = [x.strip() for x in value.split(";") if x.strip()]
    return list(dict.fromkeys(items))


def load_recovered_sites(site_files):
    site_keys = set()

    for file in site_files:
        path = Path(file)

        if not path.is_file():
            raise FileNotFoundError(f"Site file not found: {path}")

        count = 0

        with open(path, "r", encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")

            if not reader.fieldnames:
                raise ValueError(f"Empty site file: {path}")

            required = {"transcript", "cds_pos"}
            missing = required - set(reader.fieldnames)

            if missing:
                raise ValueError(
                    f"Missing columns in {path}: {sorted(missing)}"
                )

            for row in reader:
                key = make_site_key(
                    row["transcript"],
                    row["cds_pos"]
                )

                if key:
                    site_keys.add(key)
                    count += 1

        print(f"[INFO] Loaded {count} records: {path}")

    return site_keys


def main():

    parser = argparse.ArgumentParser(
        description="Add newly recovered editing sites to ID2_group."
    )

    parser.add_argument(
        "--annotated",
        required=True,
        help="Original merged_edit_sites_codon_table_annotated.csv"
    )

    parser.add_argument(
        "--sites",
        required=True,
        action="append",
        help="Novel supported sites TSV. Can be specified multiple times."
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory"
    )

    args = parser.parse_args()

    annotated_file = Path(args.annotated)
    outdir = Path(args.outdir)

    if not annotated_file.is_file():
        raise FileNotFoundError(
            f"Annotated table not found: {annotated_file}"
        )

    outdir.mkdir(parents=True, exist_ok=True)

    output_file = (
        outdir /
        "merged_edit_sites_codon_table_annotated_with_new_sites.csv"
    )

    print("[INFO] Loading recovered sites...")

    site_keys = load_recovered_sites(args.sites)

    print(f"[INFO] Unique recovered sites: {len(site_keys)}")

    print("[INFO] Loading original annotation table...")

    with open(
        annotated_file, "r", encoding="utf-8-sig", newline=""
    ) as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames

        if not fieldnames:
            raise ValueError("Annotated table has no header.")

        if "ID2_group" not in fieldnames:
            raise ValueError("ID2_group column not found.")

        cds_columns = [
            col for col in fieldnames
            if col.endswith("_cds_pos")
        ]

        if not cds_columns:
            raise ValueError("No species *_cds_pos columns found.")

        rows = list(reader)

    print(f"[INFO] Original rows: {len(rows)}")
    print(f"[INFO] Species CDS position columns: {len(cds_columns)}")

    matched_site_keys = set()
    updated_rows = 0
    added_ids = 0

    for row in rows:

        original_group = parse_group(row.get("ID2_group", ""))
        original_set = set(original_group)

        matched_keys = []

        for col in cds_columns:

            key = clean(row.get(col, ""))

            if key in site_keys:
                matched_keys.append(key)
                matched_site_keys.add(key)

        matched_keys = list(dict.fromkeys(matched_keys))

        # Only add IDs that are not already in the original group.
        new_keys = [
            key for key in matched_keys
            if key not in original_set
        ]

        new_group = original_group + new_keys

        row["new_id2_group"] = ";".join(new_group)
        row["recovered_only_group"] = ";".join(new_keys)

        if new_keys:
            updated_rows += 1
            added_ids += len(new_keys)

    unmatched = sorted(site_keys - matched_site_keys)

    print("[INFO] Writing updated annotation table...")

    # Preserve all original columns and append the two new columns.
    output_fields = [
        col for col in fieldnames
        if col not in {"new_id2_group", "recovered_only_group"}
    ]

    output_fields.extend([
        "new_id2_group",
        "recovered_only_group"
    ])

    with open(
        output_file, "w", encoding="utf-8-sig", newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=output_fields,
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)

    print("\n========== SUMMARY ==========")
    print(f"Original rows: {len(rows)}")
    print(f"Unique recovered sites: {len(site_keys)}")
    print(f"Matched recovered sites: {len(matched_site_keys)}")
    print(f"Unmatched recovered sites: {len(unmatched)}")
    print(f"Rows with newly added IDs: {updated_rows}")
    print(f"Total IDs added across rows: {added_ids}")
    print(f"Output: {output_file}")

    if unmatched:
        print("\n[WARNING] Unmatched site examples:")
        for key in unmatched[:20]:
            print(f"  {key}")

    print("\n[INFO] Finished successfully.")


if __name__ == "__main__":
    main()
