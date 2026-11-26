
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import glob
import os
import re
import sys
from pathlib import Path


BASES = ("A", "C", "G", "T")

STRAND_LABELS = {
    "+",
    "-",
    "unknown_forward_reference"
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Merge REDItoolKnown A/C/G/T counts from all samples "
            "of one species into a single wide TSV table."
        )
    )

    parser.add_argument("--sample-table", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--output", required=True)

    return parser.parse_args()


def read_samples(sample_table, species):
    samples = []
    seen = set()

    with open(
        sample_table,
        encoding="utf-8-sig",
        newline=""
    ) as handle:

        reader = csv.DictReader(handle, delimiter="\t")

        if not reader.fieldnames:
            raise ValueError(
                f"Sample table has no header: {sample_table}"
            )

        required = {"sample_id", "species"}
        missing = required - set(reader.fieldnames)

        if missing:
            raise ValueError(
                f"Sample table missing columns: {sorted(missing)}"
            )

        for line_no, row in enumerate(reader, start=2):

            sample_id = (row.get("sample_id") or "").strip()
            row_species = (row.get("species") or "").strip()

            if not sample_id:
                continue

            if row_species != species:
                continue

            if not re.fullmatch(
                r"[A-Za-z0-9][A-Za-z0-9_.-]*",
                sample_id
            ):
                raise ValueError(
                    f"Invalid sample_id at line {line_no}: {sample_id}"
                )

            if sample_id in seen:
                raise ValueError(
                    f"Duplicate sample_id in species {species}: {sample_id}"
                )

            if "type" in row and (row.get("type") or "").strip():
                if row["type"].strip().upper() != "RNA":
                    raise ValueError(
                        f"Non-RNA sample found: {sample_id}"
                    )

            seen.add(sample_id)
            samples.append(sample_id)

    if not samples:
        raise ValueError(
            f"No samples found for species: {species}"
        )

    return samples


def read_site_metadata(metadata_path):
    sites = {}
    site_order = []

    with open(
        metadata_path,
        encoding="utf-8-sig",
        newline=""
    ) as handle:

        reader = csv.DictReader(handle, delimiter="\t")

        if not reader.fieldnames:
            raise ValueError(
                f"Site metadata has no header: {metadata_path}"
            )

        required = {
            "_chromosome",
            "_position",
            "_strand_label"
        }

        missing = required - set(reader.fieldnames)

        if missing:
            raise ValueError(
                f"Site metadata missing columns: {sorted(missing)}"
            )

        for line_no, row in enumerate(reader, start=2):

            chrom = (row.get("_chromosome") or "").strip()
            pos_text = (row.get("_position") or "").strip()
            strand = (row.get("_strand_label") or "").strip()

            if not chrom or not pos_text.isdigit():
                raise ValueError(
                    f"Invalid chromosome or position at metadata line {line_no}."
                )

            if strand not in STRAND_LABELS:
                raise ValueError(
                    f"Invalid strand label '{strand}' "
                    f"at metadata line {line_no}."
                )

            key = (chrom, int(pos_text))

            if key in sites:
                if sites[key] != strand:
                    raise ValueError(
                        f"Conflicting strand assignments for "
                        f"{chrom}:{pos_text}."
                    )
                continue

            sites[key] = strand
            site_order.append(key)

    if not sites:
        raise ValueError(
            f"No sites found in metadata: {metadata_path}"
        )

    return sites, site_order


def parse_base_counts(value, path, line_no):
    numbers = re.findall(r"\d+", value)

    if len(numbers) != 4:
        raise ValueError(
            f"Could not parse four A/C/G/T counts from "
            f"'{value}' at {path}:{line_no}."
        )

    return {
        base: int(number)
        for base, number in zip(BASES, numbers)
    }


def load_sample_counts(sample_dir, sample_id):
    if not os.path.isdir(sample_dir):
        raise FileNotFoundError(
            f"REDItoolKnown directory not found for {sample_id}: "
            f"{sample_dir}"
        )

    tables = sorted(
        glob.glob(
            os.path.join(
                sample_dir,
                "**",
                "outTable_*"
            ),
            recursive=True
        )
    )

    tables = [
        path for path in tables
        if os.path.isfile(path) and os.path.getsize(path) > 0
    ]

    if not tables:
        raise FileNotFoundError(
            f"No non-empty REDItoolKnown outTable files found "
            f"for {sample_id}: {sample_dir}"
        )

    counts = {}
    parsed_rows = 0

    for path in tables:

        with open(
            path,
            encoding="utf-8",
            errors="replace"
        ) as handle:

            reader = csv.reader(handle, delimiter="\t")

            for line_no, fields in enumerate(reader, start=1):

                if len(fields) < 7:
                    continue

                if not fields[1].strip().isdigit():
                    continue

                chrom = fields[0].strip()
                pos = int(fields[1])

                if not chrom or pos < 1:
                    raise ValueError(
                        f"Invalid genomic coordinate at {path}:{line_no}."
                    )

                key = (chrom, pos)

                base_counts = parse_base_counts(
                    fields[6],
                    path,
                    line_no
                )

                if key in counts:
                    if counts[key] != base_counts:
                        raise ValueError(
                            f"Conflicting base counts for "
                            f"{chrom}:{pos} in sample {sample_id}."
                        )
                else:
                    counts[key] = base_counts

                parsed_rows += 1

    return counts, len(tables), parsed_rows


def main():
    args = parse_args()

    samples = read_samples(
        args.sample_table,
        args.species
    )

    site_strands, site_order = read_site_metadata(
        args.metadata
    )

    sample_counts = {}
    sample_stats = {}

    # Load REDItoolKnown counts directly for each sample.
    for sample_id in samples:

        sample_dir = os.path.join(
            args.input_dir,
            sample_id
        )

        counts, n_tables, n_rows = load_sample_counts(
            sample_dir,
            sample_id
        )

        sample_counts[sample_id] = counts

        sample_stats[sample_id] = {
            "tables": n_tables,
            "rows": n_rows,
            "missing": 0
        }

    # Construct the final wide table.
    output_fields = [
        "chromosome",
        "position",
        "strand"
    ]

    for sample_id in samples:
        for base in BASES:
            output_fields.append(
                f"{sample_id}_{base}_count"
            )

    output = Path(args.output)
    output.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    tmp = output.with_name(
        output.name + f".tmp.{os.getpid()}"
    )

    try:
        with open(
            tmp,
            "w",
            encoding="utf-8",
            newline=""
        ) as handle:

            writer = csv.DictWriter(
                handle,
                fieldnames=output_fields,
                delimiter="\t",
                lineterminator="\n"
            )

            writer.writeheader()

            for chrom, pos in site_order:

                key = (chrom, pos)

                row = {
                    "chromosome": chrom,
                    "position": pos,
                    "strand": site_strands[key]
                }

                for sample_id in samples:

                    record = sample_counts[sample_id].get(key)

                    if record is None:
                        record = {
                            base: 0 for base in BASES
                        }

                        sample_stats[sample_id]["missing"] += 1

                    for base in BASES:
                        row[
                            f"{sample_id}_{base}_count"
                        ] = record[base]

                writer.writerow(row)

        os.replace(tmp, output)

    finally:
        if tmp.exists():
            tmp.unlink()

    # Print summary.
    print("=" * 60)
    print(f"Species: {args.species}")
    print(f"Number of samples: {len(samples)}")
    print(f"Unique genomic sites: {len(site_order)}")
    print(f"Output columns: {len(output_fields)}")
    print("-" * 60)

    for sample_id in samples:
        stats = sample_stats[sample_id]

        print(
            f"{sample_id}: "
            f"REDItoolKnown tables={stats['tables']}; "
            f"parsed rows={stats['rows']}; "
            f"missing sites filled with zero={stats['missing']}"
        )

    print("-" * 60)
    print(f"Output: {output}")
    print("=" * 60)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
