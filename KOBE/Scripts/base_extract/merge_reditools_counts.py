#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import re
import sys


BASES = ("A", "C", "G", "T")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge REDItoolKnown base counts back to the complete input site list."
    )
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--redi-dir", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def parse_base_counts(value, path, line_number):
    numbers = re.findall(r"\d+", value)
    if len(numbers) != 4:
        raise ValueError(
            f"Could not parse four A/C/G/T counts from '{value}' at {path}:{line_number}."
        )
    return dict(zip(BASES, map(int, numbers)))


def load_redi_tables(redi_dir):
    tables = sorted(
        glob.glob(os.path.join(redi_dir, "**", "outTable_*"), recursive=True)
    )
    tables = [path for path in tables if os.path.isfile(path) and os.path.getsize(path) > 0]
    if not tables:
        raise FileNotFoundError(f"No non-empty REDItoolKnown outTable files found under: {redi_dir}")

    counts = {}
    for path in tables:
        with open(path, encoding="utf-8", errors="replace") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for line_number, fields in enumerate(reader, start=1):
                if not fields or len(fields) < 7:
                    continue
                if len(fields) >= 2 and not fields[1].strip().isdigit():
                    continue

                chrom = fields[0].strip()
                pos = int(fields[1])
                key = (chrom, pos)
                base_counts = parse_base_counts(fields[6], path, line_number)

                coverage = fields[4].strip() if len(fields) > 4 else "NA"
                mean_q = fields[5].strip() if len(fields) > 5 else "NA"
                all_subs = fields[7].strip() if len(fields) > 7 else "NA"
                frequency = fields[8].strip() if len(fields) > 8 else "NA"

                record = {
                    **base_counts,
                    "coverage_q25": coverage,
                    "meanQ": mean_q,
                    "all_substitutions": all_subs,
                    "frequency": frequency,
                }

                if key in counts and counts[key] != record:
                    raise ValueError(f"Conflicting REDItoolKnown records for {chrom}:{pos}.")
                counts[key] = record

    return counts, len(tables)


def main():
    args = parse_args()
    if not os.path.isfile(args.metadata):
        raise FileNotFoundError(f"Site metadata not found: {args.metadata}")

    redi_counts, n_tables = load_redi_tables(args.redi_dir)

    with open(args.metadata, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Site metadata has no header.")
        input_fields = list(reader.fieldnames)
        required = {
            "_input_index", "_chromosome", "_position", "_strand_code",
            "_strand_provided", "_strand_label"
        }
        missing = required - set(input_fields)
        if missing:
            raise ValueError(f"Site metadata is missing internal columns: {sorted(missing)}")
        rows = list(reader)

    internal_fields = {
        "_input_index", "_chromosome", "_position", "_strand_code",
        "_strand_provided", "_strand_label"
    }
    original_fields = [name for name in input_fields if name not in internal_fields]
    output_fields = original_fields + [
        "A_count", "C_count", "G_count", "T_count",
        "coverage_q25", "meanQ", "all_substitutions", "frequency",
        "editing_level", "count_orientation"
    ]

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    n_zero = 0
    n_editable = 0

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()

        for row in rows:
            chrom = row["_chromosome"]
            pos = int(row["_position"])
            record = redi_counts.get((chrom, pos))

            if record is None:
                n_zero += 1
                record = {
                    "A": 0, "C": 0, "G": 0, "T": 0,
                    "coverage_q25": "0", "meanQ": "NA",
                    "all_substitutions": "-", "frequency": "NA"
                }

            a_count = int(record["A"])
            c_count = int(record["C"])
            g_count = int(record["G"])
            t_count = int(record["T"])
            strand_provided = row["_strand_provided"] == "1"
            strand_label = row["_strand_label"]

            if strand_provided and strand_label in ("+", "-") and (a_count + g_count) > 0:
                editing_level = f"{g_count / (a_count + g_count):.8f}"
                n_editable += 1
            else:
                editing_level = "NA"

            if strand_label in ("+", "-"):
                orientation = "transcript"
            else:
                orientation = "genomic_reference"

            output_row = {name: row[name] for name in original_fields}
            output_row.update({
                "A_count": a_count,
                "C_count": c_count,
                "G_count": g_count,
                "T_count": t_count,
                "coverage_q25": record["coverage_q25"],
                "meanQ": record["meanQ"],
                "all_substitutions": record["all_substitutions"],
                "frequency": record["frequency"],
                "editing_level": editing_level,
                "count_orientation": orientation,
            })
            writer.writerow(output_row)

    print(f"REDItoolKnown tables parsed: {n_tables}")
    print(f"Input site rows retained: {len(rows)}")
    print(f"Rows with no REDItoolKnown record (zero-filled): {n_zero}")
    print(f"Rows with editing level calculated: {n_editable}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
