#!/usr/bin/env python3
"""Convert a supported_sites_ge1 table into an A/G count matrix."""

import argparse
import csv
import os
import re
import sys


def parse_count(value, column, line_number):
    value = str(value).strip()
    if value in ("", "NA", "NaN", "nan", ".", "-"):
        return 0

    try:
        number = float(value)
    except ValueError:
        raise ValueError(
            f"Invalid read count in column '{column}' at line "
            f"{line_number}: {value}"
        )

    if number < 0 or not number.is_integer():
        raise ValueError(
            f"Read count must be a non-negative integer in column "
            f"'{column}' at line {line_number}: {value}"
        )

    return int(number)


def main():
    parser = argparse.ArgumentParser(
        prog="supported_sites_to_editing_counts",
        description=(
            "Extract transcript:cds_pos IDs and per-sample A/G read "
            "counts from a supported_sites_ge1 TSV."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n  KOBE supported_sites_to_editing_counts \\\n    --input /project/results/invicta_supported_sites_ge1.tsv \\\n    --species invicta \\\n    --work_dir /project/results/diff_edit_site_data \\\n    --sif /project/rna_editing.sif --bind /project",
    )
    parser.add_argument(
        "--input", required=True, metavar="FILE",
        help="Input supported_sites_ge1 TSV file."
    )
    parser.add_argument(
        "--species", required=True, metavar="STRING",
        help="Species name used for the output filename."
    )
    parser.add_argument(
        "--work_dir", required=True, metavar="DIR",
        help="Output working directory."
    )
    parser.add_argument(
        "--sif",
        metavar="FILE",
        help="Apptainer SIF image (or set RNA_EDITING_SIF). Consumed by the command launcher."
    )
    parser.add_argument(
        "--bind", "-B",
        action="append",
        metavar="DIR",
        help="Additional Apptainer bind specification; may be repeated. Consumed by the command launcher."
    )

    args = parser.parse_args()

    input_file = os.path.abspath(args.input)
    work_dir = os.path.abspath(args.work_dir)
    output_file = os.path.join(
        work_dir, f"{args.species}_editing_counts.tsv"
    )

    if not os.path.isfile(input_file):
        sys.exit(f"ERROR: Input file not found: {input_file}")

    os.makedirs(work_dir, exist_ok=True)

    tmp_file = output_file + ".tmp"

    try:
        with open(
            input_file, "r", encoding="utf-8-sig", newline=""
        ) as fin:
            reader = csv.DictReader(fin, delimiter="\t")

            if not reader.fieldnames:
                raise ValueError("Input TSV has no header.")

            reader.fieldnames = [
                name.strip() for name in reader.fieldnames
            ]

            required = {"transcript", "cds_pos"}
            missing = required.difference(reader.fieldnames)
            if missing:
                raise ValueError(
                    f"Missing required columns: {sorted(missing)}"
                )

            sample_columns = []
            seen_samples = set()

            for column in reader.fieldnames:
                match = re.fullmatch(r"(.+)_A_reads", column)
                if not match:
                    continue

                sample = match.group(1)
                g_column = f"{sample}_G_reads"

                if g_column not in reader.fieldnames:
                    raise ValueError(
                        f"Found '{column}' but missing '{g_column}'."
                    )

                if sample in seen_samples:
                    raise ValueError(
                        f"Duplicate sample columns detected: {sample}"
                    )

                sample_columns.append((sample, column, g_column))
                seen_samples.add(sample)

            if not sample_columns:
                raise ValueError(
                    "No sample columns matching '*_A_reads' and "
                    "'*_G_reads' were found."
                )

            output_header = ["ID"]
            for sample, a_column, g_column in sample_columns:
                output_header.extend([a_column, g_column])

            seen_ids = set()
            n_rows = 0

            with open(
                tmp_file, "w", encoding="utf-8", newline=""
            ) as fout:
                writer = csv.writer(
                    fout, delimiter="\t", lineterminator="\n"
                )
                writer.writerow(output_header)

                for line_number, row in enumerate(reader, start=2):
                    transcript = (row.get("transcript") or "").strip()
                    cds_pos = (row.get("cds_pos") or "").strip()

                    if not transcript or not cds_pos:
                        raise ValueError(
                            f"Empty transcript or cds_pos at line "
                            f"{line_number}."
                        )

                    site_id = f"{transcript}:{cds_pos}"

                    if site_id in seen_ids:
                        raise ValueError(
                            f"Duplicate physical-site ID at line "
                            f"{line_number}: {site_id}"
                        )
                    seen_ids.add(site_id)

                    output_row = [site_id]

                    for sample, a_column, g_column in sample_columns:
                        output_row.append(
                            parse_count(
                                row.get(a_column), a_column, line_number
                            )
                        )
                        output_row.append(
                            parse_count(
                                row.get(g_column), g_column, line_number
                            )
                        )

                    writer.writerow(output_row)
                    n_rows += 1

        os.replace(tmp_file, output_file)

    except Exception as exc:
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
        sys.exit(f"ERROR: {exc}")

    print("Conversion completed.")
    print(f"Input: {input_file}")
    print(f"Species: {args.species}")
    print(f"Samples detected: {len(sample_columns)}")
    print(f"Data rows: {n_rows}")
    print(f"Output columns: {1 + 2 * len(sample_columns)}")
    print(f"Output: {output_file}")


if __name__ == "__main__":
    main()
