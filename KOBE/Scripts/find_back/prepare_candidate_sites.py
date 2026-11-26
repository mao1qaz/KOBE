
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import re
import shutil
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate REDItoolKnown site files for all species in the sample table."
    )

    parser.add_argument(
        "--sample_table",
        required=True,
        help="Input samples.tsv"
    )

    parser.add_argument(
        "--cds_pos_dir",
        required=True,
        help="Directory containing species-specific *_cds_pos.list files"
    )

    parser.add_argument(
        "--out_dir",
        required=True,
        help="find_back_result directory"
    )

    return parser.parse_args()


def read_samples(table):
    required = ["sample_id", "species"]

    with open(table, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")

        if not reader.fieldnames:
            raise ValueError("Empty sample table")

        reader.fieldnames = [x.strip() for x in reader.fieldnames]

        missing = [x for x in required if x not in reader.fieldnames]

        if missing:
            raise ValueError(
                "Missing columns: " + ", ".join(missing)
            )

        rows = []

        for line_number, row in enumerate(reader, start=2):
            row = {
                k.strip(): (v or "").strip()
                for k, v in row.items()
                if k is not None
            }

            sample = row["sample_id"]
            species = row["species"]

            if not sample or not species:
                raise ValueError(
                    f"Missing sample_id or species at line {line_number}"
                )

            if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample):
                raise ValueError(f"Invalid sample_id: {sample}")

            if not re.fullmatch(r"[A-Za-z0-9_.-]+", species):
                raise ValueError(f"Invalid species: {species}")

            rows.append((sample, species))

    if not rows:
        raise ValueError("No samples found")

    return rows


def generate_known_sites(input_file, output_file):
    records = []

    with open(input_file, encoding="utf-8-sig") as f:

        for line_number, line in enumerate(f, start=1):

            line = line.strip()

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) < 2:
                raise ValueError(
                    f"Invalid row in {input_file}, line {line_number}"
                )

            transcript = fields[0].strip()
            position = fields[1].replace("$", "").strip()

            if not transcript or not position.isdigit():
                raise ValueError(
                    f"Invalid transcript/CDS position in "
                    f"{input_file}, line {line_number}: {line}"
                )

            records.append((transcript, int(position)))

    # Keep duplicates, matching the original Bash workflow.
    records.sort(key=lambda x: (x[0], x[1]))

    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as f:
        for transcript, position in records:
            f.write(f"{transcript}\t{position}\t2\n")

    return len(records)


def main():

    args = parse_args()

    sample_table = Path(args.sample_table).resolve()
    cds_pos_dir = Path(args.cds_pos_dir).resolve()
    out_dir = Path(args.out_dir).resolve()

    mapping_dir = out_dir / "Mapping_BWA_to_CDS"
    known_dir = out_dir / "known_sites"

    if not sample_table.is_file():
        raise FileNotFoundError(sample_table)

    if not cds_pos_dir.is_dir():
        raise NotADirectoryError(cds_pos_dir)

    samples = read_samples(sample_table)

    # Group samples by species, preserving sample order.
    species_samples = {}

    for sample, species in samples:
        species_samples.setdefault(species, []).append(sample)

    print("===================================================")
    print("Preparing REDItoolKnown input files")
    print("===================================================")
    print(f"Sample table: {sample_table}")
    print(f"CDS position directory: {cds_pos_dir}")
    print(f"Output directory: {out_dir}")
    print(f"Species: {len(species_samples)}")
    print(f"Samples: {len(samples)}")

    # Check all required species files before processing.
    species_files = {}

    for species in species_samples:

        input_file = cds_pos_dir / f"{species}_cds_pos.list"

        if not input_file.is_file():
            raise FileNotFoundError(
                f"Missing CDS position list for {species}: {input_file}"
            )

        species_files[species] = input_file

    # Generate one knownEditingSites file per species.
    for species, sample_list in species_samples.items():

        print("\n===================================================")
        print(f"Processing species: {species}")
        print(f"Number of samples: {len(sample_list)}")
        print("===================================================")

        input_file = species_files[species]

        species_output = (
            known_dir / f"{species}_knownEditingSites.tab"
        )

        n_sites = generate_known_sites(
            input_file,
            species_output
        )

        print(f"Input: {input_file}")
        print(f"Sites: {n_sites}")
        print(f"Generated: {species_output}")

        # Copy the species-specific file to each sample mapping directory.
        for sample in sample_list:

            sample_dir = mapping_dir / species / sample
            sample_dir.mkdir(parents=True, exist_ok=True)

            target = sample_dir / "knownEditingSites.tab"

            shutil.copy2(species_output, target)

            print(f"Copied: {species}/{sample}")

    print("\n===================================================")
    print("All knownEditingSites files generated")
    print("===================================================")

    print(f"Species processed: {len(species_samples)}")
    print(f"Samples processed: {len(samples)}")
    print(f"Shared files: {known_dir}")
    print(f"Sample files: {mapping_dir}")


if __name__ == "__main__":
    main()
