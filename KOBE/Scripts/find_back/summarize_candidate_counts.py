
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import argparse
from collections import OrderedDict


ANNOTATION_COLS = [
    "transcript",
    "cds_pos",
    "anno",
    "codon",
    "AA_before",
    "codon_pos",
    "edited_codon",
    "AA_after",
    "handle"
]

COUNT_COLS = [
    "A_reads",
    "G_reads",
    "editing_level"
]


def read_samples(sample_table):

    species_samples = OrderedDict()

    with open(sample_table, encoding="utf-8-sig", newline="") as f:

        reader = csv.DictReader(f, delimiter="\t")

        if not reader.fieldnames:
            raise ValueError("Empty sample table")

        reader.fieldnames = [x.strip() for x in reader.fieldnames]

        def find_col(names):
            for name in names:
                if name in reader.fieldnames:
                    return name
            return None

        sample_col = find_col(["sample_id", "sample", "sample_name"])
        species_col = find_col(["species"])
        type_col = find_col(["type"])

        if not sample_col or not species_col:
            raise ValueError(
                "Sample table requires sample_id and species columns"
            )

        for row in reader:

            sample = (row.get(sample_col) or "").strip()
            species = (row.get(species_col) or "").strip()

            if type_col:
                sample_type = (row.get(type_col) or "").strip().upper()
                if sample_type and sample_type != "RNA":
                    continue

            if not sample or not species:
                raise ValueError("Empty sample_id or species")

            species_samples.setdefault(species, [])

            if sample in species_samples[species]:
                raise ValueError(
                    f"Duplicate sample: {species}/{sample}"
                )

            species_samples[species].append(sample)

    return species_samples


def parse_count(value):

    value = str(value).strip()

    if value in ("", "NA", "NaN", "-", "."):
        return 0

    number = float(value)

    if number < 0 or not number.is_integer():
        raise ValueError(f"Invalid read count: {value}")

    return int(number)


def parse_level(value):

    value = str(value).strip()

    if value in ("", "NA", "NaN", "-", "."):
        return None

    level = float(value)

    if not 0 <= level <= 1:
        raise ValueError(f"Invalid editing level: {value}")

    return level


def summarize_species(species, samples, count_dir):

    print("\n========================================")
    print(f"Processing species: {species}")
    print(f"Samples: {len(samples)}")
    print("========================================")

    sites = OrderedDict()

    # --------------------------------------------------------
    # Read per-sample results
    # --------------------------------------------------------

    for sample in samples:

        input_file = os.path.join(
            count_dir,
            species,
            f"{sample}.candidate_counts.tsv"
        )

        if not os.path.isfile(input_file):
            raise FileNotFoundError(
                f"Missing sample result: {input_file}"
            )

        print(f"Reading: {sample}")

        sample_sites = set()

        with open(input_file, encoding="utf-8-sig", newline="") as f:

            reader = csv.DictReader(f, delimiter="\t")

            required = ANNOTATION_COLS + COUNT_COLS

            missing = [c for c in required if c not in (reader.fieldnames or [])]

            if missing:
                raise ValueError(
                    f"{input_file}: missing columns {missing}"
                )

            for line_number, row in enumerate(reader, start=2):

                transcript = row["transcript"]
                cds_pos = row["cds_pos"]

                key = (transcript, cds_pos)

                if key in sample_sites:
                    raise ValueError(
                        f"Duplicate physical site in {sample}: {key}"
                    )

                sample_sites.add(key)

                annotation = [row[c] for c in ANNOTATION_COLS]

                a_reads = parse_count(row["A_reads"])
                g_reads = parse_count(row["G_reads"])
                level = parse_level(row["editing_level"])

                coverage = a_reads + g_reads

                if coverage == 0 and level is not None:
                    raise ValueError(
                        f"Invalid editing level at {key} in {sample}: "
                        "zero coverage but non-NA level"
                    )

                if coverage > 0:
                    calculated_level = g_reads / coverage

                    if level is None or abs(level - calculated_level) > 1e-7:
                        raise ValueError(
                            f"Editing level inconsistent with A/G reads: "
                            f"{key} in {sample}"
                        )

                if key not in sites:

                    sites[key] = {
                        "annotation": annotation,
                        "samples": {}
                    }

                else:

                    if sites[key]["annotation"] != annotation:
                        raise ValueError(
                            f"Conflicting annotations for physical site: {key}"
                        )

                sites[key]["samples"][sample] = {
                    "A_reads": a_reads,
                    "G_reads": g_reads,
                    "editing_level": level
                }

        print(f"  Sites: {len(sample_sites)}")

    if not sites:
        raise ValueError(f"No sites found for {species}")

    # --------------------------------------------------------
    # Prepare output columns
    # --------------------------------------------------------

    header = ANNOTATION_COLS.copy()

    for sample in samples:

        header.extend([
            f"{sample}_A_reads",
            f"{sample}_G_reads",
            f"{sample}_editing_level"
        ])

    header.append("n_samples_editing_level_gt0")

    # --------------------------------------------------------
    # Write summary
    # --------------------------------------------------------

    output_file = os.path.join(
        count_dir,
        species,
        f"{species}_all_samples_summary.tsv"
    )

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    total_positive = 0

    with open(output_file, "w", encoding="utf-8", newline="") as f:

        writer = csv.writer(
            f,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow(header)

        for key, record in sites.items():

            output_row = record["annotation"]

            n_positive = 0

            for sample in samples:

                sample_data = record["samples"].get(sample)

                if sample_data is None:

                    a_reads = 0
                    g_reads = 0
                    level = "NA"

                else:

                    a_reads = sample_data["A_reads"]
                    g_reads = sample_data["G_reads"]

                    level = sample_data["editing_level"]

                    if level is None:
                        level = "NA"
                    else:
                        if level > 0:
                            n_positive += 1

                        level = f"{level:.8f}"

                output_row.extend([
                    a_reads,
                    g_reads,
                    level
                ])

            output_row.append(n_positive)

            if n_positive > 0:
                total_positive += 1

            writer.writerow(output_row)

    print("----------------------------------------")
    print(f"Species: {species}")
    print(f"Samples: {len(samples)}")
    print(f"Unique physical sites: {len(sites)}")
    print(f"Sites edited in at least one sample: {total_positive}")
    print(f"Output: {output_file}")

    return output_file


def main():

    parser = argparse.ArgumentParser(
        description="Summarize per-sample RNA editing counts by species"
    )

    parser.add_argument("--sample_table", required=True)
    parser.add_argument("--count_dir", required=True)

    args = parser.parse_args()

    species_samples = read_samples(args.sample_table)

    print("===================================================")
    print("Multi-sample RNA editing summary")
    print("===================================================")

    print(f"Species: {len(species_samples)}")

    for species, samples in species_samples.items():

        summarize_species(
            species,
            samples,
            args.count_dir
        )

    print("\n===================================================")
    print("All species summaries completed")
    print("===================================================")


if __name__ == "__main__":
    main()
