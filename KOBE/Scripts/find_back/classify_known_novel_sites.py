#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import re
import argparse
from collections import OrderedDict


SUMMARY_ANNOTATION_COLS = [
    "transcript",
    "cds_pos",
    "anno",
    "codon",
    "AA_before",
    "codon_pos",
    "edited_codon",
    "AA_after",
    "handle",
]

TRANSCRIPT_ALIASES = {
    "transcript",
    "transcript_id",
    "transcriptid",
    "transcript_name",
    "cds",
    "cds_id",
    "seqid",
    "seq_id",
}

POSITION_ALIASES = {
    "cds_pos",
    "cds_position",
    "cdspos",
    "position",
    "pos",
    "site_pos",
}


def normalize_header(value):
    value = value.strip().lstrip("#").strip().lower()
    value = re.sub(r"[\s\-]+", "_", value)
    return value


def normalize_position(value, source):
    value = str(value).strip()

    if value.endswith(".0"):
        value = value[:-2]

    if not value.isdigit():
        raise ValueError(
            f"Invalid CDS position '{value}' in {source}"
        )

    return str(int(value))


def split_fields(line, delimiter):
    if delimiter is None:
        return line.strip().split()

    return next(csv.reader([line], delimiter=delimiter))


def detect_delimiter(line):
    if "\t" in line:
        return "\t"
    if "," in line:
        return ","
    return None


def find_known_file(edit_txt_dir, species):
    species_key = species.strip().lower()

    matches = []

    for name in os.listdir(edit_txt_dir):
        path = os.path.join(edit_txt_dir, name)

        if not os.path.isfile(path):
            continue

        lower = name.lower()

        if (
            species_key in lower
            and "candi" in lower
            and "cds" in lower
            and "sites" in lower
            and lower.endswith((".txt", ".tsv", ".tab"))
        ):
            matches.append(path)

    if len(matches) == 0:
        raise FileNotFoundError(
            f"No original known-site file matched species '{species}' "
            f"in {edit_txt_dir}. Expected a filename containing the "
            "species name and 'candi_cds_sites'."
        )

    if len(matches) > 1:
        raise RuntimeError(
            f"Multiple known-site files matched species '{species}':\n"
            + "\n".join(sorted(matches))
        )

    return matches[0]


def read_known_sites(known_file):
    with open(known_file, encoding="utf-8-sig") as f:
        lines = [line.rstrip("\r\n") for line in f if line.strip()]

    if not lines:
        raise ValueError(f"Known-site file is empty: {known_file}")

    delimiter = detect_delimiter(lines[0])
    first_fields = split_fields(lines[0], delimiter)
    normalized = [normalize_header(x) for x in first_fields]

    transcript_idx = None
    position_idx = None
    has_header = False

    for i, name in enumerate(normalized):
        if name in TRANSCRIPT_ALIASES:
            transcript_idx = i
        if name in POSITION_ALIASES:
            position_idx = i

    if transcript_idx is not None and position_idx is not None:
        has_header = True
        data_lines = lines[1:]
    else:
        transcript_idx = 0
        position_idx = 1 if len(first_fields) > 1 else None
        data_lines = lines

    known_sites = set()

    for line_number, line in enumerate(
        data_lines,
        start=2 if has_header else 1
    ):
        fields = split_fields(line, delimiter)

        if not fields:
            continue

        if position_idx is None:
            if len(fields) == 1 and ":" in fields[0]:
                transcript, pos = fields[0].rsplit(":", 1)
            else:
                raise ValueError(
                    f"Cannot identify transcript and CDS position "
                    f"in {known_file}, line {line_number}: {line}"
                )
        else:
            if len(fields) <= max(transcript_idx, position_idx):
                raise ValueError(
                    f"Insufficient columns in {known_file}, "
                    f"line {line_number}: {line}"
                )

            transcript = fields[transcript_idx].strip()
            pos = fields[position_idx].strip()

        if not transcript or not pos:
            continue

        pos = normalize_position(pos, known_file)
        known_sites.add((transcript, pos))

    if not known_sites:
        raise ValueError(
            f"No valid transcript:CDS-position records found in {known_file}"
        )

    return known_sites


def read_sample_species(sample_table):
    species_samples = OrderedDict()

    with open(sample_table, encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        if not reader.fieldnames:
            raise ValueError(f"Empty sample table: {sample_table}")

        reader.fieldnames = [
            x.strip().lstrip("\ufeff") for x in reader.fieldnames
        ]

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
                raise ValueError(
                    "Empty sample_id or species in sample table"
                )

            species_samples.setdefault(species, [])

            if sample in species_samples[species]:
                raise ValueError(
                    f"Duplicate sample in sample table: {species}/{sample}"
                )

            species_samples[species].append(sample)

    if not species_samples:
        raise ValueError("No RNA samples found in sample table")

    return species_samples


def read_summary(summary_file, min_samples):
    with open(summary_file, encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        if not reader.fieldnames:
            raise ValueError(f"Empty summary file: {summary_file}")

        reader.fieldnames = [
            x.strip().lstrip("\ufeff") for x in reader.fieldnames
        ]

        required = [
            "transcript",
            "cds_pos",
            "n_samples_editing_level_gt0",
        ]

        missing = [x for x in required if x not in reader.fieldnames]

        if missing:
            raise ValueError(
                f"Missing required columns in {summary_file}: {missing}"
            )

        rows = []
        seen = set()

        for line_number, row in enumerate(reader, start=2):
            transcript = (row.get("transcript") or "").strip()
            pos = (row.get("cds_pos") or "").strip()

            if not transcript or not pos:
                raise ValueError(
                    f"Empty transcript/CDS position in "
                    f"{summary_file}, line {line_number}"
                )

            pos = normalize_position(pos, summary_file)
            key = (transcript, pos)

            if key in seen:
                raise ValueError(
                    f"Duplicate physical site {transcript}:{pos} "
                    f"in {summary_file}"
                )

            seen.add(key)

            try:
                n_positive = int(
                    float(row["n_samples_editing_level_gt0"])
                )
            except (TypeError, ValueError):
                raise ValueError(
                    f"Invalid n_samples_editing_level_gt0 at "
                    f"{summary_file}, line {line_number}"
                )

            if n_positive < 0:
                raise ValueError(
                    f"Negative sample count at {summary_file}, "
                    f"line {line_number}"
                )

            row["transcript"] = transcript
            row["cds_pos"] = pos
            row["n_samples_editing_level_gt0"] = str(n_positive)

            rows.append((key, row, n_positive))

    if not rows:
        raise ValueError(f"No site records found in {summary_file}")

    supported = [
        (key, row, n_positive)
        for key, row, n_positive in rows
        if n_positive >= min_samples
    ]

    return rows, supported


def write_tsv(output_file, fieldnames, rows):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def process_species(
    species,
    samples,
    summary_dir,
    edit_txt_dir,
    out_dir,
    min_samples,
):
    print("\n==============================================")
    print(f"Processing species: {species}")
    print(f"RNA samples: {len(samples)}")
    print("==============================================")

    summary_file = os.path.join(
        summary_dir,
        species,
        f"{species}_all_samples_summary.tsv",
    )

    if not os.path.isfile(summary_file):
        raise FileNotFoundError(
            f"Multi-sample summary not found: {summary_file}"
        )

    known_file = find_known_file(edit_txt_dir, species)
    print(f"Known-site reference: {known_file}")
    print(f"Multi-sample summary: {summary_file}")

    known_sites = read_known_sites(known_file)
    all_rows, supported_rows = read_summary(
        summary_file,
        min_samples,
    )

    # Count support among the original known-site set.
    # Known sites absent from the merged summary are assigned zero support.
    summary_support = {
        key: n_positive
        for key, row, n_positive in all_rows
    }

    known_support_counts = [
        summary_support.get(key, 0)
        for key in known_sites
    ]

    n_known_zero_support = sum(
        count == 0 for count in known_support_counts
    )
    n_known_exactly_one = sum(
        count == 1 for count in known_support_counts
    )
    n_known_at_least_two = sum(
        count >= 2 for count in known_support_counts
    )

    # --------------------------------------------------------
    # Classify every row in the complete multi-sample summary.
    # This classification is independent of the support threshold.
    # --------------------------------------------------------

    all_classified_rows = []

    for key, row, n_positive in all_rows:
        output_row = dict(row)
        output_row["known_status"] = (
            "Known" if key in known_sites else "Novel"
        )
        all_classified_rows.append(output_row)

    # Preserve the existing supported-site output, which contains
    # both Known and Novel sites meeting the support threshold.
    supported_output_rows = []
    novel_supported_output_rows = []

    n_known_supported = 0
    n_novel_supported = 0

    for key, row, n_positive in supported_rows:
        status = "Known" if key in known_sites else "Novel"

        output_row = dict(row)
        output_row["known_status"] = status
        supported_output_rows.append(output_row)

        if status == "Known":
            n_known_supported += 1
        else:
            n_novel_supported += 1
            novel_supported_output_rows.append(output_row)

    if n_known_supported + n_novel_supported != len(supported_rows):
        raise RuntimeError(
            "Known/Novel counts do not sum to supported-site total"
        )

    stats = {
        "species": species,
        "n_RNA_samples": len(samples),
        "min_samples_required": min_samples,
        "n_summary_unique_sites": len(all_rows),
        "n_sites_supported_ge_min_samples": len(supported_rows),
        "n_known_supported_sites_ge_threshold": n_known_supported,
        "n_novel_supported_sites_ge_threshold": n_novel_supported,
        "n_original_known_unique_sites": len(known_sites),
        "n_original_known_sites_zero_support": n_known_zero_support,
        "n_original_known_sites_exactly_one_sample": n_known_exactly_one,
        "n_original_known_sites_at_least_two_samples": n_known_at_least_two,
        "n_original_known_sites_below_threshold": sum(
            count < min_samples for count in known_support_counts
        ),
        "known_site_file": known_file,
        "summary_file": summary_file,
    }

    species_out_dir = os.path.join(out_dir, species)
    os.makedirs(species_out_dir, exist_ok=True)

    stats_file = os.path.join(
        species_out_dir,
        f"{species}_editing_site_statistics_ge{min_samples}.tsv",
    )

    supported_file = os.path.join(
        species_out_dir,
        f"{species}_supported_sites_ge{min_samples}.tsv",
    )

    write_tsv(stats_file, list(stats.keys()), [stats])

    summary_fieldnames = list(all_rows[0][1].keys())
    # Recalculate this column on every run rather than retaining stale values.
    summary_fieldnames = [
        col for col in summary_fieldnames if col != "known_status"
    ]
    summary_fieldnames.append("known_status")

    # Add known_status to the complete species-level summary.
    # Write to a temporary file first, then replace the original atomically.
    summary_tmp = summary_file + ".known_status.tmp"
    write_tsv(
        summary_tmp,
        summary_fieldnames,
        all_classified_rows,
    )
    os.replace(summary_tmp, summary_file)

    # Keep the existing table of all sites meeting the support threshold.
    write_tsv(
        supported_file,
        summary_fieldnames,
        supported_output_rows,
    )

    # New output: only Novel sites meeting the support threshold.
    novel_supported_file = os.path.join(
        species_out_dir,
        f"{species}_novel_supported_sites_ge{min_samples}.tsv",
    )
    write_tsv(
        novel_supported_file,
        summary_fieldnames,
        novel_supported_output_rows,
    )

    print("----------------------------------------------")
    print(f"Original known unique sites: {len(known_sites)}")
    print(
        f"Sites supported in >= {min_samples} samples: "
        f"{len(supported_rows)}"
    )
    print(f"Known supported sites: {n_known_supported}")
    print(f"Novel supported sites: {n_novel_supported}")
    print(f"Classified all-site summary: {summary_file}")
    print(f"Supported sites (Known + Novel): {supported_file}")
    print(f"Novel supported sites only: {novel_supported_file}")

    return stats


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Classify multi-sample-supported RNA editing candidates "
            "as known or novel using species-specific original site sets."
        )
    )

    parser.add_argument("--sample_table", required=True)
    parser.add_argument("--summary_dir", required=True)
    parser.add_argument("--edit_txt_dir", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--min_samples", type=int, default=2)

    args = parser.parse_args()

    if args.min_samples < 1:
        parser.error("--min_samples must be >= 1")

    for path in [
        args.sample_table,
        args.summary_dir,
        args.edit_txt_dir,
    ]:
        if not os.path.exists(path):
            parser.error(f"Path does not exist: {path}")

    species_samples = read_sample_species(args.sample_table)

    print("===================================================")
    print("RNA editing site support and known/novel statistics")
    print("===================================================")
    print(f"Species: {len(species_samples)}")
    print(f"Minimum supporting samples: {args.min_samples}")

    all_stats = []

    for species, samples in species_samples.items():
        stats = process_species(
            species=species,
            samples=samples,
            summary_dir=args.summary_dir,
            edit_txt_dir=args.edit_txt_dir,
            out_dir=args.out_dir,
            min_samples=args.min_samples,
        )
        all_stats.append(stats)

    master_file = os.path.join(
        args.out_dir,
        f"editing_site_statistics_ge{args.min_samples}_all_species.tsv",
    )

    write_tsv(
        master_file,
        list(all_stats[0].keys()),
        all_stats,
    )

    print("\n===================================================")
    print("All species completed")
    print(f"Combined statistics: {master_file}")
    print("===================================================")


if __name__ == "__main__":
    main()