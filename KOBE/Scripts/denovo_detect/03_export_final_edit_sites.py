#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Export final RNA editing sites from Module 02 selected results.

Input:
    <input_root>/<species>/genic.selected.tsv
    <input_root>/<species>/intergenic.selected.tsv
    samples.tsv

Output:
    <output_root>/<species>/genic.final_sites.tsv
    <output_root>/<species>/intergenic.final_sites.tsv
    <output_root>/<species>/all.final_sites.tsv

For A-to-I editing:
    Genomic A>G (+): transcript A reads = REF; G reads = ALT
    Genomic T>C (-): transcript A reads = REF; G reads = ALT

Editing level:
    G / (A + G)

All RNA libraries are retained independently.
"""

import argparse
import csv
import os
import re
import sys
import tempfile

from pathlib import Path


BASE_COLUMNS = [
    "chr",
    "pos",
    "strand",
    "ref",
    "alt",
    "genomic_mut",
    "transcript_mut",
    "region",
    "anno",
    "gene_id",
    "transcript",
    "cds_pos",
    "AA_pos",
    "AA_change",
    "DNA_pass_sample_names",
    "RNA_pass_samples",
    "n_RNA_pass",
    "edit_level_pass_samples",
]


def natural_key(value):
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", str(value))
    ]


def read_sample_table(path):
    if not path.is_file():
        raise FileNotFoundError(f"Sample table not found: {path}")

    by_species = {}
    seen_global = set()

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"sample_id", "sample_name", "species", "type"}

        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(
                "samples.tsv must contain columns: "
                "sample_id, sample_name, species, type"
            )

        for line_no, row in enumerate(reader, start=2):

            sid = (row.get("sample_id") or "").strip()
            sample_name = (row.get("sample_name") or "").strip()
            species = (row.get("species") or "").strip()
            sample_type = (row.get("type") or "").strip().upper()

            if (
                not sid
                or not sample_name
                or not species
                or sample_type not in {"DNA", "RNA"}
            ):
                raise ValueError(
                    f"Invalid sample metadata at line {line_no}: {row}"
                )

            if sid in seen_global:
                raise ValueError(f"Duplicate sample_id: {sid}")

            seen_global.add(sid)

            by_species.setdefault(species, []).append({
                "sample_id": sid,
                "sample_name": sample_name,
                "species": species,
                "type": sample_type,
            })

    for species, rows in by_species.items():
        if not any(row["type"] == "RNA" for row in rows):
            raise ValueError(f"No RNA libraries found for species: {species}")

    return by_species


def parse_int(value, label, default_zero=False):

    if value is None or str(value).strip() in {"", "NA", "."}:
        if default_zero:
            return 0
        raise ValueError(f"Missing integer value for {label}")

    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Invalid integer value for {label}: {value!r}"
        ) from exc

    if result < 0:
        raise ValueError(f"Negative read count for {label}: {result}")

    return result


def infer_orientation(row, expected_region):

    """
    Infer transcript strand from genomic REF/ALT.

    A>G: positive strand.
    T>C: negative strand.

    For both orientations, transcript mutation is A>G.
    """

    ref = (row.get("ref") or "").strip().upper()
    alt = (row.get("alt") or "").strip().upper()
    mut = (row.get("mut") or "").strip().upper()

    genomic_mut = f"{ref}>{alt}"

    if (ref, alt) == ("A", "G"):
        strand = "+"

    elif (ref, alt) == ("T", "C"):
        strand = "-"

    else:
        raise ValueError(
            f"Unexpected genomic REF/ALT at "
            f"{row.get('chr')}:{row.get('pos')}: "
            f"{ref}>{alt}; expected A>G or T>C"
        )

    if expected_region == "genic":

        # SnpEff HGVS.c is transcript-oriented.
        if mut != "A>G":
            raise ValueError(
                f"Genic site is not transcript-oriented A>G at "
                f"{row.get('chr')}:{row.get('pos')}: mut={mut}"
            )

    elif expected_region == "intergenic":

        # Intergenic mutation is genomic-oriented.
        if mut != genomic_mut:
            raise ValueError(
                f"Intergenic mutation mismatch at "
                f"{row.get('chr')}:{row.get('pos')}: "
                f"mut={mut}, genomic={genomic_mut}"
            )

    return strand, genomic_mut, "A>G"


def read_selected_table(path, region, rna_samples):

    if not path.is_file():
        raise FileNotFoundError(
            f"Selected-site table not found: {path}"
        )

    rows = []

    with path.open("r", encoding="utf-8-sig", newline="") as handle:

        reader = csv.DictReader(handle, delimiter="\t")

        if not reader.fieldnames:
            raise ValueError(f"Missing header: {path}")

        fields = set(reader.fieldnames)

        required = {
            "chr",
            "pos",
            "ref",
            "alt",
            "mut",
            "region",
            "anno",
            "gene_id",
            "transcript",
            "cds_pos",
            "AA_pos",
            "AA_change",
            "DNA_pass_sample_names",
            "RNA_pass_samples",
            "n_RNA_pass",
            "edit_level_pass_samples",
        }

        missing = required - fields

        if missing:
            raise ValueError(
                f"{path} is missing columns: {sorted(missing)}"
            )

        # Check that every RNA library has its own REF/ALT counts.
        for sample in rna_samples:

            sid = sample["sample_id"]

            for col in (f"R_{sid}", f"A_{sid}"):

                if col not in fields:
                    raise ValueError(
                        f"{path} is missing read-count column: {col}"
                    )

        for line_no, row in enumerate(reader, start=2):

            if not row.get("chr") or not row.get("pos"):
                raise ValueError(
                    f"Missing chromosome/position at line {line_no}"
                )

            if row.get("region") != region:
                raise ValueError(
                    f"Wrong region at line {line_no}: "
                    f"expected {region}, got {row.get('region')}"
                )

            try:
                pos = int(row["pos"])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid position at line {line_no}: {row['pos']}"
                ) from exc

            if pos < 1:
                raise ValueError(
                    f"Position must be >= 1 at line {line_no}"
                )

            strand, genomic_mut, transcript_mut = infer_orientation(
                row, region
            )

            n_pass = parse_int(
                row.get("n_RNA_pass"),
                f"n_RNA_pass at {row['chr']}:{pos}"
            )

            passed_ids = [
                x for x in (row.get("RNA_pass_samples") or "").split(";")
                if x
            ]

            if len(set(passed_ids)) != len(passed_ids):
                raise ValueError(
                    f"Duplicate RNA sample_id in RNA_pass_samples at "
                    f"{row['chr']}:{pos}"
                )

            if len(passed_ids) != n_pass:
                raise ValueError(
                    f"n_RNA_pass ({n_pass}) does not match "
                    f"RNA_pass_samples ({len(passed_ids)}) at "
                    f"{row['chr']}:{pos}"
                )

            out = {
                "chr": row["chr"],
                "pos": pos,
                "strand": strand,
                "ref": row["ref"].upper(),
                "alt": row["alt"].upper(),
                "genomic_mut": genomic_mut,
                "transcript_mut": transcript_mut,
                "region": region,
                "anno": row.get("anno", ""),
                "gene_id": row.get("gene_id", ""),
                "transcript": row.get("transcript", ""),
                "cds_pos": row.get("cds_pos", ""),
                "AA_pos": row.get("AA_pos", ""),
                "AA_change": row.get("AA_change", ""),
                "DNA_pass_sample_names": row.get(
                    "DNA_pass_sample_names", ""
                ),
                "RNA_pass_samples": row.get(
                    "RNA_pass_samples", ""
                ),
                "n_RNA_pass": n_pass,
                "edit_level_pass_samples": row.get(
                    "edit_level_pass_samples", ""
                ),
            }

            # Convert genomic REF/ALT counts into transcript-oriented
            # A/G counts.
            #
            # Positive strand A>G:
            #   A = genomic REF
            #   G = genomic ALT
            #
            # Negative strand T>C:
            #   A = genomic REF (T on genomic negative strand)
            #   G = genomic ALT (C on genomic negative strand)
            #
            # Both therefore use the same REF/ALT count mapping.

            for sample in rna_samples:

                sid = sample["sample_id"]

                ref_count = parse_int(
                    row.get(f"R_{sid}"),
                    f"R_{sid} at {row['chr']}:{pos}",
                    default_zero=True
                )

                alt_count = parse_int(
                    row.get(f"A_{sid}"),
                    f"A_{sid} at {row['chr']}:{pos}",
                    default_zero=True
                )

                a_count = ref_count
                g_count = alt_count

                coverage = a_count + g_count

                edit_level = (
                    f"{g_count / coverage:.6f}"
                    if coverage > 0
                    else "NA"
                )

                out[f"A_{sid}"] = a_count
                out[f"G_{sid}"] = g_count
                out[f"editing_level_{sid}"] = edit_level
                out[f"coverage_{sid}"] = coverage

            rows.append(out)

    return rows


def validate_unique_sites(rows, label):

    seen = {}

    for row in rows:

        key = (row["chr"], row["pos"])

        if key in seen:

            previous = seen[key]

            raise ValueError(
                f"Duplicate genomic site in {label}: "
                f"{key[0]}:{key[1]} "
                f"({previous['region']} and {row['region']}). "
                f"Refusing to silently collapse records."
            )

        seen[key] = row


def write_tsv(path, rows, columns):

    path.parent.mkdir(parents=True, exist_ok=True)

    fd, temp_name = tempfile.mkstemp(
        prefix=path.name + ".",
        suffix=".tmp",
        dir=path.parent
    )

    try:

        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:

            writer = csv.DictWriter(
                handle,
                fieldnames=columns,
                delimiter="\t",
                lineterminator="\n",
                extrasaction="ignore"
            )

            writer.writeheader()

            for row in rows:
                writer.writerow({
                    col: row.get(col, "")
                    for col in columns
                })

        os.replace(temp_name, path)

    except Exception:

        try:
            os.unlink(temp_name)
        except OSError:
            pass

        raise


def export_species(species, input_dir, output_root, samples):

    rna_samples = sorted(
        (row for row in samples if row["type"] == "RNA"),
        key=lambda x: x["sample_id"]
    )

    if not rna_samples:
        raise ValueError(
            f"No RNA samples for species: {species}"
        )

    genic_path = input_dir / "genic.selected.tsv"
    intergenic_path = input_dir / "intergenic.selected.tsv"

    genic_rows = read_selected_table(
        genic_path, "genic", rna_samples
    )

    intergenic_rows = read_selected_table(
        intergenic_path, "intergenic", rna_samples
    )

    validate_unique_sites(
        genic_rows, f"{species}/genic"
    )

    validate_unique_sites(
        intergenic_rows, f"{species}/intergenic"
    )

    all_rows = genic_rows + intergenic_rows

    validate_unique_sites(
        all_rows, f"{species}/all"
    )

    genic_rows.sort(
        key=lambda x: (natural_key(x["chr"]), x["pos"])
    )

    intergenic_rows.sort(
        key=lambda x: (natural_key(x["chr"]), x["pos"])
    )

    all_rows.sort(
        key=lambda x: (natural_key(x["chr"]), x["pos"])
    )

    columns = list(BASE_COLUMNS)

    for sample in rna_samples:

        sid = sample["sample_id"]

        columns.extend([
            f"A_{sid}",
            f"G_{sid}",
            f"editing_level_{sid}",
            f"coverage_{sid}",
        ])

    out_dir = output_root / species

    write_tsv(
        out_dir / "genic.final_sites.tsv",
        genic_rows,
        columns
    )

    write_tsv(
        out_dir / "intergenic.final_sites.tsv",
        intergenic_rows,
        columns
    )

    write_tsv(
        out_dir / "all.final_sites.tsv",
        all_rows,
        columns
    )

    return (
        len(genic_rows),
        len(intergenic_rows),
        len(all_rows)
    )


def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument(
        "--input-root",
        type=Path,
        required=True
    )

    parser.add_argument(
        "--samples",
        type=Path,
        required=True
    )

    parser.add_argument(
        "--output-root",
        type=Path,
        required=True
    )

    parser.add_argument(
        "--species",
        nargs="+",
        default=None
    )

    args = parser.parse_args()

    input_root = args.input_root.resolve()
    output_root = args.output_root.resolve()

    sample_map = read_sample_table(
        args.samples.resolve()
    )

    if args.species:

        species_list = args.species

    else:

        species_list = sorted(
            p.name
            for p in input_root.iterdir()
            if (
                p.is_dir()
                and (p / "genic.selected.tsv").is_file()
                and (p / "intergenic.selected.tsv").is_file()
            )
        ) if input_root.is_dir() else []

    if not species_list:
        raise FileNotFoundError(
            f"No species folders containing genic.selected.tsv "
            f"and intergenic.selected.tsv found under {input_root}"
        )

    total_genic = 0
    total_intergenic = 0
    total_all = 0

    for species in species_list:

        if species not in sample_map:
            raise ValueError(
                f"Species {species} is absent from samples.tsv"
            )

        input_dir = input_root / species

        n_genic, n_intergenic, n_all = export_species(
            species,
            input_dir,
            output_root,
            sample_map[species]
        )

        total_genic += n_genic
        total_intergenic += n_intergenic
        total_all += n_all

        print(
            f"[{species}] "
            f"genic={n_genic}; "
            f"intergenic={n_intergenic}; "
            f"all={n_all}"
        )

        print(
            f"[{species}] output: {output_root / species}"
        )

    print(
        f"TOTAL: genic={total_genic}; "
        f"intergenic={total_intergenic}; "
        f"all={total_all}"
    )

    print("Export completed successfully.")


if __name__ == "__main__":

    try:
        main()

    except (OSError, ValueError, csv.Error) as exc:
        sys.exit(f"ERROR: {exc}")