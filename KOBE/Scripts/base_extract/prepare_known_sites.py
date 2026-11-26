
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import os
import re
import sys
from pathlib import Path


STRAND_MAP = {
    "+": ("1", True, "+"),
    "1": ("1", True, "+"),
    "-": ("0", True, "-"),
    "0": ("0", True, "-"),
    "2": ("2", False, "unknown_forward_reference"),
    ".": ("2", False, "unknown_forward_reference"),
    "": ("2", False, "unknown_forward_reference"),
    "NA": ("2", False, "unknown_forward_reference"),
}

RESERVED_COLUMNS = [
    "_input_index",
    "_chromosome",
    "_position",
    "_strand_code",
    "_strand_provided",
    "_strand_label",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Validate a headerless genomic known-site table and create an "
            "UNCOMPRESSED three-column REDItoolKnown site list plus metadata."
        )
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--redi-output", required=True)
    parser.add_argument("--metadata-output", required=True)
    parser.add_argument("--fai", required=True)
    return parser.parse_args()


def read_fai(fai_path):
    if not os.path.isfile(fai_path) or os.path.getsize(fai_path) == 0:
        raise SystemExit(
            f"ERROR: Genome FASTA index is missing or empty: {fai_path}"
        )

    contig_lengths = {}
    with open(fai_path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 2:
                raise SystemExit(
                    f"ERROR: Malformed FASTA index row at {fai_path}:{line_number}."
                )
            try:
                length = int(fields[1])
            except ValueError:
                raise SystemExit(
                    f"ERROR: Invalid contig length at {fai_path}:{line_number}."
                )
            if length < 1:
                raise SystemExit(
                    f"ERROR: Non-positive contig length at {fai_path}:{line_number}."
                )
            contig_lengths[fields[0]] = length

    if not contig_lengths:
        raise SystemExit(
            f"ERROR: No contig information found in FASTA index: {fai_path}"
        )

    return contig_lengths


def read_headerless_sites(input_path):
    if not os.path.isfile(input_path) or os.path.getsize(input_path) == 0:
        raise SystemExit(
            f"ERROR: Input site table is missing or empty: {input_path}"
        )

    rows = []
    max_columns = 0

    with open(input_path, encoding="utf-8-sig", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.rstrip("\r\n")
            if not stripped.strip():
                continue
            if stripped.startswith("##"):
                continue

            fields = stripped.split("\t")
            if len(fields) < 2:
                raise SystemExit(
                    f"ERROR: Expected at least 2 tab-separated columns at line {line_number}; "
                    "the input must be headerless."
                )

            rows.append((line_number, fields))
            max_columns = max(max_columns, len(fields))

    if not rows:
        raise SystemExit("ERROR: No site rows found in the headerless input table.")

    original_header = ["chromosome", "position"]
    if max_columns >= 3:
        original_header.append("strand")

    for column_index in range(len(original_header) + 1, max_columns + 1):
        original_header.append(f"column_{column_index}")

    return rows, original_header


def remove_stale_compressed_outputs(redi_output):
    for suffix in (".gz", ".gz.tbi", ".tbi", ".csi"):
        stale_path = Path(str(redi_output) + suffix)
        if stale_path.exists() or stale_path.is_symlink():
            stale_path.unlink()


def assert_plain_text(path):
    with open(path, "rb") as handle:
        magic = handle.read(2)

    if magic == b"\x1f\x8b":
        raise RuntimeError(
            f"Expected an uncompressed TSV, but gzip magic bytes were found: {path}"
        )


def main():
    args = parse_args()

    redi_output = Path(args.redi_output).expanduser().absolute()
    metadata_output = Path(args.metadata_output).expanduser().absolute()

    if str(redi_output).lower().endswith(".gz"):
        raise SystemExit(
            "ERROR: --redi-output must be an uncompressed file ending in .tab, not .gz."
        )

    if redi_output == metadata_output:
        raise SystemExit(
            "ERROR: --redi-output and --metadata-output must be different files."
        )

    contig_lengths = read_fai(args.fai)
    data_rows, original_header = read_headerless_sites(args.input)

    for reserved in RESERVED_COLUMNS:
        if reserved in original_header:
            raise SystemExit(
                f"ERROR: Reserved metadata column name encountered: {reserved}"
            )

    metadata_header = original_header + RESERVED_COLUMNS
    normalized_sites = {}
    metadata_rows = []

    for input_index, (line_number, fields) in enumerate(data_rows, start=1):

        if len(fields) < len(original_header):
            fields = fields + [""] * (len(original_header) - len(fields))
        elif len(fields) > len(original_header):
            raise SystemExit(
                f"ERROR: Found {len(fields)} columns but expected at most "
                f"{len(original_header)} at line {line_number}."
            )

        chrom = fields[0].strip()
        pos_text = fields[1].strip()

        if not chrom:
            raise SystemExit(
                f"ERROR: Empty chromosome at line {line_number}."
            )

        if not re.fullmatch(r"[0-9]+", pos_text) or int(pos_text) < 1:
            raise SystemExit(
                f"ERROR: Invalid 1-based genomic position '{pos_text}' at line {line_number}. "
                "The input must not contain a header."
            )

        pos = int(pos_text)

        if chrom not in contig_lengths:
            raise SystemExit(
                f"ERROR: Chromosome '{chrom}' at line {line_number} "
                "is absent from the genome FASTA index."
            )

        if pos > contig_lengths[chrom]:
            raise SystemExit(
                f"ERROR: Position {chrom}:{pos} exceeds contig length "
                f"{contig_lengths[chrom]}."
            )

        strand_text = fields[2].strip().upper() if len(fields) >= 3 else ""

        if strand_text not in STRAND_MAP:
            raise SystemExit(
                f"ERROR: Third column must contain strand (+, -, 1, 0, 2, NA) "
                f"or be empty; found '{fields[2]}' at line {line_number}."
            )

        strand_code, strand_provided, strand_label = STRAND_MAP[strand_text]

        key = (chrom, pos)

        if key in normalized_sites and normalized_sites[key] != strand_code:
            raise SystemExit(
                f"ERROR: Conflicting strand values for {chrom}:{pos} "
                "in the input site table."
            )

        normalized_sites[key] = strand_code

        metadata_rows.append(
            fields + [
                str(input_index),
                chrom,
                str(pos),
                strand_code,
                "1" if strand_provided else "0",
                strand_label,
            ]
        )

    redi_output.parent.mkdir(parents=True, exist_ok=True)
    metadata_output.parent.mkdir(parents=True, exist_ok=True)

    redi_tmp = Path(str(redi_output) + f".tmp.{os.getpid()}")
    metadata_tmp = Path(str(metadata_output) + f".tmp.{os.getpid()}")

    try:
        with open(redi_tmp, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")

            for (chrom, pos), strand_code in sorted(
                normalized_sites.items(),
                key=lambda item: (item[0][0], item[0][1])
            ):
                writer.writerow([chrom, pos, strand_code])

        with open(metadata_tmp, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(metadata_header)
            writer.writerows(metadata_rows)

        assert_plain_text(redi_tmp)
        assert_plain_text(metadata_tmp)

        os.replace(redi_tmp, redi_output)
        os.replace(metadata_tmp, metadata_output)

        remove_stale_compressed_outputs(redi_output)

    finally:
        for tmp_path in (redi_tmp, metadata_tmp):
            if tmp_path.exists():
                tmp_path.unlink()

    if not redi_output.is_file() or redi_output.stat().st_size == 0:
        raise RuntimeError(
            f"REDItoolKnown site list was not created: {redi_output}"
        )

    if not metadata_output.is_file() or metadata_output.stat().st_size == 0:
        raise RuntimeError(
            f"Site metadata was not created: {metadata_output}"
        )

    assert_plain_text(redi_output)

    n_unknown = sum(row[-2] == "0" for row in metadata_rows)

    print(f"Input rows: {len(metadata_rows)}")
    print(f"Unique genomic sites: {len(normalized_sites)}")
    print(f"Rows without an explicit strand: {n_unknown}")
    print(f"REDItoolKnown site list (UNCOMPRESSED): {redi_output}")
    print(f"Site metadata: {metadata_output}")
    print("Compression: disabled")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
