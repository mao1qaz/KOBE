
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import re
from pathlib import Path
from collections import Counter


KEEP_ANNOTATIONS = {
    "synonymous_variant",
    "missense_variant"
}

OUTPUT_COLUMNS = [
    "transcript",
    "cds_pos",
    "anno",
    "AA_pos",
    "AA_change",
    "gene_id"
]


def parse_args():

    parser = argparse.ArgumentParser(
        prog="extract_cds_candidate_sites",
        description="Extract synonymous and missense CDS editing sites.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n  KOBE extract_cds_candidate_sites \\\n    --input /project/results/invicta/genic.final_sites.tsv \\\n    --output /project/results/orthology/invicta_candi_cds_sites.txt \\\n    --sif /project/rna_editing.sif --bind /project",
    )

    parser.add_argument(
        "--input",
        required=True,
        metavar="FILE",
        help="Input annotated editing sites TSV."
    )

    parser.add_argument(
        "--output",
        required=True,
        metavar="FILE",
        help="Output six-column CDS candidate sites TXT."
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

    return parser.parse_args()


def main():

    args = parse_args()

    input_file = Path(args.input).resolve()
    output_file = Path(args.output).resolve()

    if not input_file.is_file():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    if input_file == output_file:
        raise ValueError("Input and output files must be different.")

    output_file.parent.mkdir(parents=True, exist_ok=True)

    counts = Counter()
    total_rows = 0
    invalid_cds_pos = 0

    with open(
        input_file, "r", encoding="utf-8-sig", newline=""
    ) as fin:

        reader = csv.DictReader(fin, delimiter="\t")

        if not reader.fieldnames:
            raise ValueError("Input TSV has no header.")

        missing = set(OUTPUT_COLUMNS) - set(reader.fieldnames)

        if missing:
            raise ValueError(
                f"Missing required columns: {sorted(missing)}"
            )

        with open(
            output_file, "w", encoding="utf-8", newline=""
        ) as fout:

            writer = csv.DictWriter(
                fout,
                fieldnames=OUTPUT_COLUMNS,
                delimiter="\t",
                lineterminator="\n"
            )

            writer.writeheader()

            for row in reader:

                total_rows += 1

                anno = (row.get("anno") or "").strip()

                if anno not in KEEP_ANNOTATIONS:
                    counts["excluded_annotation"] += 1
                    continue

                cds_pos = (row.get("cds_pos") or "").strip()

                if not re.fullmatch(r"[1-9]\d*", cds_pos):
                    invalid_cds_pos += 1
                    continue

                writer.writerow({
                    col: (row.get(col) or "").strip()
                    for col in OUTPUT_COLUMNS
                })

                counts[anno] += 1

    total_output = (
        counts["synonymous_variant"]
        + counts["missense_variant"]
    )

    print("\n========== SUMMARY ==========")
    print(f"Input: {input_file}")
    print(f"Output: {output_file}")
    print(f"Total input rows: {total_rows}")
    print(f"Synonymous sites: {counts['synonymous_variant']}")
    print(f"Missense sites: {counts['missense_variant']}")
    print(f"Invalid CDS positions skipped: {invalid_cds_pos}")
    print(f"Excluded other annotations: {counts['excluded_annotation']}")
    print(f"Total output rows: {total_output}")
    print("\nFinished successfully.")


if __name__ == "__main__":
    main()
