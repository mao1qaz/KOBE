
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import re
from pathlib import Path
from collections import defaultdict, Counter


MISSING = {
    "", "-", "---", "NA", "NAN", "NONE", "NULL", "."
}

NOT_FOUND = {
    "NOT_FOUND", "NOT FOUND"
}


def parse_args():

    parser = argparse.ArgumentParser(
        prog="add_cds_sites_genomic_coordinates",
        description="Map homologous CDS editing sites to genomic coordinates using GTF.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n  KOBE add_cds_sites_genomic_coordinates \\\n    --input /project/data/merged_edit_sites_codon_table_annotated.csv \\\n    --gtf invicta=/project/ref/invicta.gtf \\\n    --gtf mellifera=/project/ref/mellifera.gtf \\\n    --output /project/results/merged_sites_with_genomic_coordinates.tsv \\\n    --sif /project/rna_editing.sif --bind /project",
    )

    parser.add_argument(
        "--input",
        required=True,
        metavar="FILE",
        help="Input annotated CSV file."
    )

    parser.add_argument(
        "--gtf",
        action="append",
        required=True,
        metavar="STRING",
        help="Species name and GTF path. Can be specified multiple times."
    )

    parser.add_argument(
        "--output",
        required=True,
        metavar="FILE",
        help="Output TSV path."
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


def parse_gtf_arguments(gtf_args):

    result = {}

    for item in gtf_args:

        if "=" in item:
            species, path = item.split("=", 1)
        elif ":" in item:
            # Backward compatibility with earlier SPECIES:GTF syntax.
            species, path = item.split(":", 1)
        else:
            raise ValueError(
                f"Invalid GTF argument: {item}. Expected SPECIES=GTF."
            )

        species = species.strip().lower()
        path = Path(path.strip())

        if not species:
            raise ValueError(f"Empty species name: {item}")

        if species in result:
            raise ValueError(f"Duplicate species: {species}")

        if not path.is_file():
            raise FileNotFoundError(f"GTF not found: {path}")

        result[species] = path

    return result


def extract_transcript_id(attributes):

    match = re.search(
        r'(?:^|;\s*)transcript_id\s+"([^"]+)"',
        attributes
    )

    if match:
        return match.group(1)

    match = re.search(
        r'(?:^|[;\s])transcript_id=([^;]+)',
        attributes
    )

    if match:
        return match.group(1).strip().strip('"')

    return None


def load_gtf(gtf_file):

    transcript_segments = defaultdict(list)

    with open(gtf_file, "r", encoding="utf-8") as handle:

        for line in handle:

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            if fields[2] != "CDS":
                continue

            chrom = fields[0]

            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue

            strand = fields[6]
            transcript_id = extract_transcript_id(fields[8])

            if transcript_id is None:
                continue

            if strand not in {"+", "-"}:
                continue

            if start < 1 or end < start:
                continue

            transcript_segments[transcript_id].append(
                (chrom, start, end, strand)
            )

    models = {}
    invalid_models = 0

    for transcript_id, segments in transcript_segments.items():

        segments = list(dict.fromkeys(segments))

        chroms = {x[0] for x in segments}
        strands = {x[3] for x in segments}

        if len(chroms) != 1 or len(strands) != 1:
            models[transcript_id] = None
            invalid_models += 1
            continue

        strand = next(iter(strands))

        if strand == "+":
            segments.sort(key=lambda x: (x[1], x[2]))
        else:
            segments.sort(
                key=lambda x: (x[1], x[2]),
                reverse=True
            )

        cds_segments = []
        cumulative_length = 0

        for chrom, start, end, strand in segments:

            segment_length = end - start + 1

            cds_start = cumulative_length + 1
            cds_end = cumulative_length + segment_length

            cds_segments.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "cds_start": cds_start,
                "cds_end": cds_end
            })

            cumulative_length += segment_length

        models[transcript_id] = {
            "chrom": next(iter(chroms)),
            "strand": strand,
            "segments": cds_segments,
            "cds_length": cumulative_length
        }

    return models, invalid_models


def parse_cds_entry(value):

    """
    Supported examples:

    Sinv005371.1:8
    XM_016911957.2:35
    FBtr0079946:---

    Not applicable examples:

    ---
    -
    FBtr0079946:---
    NOT_FOUND
    NOT_FOUND:---
    NOT_FOUND:35
    """

    value = str(value or "").strip()

    upper_value = value.upper()

    # Entire field is absent.
    if upper_value in MISSING:
        return None, None, "not_applicable_no_position"

    # Entire field explicitly indicates no transcript.
    if upper_value in NOT_FOUND:
        return None, None, "not_applicable_transcript_not_found"

    if ":" not in value:
        return None, None, "invalid_cds_format"

    transcript_id, position = value.rsplit(":", 1)

    transcript_id = transcript_id.strip()
    position = position.strip()

    upper_transcript = transcript_id.upper()
    upper_position = position.upper()

    # Explicitly missing transcript.
    if upper_transcript in NOT_FOUND:
        return None, None, "not_applicable_transcript_not_found"

    # Missing transcript ID.
    if upper_transcript in MISSING:
        return None, None, "not_applicable_transcript_not_found"

    # Transcript exists in the original table, but no CDS position.
    if upper_position in MISSING:
        return None, None, "not_applicable_no_position"

    # Explicitly unavailable CDS position.
    if upper_position in NOT_FOUND:
        return None, None, "not_applicable_no_position"

    try:
        cds_pos = int(position)
    except ValueError:
        return transcript_id, None, "invalid_cds_position"

    if cds_pos < 1:
        return transcript_id, None, "invalid_cds_position"

    return transcript_id, cds_pos, "valid"


def map_cds_to_genome(transcript_id, cds_pos, models):

    if transcript_id not in models:
        return None, "transcript_not_found"

    model = models[transcript_id]

    if model is None:
        return None, "invalid_gtf_model"

    if cds_pos > model["cds_length"]:
        return None, "cds_position_out_of_range"

    for segment in model["segments"]:

        if not (
            segment["cds_start"] <= cds_pos <= segment["cds_end"]
        ):
            continue

        offset = cds_pos - segment["cds_start"]

        if segment["strand"] == "+":
            genome_pos = segment["start"] + offset
        else:
            genome_pos = segment["end"] - offset

        coordinate = (
            f'{segment["chrom"]}:'
            f'{genome_pos}:'
            f'{segment["strand"]}'
        )

        return coordinate, "mapped"

    return None, "mapping_failed"


def main():

    args = parse_args()

    input_file = Path(args.input)

    if not input_file.is_file():
        raise FileNotFoundError(
            f"Input CSV not found: {input_file}"
        )

    gtf_files = parse_gtf_arguments(args.gtf)

    output_file = Path(args.output)

    if output_file.resolve() == input_file.resolve():
        raise ValueError("Output file must differ from input CSV.")

    print("========== INPUT ==========")
    print(f"Input CSV: {input_file}")
    print(f"Output TSV: {output_file}")

    with open(
        input_file,
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as handle:

        reader = csv.DictReader(handle)

        if not reader.fieldnames:
            raise ValueError("Input CSV has no header.")

        fieldnames = [x.strip() for x in reader.fieldnames]

        rows = []

        for raw_row in reader:

            row = {
                (key.strip() if key else key): value
                for key, value in raw_row.items()
                if key is not None
            }

            rows.append(row)

    print(f"Original rows: {len(rows)}")

    # Validate required CDS columns.
    for species in gtf_files:

        cds_column = f"{species}_cds_pos"

        if cds_column not in fieldnames:
            raise ValueError(
                f"Missing required column: {cds_column}"
            )

    # Load GTF models.
    print("\n========== LOADING GTF ==========")

    models_by_species = {}

    for species, gtf_file in gtf_files.items():

        print(f"\nSpecies: {species}")
        print(f"GTF: {gtf_file}")

        models, invalid = load_gtf(gtf_file)

        models_by_species[species] = models

        print(f"CDS transcripts: {len(models)}")
        print(f"Invalid transcript models: {invalid}")

    # Map all species.
    print("\n========== MAPPING CDS SITES ==========")

    for species, models in models_by_species.items():

        cds_column = f"{species}_cds_pos"
        genomic_column = f"{species}_genome_pos"

        if genomic_column not in fieldnames:
            fieldnames.append(genomic_column)

        stats = Counter()

        for row in rows:

            value = row.get(cds_column, "")

            transcript_id, cds_pos, parse_status = parse_cds_entry(
                value
            )

            if parse_status != "valid":

                row[genomic_column] = "NA"
                stats[parse_status] += 1
                continue

            coordinate, status = map_cds_to_genome(
                transcript_id,
                cds_pos,
                models
            )

            if status == "mapped":
                row[genomic_column] = coordinate
            else:
                row[genomic_column] = "NA"

            stats[status] += 1

        print(f"\nSpecies: {species}")
        print(f"CDS column: {cds_column}")
        print(f"Output column: {genomic_column}")

        # Explicitly report normal missing entries separately.
        print(
            f"  mapped: {stats['mapped']}"
        )

        print(
            f"  not_applicable_no_position: "
            f"{stats['not_applicable_no_position']}"
        )

        print(
            f"  not_applicable_transcript_not_found: "
            f"{stats['not_applicable_transcript_not_found']}"
        )

        # Actual mapping problems.
        for status in [
            "transcript_not_found",
            "invalid_gtf_model",
            "cds_position_out_of_range",
            "invalid_cds_position",
            "invalid_cds_format",
            "mapping_failed"
        ]:

            print(f"  {status}: {stats[status]}")

        print(f"  Total rows: {sum(stats.values())}")

    # Write output.
    output_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with open(
        output_file,
        "w",
        encoding="utf-8",
        newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n"
        )

        writer.writeheader()
        writer.writerows(rows)

    print("\n========== FINISHED ==========")
    print(f"Original rows preserved: {len(rows)}")
    print(f"Species processed: {len(gtf_files)}")
    print(f"Output: {output_file}")


if __name__ == "__main__":
    main()
