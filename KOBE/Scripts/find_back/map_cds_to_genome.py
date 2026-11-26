#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Append GTF-derived genomic coordinates to existing find_back result tables."""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

COORD_COLUMNS = ("chrom", "genome_pos", "genome_strand")


def parse_args():
    p = argparse.ArgumentParser(
        description="Append chrom, genome_pos and genome_strand to existing find_back tables."
    )
    p.add_argument("--sample_table", required=True)
    p.add_argument("--summary_dir", required=True)
    p.add_argument("--statistics_dir", required=True)
    p.add_argument(
        "--gtf", action="append", default=[], metavar="SPECIES=FILE",
        help="Optional species-specific GTF; repeat for multiple species."
    )
    p.add_argument("--min_samples", type=int, default=1)
    p.add_argument("--species", nargs="+", default=None,
                   help="Process only the specified species, e.g. --species invicta.")
    p.add_argument("--out_dir", required=True,
                   help="Output root. Original input tables are not overwritten.")
    return p.parse_args()


def read_species(sample_table):
    species_list = []
    with open(sample_table, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Sample table has no header: {sample_table}")
        reader.fieldnames = [x.strip().lstrip("\ufeff") for x in reader.fieldnames]
        if "species" not in reader.fieldnames:
            raise ValueError("Sample table must contain a 'species' column.")
        for row in reader:
            species = (row.get("species") or "").strip()
            if not species:
                continue
            sample_type = (row.get("type") or "").strip().upper()
            if sample_type and sample_type != "RNA":
                continue
            if species not in species_list:
                species_list.append(species)
    if not species_list:
        raise ValueError(f"No RNA species found in sample table: {sample_table}")
    return species_list


def parse_gtf_attributes(text):
    attributes = {}
    for key, value in re.findall(r'([^;\s]+)\s+"([^"]*)"', text):
        attributes[key] = value
    if not attributes:
        for item in text.strip().strip(";").split(";"):
            item = item.strip()
            if "=" in item:
                key, value = item.split("=", 1)
                attributes[key.strip()] = value.strip().strip('"')
    transcript_id = attributes.get("transcript_id")
    if not transcript_id:
        parent = attributes.get("Parent")
        if parent:
            transcript_id = parent.split(",")[0].strip()
            transcript_id = re.sub(r"^(transcript:|rna-)", "", transcript_id)
    return transcript_id


def load_cds_models(gtf_path):
    raw_models = defaultdict(list)
    with open(gtf_path, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"Malformed GTF row at {gtf_path}:{line_number}; expected 9 columns."
                )
            if fields[2] != "CDS":
                continue
            chrom = fields[0]
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                raise ValueError(f"Invalid CDS coordinates at {gtf_path}:{line_number}.")
            strand = fields[6]
            if start < 1 or end < start or strand not in ("+", "-"):
                raise ValueError(f"Invalid CDS feature at {gtf_path}:{line_number}.")
            transcript_id = parse_gtf_attributes(fields[8])
            if not transcript_id:
                raise ValueError(
                    f"No transcript_id/Parent found at {gtf_path}:{line_number}."
                )
            raw_models[transcript_id].append((chrom, start, end, strand))

    models = {}
    for transcript_id, features in raw_models.items():
        chroms = {x[0] for x in features}
        strands = {x[3] for x in features}
        if len(chroms) != 1 or len(strands) != 1:
            raise ValueError(
                f"Transcript {transcript_id} has CDS features on multiple chromosomes or strands: {gtf_path}"
            )
        chrom, strand = next(iter(chroms)), next(iter(strands))
        intervals = sorted(set((x[1], x[2]) for x in features), reverse=(strand == "-"))
        segments = []
        cds_cursor = 1
        for start, end in intervals:
            seg_len = end - start + 1
            segments.append({
                "start": start,
                "end": end,
                "cds_start": cds_cursor,
                "cds_end": cds_cursor + seg_len - 1,
            })
            cds_cursor += seg_len
        models[transcript_id] = {
            "chrom": chrom,
            "strand": strand,
            "segments": segments,
        }
    return models


def map_cds_position(model, cds_pos):
    for seg in model["segments"]:
        if seg["cds_start"] <= cds_pos <= seg["cds_end"]:
            offset = cds_pos - seg["cds_start"]
            genome_pos = seg["start"] + offset if model["strand"] == "+" else seg["end"] - offset
            return model["chrom"], genome_pos, model["strand"]
    return None


def output_path_for(out_dir, category, species, input_path):
    filename = os.path.splitext(os.path.basename(input_path))[0] + ".genomic.tsv"
    return os.path.join(out_dir, category, species, filename)


def map_table(input_path, output_path, models, species):
    if not os.path.isfile(input_path):
        raise FileNotFoundError(f"Required result table not found: {input_path}")
    with open(input_path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Result table has no header: {input_path}")
        fieldnames = [x.strip().lstrip("\ufeff") for x in reader.fieldnames]
        missing = {"transcript", "cds_pos"} - set(fieldnames)
        if missing:
            raise ValueError(f"{input_path} missing required columns: {sorted(missing)}")
        rows = list(reader)

    output_fields = [x for x in fieldnames if x not in COORD_COLUMNS] + list(COORD_COLUMNS)
    mapped_count = 0
    for row in rows:
        row["chrom"] = row["genome_pos"] = row["genome_strand"] = "NA"
        transcript = (row.get("transcript") or "").strip()
        pos_text = (row.get("cds_pos") or "").strip()
        if not pos_text.isdigit() or int(pos_text) < 1:
            continue
        model = models.get(transcript)
        if model is None:
            continue
        result = map_cds_position(model, int(pos_text))
        if result is None:
            continue
        row["chrom"], genome_pos, row["genome_strand"] = result
        row["genome_pos"] = str(genome_pos)
        mapped_count += 1

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=output_fields, delimiter="\t",
            lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"{species}\t{os.path.basename(input_path)}\trows={len(rows)}"
        f"\tmapped={mapped_count}\toutput={output_path}"
    )


def main():
    args = parse_args()
    if args.min_samples < 1:
        raise ValueError("--min_samples must be >= 1")
    for path, label in (
        (args.sample_table, "sample table"),
        (args.summary_dir, "summary directory"),
        (args.statistics_dir, "statistics directory"),
    ):
        if not os.path.exists(path):
            raise FileNotFoundError(f"{label} not found: {path}")
    os.makedirs(args.out_dir, exist_ok=True)

    all_species = read_species(args.sample_table)
    species_list = args.species if args.species else all_species
    gtf_files = {}
    for value in args.gtf:
        if "=" not in value:
            raise ValueError(f"Invalid --gtf value {value!r}; expected SPECIES=FILE")
        species, path = value.split("=", 1)
        if not species or not path:
            raise ValueError(f"Invalid --gtf value {value!r}; expected SPECIES=FILE")
        if species in gtf_files:
            raise ValueError(f"GTF supplied more than once for species: {species}")
        if species not in all_species:
            raise ValueError(f"GTF species is not present in sample table: {species}")
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise FileNotFoundError(f"Species GTF missing or empty: {path}")
        gtf_files[species] = path
    unknown = sorted(set(species_list) - set(all_species))
    if unknown:
        raise ValueError(f"Species not present in sample table: {unknown}")

    print("============================================================")
    print("CDS-to-genome coordinate mapping")
    print("Added columns: chrom, genome_pos, genome_strand")
    print(f"Species: {', '.join(species_list)}")
    print(f"Output root: {args.out_dir}")
    print("============================================================")

    for species in species_list:
        gtf_path = gtf_files.get(species)
        if gtf_path is None:
            print(f"WARNING: No GTF supplied for {species}; skipping genomic-coordinate mapping.")
            continue
        models = load_cds_models(gtf_path)
        print(f"\nSpecies: {species}; GTF transcripts with CDS: {len(models)}")

        table_specs = (
            ("candidate_counts", os.path.join(
                args.summary_dir, species, f"{species}_all_samples_summary.tsv"
            )),
            ("statistics", os.path.join(
                args.statistics_dir, species,
                f"{species}_supported_sites_ge{args.min_samples}.tsv"
            )),
            ("statistics", os.path.join(
                args.statistics_dir, species,
                f"{species}_novel_supported_sites_ge{args.min_samples}.tsv"
            )),
        )
        for category, input_path in table_specs:
            output_path = output_path_for(args.out_dir, category, species, input_path)
            map_table(input_path, output_path, models, species)

    print("\nCoordinate mapping completed for supplied GTF species.")
    skipped = [sp for sp in species_list if sp not in gtf_files]
    if skipped:
        print("Skipped species without GTF: " + ", ".join(skipped))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
