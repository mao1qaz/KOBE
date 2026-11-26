#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

parser = argparse.ArgumentParser(
    description="Export conserved editing sites with CDS and optional genome coordinates."
)
parser.add_argument("--out_dir", required=True)
parser.add_argument("--min_species", type=int, required=True)
parser.add_argument("--species_order", required=True)
parser.add_argument("--site_file", required=True)
parser.add_argument("--annotated", required=True)
parser.add_argument(
    "--gtf", action="append", default=[], metavar="SPECIES=GTF_PATH",
    help="Explicit GTF path for a species. Repeat for each species to be mapped; no default GTF paths are used."
)
args = parser.parse_args()

if args.min_species < 1:
    parser.error("--min_species must be >= 1")

out_dir = Path(args.out_dir).resolve()
species_order = [x.strip() for x in args.species_order.split(",") if x.strip()]

if not species_order or len(species_order) != len(set(species_order)):
    parser.error("species_order must contain unique species names")

site_file = Path(args.site_file)
anno_file = Path(args.annotated)
out_file = out_dir / f"conserved_sites_ge{args.min_species}species_full_coordinates.tsv"
qc_file = out_dir / f"conserved_sites_ge{args.min_species}species_genome_mapping_QC.tsv"
out_dir.mkdir(parents=True, exist_ok=True)

GENETIC_CODE = {
    "TTT":"Phe","TTC":"Phe","TTA":"Leu","TTG":"Leu",
    "TCT":"Ser","TCC":"Ser","TCA":"Ser","TCG":"Ser",
    "TAT":"Tyr","TAC":"Tyr","TAA":"Stop","TAG":"Stop",
    "TGT":"Cys","TGC":"Cys","TGA":"Stop","TGG":"Trp",
    "CTT":"Leu","CTC":"Leu","CTA":"Leu","CTG":"Leu",
    "CCT":"Pro","CCC":"Pro","CCA":"Pro","CCG":"Pro",
    "CAT":"His","CAC":"His","CAA":"Gln","CAG":"Gln",
    "CGT":"Arg","CGC":"Arg","CGA":"Arg","CGG":"Arg",
    "ATT":"Ile","ATC":"Ile","ATA":"Ile","ATG":"Met",
    "ACT":"Thr","ACC":"Thr","ACA":"Thr","ACG":"Thr",
    "AAT":"Asn","AAC":"Asn","AAA":"Lys","AAG":"Lys",
    "AGT":"Ser","AGC":"Ser","AGA":"Arg","AGG":"Arg",
    "GTT":"Val","GTC":"Val","GTA":"Val","GTG":"Val",
    "GCT":"Ala","GCC":"Ala","GCA":"Ala","GCG":"Ala",
    "GAT":"Asp","GAC":"Asp","GAA":"Glu","GAG":"Glu",
    "GGT":"Gly","GGC":"Gly","GGA":"Gly","GGG":"Gly"
}

def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8")

def clean(value):
    value = str(value).strip()
    return "NA" if value.upper() in {"", "NA", "NAN", "NONE", "-", "---", "NULL"} else value

def strip_version(tx):
    return re.sub(r"\.\d+$", "", tx)

def parse_cds(value):
    value = str(value).strip()
    if value.upper().startswith("NOT_FOUND"):
        return None, None, "NO_ORTHOLOG"
    if not value or value.upper() in {"NA", "NAN", "-", "---"} or ":" not in value:
        return None, None, "no_valid_CDS_coordinate"

    tx, pos_text = value.rsplit(":", 1)
    if not tx:
        return None, None, "invalid_CDS_position"
    if pos_text.strip() in {"---", "-"}:
        return tx, None, "GAP"

    try:
        pos = int(pos_text)
    except ValueError:
        return None, None, "invalid_CDS_position"
    if pos < 1:
        return None, None, "invalid_CDS_position"
    return tx, pos, "PASS"

def parse_phase(value):
    try:
        phase = int(float(str(value)))
    except (ValueError, TypeError):
        return "NA"
    return str(phase) if phase in (1, 2, 3) else "NA"

def translate(codon):
    codon = clean(codon).upper().replace("U", "T")
    if codon == "NA" or len(codon) != 3:
        return "NA"
    return GENETIC_CODE.get(codon, "NA")

def parse_gtf_attributes(text):
    attrs = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item and " " not in item.split("=", 1)[0]:
            key, value = item.split("=", 1)
        else:
            parts = item.split(None, 1)
            if len(parts) != 2:
                continue
            key, value = parts
        attrs[key.strip()] = value.strip().strip('"')
    return attrs

def load_gtf(gtf_file):
    gtf_file = Path(gtf_file)
    if not gtf_file.is_file():
        raise FileNotFoundError(f"GTF not found: {gtf_file}")

    raw = defaultdict(set)
    with open_text(gtf_file) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = parse_gtf_attributes(fields[8])
            tx_value = attrs.get("transcript_id") or attrs.get("Parent")
            if not tx_value:
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise ValueError(f"Invalid GTF coordinate at line {line_no}") from exc
            chrom, strand = fields[0], fields[6]
            if start < 1 or end < start:
                raise ValueError(f"Invalid CDS interval at line {line_no}")
            for tx in tx_value.split(","):
                tx = tx.strip()
                if tx.startswith("transcript:"):
                    tx = tx.split(":", 1)[1]
                if tx:
                    raw[tx].add((chrom, start, end, strand))

    if not raw:
        raise ValueError(f"No CDS features with transcript IDs found in {gtf_file}")

    tx_info, bare_map = {}, defaultdict(list)
    for tx, intervals_set in raw.items():
        intervals = list(intervals_set)
        chroms, strands = {x[0] for x in intervals}, {x[3] for x in intervals}
        bare_map[strip_version(tx)].append(tx)

        if len(chroms) != 1 or len(strands) != 1:
            tx_info[tx] = {"error": "multiple_chromosomes_or_strands"}
            continue

        chrom, strand = next(iter(chroms)), next(iter(strands))
        if strand not in {"+", "-"}:
            tx_info[tx] = {"error": "invalid_strand"}
            continue

        genomic_order = sorted(intervals, key=lambda x: x[1])
        if any(genomic_order[i][1] <= genomic_order[i-1][2] for i in range(1, len(genomic_order))):
            tx_info[tx] = {"error": "overlapping_CDS_intervals"}
            continue

        ordered = sorted(intervals, key=lambda x: x[1], reverse=(strand == "-"))
        segments, cumulative = [], 0
        for seg_chrom, start, end, seg_strand in ordered:
            length = end - start + 1
            segments.append({
                "cds_start": cumulative + 1,
                "cds_end": cumulative + length,
                "chrom": seg_chrom, "start": start, "end": end,
                "strand": seg_strand
            })
            cumulative += length
        tx_info[tx] = {"segments": segments, "length": cumulative, "error": None}

    print(f"Loaded GTF: {gtf_file.name}; transcripts with CDS: {len(tx_info)}")
    return tx_info, bare_map

def map_cds_to_genome(tx, cds_pos, tx_info, bare_map):
    if tx in tx_info:
        resolved_tx = tx
    else:
        candidates = bare_map.get(strip_version(tx), [])
        if not candidates:
            return "NA", "transcript_not_found_in_GTF"
        if len(candidates) > 1:
            return "NA", "ambiguous_transcript_version"
        resolved_tx = candidates[0]

    info = tx_info[resolved_tx]
    if info.get("error"):
        return "NA", info["error"]
    if cds_pos > info["length"]:
        return "NA", "CDS_position_out_of_range"

    for seg in info["segments"]:
        if seg["cds_start"] <= cds_pos <= seg["cds_end"]:
            offset = cds_pos - seg["cds_start"]
            genomic_pos = seg["start"] + offset if seg["strand"] == "+" else seg["end"] - offset
            return f"{seg['chrom']}:{genomic_pos}:{seg['strand']}", "PASS"
    return "NA", "CDS_position_not_mapped"


# Only explicitly supplied GTF paths are used; unspecified species are skipped.
gtf_paths = {}
for item in args.gtf:
    if "=" not in item:
        parser.error(f"Invalid --gtf argument: {item}. Use species=/path/to/file.gtf")
    sp, path = item.split("=", 1)
    sp, path = sp.strip(), path.strip()
    if sp not in species_order:
        parser.error(f"Unknown species in --gtf: {sp}")
    if sp in gtf_paths:
        parser.error(f"Duplicated GTF argument for {sp}")
    gtf_paths[sp] = Path(path).expanduser().absolute()

gtf_databases = {}
for sp in species_order:
    if sp in gtf_paths:
        gtf_databases[sp] = load_gtf(gtf_paths[sp])
    else:
        gtf_databases[sp] = None
        print(f"GTF not provided for {sp}; genome mapping skipped")


# Read input tables.
if not site_file.is_file():
    raise FileNotFoundError(f"Filtered site table not found: {site_file}")
if not anno_file.is_file():
    raise FileNotFoundError(f"Annotated table not found: {anno_file}")

sites = pd.read_csv(site_file, sep="\t", dtype=str, keep_default_na=False)
try:
    anno = pd.read_csv(anno_file, dtype=str, keep_default_na=False, encoding="utf-8-sig")
except UnicodeDecodeError:
    anno = pd.read_csv(anno_file, dtype=str, keep_default_na=False, encoding="gbk")

required_site = ["ID2", "edited_species", "n_edited_species"]
required_anno = ["ID2"]
for sp in species_order:
    required_anno.extend([f"{sp}_cds_pos", f"{sp}_codon", f"{sp}_codon_pos"])

for table, required, name in [
    (sites, required_site, "Filtered site table"),
    (anno, required_anno, "Annotated table")
]:
    missing = [col for col in required if col not in table.columns]
    if missing:
        raise ValueError(f"{name} missing columns: " + ", ".join(missing))

sites["ID2"] = sites["ID2"].str.strip()
anno["ID2"] = anno["ID2"].str.strip()
if sites["ID2"].duplicated().any():
    raise ValueError("Duplicate ID2 in filtered site table")
if anno["ID2"].duplicated().any():
    raise ValueError("Duplicate ID2 in annotated table")

missing_ids = set(sites["ID2"]) - set(anno["ID2"])
if missing_ids:
    raise ValueError(f"{len(missing_ids)} filtered sites cannot be found in the annotated table")

anno_cols = ["ID2"] + required_anno[1:]
df = sites[required_site].merge(
    anno[anno_cols], on="ID2", how="left",
    validate="one_to_one", sort=False
)
df["n_edited_species"] = pd.to_numeric(df["n_edited_species"], errors="raise").astype(int)
if (df["n_edited_species"] < args.min_species).any():
    raise ValueError("Input contains sites below the requested minimum species count")

print("Input conserved sites:", len(df))
print("Species order:", ", ".join(species_order))


# Build one output row per conserved site.
output_rows, qc_rows = [], []
mapped_counts = defaultdict(int)
no_ortholog_counts = defaultdict(int)
gap_counts = defaultdict(int)

for _, row in df.iterrows():
    result = {"ID2": row["ID2"]}

    for sp in species_order:
        raw_cds = row[f"{sp}_cds_pos"]
        tx, cds_pos, parse_status = parse_cds(raw_cds)
        cds_coordinate = f"{tx}:{cds_pos}" if parse_status == "PASS" else "NA"
        genome_coordinate = "NA"

        # Missing orthologs and alignment gaps are expected, not mapping failures.
        if parse_status == "NO_ORTHOLOG":
            no_ortholog_counts[sp] += 1
        elif parse_status == "GAP":
            gap_counts[sp] += 1

        # Map only when a GTF was supplied for this species.
        if gtf_databases[sp] is not None:
            mapping_status = parse_status
            if parse_status == "PASS":
                tx_info, bare_map = gtf_databases[sp]
                genome_coordinate, mapping_status = map_cds_to_genome(
                    tx, cds_pos, tx_info, bare_map
                )
                if mapping_status == "PASS":
                    mapped_counts[sp] += 1
                else:
                    qc_rows.append({
                        "ID2": row["ID2"], "species": sp,
                        "transcript": tx, "CDS_position": cds_pos,
                        "GTF": str(gtf_paths[sp]), "reason": mapping_status
                    })
            elif parse_status not in {"NO_ORTHOLOG", "GAP"} and clean(raw_cds) != "NA":
                qc_rows.append({
                    "ID2": row["ID2"], "species": sp,
                    "transcript": "NA", "CDS_position": "NA",
                    "GTF": str(gtf_paths[sp]), "reason": parse_status
                })

        codon = clean(row[f"{sp}_codon"]).upper().replace("U", "T")
        if codon != "NA" and (len(codon) != 3 or any(base not in "ACGT" for base in codon)):
            codon = "NA"

        phase = parse_phase(row[f"{sp}_codon_pos"])
        pre_edit_aa = translate(codon)

        result[f"{sp}_CDS_coordinate"] = cds_coordinate
        result[f"{sp}_genome_coordinate"] = genome_coordinate
        result[f"{sp}_codon"] = codon
        result[f"{sp}_codon_phase"] = phase
        result[f"{sp}_pre_edit_AA"] = pre_edit_aa

    result["edited_species"] = row["edited_species"]
    result["n_edited_species"] = row["n_edited_species"]
    output_rows.append(result)


# Save output and mapping QC.
output_columns = ["ID2"]
for sp in species_order:
    output_columns.extend([
        f"{sp}_CDS_coordinate",
        f"{sp}_genome_coordinate",
        f"{sp}_codon",
        f"{sp}_codon_phase",
        f"{sp}_pre_edit_AA"
    ])
output_columns.extend(["edited_species", "n_edited_species"])

out_df = pd.DataFrame(output_rows, columns=output_columns)
out_df.to_csv(out_file, sep="\t", index=False, na_rep="NA")

qc_df = pd.DataFrame(
    qc_rows,
    columns=["ID2", "species", "transcript", "CDS_position", "GTF", "reason"]
)
qc_df.to_csv(qc_file, sep="\t", index=False)


# Report.
print("\n==============================================")
print("Conserved site coordinate extraction")
print("==============================================")
print("Conserved sites:", len(out_df))
print("Minimum edited species:", args.min_species)
print("\nGenome coordinate mapping:")

for sp in species_order:
    if sp not in gtf_paths:
        print(f"{sp}: SKIPPED (GTF not provided)")
    else:
        print(f"{sp}: {mapped_counts[sp]} / {len(out_df)} mapped")

print("\nExpected missing orthologs (excluded from mapping failures):")
for sp in species_order:
    print(f"{sp}: {no_ortholog_counts[sp]}")

print("\nAlignment gaps (excluded from mapping failures):")
for sp in species_order:
    print(f"{sp}: {gap_counts[sp]}")

print("\nMapping failures:", len(qc_df))
if not qc_df.empty:
    print("\nFailure reasons:")
    print(qc_df["reason"].value_counts().to_string())

print("\nOutput:")
print(out_file)
print("\nQC:")
print(qc_file)
print("\nDONE")
