#!/usr/bin/env python3
"""Normalize VCF records, extract annotated genic/intergenic counts, and merge results."""
import argparse
import csv
import gzip
import os
import re
import sys
from urllib.parse import unquote

BASES = set("ACGT")
GENIC_C_RE = re.compile(r"^c\.\d")
SUBSTITUTION_RE = re.compile(r"([ACGT]>[ACGT])$")
PROTEIN_RE = re.compile(r"^p\.\(?([A-Za-z*]{1,5})(\d+)([A-Za-z*=]{1,5})\)?$")
CORE_COLUMNS = ["chr", "pos", "ref", "alt"]
ANNOTATION_COLUMNS = ["anno", "gene_id", "transcript", "cds_pos", "mut", "AA_pos", "AA_change", "region"]


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, "r", encoding="utf-8")


def samples(value):
    return [x.strip() for x in value.split(",") if x.strip()]


def safe_output_file(path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    return open(path, "w", encoding="utf-8", newline="")


def split_alleles(field):
    if field in {".", ""}:
        return [0, 0]
    arr = field.split(",")
    if len(arr) < 2:
        raise ValueError(f"Expected reference and ALT counts: {field}")
    try:
        vals = [int(arr[0]), int(arr[1])]
    except ValueError as err:
        raise ValueError(f"Invalid allelic counts: {field}") from err
    if min(vals) < 0:
        raise ValueError(f"Negative allelic counts: {field}")
    return vals


def format_counts(fmt, call, line_no):
    keys = fmt.split(":")
    if "ADF" not in keys or "ADR" not in keys:
        raise ValueError(f"VCF line {line_no}: missing FORMAT/ADF or FORMAT/ADR; FORMAT={fmt}")
    values = call.split(":")
    if len(values) < len(keys):
        values += ["."] * (len(keys) - len(values))
    return split_alleles(values[keys.index("ADF")]), split_alleles(values[keys.index("ADR")])


def accept_alt(ref, alt, accept_plain):
    if len(ref) != 1 or ref not in BASES:
        return None
    if len(alt) == 1 and alt in BASES and alt != ref:
        return alt if accept_plain else None
    parts = alt.split(",")
    if len(parts) == 2 and len(parts[0]) == 1 and parts[0] in BASES and parts[0] != ref and re.fullmatch(r"<[^,<>]+>", parts[1]):
        return parts[0]
    return None


def normalize(args):
    expected = samples(args.dna) + samples(args.rna)
    if len(set(expected)) != len(expected):
        raise ValueError("Duplicate expected BAM/VCF sample IDs")
    seen_header = False
    kept = read = 0
    with open_text(args.input) as handle, safe_output_file(args.output) as out:
        for line_no, raw in enumerate(handle, 1):
            if raw.startswith("##"):
                out.write(raw)
                continue
            if raw.startswith("#CHROM"):
                fields = raw.rstrip("\n").split("\t")
                observed = fields[9:]
                if observed != expected:
                    raise ValueError(f"VCF sample order mismatch: expected={expected}, observed={observed}")
                out.write(raw)
                seen_header = True
                continue
            if raw.startswith("#"):
                out.write(raw)
                continue
            if not seen_header:
                raise ValueError("VCF #CHROM header not found before records")
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 9 + len(expected):
                raise ValueError(f"VCF line {line_no}: expected {9 + len(expected)} columns, got {len(fields)}")
            read += 1
            alt = accept_alt(fields[3].upper(), fields[4].upper(), bool(args.accept_plain))
            if alt is None:
                continue
            calls = []
            for call in fields[9:]:
                adf, adr = format_counts(fields[8], call, line_no)
                calls.append(f"{adf[0]},{adf[1]}:{adr[0]},{adr[1]}")
            fields[4], fields[7], fields[8] = alt, ".", "ADF:ADR"
            out.write("\t".join(fields[:9] + calls) + "\n")
            kept += 1
        if not seen_header:
            raise ValueError("Input VCF has no #CHROM header")
    print(f"Normalized variants: {kept}/{read}; samples: {len(expected)}")


def parse_ann(info):
    for item in info.split(";"):
        if item.startswith("ANN="):
            return [entry.split("|") for entry in item[4:].split(",") if entry]
    return []


def ann_to_cols(a, ref, alt, region):
    a = (a + [""] * 16)[:16]
    hgvs_c, hgvs_p = unquote(a[9]), unquote(a[10])
    mut = f"{ref}>{alt}"
    cds_pos = "NA"
    if region == "genic" and hgvs_c.startswith(("c.", "n.")):
        seq = hgvs_c[2:]
        match = SUBSTITUTION_RE.search(seq)
        if match:
            cds_pos = seq[:match.start()] or "NA"
            mut = match.group(1)
    aa_pos, aa_change = "0", "NNN>NNN"
    protein = PROTEIN_RE.match(hgvs_p)
    if protein:
        aa_change = f"{protein.group(1)}>{protein.group(3)}"
        aa_pos = protein.group(2)
    return [a[1] or "NA", unquote(a[4]).replace(" ", "_") or "NA", unquote(a[6]).replace(" ", "_") or "NA", cds_pos, mut, aa_pos, aa_change, region]


def extract(args):
    """Extract all annotated SNV counts; all coverage/FDR filtering is done downstream per library."""
    dna, rna = samples(args.dna), samples(args.rna)
    expected = dna + rna
    if not dna or not rna:
        raise ValueError("At least one DNA and one RNA library are required")
    if len(set(expected)) != len(expected):
        raise ValueError("Duplicate DNA/RNA sample IDs")
    header = CORE_COLUMNS + [f"{tag}_{sid}" for sid in expected for tag in ("R", "A")] + ANNOTATION_COLUMNS
    tally = {"all_variants": 0, "genic": 0, "intergenic": 0, "unannotated_or_other": 0}
    saw_ann = False
    with open_text(args.input) as handle, safe_output_file(args.genic) as gh, safe_output_file(args.intergenic) as ih:
        genic_writer = csv.writer(gh, delimiter="\t")
        inter_writer = csv.writer(ih, delimiter="\t")
        genic_writer.writerow(header)
        inter_writer.writerow(header)
        got_header = False
        for line_no, raw in enumerate(handle, 1):
            if raw.startswith("#CHROM"):
                observed = raw.rstrip("\r\n").split("\t")[9:]
                if observed != expected:
                    raise ValueError(f"Annotated VCF sample order mismatch: expected={expected}, observed={observed}")
                got_header = True
                continue
            if raw.startswith("#"):
                continue
            if not got_header:
                raise ValueError("Missing annotated VCF header")
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != len(expected) + 9:
                raise ValueError(f"Annotated VCF line {line_no}: unexpected number of samples")
            tally["all_variants"] += 1
            counts = []
            for call in fields[9:]:
                adf, adr = format_counts(fields[8], call, line_no)
                ref_count = adf[0] + adr[0]
                alt_count = adf[1] + adr[1]
                counts.extend((ref_count, alt_count))
            annotations = parse_ann(fields[7])
            if annotations:
                saw_ann = True
            if args.genic_mode != "c_numeric":
                raise ValueError(f"Unsupported genic_mode: {args.genic_mode}")
            genic = next((ann for ann in annotations if len(ann) > 9 and GENIC_C_RE.match(unquote(ann[9]))), None)
            if genic is not None:
                genic_writer.writerow(fields[:2] + fields[3:5] + counts + ann_to_cols(genic, fields[3], fields[4], "genic"))
                tally["genic"] += 1
                continue
            intergenic = next((ann for ann in annotations if len(ann) > 1 and "intergenic_region" in ann[1].split("&")), None)
            if intergenic is not None:
                inter_writer.writerow(fields[:2] + fields[3:5] + counts + ann_to_cols(intergenic, fields[3], fields[4], "intergenic"))
                tally["intergenic"] += 1
            else:
                tally["unannotated_or_other"] += 1
        if not got_header:
            raise ValueError("No annotated VCF header found")
        if tally["all_variants"] and not saw_ann:
            raise ValueError("No ANN annotations observed in annotated VCF")
    print("\n".join(f"{key}: {value}" for key, value in tally.items()))
    print("Coverage and BH-FDR filtering are deferred to 02_filter_plot.R and are performed independently per sample library.")


def combine(args):
    selected_cols = [
        "chr", "pos", "AA_change", "transcript", "DNA_pass_sample_names",
        "RNA_pass_samples", "n_RNA_pass", "anno_group", "mut", "anno",
        "edit_level_pass_samples", "region"
    ]
    records = []
    for filename in (args.genic, args.intergenic):
        with open(filename, encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = set(selected_cols) - set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"Missing columns in {filename}: {sorted(missing)}")
            for row in reader:
                records.append({col: row.get(col, "") for col in selected_cols})
    records.sort(key=lambda row: (row["chr"], int(row["pos"]), row["transcript"]))
    with safe_output_file(args.output) as out:
        writer = csv.DictWriter(out, fieldnames=selected_cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    print(f"Combined selected sites: {len(records)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="operation", required=True)
    norm = sub.add_parser("normalize", help="Retain biallelic SNVs and retain ADF/ADR counts")
    for key in ("input", "output", "dna", "rna"):
        norm.add_argument(f"--{key}", required=True)
    norm.add_argument("--accept-plain", type=int, default=1, choices=(0, 1))
    norm.set_defaults(func=normalize)
    extract_parser = sub.add_parser("extract", help="Extract annotated SNV counts without applying coverage filters")
    for key in ("input", "genic", "intergenic", "dna", "rna"):
        extract_parser.add_argument(f"--{key}", required=True)
    extract_parser.add_argument("--genic-mode", default="c_numeric", choices=("c_numeric",))
    extract_parser.set_defaults(func=extract)
    merge = sub.add_parser("combine", help="Merge selected genic and intergenic sites")
    for key in ("genic", "intergenic", "output"):
        merge.add_argument(f"--{key}", required=True)
    merge.set_defaults(func=combine)
    args = parser.parse_args()
    try:
        args.func(args)
    except (OSError, ValueError, csv.Error) as exc:
        sys.exit(f"ERROR: {exc}")


if __name__ == "__main__":
    main()
