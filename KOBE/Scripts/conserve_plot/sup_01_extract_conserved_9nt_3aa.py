#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import Counter
from pathlib import Path
import re

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Extract homologous CDS windows for supplementary recovered-site analysis.")
    p.add_argument("--anno_file", required=True, help="Updated annotation CSV containing new_id2_group")
    p.add_argument("--fasta_dir", required=True, help="Directory containing *_longestcds.fasta")
    p.add_argument("--out_dir", required=True, help="Supplementary output directory")
    p.add_argument("--source_priority", default=None, help="Comma-separated representative-species priority")
    return p.parse_args()


args = parse_args()
anno_file = Path(args.anno_file)
fasta_dir = Path(args.fasta_dir)
out_dir = Path(args.out_dir)
out_dir.mkdir(parents=True, exist_ok=True)

if not anno_file.is_file():
    raise FileNotFoundError(f"Updated annotation CSV not found: {anno_file}")
if not fasta_dir.is_dir():
    raise NotADirectoryError(f"CDS FASTA directory not found: {fasta_dir}")

GENETIC_CODE = {
    "TTT":"Phe", "TTC":"Phe", "TTA":"Leu", "TTG":"Leu",
    "TCT":"Ser", "TCC":"Ser", "TCA":"Ser", "TCG":"Ser",
    "TAT":"Tyr", "TAC":"Tyr", "TAA":"Stop", "TAG":"Stop",
    "TGT":"Cys", "TGC":"Cys", "TGA":"Stop", "TGG":"Trp",
    "CTT":"Leu", "CTC":"Leu", "CTA":"Leu", "CTG":"Leu",
    "CCT":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",
    "CAT":"His", "CAC":"His", "CAA":"Gln", "CAG":"Gln",
    "CGT":"Arg", "CGC":"Arg", "CGA":"Arg", "CGG":"Arg",
    "ATT":"Ile", "ATC":"Ile", "ATA":"Ile", "ATG":"Met",
    "ACT":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr",
    "AAT":"Asn", "AAC":"Asn", "AAA":"Lys", "AAG":"Lys",
    "AGT":"Ser", "AGC":"Ser", "AGA":"Arg", "AGG":"Arg",
    "GTT":"Val", "GTC":"Val", "GTA":"Val", "GTG":"Val",
    "GCT":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala",
    "GAT":"Asp", "GAC":"Asp", "GAA":"Glu", "GAG":"Glu",
    "GGT":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly"
}


def read_fasta(path):
    seqs = {}
    name, pieces = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    if name in seqs:
                        raise ValueError(f"Duplicate transcript ID in {path}: {name}")
                    seqs[name] = "".join(pieces).upper().replace("U", "T")
                name, pieces = line[1:].split()[0], []
            else:
                pieces.append(line)
        if name is not None:
            if name in seqs:
                raise ValueError(f"Duplicate transcript ID in {path}: {name}")
            seqs[name] = "".join(pieces).upper().replace("U", "T")
    return seqs


def parse_cds_pos(value):
    value = str(value).strip()
    if not value or value.upper() in {"-", "---", "NA", "NAN"} or value.startswith("NOT_FOUND") or ":" not in value:
        return None
    tx, pos = value.rsplit(":", 1)
    if not tx.strip() or not re.fullmatch(r"\d+", pos.strip()):
        return None
    return tx.strip(), int(pos)


def valid_phase(value):
    try:
        phase = int(float(str(value)))
    except (TypeError, ValueError):
        return None
    return phase if phase in (1, 2, 3) else None


def infer_missing_type(raw_pos, raw_phase, raw_codon):
    raw_pos = str(raw_pos).strip()
    raw_phase = str(raw_phase).strip()
    raw_codon = str(raw_codon).strip().upper()
    if raw_pos.startswith("NOT_FOUND"):
        return "NO_ORTHOLOG"
    if ":" in raw_pos and raw_pos.rsplit(":", 1)[1] == "---":
        return "GAP"
    if raw_phase in {"0", "0.0"} and raw_codon in {"---", "XXX"}:
        return "GAP"
    return "OTHER"


def translate_window(seq):
    return ";".join(GENETIC_CODE.get(seq[i:i + 3], "Xaa") for i in range(0, len(seq), 3))


def group_ids_from(value):
    parts = [x.strip() for x in str(value).split(";")]
    return list(dict.fromkeys(x for x in parts if x and x.upper() not in {"NA", "NAN", "UNFOUND", "-", "---"}))


def map_group_to_species(group_ids, row, species):
    group_set = set(group_ids)
    mapped_species = []
    recognized_ids = set()
    for sp in species:
        value = str(row.get(f"{sp}_cds_pos", "")).strip()
        if value and value in group_set:
            mapped_species.append(sp)
            recognized_ids.add(value)
    unmapped_ids = [x for x in group_ids if x not in recognized_ids]
    return mapped_species, unmapped_ids


def extract_species_window(row, sp, seqs):
    raw_pos = str(row.get(f"{sp}_cds_pos", "")).strip()
    raw_phase = str(row.get(f"{sp}_codon_pos", "")).strip()
    raw_codon = str(row.get(f"{sp}_codon", "")).strip().upper()
    parsed = parse_cds_pos(raw_pos)
    phase = valid_phase(raw_phase)

    if parsed is None or phase is None or raw_codon in {"", "-", "---", "XXX", "NA", "NAN"}:
        return None, infer_missing_type(raw_pos, raw_phase, raw_codon), "invalid_position_or_codon"
    if seqs is None:
        return None, "OTHER", "fasta_not_available"

    tx, cds_pos = parsed
    if tx not in seqs:
        return None, "OTHER", "transcript_not_in_fasta"

    seq = seqs[tx]
    codon_start = cds_pos - phase + 1
    codon_end = codon_start + 2
    if codon_start < 1 or codon_end > len(seq):
        return None, "OTHER", "central_codon_out_of_range"

    has_prev = codon_start >= 4
    has_next = codon_end + 3 <= len(seq)
    window_start = codon_start - 3 if has_prev else codon_start
    window_end = codon_end + 3 if has_next else codon_end
    seq_window = seq[window_start - 1:window_end]

    if len(seq_window) not in (3, 6, 9):
        return None, "OTHER", "invalid_window_length"

    central_offset = 3 if has_prev else 0
    fasta_codon = seq_window[central_offset:central_offset + 3]
    if fasta_codon != raw_codon:
        return None, "OTHER", f"codon_mismatch:{raw_codon}!={fasta_codon}"

    edit_base_pos = central_offset + phase
    result = {
        "transcript": tx,
        "transcript_pos": cds_pos,
        "codon_phase": phase,
        "CDS": seq_window,
        "pep": translate_window(seq_window),
        "edit_base_pos": edit_base_pos,
        "edit_base": seq_window[edit_base_pos - 1],
        "annotated_codon": raw_codon,
        "fasta_codon": fasta_codon
    }
    return result, "SEQ", "PASS"


try:
    anno = pd.read_csv(anno_file, dtype=str, keep_default_na=False, encoding="utf-8-sig")
except UnicodeDecodeError:
    anno = pd.read_csv(anno_file, dtype=str, keep_default_na=False, encoding="gbk")

required = ["ID2", "ID2_group", "new_id2_group", "recovered_only_group"]
missing = [x for x in required if x not in anno.columns]
if missing:
    raise ValueError("Updated annotation CSV is missing required columns: " + ", ".join(missing))

anno["ID2"] = anno["ID2"].str.strip()
anno = anno.loc[anno["ID2"].ne("")].copy()
duplicates = anno.loc[anno["ID2"].duplicated(keep=False)]
if not duplicates.empty:
    dup_file = out_dir / "sup_duplicated_annotation_ID2.tsv"
    duplicates.to_csv(dup_file, sep="\t", index=False)
    print(f"WARNING: duplicated ID2 rows detected ({len(duplicates)} rows); first occurrence retained. See {dup_file}")
anno = anno.drop_duplicates(subset="ID2", keep="first")

species = [c[:-len("_cds_pos")] for c in anno.columns if c.endswith("_cds_pos")]
if not species:
    raise ValueError("No species *_cds_pos columns found")

for sp in species:
    needed = [f"{sp}_codon_pos", f"{sp}_codon"]
    absent = [c for c in needed if c not in anno.columns]
    if absent:
        raise ValueError(f"Missing annotation columns for {sp}: {', '.join(absent)}")

fasta_paths = {}
for path in sorted(fasta_dir.glob("*_longestcds.fasta")):
    sp = path.name[:-len("_longestcds.fasta")].split("_")[-1].lower()
    if sp in fasta_paths:
        raise ValueError(f"Multiple longest-CDS FASTA files detected for species {sp}")
    fasta_paths[sp] = path

fasta_db = {}
for sp in species:
    fasta_path = fasta_paths.get(sp.lower())
    if fasta_path:
        fasta_db[sp] = read_fasta(fasta_path)
        print(f"Loaded {sp}: {len(fasta_db[sp])} CDS sequences")
    else:
        fasta_db[sp] = None
        print(f"WARNING: no *_longestcds.fasta found for {sp}; sequence state will be OTHER")

if args.source_priority:
    requested = [x.strip() for x in args.source_priority.split(",") if x.strip()]
    lookup = {sp.lower(): sp for sp in species}
    unknown = [x for x in requested if x.lower() not in lookup]
    if unknown:
        raise ValueError("Unknown species in --source_priority: " + ", ".join(unknown))
    if len({x.lower() for x in requested}) != len(requested):
        raise ValueError("Duplicate species in --source_priority")
    source_priority = [lookup[x.lower()] for x in requested]
    source_priority += [sp for sp in species if sp not in source_priority]
else:
    source_priority = species.copy()

print(f"Annotation rows: {len(anno)}")
print(f"Species: {', '.join(species)}")
print(f"Representative source priority: {', '.join(source_priority)}")

results, qc_rows, unmatched_rows = [], [], []
original_unmapped_rows = 0
new_unmapped_rows = 0
recovery_group_mismatches = 0
not_superset_rows = 0

for _, row in anno.iterrows():
    id2 = row["ID2"]
    original_ids = group_ids_from(row["ID2_group"])
    new_ids = group_ids_from(row["new_id2_group"])
    supplied_recovered_ids = group_ids_from(row["recovered_only_group"])

    original_set = set(original_ids)
    new_set = set(new_ids)
    derived_recovered_ids = [x for x in new_ids if x not in original_set]

    if not original_set.issubset(new_set):
        not_superset_rows += 1

    if supplied_recovered_ids and set(supplied_recovered_ids) != set(derived_recovered_ids):
        recovery_group_mismatches += 1

    original_species, original_unmapped = map_group_to_species(original_ids, row, species)
    new_species, new_unmapped = map_group_to_species(new_ids, row, species)
    recovered_species, recovered_unmapped = map_group_to_species(derived_recovered_ids, row, species)

    if original_unmapped:
        original_unmapped_rows += 1
    if new_unmapped:
        new_unmapped_rows += 1

    per_species, states, reasons = {}, {}, {}
    for sp in species:
        per_species[sp], states[sp], reasons[sp] = extract_species_window(row, sp, fasta_db[sp])

    representative = next((sp for sp in source_priority if per_species[sp] is not None), None)
    source = per_species[representative] if representative else None

    head = {
        "ID2": id2,
        "ID2_group": row["ID2_group"],
        "new_id2_group": row["new_id2_group"],
        "recovered_only_group": row["recovered_only_group"],
        "edited_species": ";".join(original_species) if original_species else "NA",
        "new_edited_species": ";".join(new_species) if new_species else "NA",
        "recovered_species": ";".join(recovered_species) if recovered_species else "NA",
        "n_ID2_group_species": len(original_ids),
        "n_new_id2_group_species": len(new_ids),
        "n_ID2_group_mapped_species": len(original_species),
        "n_new_id2_group_mapped_species": len(new_species),
        "n_recovered_species": len(derived_recovered_ids),
        "source_species": representative or "NA",
        "transcript": source["transcript"] if source else "NA",
        "transcript_pos": source["transcript_pos"] if source else "NA"
    }

    main, qc = dict(head), dict(head)
    qc["unmapped_original_group_IDs"] = ";".join(original_unmapped) if original_unmapped else "NA"
    qc["unmapped_new_group_IDs"] = ";".join(new_unmapped) if new_unmapped else "NA"
    qc["unmapped_recovered_group_IDs"] = ";".join(recovered_unmapped) if recovered_unmapped else "NA"

    for sp in species:
        result = per_species[sp]
        main[f"{sp}_state"] = states[sp]
        qc[f"{sp}_state"] = states[sp]
        qc[f"{sp}_reason"] = reasons[sp]

        for key in ("CDS", "pep", "edit_base_pos", "edit_base"):
            main[f"{sp}_{key}"] = result[key] if result else "NA"

        for key in (
            "transcript", "transcript_pos", "codon_phase",
            "annotated_codon", "fasta_codon", "CDS", "pep",
            "edit_base_pos", "edit_base"
        ):
            qc[f"{sp}_{key}"] = result[key] if result else "NA"

    results.append(main)
    qc_rows.append(qc)

    if representative is None:
        unmatched_rows.append({
            "ID2": id2,
            "reason": ";".join(f"{sp}:{states[sp]}({reasons[sp]})" for sp in species)
        })

out_file = out_dir / "sup_conserved_site_all_species_9nt_3aa.tsv"
qc_file = out_dir / "sup_conserved_site_all_species_9nt_3aa_QC.tsv"
unmatched_file = out_dir / "sup_conserved_site_all_species_unmatched.tsv"

pd.DataFrame(results).to_csv(out_file, sep="\t", index=False)
pd.DataFrame(qc_rows).to_csv(qc_file, sep="\t", index=False)
pd.DataFrame(unmatched_rows, columns=["ID2", "reason"]).to_csv(unmatched_file, sep="\t", index=False)

print("\n========== SUPPLEMENTARY EXTRACTION SUMMARY ==========")
print(f"Annotation rows: {len(anno)}")
print(f"Output site rows: {len(results)}")
print(f"Rows without a representative CDS window: {len(unmatched_rows)}")
print(f"Rows with unmapped original group IDs: {original_unmapped_rows}")
print(f"Rows with unmapped new group IDs: {new_unmapped_rows}")
print(f"Rows where new_id2_group does not contain all original IDs: {not_superset_rows}")
print(f"Rows where recovered_only_group differs from new_id2_group minus ID2_group: {recovery_group_mismatches}")
print("\nGroup-size distributions:")
print("ID2_group:")
print(pd.Series([len(group_ids_from(x)) for x in anno["ID2_group"]]).value_counts().sort_index().to_string())
print("new_id2_group:")
print(pd.Series([len(group_ids_from(x)) for x in anno["new_id2_group"]]).value_counts().sort_index().to_string())
print("\nOutputs:")
print(out_file)
print(qc_file)
print(unmatched_file)
print("Done.")
