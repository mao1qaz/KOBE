#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from collections import Counter
from pathlib import Path

import pandas as pd

# ============================================================
# 1. Input / output
# ============================================================
parser = argparse.ArgumentParser(description="Extract homologous CDS windows and edited species from ID2_group.")
parser.add_argument("--annotated", required=True, help="Annotated homologous editing-site CSV")
parser.add_argument("--cds_dir", required=True, help="Directory containing *_longestcds.fasta files")
parser.add_argument("--out_dir", required=True, help="Output directory")
parser.add_argument("--source_priority", default=None, help="Comma-separated priority for representative transcripts, e.g. sapiens,musculus")
args = parser.parse_args()

anno_file = Path(args.annotated)
fasta_dir = Path(args.cds_dir)
out_dir = Path(args.out_dir)

if not anno_file.is_file():
    raise FileNotFoundError(f"Annotation file not found: {anno_file}")
if not fasta_dir.is_dir():
    raise NotADirectoryError(f"FASTA directory not found: {fasta_dir}")
out_dir.mkdir(parents=True, exist_ok=True)

# ============================================================
# 2. Standard genetic code
# ============================================================
genetic_code = {
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

# ============================================================
# 3. Helper functions
# ============================================================
def read_fasta(path):
    seqs = {}
    name, pieces = None, []
    with open(path) as f:
        for line in f:
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
    return ";".join(genetic_code.get(seq[i:i+3], "Xaa") for i in range(0, len(seq), 3))


def group_ids_from(value):
    parts = [x.strip() for x in str(value).split(";")]
    return list(dict.fromkeys(x for x in parts if x and x.upper() not in {"NA", "NAN", "UNFOUND", "-", "---"}))


# ============================================================
# 4. Extract CDS window for one species
# ============================================================
def extract_species_window(row, sp, seqs):
    raw_pos = str(row[f"{sp}_cds_pos"]).strip()
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

    # The central codon must be complete (1-based CDS coordinates).
    codon_start = cds_pos - phase + 1
    codon_end = codon_start + 2
    if codon_start < 1 or codon_end > len(seq):
        return None, "OTHER", "central_codon_out_of_range"

    # Include adjacent codons only when complete.
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


# ============================================================
# 5. Read annotation table
# ============================================================
try:
    anno = pd.read_csv(anno_file, dtype=str, keep_default_na=False, encoding="utf-8-sig")
except UnicodeDecodeError:
    anno = pd.read_csv(anno_file, dtype=str, keep_default_na=False, encoding="gbk")

for col in ("ID2", "ID2_group"):
    if col not in anno.columns:
        raise ValueError(f"Missing annotation column: {col}")

anno["ID2"] = anno["ID2"].str.strip()
anno = anno[anno["ID2"].ne("")].copy()

duplicates = anno.loc[anno["ID2"].duplicated(keep=False)]
if not duplicates.empty:
    dup_file = out_dir / "duplicated_annotation_ID2.tsv"
    duplicates.to_csv(dup_file, sep="\t", index=False)
    print(f"WARNING: {len(duplicates)} rows have duplicated ID2; using first occurrence. See {dup_file}")
anno = anno.drop_duplicates(subset="ID2", keep="first")

# ============================================================
# 6. Automatically detect species, FASTA, and priority
# ============================================================
species = [c[:-len("_cds_pos")] for c in anno.columns if c.endswith("_cds_pos")]
if not species:
    raise ValueError("No species *_cds_pos columns found")

if args.source_priority:
    requested = [x.strip() for x in args.source_priority.split(",") if x.strip()]
    unknown = [sp for sp in requested if sp not in species]
    if unknown:
        raise ValueError(f"Unknown species in --source_priority: {', '.join(unknown)}. Available: {', '.join(species)}")
    if len(requested) != len(set(requested)):
        raise ValueError("Duplicate species in --source_priority")
    source_priority = requested + [sp for sp in species if sp not in requested]
else:
    source_priority = species.copy()

fasta_paths = {}
for path in sorted(fasta_dir.glob("*_longestcds.fasta")):
    sp = path.name[:-len("_longestcds.fasta")].split("_")[-1].lower()
    if sp in fasta_paths:
        raise ValueError(f"Multiple FASTA files detected for species {sp}")
    fasta_paths[sp] = path

fasta_db = {}
for sp in species:
    if sp in fasta_paths:
        fasta_db[sp] = read_fasta(fasta_paths[sp])
        print(f"Loaded {sp}: {len(fasta_db[sp])} CDS sequences")
    else:
        fasta_db[sp] = None
        print(f"WARNING: No CDS FASTA found for {sp}; sequence state will be OTHER when coordinates are present")

print(f"Annotation rows: {len(anno)}")
print(f"Species in annotation: {', '.join(species)}")
print(f"Representative source priority: {', '.join(source_priority)}")

# ============================================================
# 7. Recover edited species from ID2_group; extract sequences
# ============================================================
results = []
qc_rows = []
unmatched_rows = []
group_mapping_issues = 0

for _, row in anno.iterrows():
    id2 = row["ID2"]
    group_ids = group_ids_from(row["ID2_group"])
    group_set = set(group_ids)

    # One valid ID2_group entry corresponds to one edited species.
    edited_species = [sp for sp in species if str(row[f"{sp}_cds_pos"]).strip() in group_set]
    recognized_ids = {str(row[f"{sp}_cds_pos"]).strip() for sp in species if str(row[f"{sp}_cds_pos"]).strip() in group_set}
    unmapped_ids = [x for x in group_ids if x not in recognized_ids]
    if unmapped_ids:
        group_mapping_issues += 1

    per_species, states, reasons = {}, {}, {}
    for sp in species:
        per_species[sp], states[sp], reasons[sp] = extract_species_window(row, sp, fasta_db[sp])

    # Choose representative coordinates using the user-specified priority.
    representative = next((sp for sp in source_priority if per_species[sp] is not None), None)
    source = per_species[representative] if representative else None

    head = {
        "ID2": id2,
        "ID2_group": row["ID2_group"],
        "edited_species": ";".join(edited_species) if edited_species else "NA",
        "n_edited_species": len(group_ids),
        "source_species": representative or "NA",
        "transcript": source["transcript"] if source else "NA",
        "transcript_pos": source["transcript_pos"] if source else "NA"
    }
    main, qc = dict(head), dict(head)
    qc["unmapped_group_IDs"] = ";".join(unmapped_ids) if unmapped_ids else "NA"

    for sp in species:
        result = per_species[sp]
        main[f"{sp}_state"] = states[sp]
        qc[f"{sp}_state"] = states[sp]
        qc[f"{sp}_reason"] = reasons[sp]

        for key in ("CDS", "pep", "edit_base_pos", "edit_base"):
            main[f"{sp}_{key}"] = result[key] if result else "NA"
        for key in ("transcript", "transcript_pos", "codon_phase", "annotated_codon", "fasta_codon", "CDS", "pep", "edit_base_pos", "edit_base"):
            qc[f"{sp}_{key}"] = result[key] if result else "NA"

    results.append(main)
    qc_rows.append(qc)
    if representative is None:
        unmatched_rows.append({
            "ID2": id2,
            "reason": ";".join(f"{sp}:{states[sp]}({reasons[sp]})" for sp in species)
        })

# ============================================================
# 8. Save results
# ============================================================
out_file = out_dir / "conserved_site_all_species_9nt_3aa.tsv"
qc_file = out_dir / "conserved_site_all_species_9nt_3aa_QC.tsv"
unmatched_file = out_dir / "conserved_site_all_species_unmatched.tsv"

pd.DataFrame(results).to_csv(out_file, sep="\t", index=False)
pd.DataFrame(qc_rows).to_csv(qc_file, sep="\t", index=False)
pd.DataFrame(unmatched_rows, columns=["ID2", "reason"]).to_csv(unmatched_file, sep="\t", index=False)

# ============================================================
# 9. Summary
# ============================================================
print(f"\nOutput ID2 rows: {len(results)}")
print(f"All-species sequence failures: {len(unmatched_rows)}")
print(f"Rows with unmatched ID2_group entries: {group_mapping_issues}")
for sp in species:
    counts = Counter(q[f"{sp}_state"] for q in qc_rows)
    print(f"[{sp}] " + ", ".join(f"{s}={counts.get(s, 0)}" for s in ("SEQ", "NO_ORTHOLOG", "GAP", "OTHER")))

print("\nOutputs:")
print(out_file)
print(qc_file)
print(unmatched_file)
print("Done.")
