#!/usr/bin/env python3
"""Validate a Picard BAM/index against the current reference FASTA dictionary."""

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Validate a Picard BAM/index and compare reference sequence names and lengths."
    )
    p.add_argument("--bam", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--samtools", required=True)
    p.add_argument("--sample-id", required=True)
    return p.parse_args()


def fasta_dictionary(path):
    result = {}
    name = None
    length = 0
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    result[name] = length
                name = line[1:].strip().split()[0]
                if not name or name in result:
                    raise ValueError(f"Invalid/duplicate FASTA sequence name in {path}: {name!r}")
                length = 0
            else:
                length += len("".join(line.split()))
    if name is not None:
        result[name] = length
    if not result:
        raise ValueError(f"No sequences found in reference FASTA: {path}")
    return result


def run(command):
    proc = subprocess.run(command, text=True, capture_output=True)
    if proc.returncode:
        raise RuntimeError(
            f"Command failed: {' '.join(command)}\n"
            f"stdout: {proc.stdout}\nstderr: {proc.stderr}"
        )
    return proc.stdout


def main():
    args = parse_args()
    bam = Path(args.bam)
    reference = Path(args.reference)
    if not bam.is_file() or bam.stat().st_size == 0:
        raise FileNotFoundError(f"BAM missing or empty: {bam}")

    index_candidates = [Path(str(bam) + ".bai"), bam.with_suffix(".bai")]
    index = next((p for p in index_candidates if p.is_file() and p.stat().st_size > 0), None)
    if index is None:
        raise FileNotFoundError(f"BAM BAI index missing or empty: {bam}(.bai)")
    if not reference.is_file() or reference.stat().st_size == 0:
        raise FileNotFoundError(f"Reference FASTA missing or empty: {reference}")

    run([args.samtools, "quickcheck", "-v", str(bam)])
    run([args.samtools, "idxstats", str(bam)])
    header = run([args.samtools, "view", "-H", str(bam)])

    reference_dict = fasta_dictionary(reference)
    bam_sq = {}
    for line in header.splitlines():
        if not line.startswith("@SQ\t"):
            continue
        fields = dict(item.split(":", 1) for item in line.split("\t")[1:] if ":" in item)
        if "SN" in fields and "LN" in fields:
            bam_sq[fields["SN"]] = int(fields["LN"])
    if not bam_sq:
        raise RuntimeError(f"BAM header has no @SQ records: {bam}")

    read_groups = []
    for line in header.splitlines():
        if line.startswith("@RG\t"):
            fields = dict(item.split(":", 1) for item in line.split("\t")[1:] if ":" in item)
            read_groups.append(fields.get("SM"))
    if not read_groups or any(sm != args.sample_id for sm in read_groups):
        raise RuntimeError(
            f"BAM read-group sample identity does not match sample_id={args.sample_id!r}: "
            f"SM values={read_groups!r}"
        )

    missing_from_bam = sorted(set(reference_dict) - set(bam_sq))
    absent_from_reference = sorted(set(bam_sq) - set(reference_dict))
    length_mismatches = [
        (name, bam_sq[name], reference_dict[name])
        for name in sorted(set(bam_sq) & set(reference_dict))
        if bam_sq[name] != reference_dict[name]
    ]
    if missing_from_bam or absent_from_reference or length_mismatches:
        details = []
        if missing_from_bam:
            details.append("reference sequences missing from BAM: " + ", ".join(missing_from_bam[:30]))
        if absent_from_reference:
            details.append("BAM sequences absent from reference: " + ", ".join(absent_from_reference[:30]))
        if length_mismatches:
            details.append("sequence length mismatches: " + "; ".join(
                f"{name}: BAM={bam_len}, FASTA={fasta_len}"
                for name, bam_len, fasta_len in length_mismatches[:30]
            ))
        raise RuntimeError("BAM is incompatible with the current reference FASTA:\n" + "\n".join(details))

    print(
        f"BAM validation PASS: {bam}; index={index}; "
        f"reference sequences checked={len(bam_sq)}; sample_id={args.sample_id}"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
