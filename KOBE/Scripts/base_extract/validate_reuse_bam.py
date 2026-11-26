#!/usr/bin/env python3
"""Validate/stage in-workdir Picard BAMs for base_extract --force_reditool.

This single helper performs the pre-Snakemake sample scan, BAM/index/reference
validation, staging, manifest creation, and Snakemake validation-marker creation.
"""

import argparse
import csv
import json
import os
import re
import subprocess
import shutil
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Automatically validate existing in-workdir Picard BAMs for --force_reditool."
    )
    parser.add_argument("--sample-table", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--samtools", required=True)
    return parser.parse_args()


def resolve(value, base):
    value = (value or "").strip().strip('"').strip("'")
    if not value:
        return ""
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    return str(path.resolve(strict=False))


def safe_name(value, label, line_no):
    value = (value or "").strip()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise ValueError(f"Invalid {label} at sample-table line {line_no}: {value!r}")
    return value


def fasta_lengths(path):
    lengths = {}
    name = None
    length = 0
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = line[1:].strip().split()[0]
                if not name or name in lengths:
                    raise ValueError(f"Invalid or duplicate FASTA sequence name in {path}: {name!r}")
                length = 0
            else:
                length += len("".join(line.split()))
    if name is not None:
        lengths[name] = length
    if not lengths:
        raise ValueError(f"No FASTA sequences found: {path}")
    return lengths


def bam_reference_matches(bam, reference, samtools, sample_id):
    result = subprocess.run(
        [samtools, "view", "-H", str(bam)],
        text=True, capture_output=True
    )
    if result.returncode:
        return False, "cannot read BAM header: " + result.stderr.strip()
    ref_lengths = fasta_lengths(reference)
    bam_sq = {}
    for line in result.stdout.splitlines():
        if line.startswith("@SQ\t"):
            fields = dict(
                item.split(":", 1) for item in line.split("\t")[1:] if ":" in item
            )
            if "SN" in fields and "LN" in fields:
                bam_sq[fields["SN"]] = int(fields["LN"])
    if not bam_sq:
        return False, "BAM header has no @SQ records"
    read_groups = []
    for line in result.stdout.splitlines():
        if line.startswith("@RG\t"):
            fields = dict(
                item.split(":", 1) for item in line.split("\t")[1:] if ":" in item
            )
            read_groups.append(fields.get("SM"))
    if not read_groups or any(sm != sample_id for sm in read_groups):
        return False, f"BAM read-group sample identity does not match {sample_id}: {read_groups!r}"
    missing_from_bam = sorted(set(ref_lengths) - set(bam_sq))
    absent_from_reference = sorted(set(bam_sq) - set(ref_lengths))
    length_mismatches = [
        (name, bam_sq[name], ref_lengths[name])
        for name in sorted(set(bam_sq) & set(ref_lengths))
        if bam_sq[name] != ref_lengths[name]
    ]
    if missing_from_bam or absent_from_reference or length_mismatches:
        details = []
        if missing_from_bam:
            details.append("reference sequences missing from BAM: " + ", ".join(missing_from_bam[:20]))
        if absent_from_reference:
            details.append("BAM sequences absent from reference: " + ", ".join(absent_from_reference[:20]))
        if length_mismatches:
            details.append("sequence length mismatches: " + "; ".join(
                f"{name}: BAM={bam_len}, FASTA={fasta_len}"
                for name, bam_len, fasta_len in length_mismatches[:20]
            ))
        return False, "reference sequence dictionary differs: " + " | ".join(details)
    return True, f"all {len(bam_sq)} BAM reference sequences match the current FASTA dictionary"


def link_if_needed(source, destination):
    source = Path(source).resolve()
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink():
        try:
            if destination.resolve() == source:
                return
        except OSError:
            pass
    if destination.exists() or destination.is_symlink():
        if destination.is_dir() and not destination.is_symlink():
            raise RuntimeError(f"Refusing to replace directory with BAM symlink: {destination}")
        destination.unlink()
    temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
    if temporary.exists() or temporary.is_symlink():
        temporary.unlink()
    temporary.symlink_to(source)
    os.replace(temporary, destination)


def main():
    args = parse_args()
    sample_table = Path(args.sample_table).resolve()
    out_dir = Path(args.out_dir).resolve()
    samtools = args.samtools

    reused, fallback, seen = [], [], set()

    # Clear only artifacts owned by this validation helper. The launcher holds
    # the output-directory flock, so no concurrent base_extract run can be
    # modifying these paths while they are rebuilt.
    stage_root = out_dir / ".reuse_bam_inputs"
    marker_root = out_dir / ".reuse_bam_validated"
    manifest = out_dir / ".reuse_bam_manifest.json"
    shutil.rmtree(stage_root, ignore_errors=True)
    shutil.rmtree(marker_root, ignore_errors=True)
    try:
        manifest.unlink()
    except FileNotFoundError:
        pass
    stage_root.mkdir(parents=True, exist_ok=True)
    marker_root.mkdir(parents=True, exist_ok=True)
    with sample_table.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Sample table has no header.")
        reader.fieldnames = [x.strip().lstrip("\ufeff") for x in reader.fieldnames]
        required = {"sample_id", "species", "fq1", "fq2"}
        missing = required - set(reader.fieldnames)
        genome_col = next(
            (x for x in ("genome_fasta", "genome", "reference_fasta", "ref")
             if x in reader.fieldnames), None
        )
        if missing or genome_col is None:
            raise ValueError(
                "Sample table must contain sample_id, species, fq1, fq2, "
                "and genome_fasta (or genome/ref)."
            )

        for line_no, row in enumerate(reader, start=2):
            if not row or not any((v or "").strip() for v in row.values()):
                continue
            sample_id = safe_name(row.get("sample_id"), "sample_id", line_no)
            species = safe_name(row.get("species"), "species", line_no)
            key = (species, sample_id)
            if key in seen:
                raise ValueError(f"Duplicate sample: {species}/{sample_id}")
            seen.add(key)

            reference = Path(resolve(row.get(genome_col), sample_table.parent))
            if not reference.is_file():
                raise FileNotFoundError(
                    f"Reference FASTA missing for {species}/{sample_id}: {reference}"
                )

            sample_dir = out_dir / species / "03_PICARD" / sample_id
            bam = sample_dir / f"{sample_id}.sort.q10.nodup.bam"
            index_candidates = [Path(str(bam) + ".bai"), bam.with_suffix(".bai")]
            index = next(
                (candidate for candidate in index_candidates
                 if candidate.is_file() and candidate.stat().st_size > 0),
                None
            )

            valid, reason = False, "BAM or index missing/empty"
            if bam.is_file() and bam.stat().st_size > 0 and index is not None:
                quick = subprocess.run(
                    [samtools, "quickcheck", "-v", str(bam)],
                    text=True, capture_output=True
                )
                idx = subprocess.run(
                    [samtools, "idxstats", str(bam)],
                    text=True, capture_output=True
                )
                if quick.returncode != 0:
                    reason = "samtools quickcheck failed: " + quick.stderr.strip()
                elif idx.returncode != 0:
                    reason = "samtools idxstats failed: " + idx.stderr.strip()
                else:
                    valid, reason = bam_reference_matches(bam, reference, samtools, sample_id)

            if valid:
                stage_dir = out_dir / ".reuse_bam_inputs" / species
                staged_bam = stage_dir / f"{sample_id}.sort.q10.nodup.bam"
                staged_index = Path(str(staged_bam) + ".bai")
                link_if_needed(bam, staged_bam)
                link_if_needed(index, staged_index)
                marker = marker_root / species / f"{sample_id}.ok"
                marker.parent.mkdir(parents=True, exist_ok=True)
                marker_payload = {
                    "bam": str(bam.resolve()),
                    "index": str(index.resolve()),
                    "reference": str(reference),
                    "sample_id": sample_id,
                    "bam_size": bam.stat().st_size,
                    "bam_mtime_ns": bam.stat().st_mtime_ns,
                    "index_size": index.stat().st_size,
                    "index_mtime_ns": index.stat().st_mtime_ns,
                    "validation": reason,
                    "status": "PASS",
                }
                marker_tmp = marker.with_name(marker.name + f".tmp.{os.getpid()}")
                with marker_tmp.open("w", encoding="utf-8") as marker_handle:
                    json.dump(marker_payload, marker_handle, indent=2)
                    marker_handle.write("\n")
                os.replace(marker_tmp, marker)
                reused.append({
                    "species": species, "sample_id": sample_id,
                    "bam": str(bam.resolve()), "index": str(index.resolve()),
                    "reference": str(reference), "marker": str(marker),
                })
                print(f"REUSE BAM: {species}/{sample_id} -> {bam}")
                print(f"  Validation: {reason}")
            else:
                # Do not let Snakemake mistake an incompatible/incomplete final
                # BAM or its upstream mapping products for valid completed outputs.
                if "reference sequence" in reason:
                    star_index = out_dir / species / "00_REF" / "STAR"
                    if star_index.is_dir():
                        shutil.rmtree(star_index)
                    reference_fai = out_dir / species / "00_REF" / "GENOME" / f"{species}.fa.fai"
                    try:
                        reference_fai.unlink()
                    except FileNotFoundError:
                        pass
                for stale in (
                    bam, *index_candidates,
                    sample_dir / f"{sample_id}.picard.metrics.txt",
                    out_dir / species / "02_Q10" / sample_id / f"{sample_id}.sort.q10.bam",
                    out_dir / species / "02_Q10" / sample_id / f"{sample_id}.sort.q10.bam.bai",
                    out_dir / species / "01_MAP" / sample_id / f"{sample_id}.sort.bam",
                ):
                    try:
                        stale.unlink()
                    except FileNotFoundError:
                        pass
                fq1 = resolve(row.get("fq1"), sample_table.parent)
                fq2 = resolve(row.get("fq2"), sample_table.parent)
                raw_pair_ok = bool(
                    fq1 and fq2
                    and Path(fq1).is_file() and Path(fq1).stat().st_size > 0
                    and Path(fq2).is_file() and Path(fq2).stat().st_size > 0
                )
                clean_r1 = out_dir / species / "00_QC" / sample_id / "qc.R1.fastq.gz"
                clean_r2 = out_dir / species / "00_QC" / sample_id / "qc.R2.fastq.gz"
                clean_pair_ok = bool(
                    clean_r1.is_file() and clean_r1.stat().st_size > 0
                    and clean_r2.is_file() and clean_r2.stat().st_size > 0
                )
                if not raw_pair_ok and not clean_pair_ok:
                    raise FileNotFoundError(
                        f"No compatible complete Picard BAM/index, raw FASTQ pair, or retained clean FASTQ pair "
                        f"for {species}/{sample_id}. Reason: {reason}; expected BAM: {bam}; "
                        f"fq1={fq1 or '[empty]'}; fq2={fq2 or '[empty]'}; "
                        f"clean_R1={clean_r1}; clean_R2={clean_r2}"
                    )
                input_source = "raw_fastq" if raw_pair_ok else "retained_clean_fastq"
                fallback.append({"species": species, "sample_id": sample_id, "reason": reason, "input_source": input_source})
                if raw_pair_ok:
                    print(f"MAP FROM RAW FASTQ: {species}/{sample_id} ({reason})")
                else:
                    print(f"MAP FROM RETAINED CLEAN FASTQ: {species}/{sample_id} ({reason})")

    temporary = manifest.with_name(manifest.name + f".tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump({"reused": reused, "fallback": fallback}, handle, indent=2)
        handle.write("\n")
    os.replace(temporary, manifest)

    print("=" * 60)
    print(f"Samples checked: {len(seen)}")
    print(f"BAMs accepted for reuse: {len(reused)}")
    print(f"Samples falling back to FASTQ mapping: {len(fallback)}")
    print(f"Validation markers written: {len(reused)}")
    print(f"Reuse manifest: {manifest}")
    print("=" * 60)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
