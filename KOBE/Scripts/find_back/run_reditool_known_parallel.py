#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Run one independent REDItoolKnown recovery process per RNA sample."""

import argparse
import csv
import os
import re
import subprocess
import sys
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run REDItoolKnown sample jobs concurrently; each REDItoolKnown uses one thread."
    )
    parser.add_argument("--sample_table", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--candi_dir", required=True)
    parser.add_argument("--reditool", required=True)
    parser.add_argument("--python", required=True, dest="python_bin")
    parser.add_argument("--intersect_pl", required=True)
    parser.add_argument("--worker_script", required=True)
    parser.add_argument("--jobs", type=int, default=10)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def read_tasks(sample_table):
    tasks = []
    seen = set()

    with open(sample_table, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Sample table is empty")

        reader.fieldnames = [name.strip() for name in reader.fieldnames]
        required = {"sample_id", "species", "ref"}
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise ValueError("Sample table missing columns: " + ", ".join(missing))

        for line_number, raw in enumerate(reader, start=2):
            row = {
                (key or "").strip(): (value or "").strip()
                for key, value in raw.items()
                if key is not None
            }

            sample = row.get("sample_id", "")
            species = row.get("species", "")
            ref = row.get("ref", "")
            sample_type = row.get("type", "RNA").upper()

            if sample_type and sample_type != "RNA":
                continue
            if not sample or not species or not ref:
                raise ValueError(
                    f"Line {line_number}: sample_id, species, or ref is empty"
                )
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample):
                raise ValueError(f"Line {line_number}: invalid sample_id: {sample}")
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", species):
                raise ValueError(f"Line {line_number}: invalid species: {species}")

            key = (species, sample)
            if key in seen:
                raise ValueError(f"Duplicate sample: {species}/{sample}")
            seen.add(key)

            tasks.append({"sample": sample, "species": species, "ref": ref})

    if not tasks:
        raise ValueError("No RNA samples found in the validated sample table")
    return tasks


def require_file(path, label):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"{label} missing or empty: {path}")
    return path


def validate_inputs(args, tasks, out_dir, candi_dir):
    require_file(args.sample_table, "Sample table")
    require_file(args.worker_script, "REDItoolKnown worker script")
    require_file(args.reditool, "REDItoolKnown.py")
    require_file(args.intersect_pl, "Perl intersection script")
    require_file(args.python_bin, "REDItools Python interpreter")
    if not os.access(args.python_bin, os.X_OK):
        raise PermissionError(f"REDItools Python is not executable: {args.python_bin}")

    if not candi_dir.is_dir():
        raise NotADirectoryError(f"Candidate directory not found: {candi_dir}")

    # Validate the REDItools Python 2 environment once before launching jobs.
    subprocess.run(
        [args.python_bin, "-c", "import pysam"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    # Compile-check the Perl helper once before launching jobs.
    perl_check = subprocess.run(
        ["perl", "-c", args.intersect_pl],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if perl_check.returncode != 0:
        raise RuntimeError(
            "Perl intersection script failed syntax validation:\n" + perl_check.stdout
        )
    print(perl_check.stdout.rstrip())

    for task in tasks:
        species = task["species"]
        sample = task["sample"]
        ref = Path(task["ref"])

        sample_dir = out_dir / "Mapping_BWA_to_CDS" / species / sample
        bam = sample_dir / f"{sample}.sort.nodup.bam"
        site_file = sample_dir / "knownEditingSites.tab"
        candidate_file = candi_dir / f"{species}_candi.txt"

        require_file(bam, f"Sorted BAM for {species}/{sample}")
        bai_candidates = [Path(str(bam) + ".bai"), bam.with_suffix(".bai")]
        if not any(path.is_file() and path.stat().st_size > 0 for path in bai_candidates):
            raise FileNotFoundError(
                f"BAM index missing for {species}/{sample}: expected "
                f"{bai_candidates[0]} or {bai_candidates[1]}"
            )
        require_file(site_file, f"Known-site file for {species}/{sample}")
        require_file(ref, f"CDS reference for {species}/{sample}")
        require_file(Path(str(ref) + ".fai"), f"FASTA index for {species}/{sample}")
        require_file(candidate_file, f"Candidate annotation for {species}")


def candidate_hash(path):
    import hashlib
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def outputs_complete(task, out_dir):
    sample, species = task["sample"], task["species"]
    sample_dir = out_dir / "Mapping_BWA_to_CDS" / species / sample
    redi_out = sample_dir / "known"
    raw_file = out_dir / species / "merge" / f"{species}_raw_{sample}.txt"
    count_file = out_dir / "candidate_counts" / species / f"{sample}.candidate_counts.tsv"
    tables = [p for p in redi_out.rglob("outTable_*") if p.is_file() and p.stat().st_size > 0] if redi_out.is_dir() else []
    return bool(tables) and raw_file.is_file() and raw_file.stat().st_size > 0 and count_file.is_file() and count_file.stat().st_size > 0


def run_one(task, args, out_dir, candi_dir, log_dir):
    sample = task["sample"]
    species = task["species"]
    tag = f"{species}__{sample}"

    stdout_path = log_dir / f"{tag}.log"
    stderr_path = log_dir / f"{tag}.error"

    command = [
        "bash", args.worker_script,
        "--sample_table", args.sample_table,
        "--out_dir", str(out_dir),
        "--candi_dir", str(candi_dir),
        "--reditool", args.reditool,
        "--python", args.python_bin,
        "--threads", str(args.threads),
        "--intersect_pl", args.intersect_pl,
        "--sample_id", sample,
        "--species", species,
    ]

    with open(stdout_path, "w", encoding="utf-8") as stdout_handle, \
         open(stderr_path, "w", encoding="utf-8") as stderr_handle:
        result = subprocess.run(
            command,
            stdout=stdout_handle,
            stderr=stderr_handle,
            text=True,
        )

    if result.returncode != 0:
        return {
            "ok": False,
            "sample": sample,
            "species": species,
            "returncode": result.returncode,
            "stdout": str(stdout_path),
            "stderr": str(stderr_path),
        }

    # Confirm the worker produced all expected per-sample outputs.
    sample_dir = out_dir / "Mapping_BWA_to_CDS" / species / sample
    redi_out = sample_dir / "known"
    raw_file = out_dir / species / "merge" / f"{species}_raw_{sample}.txt"
    count_file = out_dir / "candidate_counts" / species / f"{sample}.candidate_counts.tsv"

    out_tables = [
        path for path in redi_out.rglob("outTable_*")
        if path.is_file() and path.stat().st_size > 0
    ] if redi_out.is_dir() else []

    missing = []
    if not out_tables:
        missing.append(f"non-empty REDItoolKnown outTable_* under {redi_out}")
    if not raw_file.is_file() or raw_file.stat().st_size == 0:
        missing.append(str(raw_file))
    if not count_file.is_file() or count_file.stat().st_size == 0:
        missing.append(str(count_file))

    if missing:
        return {
            "ok": False,
            "sample": sample,
            "species": species,
            "returncode": 1,
            "stdout": str(stdout_path),
            "stderr": str(stderr_path),
            "message": "Worker returned success but outputs are missing: " + "; ".join(missing),
        }

    return {
        "ok": True,
        "sample": sample,
        "species": species,
        "returncode": 0,
        "stdout": str(stdout_path),
        "stderr": str(stderr_path),
    }


def rebuild_sample_lists(tasks, out_dir):
    species_samples = OrderedDict()
    for task in tasks:
        species_samples.setdefault(task["species"], []).append(task["sample"])

    for species, samples in species_samples.items():
        merge_dir = out_dir / species / "merge"
        merge_dir.mkdir(parents=True, exist_ok=True)
        target = merge_dir / "sample.list"
        temporary = merge_dir / f"sample.list.tmp.{os.getpid()}"

        with open(temporary, "w", encoding="utf-8", newline="") as handle:
            for sample in samples:
                handle.write(sample + "\n")

        os.replace(temporary, target)
        print(f"Updated sample list: {target} ({len(samples)} samples)")


def main():
    args = parse_args()
    if args.jobs < 1 or args.threads < 1:
        raise ValueError("--jobs and --threads must be positive integers")

    sample_table = Path(args.sample_table).resolve()
    out_dir = Path(args.out_dir).resolve()
    candi_dir = Path(args.candi_dir).resolve()
    worker_script = Path(args.worker_script).resolve()
    reditool = Path(args.reditool).resolve()
    python_bin = Path(args.python_bin).resolve()
    intersect_pl = Path(args.intersect_pl).resolve()

    args.sample_table = str(sample_table)
    args.out_dir = str(out_dir)
    args.candi_dir = str(candi_dir)
    args.worker_script = str(worker_script)
    args.reditool = str(reditool)
    args.python_bin = str(python_bin)
    args.intersect_pl = str(intersect_pl)

    if not sample_table.is_file():
        raise FileNotFoundError(f"Sample table not found: {sample_table}")
    out_dir.mkdir(parents=True, exist_ok=True)

    tasks = read_tasks(sample_table)
    state_dir = out_dir / "state"
    state_dir.mkdir(parents=True, exist_ok=True)
    hash_file = state_dir / "candidate_site_hashes.tsv"
    old_hashes = {}
    if hash_file.is_file():
        with open(hash_file, encoding="utf-8") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) == 2:
                    old_hashes[parts[0]] = parts[1]
    current_hashes = {}
    for species in OrderedDict((task["species"], None) for task in tasks):
        candidate_file = candi_dir / f"{species}_candi.txt"
        require_file(candidate_file, f"Candidate annotation for {species}")
        current_hashes[species] = candidate_hash(candidate_file)
    log_dir = out_dir / "logs" / "reditool_known"
    log_dir.mkdir(parents=True, exist_ok=True)

    validate_inputs(args, tasks, out_dir, candi_dir)

    pending_tasks = []
    skipped_tasks = []
    for task in tasks:
        species = task["species"]
        reusable = (args.resume and not args.force and
                    old_hashes.get(species) == current_hashes[species] and
                    outputs_complete(task, out_dir))
        if reusable:
            skipped_tasks.append(task)
        else:
            pending_tasks.append(task)

    worker_count = min(args.jobs, max(1, len(pending_tasks)))
    print("============================================================")
    print("Sample-level REDItoolKnown parallel processing")
    print(f"RNA samples: {len(tasks)}")
    print(f"Concurrent sample processes: {worker_count}")
    print(f"REDItoolKnown threads per process: {args.threads}")
    print(f"Reusable completed samples: {len(skipped_tasks)}")
    print(f"Samples to run: {len(pending_tasks)}")
    print(f"Per-sample logs: {log_dir}")
    print("============================================================")

    results = []
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = {
            executor.submit(run_one, task, args, out_dir, candi_dir, log_dir): task
            for task in pending_tasks
        }
        for future in as_completed(futures):
            task = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {
                    "ok": False,
                    "sample": task["sample"],
                    "species": task["species"],
                    "returncode": 1,
                    "message": str(exc),
                }
            results.append(result)
            label = f"{result['species']}/{result['sample']}"
            if result["ok"]:
                print(f"[PASS] {label}")
            else:
                print(f"[FAIL] {label}")
                if result.get("message"):
                    print(f"       {result['message']}")
                print(f"       stdout: {result.get('stdout', 'not available')}")
                print(f"       stderr: {result.get('stderr', 'not available')}")

    failures = [result for result in results if not result["ok"]]
    if failures:
        raise RuntimeError(
            f"{len(failures)} of {len(tasks)} sample jobs failed; "
            "multi-sample summary and downstream analysis were not started. "
            f"See logs in {log_dir}"
        )

    rebuild_sample_lists(tasks, out_dir)
    tmp_hash = hash_file.with_suffix(".tmp")
    with open(tmp_hash, "w", encoding="utf-8", newline="") as handle:
        for species, digest in current_hashes.items():
            handle.write(f"{species}\t{digest}\n")
    os.replace(tmp_hash, hash_file)
    print("All REDItoolKnown sample jobs completed successfully.")
    print("Multi-sample summarization will be performed by main.sh after this stage.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
