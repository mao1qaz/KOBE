#!/usr/bin/env python3
# -*- coding: gbk -*-

"""
Generate aligned CDS from protein alignment
STOP-aware column model (stable, no overwrite, NO_HIT supported)
"""

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed


# =========================
def load_cds(cds_file):
    cds_dict = {}

    if cds_file is None or not cds_file.exists():
        return cds_dict

    with open(cds_file) as f:
        seq = []
        gene_id = None

        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if gene_id:
                    cds_dict[gene_id] = "".join(seq)
                gene_id = line[1:].split('|')[-1]
                seq = []
            else:
                seq.append(line)

        if gene_id:
            cds_dict[gene_id] = "".join(seq)

    return cds_dict


# =========================
def read_alignment(pep_file):
    sequences = {}
    headers = []

    with open(pep_file) as f:
        seq_lines = []
        gene_id = None

        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if gene_id:
                    sequences[gene_id] = "".join(seq_lines)
                gene_id = line[1:]
                headers.append(gene_id)
                seq_lines = []
            else:
                seq_lines.append(line)

        if gene_id:
            sequences[gene_id] = "".join(seq_lines)

    return headers, sequences


# =========================
def process_file(pep_file, cds_dir, out_dir):

    name = pep_file.stem
    out_file = out_dir / f"{name}.aln_cds.fasta"

    print(f"[Processing] {name}")

    headers, seqs = read_alignment(pep_file)

    # =========================
    # EMPTY FILE
    # =========================
    if not seqs or not any(seqs.values()):
        print(f"[EMPTY] {name}")
        with open(out_file, "w") as f:
            pass
        return out_file

    aln_len = len(next(iter(seqs.values())))

    stop_codons = {"TAA", "TAG", "TGA"}

    # =========================
    # NO_HIT
    # =========================
    nohit_species = set()
    for h in headers:
        if ("NO_HIT" in h) or ("NOT_FOUND" in h):
            nohit_species.add(h)

    # =========================
    # CDS cache
    # =========================
    cds_cache = {}
    pos = {}

    for h in headers:
        species = h.split("|")[0]

        if species not in cds_cache:
            cds_file = next(cds_dir.glob(f"*{species}*_longestcds.fasta"), None)
            cds_cache[species] = load_cds(cds_file)

        pos[h] = 0

    # output
    out = {h: [] for h in headers}

    # =========================
    # MAIN LOOP
    # =========================
    i = 0
    while i < aln_len:

        while True:

            col_codons = {}
            has_stop = False

            # -------- pre-check --------
            for h in headers:

                if h in nohit_species:
                    col_codons[h] = None
                    continue

                aa = seqs[h][i]

                if aa == "-":
                    col_codons[h] = None
                    continue

                species = h.split("|")[0]
                gene = h.split("|")[-1]

                cds_dict = cds_cache.get(species)

                if not cds_dict or gene not in cds_dict:
                    col_codons[h] = "XXX"
                    continue

                idx = pos[h]
                codon = cds_dict[gene][idx*3:idx*3+3]

                col_codons[h] = codon

                if codon in stop_codons:
                    has_stop = True

            # =========================
            # STOP MODE
            # =========================
            if has_stop:

                for h in headers:

                    if h in nohit_species:
                        out[h].append("XXX")
                        continue

                    aa = seqs[h][i]
                    codon = col_codons[h]

                    if aa == "-":
                        out[h].append("---")
                        continue

                    if codon in stop_codons:
                        out[h].append(codon)
                        pos[h] += 1
                    else:
                        out[h].append("---")

                # ?? 关键：不推进 i
                continue

            # =========================
            # NORMAL MODE
            # =========================
            for h in headers:

                if h in nohit_species:
                    out[h].append("XXX")
                    continue

                aa = seqs[h][i]

                if aa == "-":
                    out[h].append("---")
                    continue

                species = h.split("|")[0]
                gene = h.split("|")[-1]

                cds_dict = cds_cache.get(species)

                if not cds_dict or gene not in cds_dict:
                    out[h].append("XXX")
                    continue

                idx = pos[h]
                codon = cds_dict[gene][idx*3:idx*3+3]

                out[h].append(codon)
                pos[h] += 1

            break  # NORMAL 完成 → 下一列

        i += 1

    # =========================
    # OUTPUT
    # =========================
    with open(out_file, "w") as f:
        for h in headers:
            f.write(f">{h}\n")
            f.write("".join(out[h]) + "\n")

    return out_file


# =========================
def main(work_dir, cds_dir, threads=8):

    work_dir = Path(work_dir)
    cds_dir = Path(cds_dir)

    pep_dir = work_dir / "aln_pep"
    out_dir = work_dir / "aln_cds"

    out_dir.mkdir(parents=True, exist_ok=True)

    pep_files = sorted(pep_dir.glob("*.fa*"))

    if threads < 1:
        raise ValueError("threads must be >= 1")

    print(f"Found {len(pep_files)} files")

    if not pep_files:
        print("No protein alignment files found; nothing to back-translate.")
        return

    workers = min(threads, len(pep_files))
    print(f"Using up to {workers} worker processes (1 process per task; no internal thread pool).")

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(process_file, f, cds_dir, out_dir)
            for f in pep_files
        ]

        for fu in as_completed(futures):
            print("[Done]", fu.result())


# =========================
if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--work_dir", required=True)
    parser.add_argument("--cds_dir", required=True)
    parser.add_argument("--threads", type=int, default=8)

    args = parser.parse_args()

    main(args.work_dir, args.cds_dir, args.threads)

    print("All done!")