
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import subprocess
import shutil
from pathlib import Path
from collections import defaultdict


COMPLEMENT = bytes.maketrans(
    b"ACGTRYMKBDHVNacgtrymkbdhvn",
    b"TGCAYRKMVHDBNtgcayrkmvhdbn"
)


def parse_args():

    parser = argparse.ArgumentParser(
        prog="extract_longest_cds",
        description="Extract the longest CDS transcript per gene.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n  KOBE extract_longest_cds \\\n    --gtf /project/ref/invicta.gtf \\\n    --genome /project/ref/invicta.fa \\\n    --output /project/results/invicta_longestcds.fasta \\\n    --sif /project/rna_editing.sif --bind /project",
    )

    parser.add_argument(
        "--gtf",
        required=True,
        metavar="FILE",
        help="Input GTF annotation."
    )

    parser.add_argument(
        "--genome",
        required=True,
        metavar="FILE",
        help="Input genome FASTA."
    )

    parser.add_argument(
        "--output",
        required=True,
        metavar="FILE",
        help="Output longest CDS FASTA."
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


def parse_attributes(text):

    attrs = {}

    for item in text.strip().strip(";").split(";"):

        item = item.strip()

        if not item:
            continue

        match = re.match(r'^(\S+)\s+"([^"]*)"$', item)

        if match:
            attrs[match.group(1)] = match.group(2)
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key.strip()] = value.strip().strip('"')

    return attrs


def normalize_transcript_id(tx):

    tx = tx.strip()

    if tx.startswith("rna-"):
        tx = tx[4:]

    return tx


class FastaReader:

    def __init__(self, fasta):

        self.fasta = Path(fasta)
        self.fai = Path(str(self.fasta) + ".fai")

        if not self.fai.exists():

            samtools = shutil.which("samtools")

            if not samtools:
                raise RuntimeError(
                    "FASTA index is missing and samtools is unavailable. "
                    "Please run samtools faidx first."
                )

            print("[INFO] Building FASTA index...")

            subprocess.run(
                [samtools, "faidx", str(self.fasta)],
                check=True
            )

        self.index = {}

        with open(self.fai) as f:

            for line in f:

                fields = line.rstrip("\n").split("\t")

                if len(fields) < 5:
                    raise ValueError("Invalid FASTA index.")

                name = fields[0]

                self.index[name] = {
                    "length": int(fields[1]),
                    "offset": int(fields[2]),
                    "line_bases": int(fields[3]),
                    "line_width": int(fields[4])
                }

        self.handle = open(self.fasta, "rb")

    def fetch(self, chrom, start, end):

        if chrom not in self.index:
            raise KeyError(f"Contig not found: {chrom}")

        info = self.index[chrom]

        if start < 1 or end < start or end > info["length"]:
            raise ValueError(
                f"Invalid coordinates: {chrom}:{start}-{end}"
            )

        pos = start - 1
        end0 = end - 1

        sequence = bytearray()

        while pos <= end0:

            line_index = pos // info["line_bases"]
            column = pos % info["line_bases"]

            n = min(
                end0 - pos + 1,
                info["line_bases"] - column
            )

            byte_offset = (
                info["offset"]
                + line_index * info["line_width"]
                + column
            )

            self.handle.seek(byte_offset)

            chunk = self.handle.read(n)

            if len(chunk) != n:
                raise IOError(
                    f"Failed to extract sequence: {chrom}:{start}-{end}"
                )

            sequence.extend(chunk)

            pos += n

        return bytes(sequence).upper()

    def close(self):

        self.handle.close()


def reverse_complement(sequence):

    return sequence.translate(COMPLEMENT)[::-1]


def parse_gtf(gtf_file):

    models = defaultdict(lambda: {
        "gene_id": None,
        "chrom": None,
        "strand": None,
        "segments": set(),
        "errors": []
    })

    total_cds_records = 0
    missing_attributes = 0
    removed_rna_prefix = 0

    with open(gtf_file) as f:

        for line in f:

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            if fields[2].lower() != "cds":
                continue

            total_cds_records += 1

            chrom = fields[0]
            start = fields[3]
            end = fields[4]
            strand = fields[6]
            phase = fields[7]

            attrs = parse_attributes(fields[8])

            gene_id = attrs.get("gene_id")
            tx_raw = attrs.get("transcript_id")

            if not tx_raw or not gene_id:
                missing_attributes += 1
                continue

            tx = normalize_transcript_id(tx_raw)

            if tx != tx_raw:
                removed_rna_prefix += 1

            model = models[tx]

            if model["gene_id"] is None:
                model["gene_id"] = gene_id

            elif model["gene_id"] != gene_id:
                model["errors"].append("conflicting_gene_id")

            if model["chrom"] is None:
                model["chrom"] = chrom

            elif model["chrom"] != chrom:
                model["errors"].append("conflicting_chromosome")

            if model["strand"] is None:
                model["strand"] = strand

            elif model["strand"] != strand:
                model["errors"].append("conflicting_strand")

            try:
                start = int(start)
                end = int(end)

                if start < 1 or end < start:
                    raise ValueError

            except ValueError:
                model["errors"].append("invalid_coordinates")
                continue

            if strand not in {"+", "-"}:
                model["errors"].append("invalid_strand")

            if phase not in {"0", "1", "2", "."}:
                model["errors"].append("invalid_phase")

            model["segments"].add((start, end))

    return models, total_cds_records, missing_attributes, removed_rna_prefix


def validate_models(models, fasta):

    valid_models = {}

    invalid_count = 0

    for tx, model in models.items():

        errors = list(model["errors"])

        segments = sorted(model["segments"])

        if not segments:
            errors.append("no_CDS_segments")

        if model["chrom"] not in fasta.index:
            errors.append("contig_not_found")

        else:

            chrom_length = fasta.index[model["chrom"]]["length"]

            for start, end in segments:

                if end > chrom_length:
                    errors.append("coordinate_out_of_range")
                    break

        for i in range(1, len(segments)):

            if segments[i][0] <= segments[i - 1][1]:
                errors.append("overlapping_CDS_segments")
                break

        if errors:

            invalid_count += 1
            continue

        cds_length = sum(
            end - start + 1
            for start, end in segments
        )

        model["segments"] = segments
        model["cds_length"] = cds_length

        valid_models[tx] = model

    return valid_models, invalid_count


def select_longest_per_gene(models):

    selected = {}

    for tx, model in models.items():

        gene_id = model["gene_id"]

        current = selected.get(gene_id)

        if current is None:

            selected[gene_id] = (tx, model)

        else:

            current_tx, current_model = current

            if (
                model["cds_length"] > current_model["cds_length"]
                or (
                    model["cds_length"] == current_model["cds_length"]
                    and tx < current_tx
                )
            ):

                selected[gene_id] = (tx, model)

    return selected


def write_fasta(selected, fasta, output):

    output = Path(output)

    output.parent.mkdir(parents=True, exist_ok=True)

    temp_output = Path(str(output) + ".tmp")

    written = 0
    non_multiple_of_three = 0

    with open(temp_output, "w") as fout:

        for gene_id in sorted(selected):

            tx, model = selected[gene_id]

            sequence = b"".join(
                fasta.fetch(model["chrom"], start, end)
                for start, end in model["segments"]
            )

            if model["strand"] == "-":
                sequence = reverse_complement(sequence)

            if len(sequence) != model["cds_length"]:
                raise ValueError(
                    f"Sequence length mismatch: {tx}"
                )

            if len(sequence) % 3 != 0:
                non_multiple_of_three += 1

            fout.write(f">{tx}\n")

            for i in range(0, len(sequence), 60):

                fout.write(
                    sequence[i:i + 60].decode("ascii") + "\n"
                )

            written += 1

    temp_output.replace(output)

    return written, non_multiple_of_three


def main():

    args = parse_args()

    gtf = Path(args.gtf).resolve()
    genome = Path(args.genome).resolve()
    output = Path(args.output).resolve()

    if not gtf.is_file():
        raise FileNotFoundError(f"GTF not found: {gtf}")

    if not genome.is_file():
        raise FileNotFoundError(f"Genome FASTA not found: {genome}")

    if genome.suffix == ".gz":
        raise ValueError("Compressed genome FASTA is not supported.")

    print("========== INPUT ==========")
    print(f"GTF: {gtf}")
    print(f"Genome: {genome}")
    print(f"Output: {output}")

    fasta = FastaReader(genome)

    print(f"Genome contigs: {len(fasta.index)}")

    print("\n========== PARSING GTF CDS ==========")

    models, total_records, missing, stripped = parse_gtf(gtf)

    print(f"CDS records: {total_records}")
    print(f"CDS transcript models: {len(models)}")
    print(f"Missing gene/transcript attributes: {missing}")
    print(f"Transcript rna- prefixes removed: {stripped}")

    valid_models, invalid = validate_models(models, fasta)

    print(f"Valid CDS transcripts: {len(valid_models)}")
    print(f"Invalid transcript models: {invalid}")

    if not valid_models:
        raise RuntimeError("No valid CDS transcript models found.")

    selected = select_longest_per_gene(valid_models)

    print(f"Genes with valid CDS: {len(selected)}")

    written, non_multiple = write_fasta(
        selected,
        fasta,
        output
    )

    fasta.close()

    print("\n========== SUMMARY ==========")
    print(f"Genes selected: {len(selected)}")
    print(f"FASTA sequences written: {written}")
    print(f"CDS lengths not divisible by 3: {non_multiple}")
    print(f"Output FASTA: {output}")
    print("\nFinished successfully.")


if __name__ == "__main__":
    main()
