#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Extract one longest CDS per gene and map genomic A-to-I sites to CDS coordinates."""

import argparse
import logging
import re
import sys
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path


AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Ter",
}

CODON_TABLE = {
    "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
    "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
}

OUTPUT_COLUMNS = [
    "transcript", "cds_pos", "anno",
    "AA_pos", "AA_change", "gene_id"
]

COMPLEMENT = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn",
    "TGCAYRKMVHDBNtgcayrkmvhdbn"
)


def setup_logger(log_path):
    logger = logging.getLogger("genome_site_cds_position")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        "%Y-%m-%d %H:%M:%S"
    )

    file_handler = logging.FileHandler(
        log_path, mode="w", encoding="utf-8"
    )
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger


class FastaIndex:
    """Random-access index for an uncompressed, regularly wrapped FASTA."""

    def __init__(self, fasta_path):
        self.path = Path(fasta_path)

        if not self.path.is_file():
            raise FileNotFoundError(
                f"Genome FASTA not found: {self.path}"
            )

        if self.path.suffix.lower() == ".gz":
            raise ValueError(
                "Compressed FASTA is not supported. "
                "Provide an uncompressed FASTA file."
            )

        self.records = {}
        self._build_index()

        if not self.records:
            raise ValueError(
                f"No FASTA records found: {self.path}"
            )

        self.handle = self.path.open("rb")

    def _build_index(self):
        current = None
        length = 0
        offset = None
        line_bases = None
        line_width = None
        saw_short = False

        def save_record():
            if current is None:
                return

            if offset is None or length == 0:
                raise ValueError(
                    f"Empty FASTA sequence: {current}"
                )

            self.records[current] = {
                "length": length,
                "offset": offset,
                "line_bases": line_bases,
                "line_width": line_width,
            }

        with self.path.open("rb") as handle:
            for line_number, raw in enumerate(handle, 1):

                if raw.startswith(b">"):
                    save_record()

                    header = raw[1:].strip().split()

                    if not header:
                        raise ValueError(
                            f"Empty FASTA header at line {line_number}"
                        )

                    current = header[0].decode("utf-8")

                    if current in self.records:
                        raise ValueError(
                            f"Duplicate FASTA sequence name: {current}"
                        )

                    length = 0
                    offset = None
                    line_bases = None
                    line_width = None
                    saw_short = False
                    continue

                if not raw.strip():
                    if current is not None:
                        raise ValueError(
                            f"Blank line inside FASTA sequence "
                            f"at line {line_number}"
                        )
                    continue

                if current is None:
                    raise ValueError(
                        f"Sequence data before first FASTA header "
                        f"at line {line_number}"
                    )

                seq_line = raw.rstrip(b"\r\n")

                if any(x in seq_line for x in (b" ", b"\t")):
                    raise ValueError(
                        f"Whitespace inside FASTA sequence "
                        f"at line {line_number}"
                    )

                bases = len(seq_line)

                if not bases:
                    continue

                width = len(raw)

                if offset is None:
                    offset = handle.tell() - len(raw)
                    line_bases = bases
                    line_width = width

                else:
                    if saw_short:
                        raise ValueError(
                            f"Inconsistent FASTA line wrapping in {current}"
                        )

                    if bases < line_bases:
                        saw_short = True

                    elif bases != line_bases or width != line_width:
                        raise ValueError(
                            f"Inconsistent FASTA line wrapping in {current}"
                        )

                length += bases

        save_record()

    def fetch(self, name, start0, end0):
        if name not in self.records:
            raise KeyError(
                f"Contig not found in genome FASTA: {name}"
            )

        rec = self.records[name]

        if start0 < 0 or end0 < start0 or end0 > rec["length"]:
            raise ValueError(
                f"Invalid FASTA interval: "
                f"{name}:{start0 + 1}-{end0}"
            )

        chunks = []
        pos = start0

        while pos < end0:
            line_index, column = divmod(
                pos, rec["line_bases"]
            )

            take = min(
                end0 - pos,
                rec["line_bases"] - column
            )

            byte_offset = (
                rec["offset"]
                + line_index * rec["line_width"]
                + column
            )

            self.handle.seek(byte_offset)
            chunk = self.handle.read(take)

            if len(chunk) != take:
                raise IOError(
                    f"Incomplete FASTA read: "
                    f"{name}:{start0 + 1}-{end0}"
                )

            chunks.append(chunk.decode("ascii"))
            pos += take

        return "".join(chunks).upper()

    def close(self):
        if getattr(self, "handle", None) is not None:
            self.handle.close()
            self.handle = None


def parse_gtf_attributes(raw):
    attrs = {}

    for item in raw.strip().strip(";").split(";"):
        item = item.strip()

        if not item:
            continue

        if "=" in item and (
            " " not in item or item.index("=") < item.index(" ")
        ):
            key, value = item.split("=", 1)

        else:
            parts = item.split(None, 1)

            if len(parts) != 2:
                continue

            key, value = parts

        attrs[key.strip()] = value.strip().strip('"')

    return attrs


def merge_intervals(intervals):
    intervals = sorted(set(intervals))

    if not intervals:
        return []

    merged = []
    start, end = intervals[0]

    for next_start, next_end in intervals[1:]:

        if next_start <= end + 1:
            end = max(end, next_end)

        else:
            merged.append((start, end))
            start, end = next_start, next_end

    merged.append((start, end))

    return merged


def load_longest_cds(gtf_path):
    """Select one transcript per gene by summed non-overlapping CDS length."""

    models = {}
    gene_order = {}
    next_gene_order = 0

    with Path(gtf_path).open(
        "r", encoding="utf-8-sig"
    ) as handle:

        for line_number, line in enumerate(handle, 1):

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\r\n").split("\t")

            if len(fields) != 9:
                raise ValueError(
                    f"Invalid GTF at line {line_number}: "
                    "expected 9 columns"
                )

            if fields[2].lower() != "cds":
                continue

            attrs = parse_gtf_attributes(fields[8])

            gene_id = attrs.get("gene_id")
            transcript_id = attrs.get("transcript_id")

            if not gene_id or not transcript_id:
                raise ValueError(
                    f"CDS missing gene_id/transcript_id "
                    f"at GTF line {line_number}"
                )

            try:
                start = int(fields[3])
                end = int(fields[4])

            except ValueError as error:
                raise ValueError(
                    f"Invalid CDS coordinates at GTF line {line_number}"
                ) from error

            if start < 1 or end < start:
                raise ValueError(
                    f"Invalid CDS interval at GTF line {line_number}"
                )

            strand = fields[6]

            if strand not in {"+", "-"}:
                raise ValueError(
                    f"Invalid CDS strand at GTF line {line_number}: "
                    f"{strand}"
                )

            if gene_id not in gene_order:
                gene_order[gene_id] = next_gene_order
                next_gene_order += 1

            key = (gene_id, transcript_id)

            if key not in models:
                models[key] = {
                    "gene_id": gene_id,
                    "transcript_id": transcript_id,
                    "chrom": fields[0],
                    "strand": strand,
                    "raw_intervals": [],
                }

            model = models[key]

            if (
                model["chrom"] != fields[0]
                or model["strand"] != strand
            ):
                raise ValueError(
                    f"Transcript {transcript_id} has CDS features "
                    "on multiple contigs or strands"
                )

            model["raw_intervals"].append((start, end))

    if not models:
        raise ValueError(
            "No CDS features with gene_id/transcript_id "
            "were found in the GTF"
        )

    best = {}

    for model in models.values():

        intervals = merge_intervals(
            model.pop("raw_intervals")
        )

        model["intervals"] = intervals

        model["cds_length"] = sum(
            end - start + 1
            for start, end in intervals
        )

        gene_id = model["gene_id"]
        old = best.get(gene_id)

        if (
            old is None
            or model["cds_length"] > old["cds_length"]
            or (
                model["cds_length"] == old["cds_length"]
                and model["transcript_id"] < old["transcript_id"]
            )
        ):
            best[gene_id] = model

    selected = sorted(
        best.values(),
        key=lambda x: gene_order[x["gene_id"]]
    )

    return selected, len(models)


def read_sites(path):
    """Read two- or three-column genomic site TSV."""

    sites = []
    seen = set()
    duplicate_count = 0

    with Path(path).open(
        "r", encoding="utf-8-sig"
    ) as handle:

        for line_number, line in enumerate(handle, 1):

            if not line.strip() or line.lstrip().startswith("#"):
                continue

            raw = line.rstrip("\r\n")

            if "\t" in raw:
                fields = [
                    value.strip()
                    for value in raw.split("\t")
                ]
            else:
                fields = raw.split()

            if (
                not sites
                and not seen
                and len(fields) >= 2
            ):
                try:
                    int(fields[1])

                except ValueError:
                    if (
                        fields[0].lower()
                        in {"chr", "chrom", "chromosome", "seqid", "seqname"}
                        and fields[1].lower()
                        in {"pos", "position", "genomic_position"}
                    ):
                        continue

                    raise ValueError(
                        f"Position must be an integer "
                        f"at line {line_number}"
                    )

            if len(fields) not in {2, 3}:
                raise ValueError(
                    f"Expected chromosome position [strand] "
                    f"at line {line_number}"
                )

            chrom = fields[0]

            try:
                pos = int(fields[1])

            except ValueError as error:
                raise ValueError(
                    f"Invalid genomic position at line {line_number}: "
                    f"{fields[1]}"
                ) from error

            if pos < 1:
                raise ValueError(
                    f"Genomic position must be 1-based and positive "
                    f"at line {line_number}"
                )

            strand = None

            if len(fields) == 3:
                value = fields[2].strip()

                if value in {"+", "1"}:
                    strand = "+"

                elif value in {"-", "0"}:
                    strand = "-"

                elif value.lower() in {
                    "", ".", "2", "na", "nan",
                    "unknown", "unknown_forward_reference"
                }:
                    strand = None

                else:
                    raise ValueError(
                        f"Invalid strand at line {line_number}: {value}"
                    )

            key = (chrom, pos, strand)

            if key in seen:
                duplicate_count += 1
                continue

            seen.add(key)

            sites.append({
                "chrom": chrom,
                "position": pos,
                "input_strand": strand
            })

    if not sites:
        raise ValueError(
            f"No genomic sites found in {path}"
        )

    return sites, duplicate_count


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1].upper()


def build_cds_sequences(selected_models, genome, fasta_out):
    """Extract CDS sequences and create a per-contig interval index."""

    segments_by_chrom = defaultdict(list)

    temp = fasta_out.with_name(
        fasta_out.name + ".tmp"
    )

    try:
        with temp.open(
            "w", encoding="utf-8"
        ) as out:

            for model in selected_models:

                chrom = model["chrom"]

                if chrom not in genome.records:
                    raise KeyError(
                        f"GTF contig absent from genome FASTA: {chrom}"
                    )

                ordered = (
                    model["intervals"]
                    if model["strand"] == "+"
                    else list(reversed(model["intervals"]))
                )

                parts = []
                cursor = 1
                segments = []

                for start, end in ordered:

                    seq = genome.fetch(
                        chrom, start - 1, end
                    )

                    if model["strand"] == "-":
                        seq = reverse_complement(seq)

                    parts.append(seq)

                    length = end - start + 1

                    segments.append({
                        "start": start,
                        "end": end,
                        "cds_start": cursor,
                        "cds_end": cursor + length - 1,
                    })

                    cursor += length

                cds_seq = "".join(parts)

                if len(cds_seq) != model["cds_length"]:
                    raise RuntimeError(
                        f"CDS length mismatch for "
                        f"{model['transcript_id']}"
                    )

                model["segments"] = segments
                model["cds_sequence"] = cds_seq

                out.write(
                    f">{model['transcript_id']}\n"
                )

                for i in range(0, len(cds_seq), 60):
                    out.write(
                        cds_seq[i:i + 60] + "\n"
                    )

                for seg in segments:

                    segments_by_chrom[chrom].append({
                        **seg,
                        "gene_id": model["gene_id"],
                        "transcript_id": model["transcript_id"],
                        "strand": model["strand"],
                        "model": model,
                    })

        temp.replace(fasta_out)

    except BaseException:
        temp.unlink(missing_ok=True)
        raise

    indexes = {}

    for chrom, records in segments_by_chrom.items():

        records.sort(
            key=lambda x: (
                x["start"],
                x["end"],
                x["gene_id"],
                x["transcript_id"]
            )
        )

        starts = [x["start"] for x in records]

        prefix_max_end = []
        maximum = 0

        for record in records:
            maximum = max(maximum, record["end"])
            prefix_max_end.append(maximum)

        indexes[chrom] = {
            "records": records,
            "starts": starts,
            "prefix_max_end": prefix_max_end
        }

    return indexes


def translate_codon(codon):
    return CODON_TABLE.get(codon)


def annotate_a_to_g(cds_sequence, cds_pos):
    """Annotate a transcript-oriented A>G substitution."""

    zero = cds_pos - 1
    aa_pos = zero // 3 + 1
    codon_start = (zero // 3) * 3

    ref_codon = cds_sequence[
        codon_start:codon_start + 3
    ]

    if (
        len(ref_codon) != 3
        or any(base not in "ACGT" for base in ref_codon)
    ):
        return "unknown_variant", aa_pos, "NNN>NNN"

    offset = zero % 3

    if ref_codon[offset] != "A":
        return "non_A_reference", aa_pos, "NNN>NNN"

    alt_codon = (
        ref_codon[:offset]
        + "G"
        + ref_codon[offset + 1:]
    )

    ref_aa = translate_codon(ref_codon)
    alt_aa = translate_codon(alt_codon)

    change = f"{AA3[ref_aa]}>{AA3[alt_aa]}"

    if ref_aa == alt_aa:
        anno = "synonymous_variant"

    elif codon_start == 0 and ref_aa == "M":
        anno = "start_lost"

    elif ref_aa == "*" and alt_aa != "*":
        anno = "stop_lost"

    elif ref_aa != "*" and alt_aa == "*":
        anno = "stop_gained"

    else:
        anno = "missense_variant"

    return anno, aa_pos, change


def map_site(site, indexes):
    """Map one genomic site, enforcing explicit strand."""

    index = indexes.get(site["chrom"])

    if index is None:
        return [], "contig_without_selected_cds", 0, 0

    pos = site["position"]

    records = index["records"]

    i = bisect_right(
        index["starts"], pos
    ) - 1

    overlaps = []

    while (
        i >= 0
        and index["prefix_max_end"][i] >= pos
    ):

        seg = records[i]

        if seg["start"] <= pos <= seg["end"]:
            overlaps.append(seg)

        i -= 1

    if not overlaps:
        return [], "no_longest_cds_overlap", 0, 0

    input_strand = site["input_strand"]

    if input_strand is None:

        overlapping_strands = {
            seg["strand"]
            for seg in overlaps
        }

        if len(overlapping_strands) != 1:
            return [], "ambiguous_gtf_strand", len(overlaps), 0

        target_strand = next(
            iter(overlapping_strands)
        )

    else:
        target_strand = input_strand

        strand_matches = [
            seg for seg in overlaps
            if seg["strand"] == target_strand
        ]

        if not strand_matches:
            return [], "input_gtf_strand_mismatch", len(overlaps), 0

        overlaps = strand_matches

    matches = []
    non_a = 0

    for seg in overlaps:

        if seg["strand"] == "+":
            cds_pos = (
                seg["cds_start"]
                + pos - seg["start"]
            )

        else:
            cds_pos = (
                seg["cds_start"]
                + seg["end"] - pos
            )

        model = seg["model"]

        transcript_base = model["cds_sequence"][
            cds_pos - 1
        ]

        if transcript_base != "A":
            non_a += 1
            continue

        anno, aa_pos, aa_change = annotate_a_to_g(
            model["cds_sequence"],
            cds_pos
        )

        matches.append({
            "transcript": model["transcript_id"],
            "cds_pos": cds_pos,
            "anno": anno,
            "AA_pos": aa_pos,
            "AA_change": aa_change,
            "gene_id": model["gene_id"],
        })

    matches.sort(
        key=lambda x: (
            x["gene_id"],
            x["transcript"],
            x["cds_pos"]
        )
    )

    if matches:
        return matches, None, len(overlaps), non_a

    return [], "non_A_transcript_reference", len(overlaps), non_a


def write_candidate_table(path, rows):
    temp = path.with_name(
        path.name + ".tmp"
    )

    try:
        with temp.open(
            "w", encoding="utf-8", newline=""
        ) as out:

            out.write(
                "\t".join(OUTPUT_COLUMNS) + "\n"
            )

            for row in rows:
                out.write(
                    "\t".join(
                        str(row[col])
                        for col in OUTPUT_COLUMNS
                    ) + "\n"
                )

        temp.replace(path)

    except BaseException:
        temp.unlink(missing_ok=True)
        raise


def run(args):

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fasta_out = outdir / (
        f"{args.species}_longestcds.fasta"
    )

    sites_out = outdir / (
        f"{args.species}_candi_cds_sites.txt"
    )

    log_out = outdir / (
        f"{args.species}_log.txt"
    )

    logger = setup_logger(log_out)

    logger.info("Species: %s", args.species)
    logger.info("Sites: %s", args.sites)
    logger.info("GTF: %s", args.gtf)
    logger.info("Genome: %s", args.genome)
    logger.info("Output directory: %s", outdir)

    sites, duplicate_count = read_sites(
        args.sites
    )

    logger.info(
        "Unique input genomic sites: %d",
        len(sites)
    )

    logger.info(
        "Duplicate input rows skipped: %d",
        duplicate_count
    )

    logger.info(
        "Parsing GTF and selecting longest CDS per gene..."
    )

    selected_models, transcript_count = load_longest_cds(
        args.gtf
    )

    logger.info(
        "CDS-bearing transcripts examined: %d",
        transcript_count
    )

    logger.info(
        "Genes with a CDS: %d",
        len(selected_models)
    )

    genome = FastaIndex(args.genome)

    try:

        for site in sites:

            if site["chrom"] not in genome.records:
                raise KeyError(
                    f"Input contig absent from genome FASTA: "
                    f"{site['chrom']}"
                )

            if (
                site["position"]
                > genome.records[site["chrom"]]["length"]
            ):
                raise ValueError(
                    f"Input position exceeds contig length: "
                    f"{site['chrom']}:{site['position']}"
                )

        logger.info(
            "Genome contigs indexed: %d",
            len(genome.records)
        )

        logger.info(
            "Extracting longest CDS FASTA and building coordinate index..."
        )

        indexes = build_cds_sequences(
            selected_models,
            genome,
            fasta_out
        )

    finally:
        genome.close()

    rows = []
    mapped_site_count = 0
    non_a_count = 0

    reason_counts = defaultdict(int)
    reason_examples = defaultdict(list)

    for site in sites:

        matches, reason, cds_overlaps, non_a = map_site(
            site,
            indexes
        )

        non_a_count += non_a

        if matches:

            mapped_site_count += 1
            rows.extend(
                row for row in matches
                if row["anno"] in {"synonymous_variant", "missense_variant"}
            )

        elif reason:

            reason_counts[reason] += 1

            if len(reason_examples[reason]) < 20:

                strand_label = (
                    site["input_strand"]
                    if site["input_strand"] is not None
                    else "inferred"
                )

                reason_examples[reason].append(
                    f"{site['chrom']}:{site['position']} "
                    f"(input_strand={strand_label})"
                )

    write_candidate_table(
        sites_out,
        rows
    )

    logger.info(
        "Input sites mapped to at least one longest CDS: %d / %d",
        mapped_site_count,
        len(sites)
    )

    logger.info(
        "Candidate site-to-transcript records written: %d",
        len(rows)
    )

    logger.info(
        "Sites skipped because transcript-oriented reference was not A: %d",
        reason_counts["non_A_transcript_reference"]
    )

    logger.info(
        "Sites skipped because no selected longest CDS overlaps: %d",
        reason_counts["no_longest_cds_overlap"]
    )

    logger.info(
        "Sites skipped because explicit input strand conflicts with GTF: %d",
        reason_counts["input_gtf_strand_mismatch"]
    )

    logger.info(
        "Sites skipped because unknown strand overlaps both strands: %d",
        reason_counts["ambiguous_gtf_strand"]
    )

    logger.info(
        "Sites on contigs without selected CDS models: %d",
        reason_counts["contig_without_selected_cds"]
    )

    logger.info(
        "Overlapping CDS records with non-A transcript reference skipped: %d",
        non_a_count
    )

    for reason, examples in sorted(
        reason_examples.items()
    ):

        logger.info(
            "Skip examples [%s]: %s",
            reason,
            "; ".join(examples)
        )

    logger.info(
        "Longest-CDS FASTA: %s",
        fasta_out
    )

    logger.info(
        "Candidate CDS sites: %s",
        sites_out
    )

    logger.info(
        "Log: %s",
        log_out
    )

    logger.info(
        "CDS position convention: 1-based, spliced CDS, "
        "transcript 5'-to-3' orientation"
    )

    logger.info(
        "Input strand aliases: +/1 = plus; -/0 = minus"
    )

    logger.info(
        "Blank/./NA/2 strand is treated as unspecified and inferred from GTF"
    )

    logger.info(
        "Explicit strand conflicting with selected CDS strand is skipped"
    )

    logger.info(
        "Unknown strand overlapping both plus and minus CDS is skipped"
    )

    logger.info(
        "Candidate effect assumes transcript-oriented A>G editing"
    )


def main():

    parser = argparse.ArgumentParser(
        prog="genome_site_cds_position",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n  KOBE genome_site_cds_position \\\n    --species invicta \\\n    --sites /project/data/invicta_known_sites.tsv \\\n    --gtf /project/ref/invicta.gtf \\\n    --genome /project/ref/invicta.fa \\\n    --outdir /project/results/genome_site_cds_position \\\n    --sif /project/rna_editing.sif --bind /project",
    )

    parser.add_argument(
        "--species",
        required=True,
        metavar="STRING",
        help="Species prefix used in output filenames"
    )

    parser.add_argument(
        "--sites",
        required=True,
        metavar="FILE",
        help="TSV with chromosome, 1-based position, optional strand"
    )

    parser.add_argument(
        "--gtf",
        required=True,
        metavar="FILE",
        help="Genome annotation GTF"
    )

    parser.add_argument(
        "--genome",
        required=True,
        metavar="FILE",
        help="Uncompressed reference genome FASTA"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        metavar="DIR",
        help="Output directory"
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

    args = parser.parse_args()

    run(args)


if __name__ == "__main__":
    main()