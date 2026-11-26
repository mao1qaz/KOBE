
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Prepare species-specific SnpEff references.

Coding genes:
    Retain the transcript with the longest CDS.

Noncoding genes:
    Retain the transcript with the longest spliced exon length.
    If exon features are unavailable, use transcript feature span.

Tie-breaking:
    Lexicographically smallest transcript_id.

The script also creates the output directory and links the genome FASTA.

SnpEff GTF output:
    Omit all gene feature rows.
    Retain selected transcript-associated features only.
    Normalize attributes with transcript_id first, gene_id second,
    and a terminal semicolon for every attribute.
"""

import argparse
import os
import re
import tempfile
from collections import defaultdict
from pathlib import Path


ATTR_RE = re.compile(
    r'(?:^|;)\s*([A-Za-z_][\w.:-]*)\s+(?:"([^"]*)"|([^;\s]+))'
)


def parse_attrs(raw):
    return {
        m.group(1): m.group(2) if m.group(2) is not None else m.group(3)
        for m in ATTR_RE.finditer(raw)
    }


def format_attrs(raw):
    """Normalize GTF attributes and place transcript_id/gene_id first."""

    attributes = parse_attrs(raw)

    ordered = []

    for key in ("transcript_id", "gene_id"):
        if key in attributes:
            ordered.append((key, attributes[key]))

    ordered.extend(
        (key, value)
        for key, value in attributes.items()
        if key not in {"gene_id", "transcript_id"}
    )

    return " ".join(
        f'{key} "{value}";'
        for key, value in ordered
    )


def covered_length(intervals):
    """Calculate union length of intervals, grouped by sequence."""
    by_seq = defaultdict(list)

    for seq, start, end in intervals:
        by_seq[seq].append((start, end))

    total = 0

    for seq_intervals in by_seq.values():
        seq_intervals.sort()
        start, end = seq_intervals[0]

        for a, b in seq_intervals[1:]:
            if a <= end + 1:
                end = max(end, b)
            else:
                total += end - start + 1
                start, end = a, b

        total += end - start + 1

    return total


def ensure_genome_link(source, destination, force=False):
    """Create a genome symlink without overwriting regular files."""

    source = source.resolve(strict=True)

    if not source.is_file() or source.stat().st_size == 0:
        raise ValueError(f"Invalid genome FASTA: {source}")

    if destination.is_symlink():
        if destination.resolve(strict=False) == source:
            print(f"[Genome] Existing correct link: {destination}")
            return

        if not force:
            raise FileExistsError(
                f"Genome link points elsewhere: {destination}. "
                "Use --force to replace it."
            )

        destination.unlink()

    elif destination.exists():
        raise FileExistsError(
            f"Destination already exists as a regular file: {destination}"
        )

    destination.symlink_to(source)
    print(f"[Genome] Linked: {destination} -> {source}")


def read_gtf(input_path):

    cds_intervals = defaultdict(list)
    exon_intervals = defaultdict(list)
    transcript_intervals = defaultdict(list)

    genes = set()
    transcripts_by_gene = defaultdict(set)

    with input_path.open(
        "r", encoding="utf-8-sig", newline=""
    ) as stream:

        for number, line in enumerate(stream, 1):

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\r\n").split("\t")

            if len(fields) != 9:
                raise ValueError(
                    f"Invalid GTF at line {number}: expected 9 columns"
                )

            attributes = parse_attrs(fields[8])

            gene = attributes.get("gene_id")
            tx = attributes.get("transcript_id")
            feature = fields[2].lower()

            if not gene:
                raise ValueError(
                    f"Missing gene_id at line {number}"
                )

            genes.add(gene)

            if tx:
                transcripts_by_gene[gene].add(tx)

            if feature not in {"cds", "exon", "transcript", "mrna"}:
                continue

            if not tx:
                raise ValueError(
                    f"{feature} missing transcript_id at line {number}"
                )

            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError as error:
                raise ValueError(
                    f"Invalid coordinates at line {number}"
                ) from error

            if start < 1 or end < start:
                raise ValueError(
                    f"Invalid interval at line {number}"
                )

            interval = (fields[0], start, end)
            key = (gene, tx)

            if feature == "cds":
                cds_intervals[key].append(interval)

            elif feature == "exon":
                exon_intervals[key].append(interval)

            else:
                transcript_intervals[key].append(interval)

    return (
        genes,
        transcripts_by_gene,
        cds_intervals,
        exon_intervals,
        transcript_intervals
    )


def select_transcripts(
    genes,
    transcripts_by_gene,
    cds_intervals,
    exon_intervals,
    transcript_intervals
):

    coding_genes = {
        gene for gene, tx in cds_intervals
    }

    noncoding_genes = genes - coding_genes

    selected = {}

    length_method = {
        "exon": 0,
        "transcript_span": 0
    }

    for gene in sorted(genes):

        candidates = []

        for tx in transcripts_by_gene[gene]:

            key = (gene, tx)

            if gene in coding_genes:

                if key not in cds_intervals:
                    continue

                length = covered_length(cds_intervals[key])

            else:

                if exon_intervals.get(key):

                    length = covered_length(exon_intervals[key])
                    length_method["exon"] += 1

                elif transcript_intervals.get(key):

                    length = covered_length(
                        transcript_intervals[key]
                    )
                    length_method["transcript_span"] += 1

                else:

                    raise ValueError(
                        f"Cannot determine transcript length: "
                        f"gene={gene}, transcript={tx}"
                    )

            candidates.append((length, tx))

        if candidates:

            candidates.sort(
                key=lambda x: (-x[0], x[1])
            )

            selected[gene] = candidates[0][1]

    return (
        selected,
        coding_genes,
        noncoding_genes,
        length_method
    )


def write_filtered_gtf(
    input_path,
    output_path,
    selected,
    coding_genes
):
    """Write selected transcript features without gene feature rows."""

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    kept = 0
    removed = 0
    removed_gene_rows = 0
    kept_coding = 0
    kept_noncoding = 0

    retained_transcripts = set()
    expected_transcripts = set(selected.items())

    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.",
        suffix=".tmp",
        dir=output_path.parent
    )

    temporary = Path(temporary_name)

    try:

        with input_path.open(
            "r", encoding="utf-8-sig", newline=""
        ) as src, os.fdopen(
            fd, "w", encoding="utf-8", newline=""
        ) as dst:

            for number, line in enumerate(src, 1):

                if not line.strip() or line.startswith("#"):
                    dst.write(line)
                    continue

                fields = line.rstrip("\r\n").split("\t")

                if len(fields) != 9:
                    raise ValueError(
                        f"Invalid GTF at line {number}: expected 9 columns"
                    )

                attributes = parse_attrs(fields[8])
                gene = attributes.get("gene_id")
                tx = attributes.get("transcript_id")
                feature = fields[2].lower()

                if feature == "gene":
                    removed += 1
                    removed_gene_rows += 1
                    continue

                if not gene:
                    raise ValueError(
                        f"Missing gene_id at line {number}"
                    )

                if not tx or selected.get(gene) != tx:
                    removed += 1
                    continue

                # Normalize all retained feature attributes.
                fields[8] = format_attrs(fields[8])

                dst.write("\t".join(fields) + "\n")

                kept += 1
                retained_transcripts.add((gene, tx))

                if gene in coding_genes:
                    kept_coding += 1
                else:
                    kept_noncoding += 1

        missing_transcripts = expected_transcripts - retained_transcripts

        if missing_transcripts:
            examples = sorted(missing_transcripts)[:10]
            raise ValueError(
                "Selected transcripts have no retained GTF feature rows. "
                f"Examples: {examples}"
            )

        if temporary.stat().st_size == 0:
            raise ValueError("Filtered GTF is empty")

        os.replace(temporary, output_path)

    except BaseException:
        temporary.unlink(missing_ok=True)
        raise

    return (
        kept,
        removed,
        kept_coding,
        kept_noncoding,
        removed_gene_rows,
        len(retained_transcripts)
    )


def run(reference_dir, output_dir, species, force=False):

    reference_dir = reference_dir.resolve()
    output_dir = output_dir.resolve()

    if not re.fullmatch(r"[A-Za-z0-9_.-]+", species):
        raise ValueError(f"Invalid species identifier: {species}")

    input_fa = reference_dir / f"{species}.fa"
    input_gtf = reference_dir / f"{species}.gtf"

    output_fa = output_dir / f"{species}.fa"
    output_gtf = output_dir / f"{species}_longestCDS.gtf"

    if not input_fa.is_file() or input_fa.stat().st_size == 0:
        raise FileNotFoundError(f"Genome FASTA not found: {input_fa}")

    if not input_gtf.is_file() or input_gtf.stat().st_size == 0:
        raise FileNotFoundError(f"GTF not found: {input_gtf}")

    if input_gtf.resolve() == output_gtf.resolve():
        raise ValueError("Input and output GTF must be different")

    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print(f"Preparing SnpEff references: {species}")
    print("=" * 60)

    # 1. Link genome FASTA.
    ensure_genome_link(
        input_fa,
        output_fa,
        force=force
    )

    # 2. Skip an existing non-empty GTF unless force is specified.
    if output_gtf.exists() and not output_gtf.is_symlink():

        if output_gtf.stat().st_size > 0 and not force:

            print(f"[GTF] Existing output found; skipping: {output_gtf}")
            print("[GTF] Use --force to regenerate.")
            return

        if not force:
            raise FileExistsError(
                f"Invalid or empty output GTF: {output_gtf}"
            )

    elif output_gtf.is_symlink():

        if not force:
            raise FileExistsError(
                f"Output GTF is a symlink: {output_gtf}. "
                "Use --force to replace it."
            )

    # 3. Read GTF and select representative transcripts.
    print("[GTF] Reading annotations...")

    (
        genes,
        transcripts_by_gene,
        cds_intervals,
        exon_intervals,
        transcript_intervals
    ) = read_gtf(input_gtf)

    (
        selected,
        coding_genes,
        noncoding_genes,
        length_method
    ) = select_transcripts(
        genes,
        transcripts_by_gene,
        cds_intervals,
        exon_intervals,
        transcript_intervals
    )

    # 4. Write selected transcript features without gene rows.
    (
        kept,
        removed,
        kept_coding,
        kept_noncoding,
        removed_gene_rows,
        retained_transcript_count
    ) = write_filtered_gtf(
        input_gtf,
        output_gtf,
        selected,
        coding_genes
    )

    coding_tx_count = len(cds_intervals)

    noncoding_tx_count = sum(
        len(transcripts_by_gene[g])
        for g in noncoding_genes
    )

    retained_noncoding_tx = sum(
        1 for g in noncoding_genes if g in selected
    )

    print()
    print("=" * 60)
    print(f"Completed: {species}")
    print("=" * 60)

    print(f"Genome: {output_fa}")
    print(f"Input GTF: {input_gtf}")
    print(f"Output GTF: {output_gtf}")

    print()
    print(f"Total genes: {len(genes)}")
    print(f"Coding genes: {len(coding_genes)}")
    print(f"CDS-bearing transcripts examined: {coding_tx_count}")
    print(
        f"Coding transcripts removed: "
        f"{coding_tx_count - len(coding_genes)}"
    )

    print()
    print(f"Noncoding genes: {len(noncoding_genes)}")
    print(f"Noncoding transcripts examined: {noncoding_tx_count}")
    print(
        f"Noncoding transcripts removed: "
        f"{noncoding_tx_count - retained_noncoding_tx}"
    )

    print()
    print("Noncoding transcript length method:")
    print(f"  Exon union length: {length_method['exon']}")
    print(
        f"  Transcript span fallback: "
        f"{length_method['transcript_span']}"
    )

    print()
    print(f"Unique selected transcripts in output: {retained_transcript_count}")
    print(f"Feature rows retained (gene rows excluded): {kept}")
    print(f"  Coding feature rows: {kept_coding}")
    print(f"  Noncoding feature rows: {kept_noncoding}")
    print(f"Gene feature rows removed: {removed_gene_rows}")
    print(f"Other feature rows removed: {removed - removed_gene_rows}")

    print()
    print("GTF attribute order: transcript_id, gene_id, other attributes")
    print("GTF attribute format: every attribute ends with a semicolon")
    print("Gene feature rows: excluded")
    print("Tie-breaking: lexicographically smallest transcript_id")
    print("Preparation completed successfully.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog="prepare_snpeff_gtf",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n  KOBE prepare_snpeff_gtf \\\n    --reference_dir /project/ref \\\n    --output_dir /project/snpeff_data \\\n    --species invicta --force \\\n    --sif /project/rna_editing.sif --bind /project",
    )

    parser.add_argument(
        "--reference_dir",
        type=Path,
        required=True,
        metavar="DIR",
        help="Directory containing species.fa and species.gtf"
    )

    parser.add_argument(
        "--output_dir",
        type=Path,
        required=True,
        metavar="DIR",
        help="SnpEff reference output directory"
    )

    parser.add_argument(
        "--species",
        required=True,
        metavar="STRING",
        help="Species identifier, e.g. invicta"
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Regenerate an existing representative GTF and genome link"
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

    run(
        args.reference_dir,
        args.output_dir,
        args.species,
        force=args.force
    )
