
#!/bin/bash

set -Eeuo pipefail

# ============================================================
# RNA-only homologous editing site recovery
#
# 1. REDItoolKnown
# 2. Extract A/G read counts
# 3. Intersect candidate sites
# 4. Deduplicate physical CDS sites
# 5. Calculate editing levels
# 6. Summarize all samples by species
#
# Physical site key:
#   transcript + cds_pos
#
# Multiple handles:
#   Retained and separated by semicolon.
#
# Editing level:
#   G / (A + G)
#
# No P-value or FDR calculation.
# ============================================================


# ============================================================
# 1. Default settings
# ============================================================

SAMPLE_TABLE=""
OUT_DIR=""
CANDI_DIR=""
REDITool=""

PYTHON_BIN="python2"
THREADS=1
ONLY_SAMPLE=""
ONLY_SPECIES=""
REUSE_EXISTING=0

INTERSECT_PL="${INTERSECT_PL:-}"


# ============================================================
# 2. Parse arguments
# ============================================================

while [[ $# -gt 0 ]]; do

    case "$1" in

        --sample_table)
            SAMPLE_TABLE="$2"
            shift 2
            ;;

        --out_dir)
            OUT_DIR="$2"
            shift 2
            ;;

        --candi_dir)
            CANDI_DIR="$2"
            shift 2
            ;;

        --reditool)
            REDITool="$2"
            shift 2
            ;;

        --python)
            PYTHON_BIN="$2"
            shift 2
            ;;

        --threads)
            THREADS="$2"
            shift 2
            ;;

        --sample_id)
            ONLY_SAMPLE="$2"
            shift 2
            ;;

        --species)
            ONLY_SPECIES="$2"
            shift 2
            ;;

        --reuse_existing)
            REUSE_EXISTING=1
            shift
            ;;

        --intersect_pl)
            INTERSECT_PL="$2"
            shift 2
            ;;

        -h|--help)

            echo "Usage:"
            echo "  bash run_reditool_known.sh \\"
            echo "    --sample_table samples.tsv \\"
            echo "    --out_dir output \\"
            echo "    --candi_dir candi_per_species \\"
            echo "    --reditool REDItoolKnown.py \\"
            echo "    --python python2 \\"
            echo "    --intersect_pl add.intersect.vcf.info.pl \\"
            echo "    --threads 1 \\"
            echo "    [--sample_id STRING --species STRING] \\"
            echo "    [--reuse_existing]"

            exit 0
            ;;

        *)
            echo "ERROR: Unknown option: $1"
            exit 1
            ;;

    esac

done

if [[ -n "$ONLY_SAMPLE" && -z "$ONLY_SPECIES" ]] || \
   [[ -z "$ONLY_SAMPLE" && -n "$ONLY_SPECIES" ]]; then
    echo "ERROR: --sample_id and --species must be supplied together"
    exit 1
fi


# ============================================================
# 3. Validate arguments
# ============================================================

for name in SAMPLE_TABLE OUT_DIR CANDI_DIR REDITool; do

    if [[ -z "${!name}" ]]; then
        echo "ERROR: Required argument missing: $name"
        exit 1
    fi

done


[[ -f "$SAMPLE_TABLE" ]] || {
    echo "ERROR: Sample table not found: $SAMPLE_TABLE"
    exit 1
}

[[ -d "$CANDI_DIR" ]] || {
    echo "ERROR: Candidate directory not found: $CANDI_DIR"
    exit 1
}

[[ -f "$INTERSECT_PL" ]] || {
    echo "ERROR: Perl intersection script not found: $INTERSECT_PL"
    exit 1
}

if [[ ! -f "$REDITool" ]]; then
    REDITool="$(command -v "$REDITool" || true)"
fi

[[ -n "$REDITool" && -f "$REDITool" ]] || {
    echo "ERROR: REDItoolKnown.py not found"
    exit 1
}

if [[ ! -x "$PYTHON_BIN" ]]; then
    PYTHON_BIN="$(command -v "$PYTHON_BIN" || true)"
fi

[[ -n "$PYTHON_BIN" && -x "$PYTHON_BIN" ]] || {
    echo "ERROR: Python interpreter not found"
    exit 1
}

command -v python3 >/dev/null 2>&1 || {
    echo "ERROR: python3 not found"
    exit 1
}

command -v perl >/dev/null 2>&1 || {
    echo "ERROR: perl not found"
    exit 1
}

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
    echo "ERROR: THREADS must be a positive integer"
    exit 1
}

"$PYTHON_BIN" -c 'import pysam' || {
    echo "ERROR: pysam unavailable in REDItools Python"
    exit 1
}

mkdir -p "$OUT_DIR"

OUT_DIR="$(cd "$OUT_DIR" && pwd)"
CANDI_DIR="$(cd "$CANDI_DIR" && pwd)"


# ============================================================
# 4. Read sample table
# ============================================================

MANIFEST=$(mktemp)

trap 'rm -f "$MANIFEST"' EXIT

python3 - "$SAMPLE_TABLE" "$ONLY_SAMPLE" "$ONLY_SPECIES" > "$MANIFEST" <<'PY'

import csv
import sys

table = sys.argv[1]
sample_filter = sys.argv[2]
species_filter = sys.argv[3]

with open(table, encoding="utf-8-sig", newline="") as f:

    reader = csv.DictReader(f, delimiter="\t")

    if not reader.fieldnames:
        raise SystemExit("ERROR: Empty sample table")

    reader.fieldnames = [x.strip() for x in reader.fieldnames]

    def find_column(names):

        for name in names:
            if name in reader.fieldnames:
                return name

        return None

    sample_col = find_column([
        "sample_id",
        "sample",
        "sample_name"
    ])

    species_col = find_column([
        "species"
    ])

    ref_col = find_column([
        "ref",
        "reference",
        "cds_ref",
        "reference_fasta"
    ])

    type_col = find_column([
        "type"
    ])

    if not sample_col or not species_col or not ref_col:

        raise SystemExit(
            "ERROR: Required columns: sample_id, species, ref"
        )

    seen = set()
    n = 0

    for row in reader:

        sample = (row.get(sample_col) or "").strip()
        species = (row.get(species_col) or "").strip()
        ref = (row.get(ref_col) or "").strip()

        if sample_filter and sample != sample_filter:
            continue

        if species_filter and species != species_filter:
            continue

        if type_col:

            sample_type = (
                row.get(type_col) or ""
            ).strip().upper()

            if sample_type and sample_type != "RNA":
                continue

        if not sample or not species or not ref:

            raise SystemExit(
                "ERROR: Empty sample_id/species/ref"
            )

        key = (species, sample)

        if key in seen:

            raise SystemExit(
                f"ERROR: Duplicate sample: {species}/{sample}"
            )

        seen.add(key)

        print(f"{sample}\t{species}\t{ref}")

        n += 1

    if n == 0:

        raise SystemExit(
            "ERROR: No matching RNA samples found"
        )

PY


# ============================================================
# 5. Pipeline information
# ============================================================

echo "==================================================="
echo "REDItoolKnown candidate counting pipeline"
echo "==================================================="

echo "Sample table: $SAMPLE_TABLE"
echo "Output directory: $OUT_DIR"
echo "Candidate directory: $CANDI_DIR"
echo "REDItoolKnown: $REDITool"
echo "Python: $PYTHON_BIN"
echo "Threads: $THREADS"
echo "Sample: ${ONLY_SAMPLE:-ALL}"
echo "Species: ${ONLY_SPECIES:-ALL}"
echo "Reuse existing: $REUSE_EXISTING"


# ============================================================
# 6. Process samples
# ============================================================

while IFS=$'\t' read -r SAMPLE SPECIES REF || [[ -n "${SAMPLE:-}" ]]; do

    [[ -z "${SAMPLE:-}" ]] && continue

    echo ""
    echo "==================================================="
    echo "Processing: $SPECIES / $SAMPLE"
    echo "==================================================="

    SAMPLE_DIR="$OUT_DIR/Mapping_BWA_to_CDS/$SPECIES/$SAMPLE"

    BAM="$SAMPLE_DIR/$SAMPLE.sort.nodup.bam"
    SITE_FILE="$SAMPLE_DIR/knownEditingSites.tab"

    CANDI_FILE="$CANDI_DIR/${SPECIES}_candi.txt"

    REDI_OUT="$SAMPLE_DIR/known"

    MERGE_DIR="$OUT_DIR/$SPECIES/merge"
    COUNT_DIR="$OUT_DIR/candidate_counts/$SPECIES"

    mkdir -p "$MERGE_DIR" "$COUNT_DIR"


    # --------------------------------------------------------
    # 6.1 Check candidate annotation
    # --------------------------------------------------------

    if [[ ! -s "$CANDI_FILE" ]]; then

        echo "ERROR: Candidate file missing: $CANDI_FILE"
        exit 1

    fi


    # --------------------------------------------------------
    # 6.2 Check existing REDItoolKnown output
    # --------------------------------------------------------

    EXISTING_TABLES=()

    if [[ "$REUSE_EXISTING" -eq 1 && -d "$REDI_OUT" ]]; then

        mapfile -d '' EXISTING_TABLES < <(
            find "$REDI_OUT" \
                -type f \
                -name 'outTable_*' \
                -size +0c \
                -print0 | sort -z
        )

    fi


    if [[ ${#EXISTING_TABLES[@]} -gt 0 ]]; then

        echo "Step 1: Reusing existing REDItoolKnown output"

        echo "Existing outTable files: ${#EXISTING_TABLES[@]}"

        USE_EXISTING=1

    else

        USE_EXISTING=0

    fi


    # --------------------------------------------------------
    # 6.3 Run REDItoolKnown when necessary
    # --------------------------------------------------------

    if [[ "$USE_EXISTING" -eq 0 ]]; then

        echo "Step 1: Running REDItoolKnown"

        for file in "$BAM" "$SITE_FILE" "$REF"; do

            if [[ ! -s "$file" ]]; then

                echo "ERROR: Missing or empty file: $file"
                exit 1

            fi

        done


        if [[ ! -s "${BAM}.bai" && ! -s "${BAM%.bam}.bai" ]]; then

            echo "ERROR: BAM index missing: $BAM"
            exit 1

        fi


        if [[ ! -s "${REF}.fai" ]]; then

            echo "ERROR: FASTA index missing: ${REF}.fai"
            exit 1

        fi


        rm -rf "$REDI_OUT"
        mkdir -p "$REDI_OUT"

        rm -f "${SITE_FILE}.gz" "${SITE_FILE}.gz.tbi"
        "$PYTHON_BIN" "$REDITool" \
            -i "$BAM" \
            -f "$REF" \
            -l "$SITE_FILE" \
            -t "$THREADS" \
            -o "$REDI_OUT" \
            -c 1 \
            -q 25 \
            -m 25 \
            -v 0 \
            -n 0

    fi


    # --------------------------------------------------------
    # 6.4 Find outTable files
    # --------------------------------------------------------

    mapfile -d '' OUT_TABLES < <(
        find "$REDI_OUT" \
            -type f \
            -name 'outTable_*' \
            -size +0c \
            -print0 | sort -z
    )


    if [[ ${#OUT_TABLES[@]} -eq 0 ]]; then

        echo "ERROR: No non-empty outTable files found: $REDI_OUT"
        exit 1

    fi


    echo "Found outTable files: ${#OUT_TABLES[@]}"


    # --------------------------------------------------------
    # 6.5 Extract A/G counts and intersect candidates
    # --------------------------------------------------------

    echo "Step 2: Extract A/G counts and intersect candidates"

    RAW_FILE="$MERGE_DIR/${SPECIES}_raw_${SAMPLE}.txt"
    RAW_TMP="${RAW_FILE}.tmp.$$"


    {
        for file in "${OUT_TABLES[@]}"; do
            cat "$file"
        done

    } \
    | awk -F '\t' 'NF >= 7 && $2 ~ /^[0-9]+$/ {print $1 "\t" $2 "\t" $7}' \
    | perl -wpe 's/\[([0-9]+),[0-9\s]+,\s*([0-9]+),.*\n/$1___$2\n/g' \
    | perl "$INTERSECT_PL" - "$CANDI_FILE" 3 0___0 \
    > "$RAW_TMP"


    if [[ ! -s "$RAW_TMP" ]]; then

        echo "ERROR: Candidate intersection output is empty"

        rm -f "$RAW_TMP"

        exit 1

    fi


    mv "$RAW_TMP" "$RAW_FILE"


    # --------------------------------------------------------
    # 6.6 Calculate editing levels and deduplicate
    # --------------------------------------------------------

    echo "Step 3: Calculate editing levels and deduplicate physical sites"

    FINAL_FILE="$COUNT_DIR/${SAMPLE}.candidate_counts.tsv"
    FINAL_TMP="${FINAL_FILE}.tmp.$$"


    python3 - "$RAW_FILE" "$FINAL_TMP" <<'PY'

import csv
import sys
from collections import OrderedDict

raw_file = sys.argv[1]
out_file = sys.argv[2]


header = [
    "transcript",
    "cds_pos",
    "anno",
    "codon",
    "AA_before",
    "codon_pos",
    "edited_codon",
    "AA_after",
    "handle",
    "A_reads",
    "G_reads",
    "AG_coverage",
    "editing_level"
]


def get_count(value):

    value = value.strip()

    if value in ("", "-", "NA", "NaN", "."):
        return 0

    try:
        number = float(value)

    except ValueError:

        raise ValueError(
            f"Invalid read count: {value}"
        )

    if number < 0 or not number.is_integer():

        raise ValueError(
            f"Invalid integer read count: {value}"
        )

    return int(number)


def parse_ag_counts(value):

    parts = value.strip().split("___")

    if len(parts) != 2:

        raise ValueError(
            f"Invalid A/G count field: {value}"
        )

    return get_count(parts[0]), get_count(parts[1])


sites = OrderedDict()
n_raw = 0


with open(raw_file) as fin:

    for line_number, line in enumerate(fin, start=1):

        if not line.strip():
            continue

        fields = line.rstrip("\r\n").split("\t")

        if len(fields) != 10:

            raise SystemExit(
                f"ERROR: Expected 10 columns, got {len(fields)} "
                f"at line {line_number}"
            )

        try:

            a_reads, g_reads = parse_ag_counts(fields[9])

        except ValueError as e:

            raise SystemExit(
                f"ERROR at line {line_number}: {e}"
            )


        transcript = fields[0]
        cds_pos = fields[1]

        if not transcript or not cds_pos or not cds_pos.isdigit():

            raise SystemExit(
                f"ERROR: Invalid physical site at line {line_number}"
            )


        key = (transcript, cds_pos)

        annotation = fields[:8]
        handle = fields[8]

        n_raw += 1


        if key not in sites:

            sites[key] = {
                "annotation": annotation,
                "handles": [],
                "A": a_reads,
                "G": g_reads
            }

        else:

            record = sites[key]

            if record["annotation"] != annotation:

                raise SystemExit(
                    f"ERROR: Conflicting annotations for "
                    f"{transcript}:{cds_pos}"
                )

            if (
                record["A"] != a_reads or
                record["G"] != g_reads
            ):

                raise SystemExit(
                    f"ERROR: Conflicting A/G counts for "
                    f"{transcript}:{cds_pos}"
                )


        if handle and handle not in sites[key]["handles"]:

            sites[key]["handles"].append(handle)


if not sites:

    raise SystemExit("ERROR: No candidate sites generated")


multi_handle_sites = 0


with open(out_file, "w", newline="") as fout:

    writer = csv.writer(
        fout,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow(header)

    for key, record in sites.items():

        a_reads = record["A"]
        g_reads = record["G"]

        coverage = a_reads + g_reads

        if coverage > 0:

            level = f"{g_reads / coverage:.8f}"

        else:

            level = "NA"


        handles = ";".join(record["handles"])

        if len(record["handles"]) > 1:
            multi_handle_sites += 1


        writer.writerow(
            record["annotation"] + [
                handles,
                a_reads,
                g_reads,
                coverage,
                level
            ]
        )


print(f"Raw candidate rows: {n_raw}")
print(f"Unique physical sites: {len(sites)}")
print(f"Duplicated rows merged: {n_raw - len(sites)}")
print(f"Sites with multiple homologous handles: {multi_handle_sites}")
print(f"Candidate sites written: {len(sites)}")
print(f"Output: {out_file}")

PY


    mv "$FINAL_TMP" "$FINAL_FILE"


    # --------------------------------------------------------
    # 6.7 Record processed sample
    # --------------------------------------------------------

    if [[ -z "$ONLY_SAMPLE" ]]; then
        SAMPLE_LIST="$MERGE_DIR/sample.list"
        touch "$SAMPLE_LIST"

        if ! grep -Fxq "$SAMPLE" "$SAMPLE_LIST"; then
            echo "$SAMPLE" >> "$SAMPLE_LIST"
        fi
    fi


    echo ""
    echo "DONE: $SPECIES / $SAMPLE"

    echo "Raw output: $RAW_FILE"
    echo "Final output: $FINAL_FILE"

    echo "Output rows: $(($(wc -l < "$FINAL_FILE") - 1))"


done < "$MANIFEST"


# ============================================================
# 7. Multi-sample summary by species
# ============================================================

if [[ -z "$ONLY_SAMPLE" ]]; then

    echo ""
    echo "==================================================="
    echo "Generating multi-sample species summaries"
    echo "==================================================="


    python3 - "$MANIFEST" "$OUT_DIR/candidate_counts" <<'PY'

import os
import csv
import sys
from collections import OrderedDict

manifest = sys.argv[1]
count_dir = sys.argv[2]


ANNOTATION_COLS = [
    "transcript",
    "cds_pos",
    "anno",
    "codon",
    "AA_before",
    "codon_pos",
    "edited_codon",
    "AA_after",
    "handle"
]

COUNT_COLS = [
    "A_reads",
    "G_reads",
    "editing_level"
]


# ------------------------------------------------------------
# Read species and samples from manifest
# ------------------------------------------------------------

species_samples = OrderedDict()

with open(manifest) as f:

    for line in f:

        if not line.strip():
            continue

        sample, species, ref = line.rstrip("\n").split("\t")

        species_samples.setdefault(species, []).append(sample)


def parse_count(value):

    value = str(value).strip()

    if value in ("", "NA", "NaN", "-", "."):
        return 0

    number = float(value)

    if number < 0 or not number.is_integer():
        raise ValueError(f"Invalid read count: {value}")

    return int(number)


def parse_level(value):

    value = str(value).strip()

    if value in ("", "NA", "NaN", "-", "."):
        return None

    level = float(value)

    if not 0 <= level <= 1:
        raise ValueError(f"Invalid editing level: {value}")

    return level


# ------------------------------------------------------------
# Summarize each species
# ------------------------------------------------------------

for species, samples in species_samples.items():

    print("")
    print("========================================")
    print(f"Processing species: {species}")
    print(f"Samples: {len(samples)}")
    print("========================================")

    sites = OrderedDict()


    # --------------------------------------------------------
    # Read each sample
    # --------------------------------------------------------

    for sample in samples:

        input_file = os.path.join(
            count_dir,
            species,
            f"{sample}.candidate_counts.tsv"
        )

        if not os.path.isfile(input_file):

            raise SystemExit(
                f"ERROR: Missing sample result: {input_file}"
            )

        print(f"Reading: {sample}")

        sample_sites = set()

        with open(input_file, encoding="utf-8-sig", newline="") as f:

            reader = csv.DictReader(f, delimiter="\t")

            required = ANNOTATION_COLS + COUNT_COLS

            missing = [
                col for col in required
                if col not in (reader.fieldnames or [])
            ]

            if missing:

                raise SystemExit(
                    f"ERROR: Missing columns in {input_file}: {missing}"
                )


            for row in reader:

                transcript = row["transcript"]
                cds_pos = row["cds_pos"]

                key = (transcript, cds_pos)

                if key in sample_sites:

                    raise SystemExit(
                        f"ERROR: Duplicate physical site in "
                        f"{sample}: {key}"
                    )

                sample_sites.add(key)

                annotation = [
                    row[col] for col in ANNOTATION_COLS
                ]

                a_reads = parse_count(row["A_reads"])
                g_reads = parse_count(row["G_reads"])
                level = parse_level(row["editing_level"])

                coverage = a_reads + g_reads

                if coverage == 0 and level is not None:

                    raise SystemExit(
                        f"ERROR: Non-NA editing level with zero "
                        f"coverage: {key} in {sample}"
                    )

                if coverage > 0:

                    calculated_level = g_reads / coverage

                    if (
                        level is None or
                        abs(level - calculated_level) > 1e-7
                    ):

                        raise SystemExit(
                            f"ERROR: Editing level inconsistent "
                            f"with A/G reads: {key} in {sample}"
                        )


                if key not in sites:

                    sites[key] = {
                        "annotation": annotation,
                        "samples": {}
                    }

                elif sites[key]["annotation"] != annotation:

                    raise SystemExit(
                        f"ERROR: Conflicting annotations for {key}"
                    )


                sites[key]["samples"][sample] = {
                    "A_reads": a_reads,
                    "G_reads": g_reads,
                    "editing_level": level
                }


        print(f"  Sites: {len(sample_sites)}")


    if not sites:

        raise SystemExit(
            f"ERROR: No sites found for species {species}"
        )


    # --------------------------------------------------------
    # Construct summary header
    # --------------------------------------------------------

    header = ANNOTATION_COLS.copy()

    for sample in samples:

        header.extend([
            f"{sample}_A_reads",
            f"{sample}_G_reads",
            f"{sample}_editing_level"
        ])

    header.append("n_samples_editing_level_gt0")


    # --------------------------------------------------------
    # Write summary table
    # --------------------------------------------------------

    output_file = os.path.join(
        count_dir,
        species,
        f"{species}_all_samples_summary.tsv"
    )

    total_positive = 0

    with open(output_file, "w", encoding="utf-8", newline="") as f:

        writer = csv.writer(
            f,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow(header)


        for key, record in sites.items():

            output_row = record["annotation"]

            n_positive = 0


            for sample in samples:

                sample_data = record["samples"].get(sample)

                if sample_data is None:

                    a_reads = 0
                    g_reads = 0
                    level = "NA"

                else:

                    a_reads = sample_data["A_reads"]
                    g_reads = sample_data["G_reads"]

                    level = sample_data["editing_level"]

                    if level is None:

                        level = "NA"

                    else:

                        if level > 0:
                            n_positive += 1

                        level = f"{level:.8f}"


                output_row.extend([
                    a_reads,
                    g_reads,
                    level
                ])


            output_row.append(n_positive)

            if n_positive > 0:
                total_positive += 1

            writer.writerow(output_row)


    print("----------------------------------------")
    print(f"Species: {species}")
    print(f"Samples: {len(samples)}")
    print(f"Unique physical sites: {len(sites)}")
    print(f"Sites with editing level > 0: {total_positive}")
    print(f"Output: {output_file}")


print("")
print("All species summaries completed")

PY

else

    echo ""
    echo "Single-sample run detected."
    echo "Skipping all-sample summary."

fi


# ============================================================
# 8. Final summary
# ============================================================

echo ""
echo "==================================================="
echo "REDItoolKnown pipeline completed"
echo "==================================================="

if [[ -n "$ONLY_SAMPLE" ]]; then
    echo "Completed sample: $ONLY_SPECIES / $ONLY_SAMPLE"
    echo "Raw output: $OUT_DIR/$ONLY_SPECIES/merge/${ONLY_SPECIES}_raw_${ONLY_SAMPLE}.txt"
    echo "Count output: $OUT_DIR/candidate_counts/$ONLY_SPECIES/${ONLY_SAMPLE}.candidate_counts.tsv"
else
    echo ""
    echo "Individual sample results:"
    find "$OUT_DIR/candidate_counts" -type f -name '*.candidate_counts.tsv' | sort

    echo ""
    echo "Species-level summary tables:"
    find "$OUT_DIR/candidate_counts" -type f -name '*_all_samples_summary.tsv' | sort

    echo ""
    echo "Raw intersection results:"
    find "$OUT_DIR" -type f -name '*_raw_*.txt' | sort
fi

echo ""
echo "ALL DONE"
