
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea
from matplotlib.patches import Rectangle

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42


def parse_args():
    p = argparse.ArgumentParser(description="Plot supplementary homologous RNA editing sites.")
    p.add_argument("transcript", nargs="?", default=None, help="Representative transcript ID")
    p.add_argument("--all", action="store_true", help="Plot all transcripts in one multi-page PDF")
    p.add_argument("--site_file", required=True, help="Filtered supplementary site table")
    p.add_argument("--anno_file", required=True, help="Updated annotation CSV containing AA-change columns")
    p.add_argument("--fasta_dir", required=True, help="Directory containing *_longestcds.fasta")
    p.add_argument("--species_order", required=True, help="Comma-separated species order")
    p.add_argument("--out_dir", required=True, help="PDF output directory")
    p.add_argument("--ncols", type=int, default=8)
    return p.parse_args()


args = parse_args()
if args.all and args.transcript:
    raise SystemExit("Use either --all or a transcript ID, not both")
if not args.all and not args.transcript:
    raise SystemExit("Provide a transcript ID or use --all")
if args.ncols < 1:
    raise SystemExit("--ncols must be >= 1")

site_file = Path(args.site_file)
anno_file = Path(args.anno_file)
fasta_dir = Path(args.fasta_dir)
out_dir = Path(args.out_dir)
out_dir.mkdir(parents=True, exist_ok=True)

species_order = [x.strip() for x in args.species_order.split(",") if x.strip()]
if not species_order or len(species_order) != len(set(species_order)):
    raise SystemExit("--species_order must contain unique species names")

font_dir = Path("/opt/rna_editing/fonts")
regular_path = font_dir / "Times.TTF"
bold_path = font_dir / "Timesbd.TTF"
italic_path = font_dir / "Timesi.TTF"
bold_italic_path = font_dir / "Timesbi.TTF"
for path in (regular_path, bold_path, italic_path, bold_italic_path):
    if not path.is_file():
        raise FileNotFoundError(f"Container font file not found: {path}")
    font_manager.fontManager.addfont(str(path))
FONT_REGULAR = font_manager.FontProperties(fname=str(regular_path))
FONT_BOLD = font_manager.FontProperties(fname=str(bold_path))
FONT_ITALIC = font_manager.FontProperties(fname=str(italic_path))
FONT_BOLD_ITALIC = font_manager.FontProperties(fname=str(bold_italic_path))

# Color definitions
EDIT_BG = "#f6d8df"
UNEDITED_A_BG = "#eeeeee"
EDIT_RED = "#cd2626"

# Kobe-era Los Angeles Lakers purple and gold
RECOVERED_BG = "#B39DDB"
RECOVERED_EDGE = "#D4AF37"
###e9d53b
AA_COLOR = "#b7a7a7"
CLASS_COLORS = {"Nonsyn": "#cd2626", "Syn": "#267acd", "Mixed": "#26cd7a"}

NCOLS = args.ncols
PANEL_W = 2.35
PANEL_H_BASE = 1.05
PER_SPECIES_H = 0.52
HEADER_H = 0.95
COL_GAP = 0.10
ROW_GAP = 0.18
SIDE_MARGIN = 0.16
BOTTOM_MARGIN = 0.12


def parse_species_list(value):
    return [
        x.strip() for x in str(value).split(";")
        if x.strip() and x.strip().upper() not in {"NA", "NAN", "-"}
    ]


def format_aa_change(value):
    s = str(value).strip()
    if s.upper() in {"", "NA", "NAN", "-", "NONE"} or s.count(">") != 1:
        return None
    aa1, aa2 = (x.strip() for x in s.split(">", 1))
    return f"{aa1}>{aa2}" if aa1 and aa2 else None


def effect_label(row):
    # Classification is based on the original editome evidence only.
    edited_species = list(dict.fromkeys(parse_species_list(row.get("edited_species", ""))))
    if not edited_species:
        return "NA", ""
    changes = []
    for sp in edited_species:
        change = format_aa_change(row.get(f"{sp}_AA_change", ""))
        if change is None:
            return "NA", ""
        changes.append(change)
    unique_changes = list(dict.fromkeys(changes))
    if len(unique_changes) > 1:
        return "Mixed", ";".join(unique_changes)
    change = unique_changes[0]
    aa1, aa2 = change.split(">", 1)
    return ("Syn" if aa1 == aa2 else "Nonsyn"), change


def species_label(path):
    suffix = "_longestcds.fasta"
    stem = path.name[:-len(suffix)]
    parts = stem.split("_")
    if len(parts) < 2:
        raise ValueError(f"Cannot parse species name from {path.name}")
    return f"{parts[-2][0].upper()}. {parts[-1][:3].lower()}"


def find_species_labels(directory, requested_species):
    labels = {}
    for path in sorted(directory.glob("*_longestcds.fasta")):
        sp = path.name[:-len("_longestcds.fasta")].split("_")[-1].lower()
        if sp in labels:
            raise ValueError(f"Multiple longest-CDS FASTA files detected for species: {sp}")
        labels[sp] = species_label(path)
    missing = [sp for sp in requested_species if sp.lower() not in labels]
    if missing:
        raise ValueError("Missing *_longestcds.fasta for: " + ", ".join(missing))
    return {sp: labels[sp.lower()] for sp in requested_species}


def split_seq_pep(seq, pep, base_pos):
    seq = str(seq).strip().upper()
    aa = [x.strip() for x in str(pep).split(";")]
    try:
        base_pos = int(str(base_pos).strip())
    except (TypeError, ValueError):
        return [None] * 3, [None] * 3

    if len(seq) == 9 and len(aa) == 3 and 1 <= base_pos <= 9:
        return [(seq[:3], 1), (seq[3:6], 4), (seq[6:9], 7)], aa
    if len(seq) == 6 and len(aa) == 2 and 1 <= base_pos <= 6:
        if base_pos <= 3:
            return [None, (seq[:3], 1), (seq[3:6], 4)], [None, aa[0], aa[1]]
        return [(seq[:3], 1), (seq[3:6], 4), None], [aa[0], aa[1], None]
    if len(seq) == 3 and len(aa) == 1 and 1 <= base_pos <= 3:
        return [None, (seq, 1), None], [None, aa[0], None]
    return [None] * 3, [None] * 3


def text_box(text, fp, size, color="black"):
    prop = fp.copy()
    prop.set_size(size)
    return TextArea(text, textprops={"fontproperties": prop, "color": color})


def add_centered_text_pair(holder, x, y, kind, detail, size, transform):
    parts = [text_box(kind, FONT_REGULAR, size, CLASS_COLORS.get(kind, "black"))]
    if detail:
        parts.append(text_box(" " + detail, FONT_REGULAR, size, "black"))
    packed = HPacker(children=parts, align="center", pad=0, sep=0)
    holder.add_artist(AnnotationBbox(
        packed, (x, y), xycoords=transform,
        box_alignment=(0.5, 1), frameon=False, pad=0
    ))


def draw_placeholder(ax, state, x, y, width, height):
    if state == "GAP":
        codon_gap = 0.025
        codon_w = (width - 2 * codon_gap) / 3

        for codon_idx in range(3):
            codon_center = x + codon_idx * (codon_w + codon_gap) + codon_w / 2
            ax.text(codon_center, y + height / 2, "---",
                    ha="center", va="center", fontsize=11,
                    fontproperties=FONT_REGULAR, color="#8a8a8a")
    else:
        ax.add_patch(Rectangle(
            (x, y + 0.02), width, height - 0.04,
            facecolor="white", edgecolor="#8a8a8a", linewidth=0.8
        ))
        label = "NO ORTHOLOG" if state == "NO_ORTHOLOG" else "OTHER"
        ax.text(x + width / 2, y + height / 2, label,
                ha="center", va="center", fontsize=5.8,
                fontproperties=FONT_REGULAR, color="#666666")


def draw_panel(ax, row, species, labels):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    panel_h = PANEL_H_BASE + PER_SPECIES_H * len(species)

    def from_top(inches):
        return 1 - inches / panel_h

    kind, detail = effect_label(row)
    title_length = len(kind) + (len(detail) + 1 if detail else 0)
    title_size = min(11.0, max(6.0, 11.0 * 23 / max(23, title_length)))
    add_centered_text_pair(ax, 0.5, from_top(0.15), kind, detail, title_size, ax.transAxes)

    ax.text(0.5, from_top(0.38), str(row.get("transcript_pos", "NA")),
            ha="center", va="top", fontsize=10, fontproperties=FONT_BOLD)

    edited_set = set(parse_species_list(row.get("edited_species", "")))
    recovered_set = set(parse_species_list(row.get("recovered_species", "")))

    seq_x0 = 0.28
    seq_w = 0.64
    codon_gap = 0.025
    codon_w = (seq_w - 2 * codon_gap) / 3
    base_w = codon_w / 3
    base_y_offset = 0.04 / panel_h
    aa_y_offset = -0.11 / panel_h
    box_h = 0.14 / panel_h

    for i, sp in enumerate(species):
        y = from_top(0.74 + i * PER_SPECIES_H)
        ax.text(0.04, y, labels[sp], ha="left", va="center",
                fontsize=8.6, fontproperties=FONT_ITALIC)

        state = str(row.get(f"{sp}_state", "OTHER")).strip().upper()
        if state != "SEQ":
            draw_placeholder(ax, state, seq_x0, y - 0.08 / panel_h, seq_w, 0.18 / panel_h)
            continue

        raw_pos = row.get(f"{sp}_edit_base_pos", "")
        slots, pep_slots = split_seq_pep(
            row.get(f"{sp}_CDS", ""), row.get(f"{sp}_pep", ""), raw_pos
        )
        try:
            edit_pos = int(str(raw_pos).strip())
        except (TypeError, ValueError):
            edit_pos = None

        for codon_idx, slot in enumerate(slots):
            if slot is None:
                continue
            codon, start_pos = slot
            x_slot = seq_x0 + codon_idx * (codon_w + codon_gap)

            for j, base in enumerate(codon):
                base_x = x_slot + (j + 0.5) * base_w
                color = "black"

                if edit_pos == start_pos + j and base == "A":
                    if sp in recovered_set and sp not in edited_set:
                        bg, color, edge, linewidth = RECOVERED_BG, EDIT_RED, RECOVERED_EDGE, 1.3
                    elif sp in edited_set:
                        bg, color, edge, linewidth = EDIT_BG, EDIT_RED, "none", 0
                    else:
                        bg, color, edge, linewidth = UNEDITED_A_BG, "black", "none", 0

                    ax.add_patch(Rectangle(
                        (x_slot + j * base_w + 0.002, y + base_y_offset - box_h / 2),
                        base_w - 0.004, box_h,
                        facecolor=bg, edgecolor=edge, linewidth=linewidth
                    ))

                ax.text(base_x, y + base_y_offset, base,
                        ha="center", va="center", fontsize=10.3,
                        fontproperties=FONT_REGULAR, color=color)

            aa = pep_slots[codon_idx]
            if aa and aa.upper() not in {"NA", "NAN"}:
                ax.text(x_slot + codon_w / 2, y + aa_y_offset, aa,
                        ha="center", va="center", fontsize=7,
                        fontproperties=FONT_REGULAR, color=AA_COLOR)


def read_table():
    if not site_file.is_file():
        raise FileNotFoundError(f"Filtered site table not found: {site_file}")
    if not anno_file.is_file():
        raise FileNotFoundError(f"Updated annotation CSV not found: {anno_file}")

    df = pd.read_csv(site_file, sep="\t", dtype=str, keep_default_na=False)
    required = ["ID2", "edited_species", "recovered_species", "source_species", "transcript", "transcript_pos"]
    for sp in species_order:
        required += [f"{sp}_state", f"{sp}_CDS", f"{sp}_pep", f"{sp}_edit_base_pos"]
    missing = [x for x in required if x not in df.columns]
    if missing:
        raise ValueError("Filtered site table is missing columns: " + ", ".join(missing))

    try:
        anno = pd.read_csv(anno_file, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    except UnicodeDecodeError:
        anno = pd.read_csv(anno_file, dtype=str, keep_default_na=False, encoding="gbk")
    if "ID2" not in anno.columns:
        raise ValueError("Updated annotation CSV is missing ID2")

    aa_cols = [c for c in anno.columns if c.endswith("_AA_change")]
    if aa_cols:
        anno = anno[["ID2"] + aa_cols].drop_duplicates(subset="ID2", keep="first")
        missing_aa = [c for c in aa_cols if c not in df.columns]
        if missing_aa:
            df = df.merge(
                anno[["ID2"] + missing_aa],
                on="ID2", how="left", validate="many_to_one"
            )

    df["_sort_pos"] = pd.to_numeric(df["transcript_pos"], errors="coerce")
    return df.sort_values(
        ["transcript", "_sort_pos", "ID2"], kind="stable"
    ).drop(columns="_sort_pos")


def plot_transcript(pdf, transcript, sites, labels):
    n_sites = len(sites)
    sources = list(dict.fromkeys(
        str(x).strip() for x in sites["source_species"]
        if str(x).strip().upper() not in {"", "NA", "NAN", "-"}
    ))
    source_text = "; ".join(labels.get(s, s) for s in sources) if sources else "NA"

    panel_h = PANEL_H_BASE + PER_SPECIES_H * len(species_order)
    nrows = math.ceil(n_sites / NCOLS)
    fig_w = SIDE_MARGIN * 2 + NCOLS * PANEL_W + (NCOLS - 1) * COL_GAP
    fig_h = HEADER_H + nrows * panel_h + (nrows - 1) * ROW_GAP + BOTTOM_MARGIN
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")

    header = HPacker(
        children=[text_box(source_text + " ", FONT_ITALIC, 16), text_box(transcript, FONT_BOLD, 16)],
        align="center", pad=0, sep=0
    )
    fig.add_artist(AnnotationBbox(
        header, (0.5, 1 - 0.18 / fig_h), xycoords=fig.transFigure,
        box_alignment=(0.5, 1), frameon=False, pad=0
    ))
    fig.text(0.5, 1 - 0.48 / fig_h, f"Number of sites: {n_sites}",
             ha="center", va="top", fontproperties=FONT_REGULAR, fontsize=11.5)

    for i, (_, row) in enumerate(sites.iterrows()):
        row_idx, col_idx = divmod(i, NCOLS)
        left_in = SIDE_MARGIN + col_idx * (PANEL_W + COL_GAP)
        bottom_in = fig_h - HEADER_H - panel_h - row_idx * (panel_h + ROW_GAP)
        ax = fig.add_axes([
            left_in / fig_w, bottom_in / fig_h,
            PANEL_W / fig_w, panel_h / fig_h
        ])
        draw_panel(ax, row, species_order, labels)

    pdf.savefig(fig)
    plt.close(fig)


def main():
    labels = find_species_labels(fasta_dir, species_order)
    df = read_table()

    if args.all:
        valid = ~df["transcript"].str.strip().str.upper().isin({"", "NA", "NAN", "-"})
        df = df.loc[valid].copy()
        transcripts = list(dict.fromkeys(df["transcript"]))
        if not transcripts:
            raise SystemExit("No valid transcripts found in filtered table")
        out_pdf = out_dir / f"{site_file.stem}_all_transcripts.pdf"
        with PdfPages(out_pdf) as pdf:
            for tx in transcripts:
                sub = df.loc[df["transcript"] == tx].copy()
                sub["_sort_pos"] = pd.to_numeric(sub["transcript_pos"], errors="coerce")
                sub = sub.sort_values(["_sort_pos", "ID2"], kind="stable").drop(columns="_sort_pos")
                plot_transcript(pdf, tx, sub, labels)
        print(f"Transcripts plotted: {len(transcripts)}")
    else:
        tx = args.transcript
        sub = df.loc[df["transcript"] == tx].copy()
        if sub.empty:
            raise SystemExit(f"No sites found for transcript: {tx}")
        out_pdf = out_dir / f"{tx}.pdf"
        with PdfPages(out_pdf) as pdf:
            plot_transcript(pdf, tx, sub, labels)

    print(f"PDF saved: {out_pdf}")


if __name__ == "__main__":
    main()
