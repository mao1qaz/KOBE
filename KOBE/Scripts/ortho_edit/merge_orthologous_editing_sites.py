#!/usr/bin/env python3
# -*- coding: gbk -*-
"""
Merge cross-species editing sites into codon table (aligned CDS)
Generates one row per editing site with all species columns (short names)
Adds debug info for missing/mismatched rows
"""

import pandas as pd
from pathlib import Path

def load_aln_cds(file):
    aln_dict = {}
    with open(file, encoding='utf-8') as f:
        header, seq_lines = None, []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    aln_dict[header] = "".join(seq_lines)
                header, seq_lines = line[1:], []
            else:
                seq_lines.append(line)
        if header:
            aln_dict[header] = "".join(seq_lines)
    return aln_dict

def get_codon_info(aln_seq, aln_index):
    if aln_index >= len(aln_seq):
        return 0, 0, 0
    if aln_seq[aln_index] == "-":
        return "---", 0, "---"
    ungapped_count = sum(1 for nt in aln_seq[:aln_index + 1] if nt != "-")
    codon_start = (ungapped_count - 1) // 3 * 3
    ungapped_seq = aln_seq.replace("-", "")
    codon = ungapped_seq[codon_start: codon_start + 3]
    codon_phase = ((ungapped_count - 1) % 3) + 1
    return codon, codon_phase, ungapped_count

def main(work_dir_path, edit_txt_dir_path, result_dir_path):
    work_dir = Path(work_dir_path)
    txt_dir = Path(edit_txt_dir_path)
    result_dir = Path(result_dir_path)
    result_dir.mkdir(parents=True, exist_ok=True)
    aln_dir = work_dir / "aln_cds"
    besthit_file = work_dir / "besthit_matrix.tsv"
    output_file = result_dir / "merged_edit_sites_codon_table.csv"

    # ---- 加载 besthit ----
    besthit = pd.read_csv(besthit_file, sep="\t", encoding="utf-8")
    queries = besthit.iloc[:, 0].tolist()
    short_names = [col.split("_")[-1].lower() for col in besthit.columns[1:]]
    species_list = short_names

    # ---- 遍历每个编辑位点文件 ----
    all_results = []
    miss_list = []

    for txt_file in txt_dir.glob("*_candi_cds_sites.txt"):
        sp_name = txt_file.stem.replace("_candi_cds_sites", "")
        sp_short = sp_name.split("_")[-1].lower()

        df = pd.read_csv(txt_file, sep="\t", encoding="gbk", usecols=["transcript", "cds_pos", "AA_change"])

        for _, row in df.iterrows():
            transcript = str(row["transcript"]).strip()
            cds_pos = int(row["cds_pos"])
            aa_change = row["AA_change"]

            # 匹配 query
            matched_queries = [q for q in queries if q == transcript or q.split("|")[-1] == transcript]
            if not matched_queries:
                miss_list.append(f"{sp_name}\t{transcript}\t[MISS_QUERY]")
                continue

            for query in matched_queries:
                aln_file = aln_dir / f"{query}.aln.aln_cds.fasta"
                if not aln_file.exists():
                    miss_list.append(f"{sp_name}\t{transcript}\t[NO_FILE:{aln_file.name}]")
                    continue

                aln_dict = load_aln_cds(aln_file)

                key_ref = next((h for h in aln_dict if h.startswith(sp_short)), None)
                if not key_ref:
                    miss_list.append(f"{sp_name}\t{transcript}\t[NO_SEQ_IN_ALN]")
                    continue

                aln_seq_ref = aln_dict[key_ref]
                ungapped, aln_index = 0, None
                for i, nt in enumerate(aln_seq_ref):
                    if nt != "-":
                        ungapped += 1
                    if ungapped == cds_pos:
                        aln_index = i
                        break
                if aln_index is None:
                    miss_list.append(f"{sp_name}\t{transcript}\t[CDS_POS_OUT_OF_RANGE]")
                    continue

                row_dict = {"ID": f"{query}:{aln_index+1}", "AA_change": f"{sp_short}:{aa_change}"}

                # 初始化所有 species 列
                for sp in species_list:
                    row_dict[f"{sp}_cds_pos"] = "NO_HIT"
                    row_dict[f"{sp}_codon_pos"] = 0
                    row_dict[f"{sp}_codon"] = "NO_HIT"

                for sp in species_list:
                    key_sp = next((h for h in aln_dict if h.startswith(sp)), None)
                    if not key_sp:
                        continue
                    aln_seq_sp = aln_dict[key_sp]
                    codon, codon_phase, cds_pos_sp = get_codon_info(aln_seq_sp, aln_index)
                    gene_name = key_sp.split("|")[-1]
                    row_dict[f"{sp}_cds_pos"] = f"{gene_name}:{cds_pos_sp}"
                    row_dict[f"{sp}_codon_pos"] = codon_phase
                    row_dict[f"{sp}_codon"] = codon

                all_results.append(row_dict)

    # ---- 输出 CSV ----
    df_out = pd.DataFrame(all_results)
    df_out.to_csv(output_file, index=False, encoding="utf-8")
    print(f"? Done! 输出保存至: {output_file}")
    print(f"? 总输出行数: {len(df_out)}")

    # ---- 输出未成功处理行 ----
    if miss_list:
        miss_file = result_dir / "missed_rows_debug.txt"
        with open(miss_file, "w", encoding="utf-8") as f:
            for line in miss_list:
                f.write(line + "\n")
        print(f"? 有 {len(miss_list)} 行未成功处理，详情见: {miss_file}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Merge cross-species editing sites into codon table with debug info")
    parser.add_argument("--work_dir", required=True, help="Working directory containing aln_cds/ and besthit_matrix.tsv")
    parser.add_argument("--edit_txt_dir", required=True, help="Directory containing *_candi_cds_sites.txt files")
    parser.add_argument("--result_dir", required=True, help="Directory for merged editing-site result tables")
    args = parser.parse_args()
    main(args.work_dir, args.edit_txt_dir, args.result_dir)
