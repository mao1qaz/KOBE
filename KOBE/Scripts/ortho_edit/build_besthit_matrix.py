#!/usr/bin/env python3
# -*- coding: gbk -*-
"""
自动生成 besthit_8 矩阵（增量式）

功能：
1. 自动扫描 CDS 文件夹 → 翻译蛋白质（跳过已有的）
2. 自动识别物种名（*_longestcds.fasta 前缀部分）
3. 自动收集每个物种的编辑位点（*_candi_cds_sites.txt）
4. 逐物种进行 blastp 比对，增量更新 besthit_8 矩阵
"""

import os
import argparse
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import pandas as pd


# ====================================================
# 命令行参数
# ====================================================
parser = argparse.ArgumentParser(description="增量生成 besthit_8 矩阵")
parser.add_argument("--cds_dir", required=True, help="CDS文件夹路径")
parser.add_argument("--edit_txt_dir", required=True, help="编辑位点txt文件夹路径")
parser.add_argument("--work_dir", required=True, help="工作目录（pep与输出矩阵等）")
parser.add_argument("--threads", type=int, default=8, help="Maximum DIAMOND threads")
args = parser.parse_args()
if args.threads < 1:
    parser.error("--threads must be >= 1")
threads = args.threads

cds_dir = Path(args.cds_dir)
edit_txt_dir = Path(args.edit_txt_dir)
work_dir = Path(args.work_dir)
pep_dir = work_dir / "pep"
besthit_dir = work_dir / "besthit_results"

pep_dir.mkdir(parents=True, exist_ok=True)
besthit_dir.mkdir(parents=True, exist_ok=True)

non3_log = work_dir / "non_multiple_of_3.txt"
besthit8_file = work_dir / "besthit_matrix.tsv"


# ====================================================
# Step 1. 自动获取物种列表
# ====================================================
species_full_names = sorted([
    f.stem.replace("_longestcds", "")
    for f in cds_dir.glob("*_longestcds.fasta")
])
print(f"? 检测到物种 {len(species_full_names)} 个:")
for s in species_full_names:
    print("  -", s)


# ====================================================
# Step 2. 翻译 CDS → 蛋白质
# 非 3 倍数 CDS 不再直接跳过，而是去掉末尾 1 或 2 bp 后翻译
# 同时仍然记录到 non_multiple_of_3.txt
# ====================================================
for cds_file in cds_dir.glob("*_longestcds.fasta"):
    species_full = cds_file.stem.replace("_longestcds", "")
    species_short = species_full.split("_")[-1].lower()
    pep_file = pep_dir / f"{species_short}_longestpep.fasta"

    if pep_file.exists():
        print(f"? {species_short} 已有蛋白文件，跳过翻译")
        continue

    records = []

    with open(non3_log, "a") as logf:
        for record in SeqIO.parse(cds_file, "fasta"):
            seq = str(record.seq).upper().replace("U", "T")

            orig_len = len(seq)
            remainder = orig_len % 3

            if remainder != 0:
                trimmed_len = orig_len - remainder

                logf.write(
                    f"{species_short}\t{record.id}\t{orig_len}\ttrim_last_{remainder}_bp\t{trimmed_len}\n"
                )

                seq = seq[:trimmed_len]

            # 修剪后如果太短，则无法翻译
            if len(seq) < 3:
                continue

            pep_seq = Seq(seq).translate(to_stop=False)

            record.id = f"{species_short}|{record.id}"
            record.description = ""
            record.seq = pep_seq
            records.append(record)

    if records:
        SeqIO.write(records, pep_file, "fasta")
        print(f"?? {species_short}: 翻译完成 {len(records)} 条蛋白")


# ====================================================
# Step 3. 收集编辑位点蛋白
# ====================================================
all_edit_proteins = defaultdict(set)
for txt_file in sorted(edit_txt_dir.glob("*_candi_cds_sites.txt")):
    species_short = txt_file.stem.split("_")[2].lower()
    df = pd.read_csv(txt_file, sep="\t", encoding="gbk", usecols=["transcript"])
    for pid in df["transcript"]:
        all_edit_proteins[species_short].add(str(pid))

species_order = sorted(all_edit_proteins.keys())
print(f"\n? 检测到有编辑位点的物种: {species_order}")


# ====================================================
# Step 4. blastp函数
# ====================================================
def run_blastp(query_file, db_file, out_file, threads):
    if Path(out_file).exists():
        return
    db_prefix = str(db_file).replace(".fasta", "")
    if not Path(db_prefix + ".dmnd").exists():
        subprocess.run(["diamond", "makedb", "--in", str(db_file), "-d", db_prefix, "--threads", str(threads)], check=True)
    subprocess.run([
        "diamond", "blastp",
        "-q", str(query_file),
        "-d", db_prefix,
        "-o", str(out_file),
        "-f", "6", "qseqid", "sseqid", "bitscore", "evalue",
        "--max-target-seqs", "1", "--evalue", "1e-5",
        "--threads", str(threads)
    ], check=True)


# ====================================================
# Step 5. 初始化 besthit_8 矩阵
# ====================================================
if besthit8_file.exists():
    besthit8_df = pd.read_csv(besthit8_file, sep="\t", index_col=0, encoding="gbk")
    print(f"?? 载入已有矩阵: {besthit8_file}")
else:
    besthit8_df = pd.DataFrame(columns=species_full_names)
    print(f"?? 新建空矩阵")


# ====================================================
# Step 6. 逐物种处理
# ====================================================
for sp_short in species_order:
    print(f"\n?? 处理 {sp_short} ...")
    pep_file_sp = pep_dir / f"{sp_short}_longestpep.fasta"

    # 找出还未存在于矩阵中的新蛋白
    existing = set(besthit8_df.index)
    missing_proteins = [p for p in all_edit_proteins[sp_short] if p not in existing]
    if not missing_proteins:
        print(f"? {sp_short} 所有蛋白已存在矩阵中，跳过")
        continue

    # 构建 query 文件
    query_records = [r for r in SeqIO.parse(pep_file_sp, "fasta") if r.id.split("|")[-1] in missing_proteins]
    query_file = besthit_dir / f"{sp_short}_edit_proteins.fasta"
    SeqIO.write(query_records, query_file, "fasta")

    # 针对所有物种 blastp
    for target_full in species_full_names:
        target_short = target_full.split("_")[-1].lower()
        db_file = pep_dir / f"{target_short}_longestpep.fasta"
        out_file = besthit_dir / f"{sp_short}_vs_{target_short}.tsv"
        run_blastp(query_file, db_file, out_file, threads)

        df = pd.read_csv(out_file, sep="\t", header=None, names=["query", "hit", "bitscore", "evalue"])
        hit_map = {row["query"].split("|")[-1]: row["hit"] for _, row in df.iterrows()}

        for prot in missing_proteins:
            if prot not in besthit8_df.index:
                besthit8_df.loc[prot] = ["nohit"] * len(species_full_names)
            col_idx = species_full_names.index(target_full)
            besthit8_df.iloc[besthit8_df.index.get_loc(prot), col_idx] = hit_map.get(prot, f"{target_short}|nohit")

    # 保存矩阵
    besthit8_df.to_csv(besthit8_file, sep="\t", index_label="query", encoding="gbk")
    print(f"?? 已保存矩阵至: {besthit8_file}")

print("\n?? 所有编辑物种的 besthit_8 矩阵已生成完成！")
