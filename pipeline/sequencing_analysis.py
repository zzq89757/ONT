import argparse
import re
import subprocess
import time
from os import system
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation
from docx.shared import Mm
from docxtpl import DocxTemplate, InlineImage
from pysam import index, depth, flagstat, FastqFile, VariantFile

# from other import init_progress_step, step_fail_adjust, update_progress_bar

parser = argparse.ArgumentParser()


# def init_progress_plan():
#     """
#     初始化进度条
#     """

#     progress_plan = [
#         '校验并初始化文件', '过滤接头', 'NanoFilt 过滤数据', 'NanoFilt 制作QC报告', '绘制分布图', '提取QC信息',
#         '提取GBK信息', '运行minimap2', '生成深度文件和比对率文件', '执行Clair3与Sniffles', '制作新GBK', '生成结果报告'
#     ]
#     init_progress_step(args.job_id, progress_plan)


def plot_read_length_distribution(fq_path: str, out_png: str):
    read_lengths = []
    with FastqFile(fq_path) as fq:
        for seq in fq:
            read_lengths.append(len(seq.sequence))
    fq.close()
    # 画图
    x_lim = max(read_lengths) // 2000 * 2000 + 2000
    # print(x_lim)
    sns.set_theme(style="darkgrid")
    plt.figure(figsize=(20, 6))
    sns.histplot(read_lengths, bins=100, color="darkcyan")
    plt.title("Reads Length Distribution")
    plt.xlabel("Read Length")
    # plt.xlim(0, x_lim)
    # plt.xticks(np.arange(0, x_lim, 2000))
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def extract_qc_info(qc_res_path: str) -> dict:
    # 读取nanoplot的summary文件
    sum_df = pd.read_csv(f"{qc_res_path}/NanoStats.txt", sep="\t", skiprows=1, header=None).head(8)
    # 取整数
    int_li = [x.replace(".0", "") for x in sum_df[1]]
    qc_info_dict = dict(zip(sum_df[0], int_li))
    # update_progress_bar(step_name='提取QC信息', job_id=args.job_id)

    return qc_info_dict


def quality_check(fq_path: str, output_data_path: str) -> dict:
    qc_res_path = f"{output_data_path}/qc_summary/"
    if not Path(qc_res_path).exists():
        Path(qc_res_path).mkdir(exist_ok=1, parents=1)
    # 过滤接头
    porechop_cmd = f"porechop -i {fq_path} -o {output_data_path}/adapter_trimmed.fq --threads 24 \
         --adapter_threshold 85 --discard_middle"
    # update_progress_bar(step_name='过滤接头', job_id=args.job_id)

    # 运行NanoFilt 过滤长度低于阈值以及低质量数据
    nanofilt_cmd = f"NanoFilt -q 10 -l 500 {output_data_path}/adapter_trimmed.fq > {output_data_path}/clean_reads.fq"
    # update_progress_bar(step_name='NanoFilt 过滤数据', job_id=args.job_id)

    # 运行Nanoplot 生成qc报告
    nanoplot_cmd = f"NanoPlot -t 24 --fastq {output_data_path}/clean_reads.fq -o {qc_res_path} --tsv_stats --no_static"
    if not Path(f"{output_data_path}/adapter_trimmed.fq").exists():
        system(porechop_cmd)
    if not Path(f"{output_data_path}/clean_reads.fq").exists():
        system(nanofilt_cmd)
    if not Path(f"{qc_res_path}/NanoStats.txt").exists():
        system(nanoplot_cmd)
    # update_progress_bar(step_name='NanoFilt 制作QC报告', job_id=args.job_id)

    # 绘制read len 分布图
    if not Path(f"{output_data_path}/read_length_distribution.png").exists():
        plot_read_length_distribution(f"{output_data_path}/clean_reads.fq",
                                      f"{output_data_path}/read_length_distribution.png")
    # update_progress_bar(step_name='绘制分布图', job_id=args.job_id)

    # 提取qc 信息
    return extract_qc_info(qc_res_path)


def reference_info_from_gbk(gbk_file: str, output_data_path: str) -> int:
    # 根据gbk生成fa和储存位置信息的字典
    output_fa = f"{output_data_path}/ref.fa"
    fa_handle = open(output_fa, 'w')
    ref_seq = ''
    # 读取 GenBank 文件中的记录
    for record in SeqIO.parse(gbk_file, "genbank"):
        # 写入fa文件
        print(f">vector", file=fa_handle)
        print(record.seq, file=fa_handle)
        ref_seq = str(record.seq).upper()
    fa_handle.close()
    system("samtools faidx " + output_fa)
    return ref_seq


def map2reference(output_data_path: str) -> None:
    # 运行minimap2
    aln_cmd_str = f"minimap2 -x map-ont -a -t 24 {output_data_path}/ref.fa {output_data_path}/clean_reads.fq | samtools sort -@ 24 -O BAM - > {output_data_path}/aln.bam"
    if not Path(f"{output_data_path}/aln.bam").exists():
        system(aln_cmd_str)
    # update_progress_bar(step_name='运行minimap2', job_id=args.job_id)


def plot_depth_per_base(x: list, depths: list, output_png: str) -> None:
    sns.set_theme(style="darkgrid")
    plt.figure(figsize=(16, 5))  # 宽高设置
    plt.plot(x, depths, color="darkcyan", linewidth=0.6)
    plt.fill_between(x, depths, color="darkcyan")

    plt.xlabel("POS", fontsize=12)
    plt.ylabel("Depth", fontsize=12)

    plt.title("Sequencing depth pattern map of full length", fontsize=14)

    plt.tight_layout()

    # 保存图像，去除黑边，白色背景
    plt.savefig(output_png, dpi=300, facecolor='white')
    plt.show()


def float_leave_1(num: float) -> str:
    '''生成保留一位小数的百分比格式字符串'''
    float_li = str(num).split(".")
    return f"{float_li[0]}.{float_li[1][:1]}%"


def obtain_map_result(ref_len: int, bam_path: str, out_png: str) -> dict:
    # 构建索引，生成深度文件和比对率文件
    depth_path = bam_path + ".depth"
    index(bam_path, "-@", "24")
    depth(bam_path, "-@", "24", "-o", depth_path)
    # 提取比对率
    flag_info = flagstat(bam_path, "-@", "24", "-O", "tsv")
    match = re.search(r"(\d+\.\d+%)\s+N/A\s+mapped %", flag_info)
    map_ratio_li = match.group(1).split(".")
    map_ratio = map_ratio_li[0] + "." + map_ratio_li[1][0] + "%"
    # 深度文件处理
    depths = pd.read_csv(depth_path, sep="\t", header=None)
    pos_li = depths[1].to_numpy()
    dep_li = depths[2].to_numpy()
    # 画图
    if not Path(out_png).exists():
        plot_depth_per_base(pos_li, dep_li, out_png)
    # 计算Avg depth, Median depth,Coverage, Cov 30x, Cov 100x
    avg_depth = f"{int(np.mean(dep_li)):,}"
    median_depth = f"{int(np.median(dep_li)):,}"
    cov = float_leave_1(len(dep_li) / ref_len * 100)
    cov_30x = float_leave_1((dep_li > 30).sum() / ref_len * 100)
    cov_100x = float_leave_1((dep_li > 100).sum() / ref_len * 100)

    # 存入字典
    map_info_dict = {
        "ref_len": ref_len,
        "map_ratio": map_ratio,
        "avg_depth": avg_depth,
        "median_depth": median_depth,
        "coverage": cov,
        "cov_30x": cov_30x,
        "cov_100x": cov_100x,
    }
    # update_progress_bar(step_name='生成深度文件和比对率文件', job_id=args.job_id)

    return map_info_dict


def get_basecall_model_version_id(fq_path: str) -> str:
    with open(fq_path, 'r') as fq:
        for line in fq:
            if line.startswith('@') and 'basecall_model_version_id=' in line:
                for field in line.strip().split():
                    if field.startswith('basecall_model_version_id='):
                        model_str = field.split('=')[1].replace(".", "").replace("@v", "_g")
                        medaka_model_formmat = model_str.split("_", 1)[1]
                        return medaka_model_formmat
            break  # 只检查第一条 read 的注释信息
    return "Not found"


def get_medaka_model(fq_path: str) -> list:
    # 读取fq 获取basecall_model_version_id
    basecall_model_version_id = get_basecall_model_version_id(fq_path)
    medaka_model_li = [
        "r103_fast_g507", "r103_fast_snp_g507", "r103_fast_variant_g507", "r103_hac_g507", "r103_hac_snp_g507",
        "r103_hac_variant_g507", "r103_min_high_g345", "r103_min_high_g360", "r103_prom_high_g360",
        "r103_prom_snp_g3210", "r103_prom_variant_g3210", "r103_sup_g507", "r103_sup_snp_g507", "r103_sup_variant_g507",
        "r1041_e82_260bps_fast_g632", "r1041_e82_260bps_fast_variant_g632", "r1041_e82_260bps_hac_g632",
        "r1041_e82_260bps_hac_variant_g632", "r1041_e82_260bps_sup_g632", "r1041_e82_260bps_sup_variant_g632",
        "r1041_e82_400bps_fast_g615", "r1041_e82_400bps_fast_g632", "r1041_e82_400bps_fast_variant_g615",
        "r1041_e82_400bps_fast_variant_g632", "r1041_e82_400bps_hac_g615", "r1041_e82_400bps_hac_g632",
        "r1041_e82_400bps_hac_variant_g615", "r1041_e82_400bps_hac_variant_g632", "r1041_e82_400bps_sup_g615",
        "r1041_e82_400bps_sup_variant_g615", "r104_e81_fast_g5015", "r104_e81_fast_variant_g5015", "r104_e81_hac_g5015",
        "r104_e81_hac_variant_g5015", "r104_e81_sup_g5015", "r104_e81_sup_g610", "r104_e81_sup_variant_g610",
        "r10_min_high_g303", "r10_min_high_g340", "r941_e81_fast_g514", "r941_e81_fast_variant_g514",
        "r941_e81_hac_g514", "r941_e81_hac_variant_g514", "r941_e81_sup_g514", "r941_e81_sup_variant_g514",
        "r941_min_fast_g303", "r941_min_fast_g507", "r941_min_fast_snp_g507", "r941_min_fast_variant_g507",
        "r941_min_hac_g507", "r941_min_hac_snp_g507", "r941_min_hac_variant_g507", "r941_min_high_g303",
        "r941_min_high_g330", "r941_min_high_g340_rle", "r941_min_high_g344", "r941_min_high_g351",
        "r941_min_high_g360", "r941_min_sup_g507", "r941_min_sup_snp_g507", "r941_min_sup_variant_g507",
        "r941_prom_fast_g303", "r941_prom_fast_g507", "r941_prom_fast_snp_g507", "r941_prom_fast_variant_g507",
        "r941_prom_hac_g507", "r941_prom_hac_snp_g507", "r941_prom_hac_variant_g507", "r941_prom_high_g303",
        "r941_prom_high_g330", "r941_prom_high_g344", "r941_prom_high_g360", "r941_prom_high_g4011",
        "r941_prom_snp_g303", "r941_prom_snp_g322", "r941_prom_snp_g360", "r941_prom_sup_g507",
        "r941_prom_sup_snp_g507", "r941_prom_sup_variant_g507", "r941_prom_variant_g303", "r941_prom_variant_g322",
        "r941_prom_variant_g360", "r941_sup_plant_g610", "r941_sup_plant_variant_g610"
    ]
    # 根据basecall_model_version_id选择medaka模型版本（guppy版本只能低于medaka）

    # 提取basecall_model_version_id guppy前面的部分
    bc_guppy_version = int(basecall_model_version_id.split("_g", )[-1])
    basecall_model_version_id_prefix = basecall_model_version_id.replace(f"_g{bc_guppy_version}", "")
    # 去medaka model li 寻找前缀相同 guppy版本低于basecall_model_version_id guppy的模型
    prefix_catch_li = [x for x in medaka_model_li if x.startswith(basecall_model_version_id_prefix)]
    # 若无结果 直接退出
    if not prefix_catch_li:
        exit("No medaka model found,check you data,exit!!!")
    # 若只有一个结果 直接采用
    if len(prefix_catch_li) == 1:
        medaka_model = prefix_catch_li[0]
        medaka_var_model = medaka_model.replace("_g", "_variant_g")
        return [medaka_model, medaka_var_model]
    # 若多个结果 寻找前缀相同 guppy版本低于basecall_model_version_id guppy的模型
    min_medaka_model = ""
    min_medaka_var_model = ""
    min_medaka_guppy_version = 10000
    for medaka_model in prefix_catch_li:
        medaka_var_model = medaka_model.replace("_g", "_variant_g")
        medaka_guppy_version = int(medaka_model.split("_g")[1])
        if medaka_guppy_version == bc_guppy_version:
            return [medaka_model, medaka_var_model]
        if medaka_guppy_version < min_medaka_guppy_version:
            min_medaka_guppy_version = medaka_guppy_version
            min_medaka_model = medaka_model
            min_medaka_var_model = medaka_var_model
    return [min_medaka_model, min_medaka_var_model]


def is_low_complexity(seq: str) -> int:
    # 切2-mer 统计种类
    kmers = [seq[i:i + 2] for i in range(len(seq) - 2 + 1)]
    return len(set(kmers)) < 3


def clair_mutation_classify(output_vcf: str, ref_seq: str) -> list:
    # 读取 VCF 文件中 SNP 位点（忽略 header）
    snps = []
    with VariantFile(output_vcf) as vcf:
        confidence = "High"
        for rec in vcf.fetch():
            chrom = rec.chrom
            pos = rec.pos
            ref = rec.ref
            # 过滤掉没有替代碱基的记录
            if not rec.alts: continue
            alt = ','.join(str(a) for a in rec.alts)
            gt, gq, dp, ad, af = str(rec).split("\t")[-1].split(":")
            # print(af)
            # AF 过滤
            if af.find(",") != -1:
                af_li = [float(x) for x in af.split(",")]
                if not any(x >= 0.3 for x in af_li):
                    continue
                af_li = [float_leave_1(float(x) * 100) for x in af_li]
                af = ",".join(af_li)
            else: 
                if float(af) < 0.3:
                    continue
                af = float_leave_1(float(af) * 100)
            # QUAL 过滤
            if rec.qual is not None and rec.qual < 2:
                continue

            # GQ 过滤（从 FORMAT 字段的 SAMPLE 里取）
            if gq is None or float(gq) < 3:
                continue
            # DP 过滤
            # if dp is None or int(dp) < 500:
            #     continue
            # if dp is None or int(dp) < 1000:
            #     confidence = "Low"
            # GT 过滤
            if gt == "./.":
                continue
            # type:SNP/SV
            max_alt_len = max([len(a) for a in rec.alts])
            indel_num = abs(len(ref) - max_alt_len)
            type = "SNP" if indel_num <= 50 else "SV"
            if type == "SNP":
                snps.append({"confidence": confidence, "pos": pos, "ref": ref, "alt": alt, "af": af, "type": type})
            else:
                sv_type = "Duplication" if len(ref) < max_alt_len else "Deletion"
                snps.append({"confidence": confidence, "pos": pos, "ref": ref, "alt": alt, "sv_type": sv_type,
                             "sv_len": indel_num, "af": af, "type": type, })

    # 检查SNP上下游是否为polyA和polyT
    for i, site in enumerate(snps):
        pos = site["pos"]
        up_stream_pos = (pos - 5, pos)
        down_stream_pos = (pos, pos + 5)
        up_stream_seq = ref_seq[up_stream_pos[0]:up_stream_pos[1]]
        down_stream_seq = ref_seq[down_stream_pos[0]:down_stream_pos[1]]
        if is_low_complexity(up_stream_seq) or is_low_complexity(down_stream_seq):
            snps[i]["confidence"] = "Low"
    return snps


def sniffles_mutation_classify(output_vcf: str, ref_seq: str) -> list:
    # 读取 VCF 文件中 SNP 位点（忽略 header）
    snps = []
    with VariantFile(output_vcf) as vcf:
        confidence = "High"
        for rec in vcf.fetch():
            chrom = rec.chrom
            pos = rec.pos
            ref = rec.ref
            alt = ','.join(str(a) for a in rec.alts)
            af = rec.info.get('AF', 0.2)
            precise = rec.info.get('PRECISE', 0)
            sv_len = abs(rec.info.get('SVLEN', 0))
            sv_type = rec.info.get('SVTYPE', '')
            gt, gq, dr, dv = str(rec).split("\t")[-1].split(":")

            # precise 过滤
            if not precise:
                continue
            # AF 过滤
            if af < 0.3:
                continue
            af = float_leave_1(af * 100)
            # QUAL 过滤
            if rec.qual is not None and rec.qual < 2:
                continue

            # DP 过滤
            dp = rec.info.get("SUPPORT", 0)
            # if dp is None or dp < 1000:
            #     continue

            # STDEV_LEN 过滤
            std_len = rec.info.get("STDEV_LEN", 0)
            if not std_len:
                continue
            sample = rec.samples[0]  # 假设只有一个样本
            gq = sample.get("GQ")
            if gq is None or gq < 3:
                continue
            # GT 过滤
            if gt == "./.":
                continue
            # type:SNP/SV
            type = "SNP" if sv_len <= 50 else "SV"
            if type == "SNP":
                snps.append({"confidence": confidence, "pos": pos, "ref": ref, "alt": alt, "af": af, "type": type})
            else:
                sv_type = "Duplication" if sv_type == "INS" else "Deletion"
                snps.append(
                    {"confidence": confidence, "pos": pos, "ref": ref, "alt": alt, "sv_type": sv_type, "sv_len": sv_len,
                     "af": af, "type": type, })

    # 检查SNP上下游是否为polyA和polyT
    for i, site in enumerate(snps):
        pos = site["pos"]
        up_stream_pos = (pos - 5, pos)
        down_stream_pos = (pos, pos + 5)
        up_stream_seq = ref_seq[up_stream_pos[0]:up_stream_pos[1]]
        down_stream_seq = ref_seq[down_stream_pos[0]:down_stream_pos[1]]
        if is_low_complexity(up_stream_seq) or is_low_complexity(down_stream_seq):
            snps[i]["confidence"] = "Low"
    return snps


def process_mutation(ref_seq: str, output_dir: str) -> dict:
    # 使用clair3和sniffles call SNP
    clair3_cmd = f"bash /mnt/ntc_data/wayne/Repositories/Clair3/run_clair3.sh --include_all_ctgs --no_phasing_for_fa --snp_min_af=0.03\
                    -b {output_dir}/aln.bam -f {output_dir}/ref.fa -t 24 -p ont \
                    -m /mnt/ntc_data/wayne/Repositories/Clair3/models/r1041_e82_400bps_sup_v420/ \
                    -o {output_dir}/clair"
    sniffles_cmd = f"sniffles --input {output_dir}/aln.bam --threads 24 --vcf {output_dir}/snf.vcf"
    if not Path(f"{output_dir}/clair/merge_output.vcf").exists():
        system(clair3_cmd)
        system(f"gzip -d {output_dir}/clair/merge_output.vcf.gz")
    if not Path(f"{output_dir}/snf.vcf").exists():
        system(sniffles_cmd)

    # 将mutation分类至high confidence 和 low confidence并存入字典
    variant_li = sniffles_mutation_classify(f"{output_dir}/snf.vcf", ref_seq) + clair_mutation_classify(
        f"{output_dir}/clair/merge_output.vcf", ref_seq)
    # update_progress_bar(step_name='执行Clair3与Sniffles', job_id=args.job_id)

    return variant_li


def annotate_mutation_to_gbk(gbk_file: str, snp_li: list, output_file: str) -> None:
    # 读取 gbk 文件
    record = SeqIO.read(gbk_file, "genbank")
    # 添加 feature 注释
    for snp in snp_li:
        type = snp["type"]
        if type == "SNP":
            pos = snp["pos"]
            ref = snp["ref"]
            alt = snp["alt"]
            af = snp["af"]
            confidence = snp["confidence"]
            feature = SeqFeature(
                location=SimpleLocation(pos - 1, pos),
                type="SNP",
                qualifiers={
                    "note": [f"{type}: {ref}>{alt},AF:{af}({confidence} Confidence)"],
                    "ref": ref,
                    "alt": alt
                }
            )
            record.features.append(feature)
        else:
            pos = snp["pos"]
            len = snp["sv_len"]
            sv_type = snp["sv_type"]
            af = snp["af"]
            confidence = snp["confidence"]
            feature = SeqFeature(
                location=SimpleLocation(pos - 1, pos),
                type="SV",
                qualifiers={
                    "note": [f"{sv_type}:{len}bp,AF:{af}({confidence} Confidence)"]
                }
            )
            record.features.append(feature)
    # 写出新的 gbk 文件
    with open(output_file, "w") as out_handle:
        SeqIO.write(record, out_handle, "genbank")

    # update_progress_bar(step_name='制作新GBK', job_id=args.job_id)


def image_insert(doc: DocxTemplate, data_dict: dict, output_dir: str) -> None:
    read_len_distribution_image_path = f"{output_dir}/read_length_distribution.png"
    depth_distribution_image_path = f"{output_dir}/depth_per_base.png"
    data_dict['image1'] = InlineImage(doc, read_len_distribution_image_path, width=Mm(170))
    data_dict['image2'] = InlineImage(doc, depth_distribution_image_path, width=Mm(170))


def convert_to_pdf(docx_path: str, output_dir: str = "."):
    subprocess.run([
        "libreoffice", "--headless", "--convert-to", "pdf", "--outdir",
        output_dir, docx_path
    ])
    

def update_data_dict(user_info_dict: dict, qc_info_dict: dict, map_info_dict: dict, variant_li: list) -> dict:
    data_dict = {}
    # user input message
    data_dict.update(user_info_dict)

    # 将qc map snp 信息存入同一字典
    data_dict.update(qc_info_dict)
    # 格式化qc整数数字
    number_of_bases = data_dict["number_of_bases"]
    number_of_reads = data_dict["number_of_reads"]
    median_read_length = float(data_dict["median_read_length"])
    mean_read_length = float(data_dict["mean_read_length"])
    read_length_stdev = float(data_dict["read_length_stdev"])
    n50 = float(data_dict["n50"])
    data_dict["number_of_bases"] = f"{int(number_of_bases):,}"
    data_dict["number_of_reads"] = f"{int(number_of_reads):,}"
    data_dict["median_read_length"] = f"{int(median_read_length):,}"
    data_dict["mean_read_length"] = f"{int(mean_read_length):,}"
    data_dict["read_length_stdev"] = f"{int(read_length_stdev):,}"
    data_dict["n50"] = f"{int(n50):,}"
    data_dict.update(map_info_dict)
    # 处理overview信息
    ref_len = data_dict["ref_len"]
    data_dict["ref_len"] = f"{int(ref_len):,}"
    data_dict["QC"] = "PASS" if float(qc_info_dict["median_qual"]) > 15 else "FAILED"
    data_dict["Mapping"] = "PASS" if float(map_info_dict["map_ratio"].replace("%", "")) > 80 and float(
        map_info_dict["coverage"].replace("%", "")) > 80 else "FAILED"
    # 拆分variant_li为snp_list和sv_list并计数
    data_dict["snp_list"] = [] if not variant_li else [x for x in variant_li if x["type"] == "SNP"]
    data_dict["mutation_count"] = len(data_dict["snp_list"])
    data_dict["sv_list"] = [] if not variant_li else [x for x in variant_li if x["type"] == "SV"]
    data_dict["sv_count"] = len(data_dict["sv_list"])
    print(data_dict)
    return data_dict


def report_generate(user_info_dict: dict, qc_info_dict: dict, map_info_dict: dict, variant_li: list,
                    output_dir: str) -> None:
    # 更新用于填充模版的字典信息
    data_dict = update_data_dict(user_info_dict, qc_info_dict, map_info_dict, variant_li)
    # 读取模版
    doc = DocxTemplate(r"./nanopore sequencing report-final.docx")
    # 插入图片
    image_insert(doc, data_dict, output_dir)
    # 渲染报告
    order = user_info_dict["order"]
    doc.render(data_dict)
    # 存为docx文档
    doc.save(f"{output_dir}/{order}.docx")
    # 转为PDF
    # convert(f"{output_dir}/report.docx")
    output_pdf = f"{output_dir}/{order}.pdf"
    convert_to_pdf(f"{output_dir}/{order}.docx", output_dir)

    # update_progress_bar(step_name='生成结果报告', job_id=args.job_id)


def main() -> None:
    t1 = time.time()
    user_info_dict = {
        "sample_id": args.sample_id,
        "order": args.order,
        "proposal": args.proposal,
        "proposal_name": args.proposal_name
    }
    # update_progress_bar(step_name='校验并初始化文件', job_id=args.job_id)

    # QC and record info and plot plot_read_length_distribution png
    qc_info_dict = quality_check(input_fq_file, output_dir)
    # extract info and seq from gbk
    ref_seq = reference_info_from_gbk(gbk_file, output_dir)
    ref_len = len(ref_seq)
    # map to referrence by minimap2
    map2reference(output_dir)
    # extract map info from bam file and plot depth per base png
    map_info_dict = obtain_map_result(ref_len, f"{output_dir}/aln.bam", f"{output_dir}/depth_per_base.png")
    # call varaition by sniffes(SV) and clair3(SNP)
    variant_li = process_mutation(ref_seq, output_dir)
    annotated_gbk = output_dir + "/" + Path(gbk_file).name.replace(".gb", "_annotated.gb")
    # generate annotated gbk file by mutation result
    annotate_mutation_to_gbk(gbk_file, variant_li, annotated_gbk)
    # generate docx report and convert to pdf
    report_generate(user_info_dict, qc_info_dict, map_info_dict, variant_li, output_dir)
    print(f"pipeline time cost:{round(time.time() - t1, 2)} s")


if __name__ == "__main__":
    """
    Nanopore 测序工具
    """
    parser.add_argument('-f1', dest='input_file1', required=True)  # 输入gbk文件  ./data03/99/gbk_file.gbk
    parser.add_argument('-f2', dest='input_file2', required=True)  # 输入fq文件  ./data03/99/fq_file.fastq.gz
    parser.add_argument('-dir', dest='output_dir', required=True)  # 输出文件夹  ./data03/99/result/
    parser.add_argument('-sample_id', dest='sample_id', required=True)
    parser.add_argument('-order', dest='order', required=True)
    parser.add_argument('-proposal', dest='proposal', default="-")
    parser.add_argument('-proposal_name', dest='proposal_name', default="-")
    parser.add_argument('-id', dest='job_id', default=0)  # 任务ID
    args = parser.parse_args()

    output_dir = args.output_dir
    gbk_file = args.input_file1
    input_fq_file = args.input_file2

    try:
        main()
    except Exception as e:
        # step_fail_adjust(args.job_id, e)
        print(e)
