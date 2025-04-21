from collections import defaultdict
from os import system
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import numpy as np
from pysam import AlignmentFile, AlignedSegment, index, depth, flagstat, FastqFile
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import seaborn as sns
import re


def plot_read_length_distribution(fq_path: str, out_png: str):
    read_lengths = []
    with FastqFile(fq_path) as fq:
        for seq in fq:
            read_lengths.append(len(seq.sequence))
    fq.close()
    # 画图
    x_lim = max(read_lengths) // 2000 * 2000 + 2000
    print(x_lim)
    sns.set_theme(style="darkgrid")
    plt.figure(figsize=(20, 6))
    sns.histplot(read_lengths, bins=100, color="mediumseagreen")
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
    sum_df = pd.read_csv(f"{qc_res_path}/NanoStats.txt",sep="\t",skiprows=1,header=None).head(8)
    # 取整数
    int_li = [x.replace(".0","") for x in sum_df[1]]
    qc_info_dict = dict(zip(sum_df[0],int_li))
    return qc_info_dict


def quality_check(fq_path: str, output_data_path: str) -> dict:
    qc_res_path = f"{output_data_path}/qc_summary/"
    if not Path(qc_res_path).exists():
        Path(qc_res_path).mkdir(exist_ok=1,parents=1)
    # 过滤接头
    porechop_cmd = f"porechop -i {fq_path} -o {output_data_path}/adapter_trimmed.fq --threads 24 \
         --adapter_threshold 85 --discard_middle"
    # 运行NanoFilt 过滤长度低于阈值以及低质量数据
    nanofilt_cmd = f"NanoFilt -q 10 -l 500 {output_data_path}/adapter_trimmed.fq > {output_data_path}/clean_reads.fq"
    # 运行Nanoplot 生成qc报告
    nanoplot_cmd = f"NanoPlot -t 24 --fastq {output_data_path}/clean_reads.fq -o {qc_res_path} \
        --tsv_stats"
    # system(porechop_cmd)
    # system(nanofilt_cmd)
    # system(nanoplot_cmd)
    # 绘制read len 分布图
    # plot_read_length_distribution(f"{output_data_path}/clean_reads.fq",f"{output_data_path}/read_length_distribution.png")
    # 提取qc 信息
    return extract_qc_info(qc_res_path)
     

def reference_info_from_gbk(gbk_file: str, output_data_path: str) -> int:
    # 根据gbk生成fa和储存位置信息的字典
    output_fa = f"{output_data_path}/ref.fa"
    fa_handle=open(output_fa,'w')
    ref_len = 0
    # 读取 GenBank 文件中的记录
    for record in SeqIO.parse(gbk_file, "genbank"):
        # 写入fa文件
        print(f">vector",file=fa_handle)
        print(record.seq,file=fa_handle)
        ref_len = len(record.seq)
        continue
        print(f"记录ID: {record.id}")
        print(f"序列长度: {len(record.seq)} bp")
        print(f"描述: {record.description}\n")

        # 遍历所有注释特征（feature）
        for feature in record.features:
            # if feature.type in ["CDS", "gene", "rRNA", "tRNA"]:
            if feature.type:
                # 区域位置
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand  # 1: 正链，-1: 负链
                
                # 尝试获取 gene 名字
                gene = feature.qualifiers.get("gene", ["未知"])[0]
                product = feature.qualifiers.get("product", ["未知"])[0]

                # 提取对应的序列
                seq = feature.extract(record.seq)

                print(f"类型: {feature.type}")
                print(f"位置: {start} - {end} (strand: {strand})")
                print(f"基因名: {gene}")
                print(f"功能: {product}")
                print(f"序列片段（前50bp）: {seq[:50]}...\n")
    return ref_len


def map2reference(output_data_path: str) -> None:
    # 运行minimap2
    aln_cmd_str = f"minimap2 -x map-ont -a -t 24 {output_data_path}/ref.fa {output_data_path}/clean_reads.fq | samtools sort -@ 24 -O BAM - > {output_data_path}/aln.bam"
    system(aln_cmd_str)


def plot_depth_per_base(x:list, depths: list, output_png: str) -> None:
    sns.set_theme(style="darkgrid")
    plt.figure(figsize=(16, 5))  # 宽高设置
    plt.plot(x, depths, color="dimgrey", linewidth=0.6)
    plt.fill_between(x, depths, color="dimgrey")

    plt.xlabel("POS", fontsize=12)
    plt.ylabel("Depth", fontsize=12)

    plt.title("Sequencing depth pattern map of full length", fontsize=14)

    plt.tight_layout()

    # 保存图像，去除黑边，白色背景
    plt.savefig(output_png, dpi=300, facecolor='white')
    plt.show()


def float_leave_1(num: float) -> str:
    float_li = str(num).split(".")
    return f"{float_li[0]}.{float_li[1][:1]}%"

    
def obtain_map_result(ref_len: int, bam_path: str, out_png: str) -> dict:
    # 构建索引，生成深度文件和比对率文件
    depth_path = bam_path + ".depth"
    index(bam_path,"-@","24")
    depth(bam_path,"-@","24","-o",depth_path)
    # 提取比对率
    flag_info = flagstat(bam_path,"-@","24","-O","tsv")
    match = re.search(r"(\d+\.\d+%)\s+N/A\s+mapped %", flag_info)
    map_ratio = match.group(1)
    # 深度文件处理
    depths = pd.read_csv(depth_path, sep="\t", header=None)
    pos_li = depths[1].to_numpy()
    dep_li = depths[2].to_numpy()
    # 画图
    plot_depth_per_base(pos_li, dep_li, out_png)
    # 计算Avg depth, Median depth,Coverage, Cov 30x, Cov 100x
    avg_depth = int(np.mean(dep_li))
    median_depth = int(np.median(dep_li))
    cov = float_leave_1(len(dep_li) / ref_len * 100)
    cov_30x = float_leave_1(len(dep_li > 30) / ref_len * 100)
    cov_100x = float_leave_1(len(dep_li > 100) / ref_len * 100)
    
    # 存入字典
    map_info_dict = {
        "ref_len":ref_len,
        "map_ratio":map_ratio,
        "avg_depth":avg_depth,
        "median_depth":median_depth,
        "coverage":cov,
        "cov_30x":cov_30x,
        "cov_100x":cov_100x,
    }

    return map_info_dict
    
    
def process_snp(output_dir: str) -> dict:
    snp_info_dict = {}
    # 使用medaka call SNP
    output_vcf = f"{output_dir}/var.vcf"
    snp_call_cmd = f"medaka_variant \
        -i {output_dir}/aln.bam \
        -f {output_dir}/ref.fa \
        -m r941_min_sup_g507 \
        -o {output_dir} "
    
    system(snp_call_cmd)
    return snp_info_dict
    

def annotate_snp_to_gbk(gbk_file: str, vcf_file: str, output_file: str) -> None:
    # 读取 gbk 文件
    record = SeqIO.read(gbk_file, "genbank")

    # 读取 VCF 文件中 SNP 位点（忽略 header）
    snps = []
    with open(vcf_file) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            pos = int(parts[1]) - 1  # VCF 是1-based，Biopython是0-based
            ref = parts[3]
            alt = parts[4]
            snps.append((pos, ref, alt))

    # 添加 feature 注释
    for pos, ref, alt in snps:
        feature = SeqFeature(
            location=FeatureLocation(pos, pos + 1),
            type="variation",
            qualifiers={
                "note": [f"SNP: {ref}>{alt}"],
                "ref": ref,
                "alt": alt
            }
        )
        record.features.append(feature)

    # 写出新的 gbk 文件
    with open(output_file, "w") as out_handle:
        SeqIO.write(record, out_handle, "genbank")
    

def report_generate(qc_info_dict: dict, map_info_dict: dict, snp_info_dict: dict, output_dir: str) -> None:
    # 将qc map snp 信息存入同一字典
    data_dict = {}
    data_dict["qc_info"] = qc_info_dict
    data_dict["map_info"] = map_info_dict
    data_dict["snp_info"] = snp_info_dict
    # 
    

def main() -> None:
    gbk_file = "../painted_fq/C2931XKUG0-1_c-ps232691-1_vars_pLann.gbk"
    input_fq_file = "../painted_fq/C2931XKUG0-1_c-ps232691-1.fastq"
    output_dir = "./test/"
    qc_info_dict = quality_check(input_fq_file, output_dir)
    ref_len = reference_info_from_gbk(gbk_file, output_dir)
    # map2reference(output_dir)
    map_info_dict = obtain_map_result(ref_len, f"{output_dir}/aln.bam",f"{output_dir}/depth_per_base.png")
    snp_info_dict = process_snp(output_dir)
    annotated_gbk = gbk_file.replace(".gbk","_annotated.gbk")
    annotate_snp_to_gbk(gbk_file,f"{output_dir}/var.vcf",annotated_gbk)
    report_generate(qc_info_dict, map_info_dict, snp_info_dict, output_dir)
    
    
if __name__ == "__main__":
    main()