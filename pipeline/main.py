from collections import defaultdict
from os import system
from pathlib import Path
from Bio import SeqIO
from pysam import AlignmentFile, AlignedSegment, index, FastqFile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_read_length_distribution(fq_path: str, out_png: str):
    read_lengths = []
    with FastqFile(fq_path) as fq:
        for seq in fq:
            read_lengths.append(len(seq.sequence))
    fq.close()
    # 画图
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.histplot(read_lengths, bins=100, kde=True)
    plt.title("Reads Length Distribution")
    plt.xlabel("Read Length")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()
    

def plot_depth_per_base(bam_path: str, out_png: str):
    index(bam_path)
    bamfile = AlignmentFile(bam_path, "r")
    
    depths = []

    # 遍历所有参考序列
    for ref in bamfile.references:
        for pileupcolumn in bamfile.pileup(ref, truncate=True):
            depths.append(pileupcolumn.nsegments)

    # 画图
    plt.figure(figsize=(12, 5))
    plt.plot(depths, color="dodgerblue", linewidth=0.8)
    plt.title("Depth per Base")
    plt.xlabel("Genomic Position")
    plt.ylabel("Depth")
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
    plot_read_length_distribution(f"{output_data_path}/clean_reads.fq",f"{output_data_path}/read_length_distribution.png")
    plot_depth_per_base("/mnt/ntc_data/wayne/Project/NTC/ONT/pipeline/test/aln.sam",f"{output_data_path}/depth_per_base.png")
    # 提取qc 信息
    return extract_qc_info(qc_res_path)
     
    
def coverage_check() -> None:
    ...


def mismatch_check() -> None:
    ...
    

def reference_info_from_gbk(gbk_file: str, output_data_path: str) -> dict:
    # 根据gbk生成fa和储存位置信息的字典
    output_fa = f"{output_data_path}/ref.fa"
    fa_handle=open(output_fa,'w')
    # 读取 GenBank 文件中的记录
    for record in SeqIO.parse(gbk_file, "genbank"):
        # 写入fa文件
        print(f">vector",file=fa_handle)
        print(record.seq,file=fa_handle)
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


def map2reference(output_data_path: str) -> None:
    # 运行minimap2
    aln_cmd_str = f"minimap2 -x map-ont -a -t 24 {output_data_path}/ref.fa {output_data_path}/clean_reads.fq | samtools sort -@ 24 -O BAM - > {output_data_path}/aln.sam"
    system(aln_cmd_str)
    
    
def parse_alignment_result(bam_file_path: str, res_table_path: str) -> None:
    # 处理比对结果 统计覆盖率 SNP
    ...
    

def main() -> None:
    gbk_file = "../painted_fq/C2931XKUG0-1_c-ps232691-1_vars_pLann.gbk"
    input_fq_file = "../painted_fq/C2931XKUG0-1_c-ps232691-1.fastq"
    output_dir = "./test/"
    qc_info_dict = quality_check(input_fq_file, output_dir)
    # reference_info_from_gbk(gbk_file, output_dir)
    # map2reference(output_dir)
    
    
if __name__ == "__main__":
    main()