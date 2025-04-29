from os import system
from pathlib import Path
import time
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import numpy as np
from pysam import index, depth, flagstat, FastqFile, VariantFile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
from docxtpl import DocxTemplate, InlineImage
from docx.shared import Mm
import subprocess


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
    nanoplot_cmd = f"/mnt/ntc_data/wayne/Software/miniconda3/envs/NTC/bin/NanoPlot -t 24 --fastq {output_data_path}/clean_reads.fq -o {qc_res_path} \
        --tsv_stats --no_static"
    if not Path(f"{output_data_path}/adapter_trimmed.fq").exists():
        system(porechop_cmd)
    if not Path(f"{output_data_path}/clean_reads.fq").exists():
        system(nanofilt_cmd)
    if not Path(f"{qc_res_path}/NanoStats.txt").exists():
        system(nanoplot_cmd)
    # 绘制read len 分布图
    if not Path(f"{output_data_path}/read_length_distribution.png").exists():
        plot_read_length_distribution(f"{output_data_path}/clean_reads.fq",f"{output_data_path}/read_length_distribution.png")
    # 提取qc 信息
    return extract_qc_info(qc_res_path)
     

def reference_info_from_gbk(gbk_file: str, output_data_path: str) -> int:
    # 根据gbk生成fa和储存位置信息的字典
    output_fa = f"{output_data_path}/ref.fa"
    fa_handle=open(output_fa,'w')
    ref_seq = ''
    # 读取 GenBank 文件中的记录
    for record in SeqIO.parse(gbk_file, "genbank"):
        # 写入fa文件
        print(f">vector",file=fa_handle)
        print(record.seq,file=fa_handle)
        ref_seq = str(record.seq).upper()
    return ref_seq


def map2reference(output_data_path: str) -> None:
    # 运行minimap2
    aln_cmd_str = f"minimap2 -x map-ont -a -t 24 {output_data_path}/ref.fa {output_data_path}/clean_reads.fq | samtools sort -@ 24 -O BAM - > {output_data_path}/aln.bam"
    if not Path(f"{output_data_path}/aln.bam").exists():
        system(aln_cmd_str)


def plot_depth_per_base(x:list, depths: list, output_png: str) -> None:
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
    if not Path(out_png).exists():
        plot_depth_per_base(pos_li, dep_li, out_png)
    # 计算Avg depth, Median depth,Coverage, Cov 30x, Cov 100x
    avg_depth = f"{int(np.mean(dep_li)):,}"
    median_depth = f"{int(np.median(dep_li)):,}"
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


def get_basecall_model_version_id(fq_path: str) -> str:
    with open(fq_path, 'r') as fq:
        for line in fq:
            if line.startswith('@') and 'basecall_model_version_id=' in line:
                for field in line.strip().split():
                    if field.startswith('basecall_model_version_id='):
                        model_str = field.split('=')[1].replace(".","").replace("@v","_g")
                        medaka_model_formmat = model_str.split("_",1)[1]
                        return medaka_model_formmat
            break  # 只检查第一条 read 的注释信息
    return "Not found"


def get_medaka_model(fq_path: str) -> list:
    # 读取fq 获取basecall_model_version_id
    basecall_model_version_id = get_basecall_model_version_id(fq_path)
    medaka_model_li = [
        "r103_fast_g507", "r103_fast_snp_g507", "r103_fast_variant_g507", "r103_hac_g507", "r103_hac_snp_g507", "r103_hac_variant_g507", "r103_min_high_g345", "r103_min_high_g360", "r103_prom_high_g360", "r103_prom_snp_g3210", "r103_prom_variant_g3210", "r103_sup_g507", "r103_sup_snp_g507", "r103_sup_variant_g507", "r1041_e82_260bps_fast_g632", "r1041_e82_260bps_fast_variant_g632", "r1041_e82_260bps_hac_g632", "r1041_e82_260bps_hac_variant_g632", "r1041_e82_260bps_sup_g632", "r1041_e82_260bps_sup_variant_g632", "r1041_e82_400bps_fast_g615", "r1041_e82_400bps_fast_g632", "r1041_e82_400bps_fast_variant_g615", "r1041_e82_400bps_fast_variant_g632", "r1041_e82_400bps_hac_g615", "r1041_e82_400bps_hac_g632", "r1041_e82_400bps_hac_variant_g615", "r1041_e82_400bps_hac_variant_g632", "r1041_e82_400bps_sup_g615", "r1041_e82_400bps_sup_variant_g615", "r104_e81_fast_g5015", "r104_e81_fast_variant_g5015", "r104_e81_hac_g5015", "r104_e81_hac_variant_g5015", "r104_e81_sup_g5015", "r104_e81_sup_g610", "r104_e81_sup_variant_g610", "r10_min_high_g303", "r10_min_high_g340", "r941_e81_fast_g514", "r941_e81_fast_variant_g514", "r941_e81_hac_g514", "r941_e81_hac_variant_g514", "r941_e81_sup_g514", "r941_e81_sup_variant_g514", "r941_min_fast_g303", "r941_min_fast_g507", "r941_min_fast_snp_g507", "r941_min_fast_variant_g507", "r941_min_hac_g507", "r941_min_hac_snp_g507", "r941_min_hac_variant_g507", "r941_min_high_g303", "r941_min_high_g330", "r941_min_high_g340_rle", "r941_min_high_g344", "r941_min_high_g351", "r941_min_high_g360", "r941_min_sup_g507", "r941_min_sup_snp_g507", "r941_min_sup_variant_g507", "r941_prom_fast_g303", "r941_prom_fast_g507", "r941_prom_fast_snp_g507", "r941_prom_fast_variant_g507", "r941_prom_hac_g507", "r941_prom_hac_snp_g507", "r941_prom_hac_variant_g507", "r941_prom_high_g303", "r941_prom_high_g330", "r941_prom_high_g344", "r941_prom_high_g360", "r941_prom_high_g4011", "r941_prom_snp_g303", "r941_prom_snp_g322", "r941_prom_snp_g360", "r941_prom_sup_g507", "r941_prom_sup_snp_g507", "r941_prom_sup_variant_g507", "r941_prom_variant_g303", "r941_prom_variant_g322", "r941_prom_variant_g360", "r941_sup_plant_g610", "r941_sup_plant_variant_g610"
    ]
    # 根据basecall_model_version_id选择medaka模型版本（guppy版本只能低于medaka）

    # 提取basecall_model_version_id guppy前面的部分
    bc_guppy_version = int(basecall_model_version_id.split("_g",)[-1])
    basecall_model_version_id_prefix = basecall_model_version_id.replace(f"_g{bc_guppy_version}","")
    # 去medaka model li 寻找前缀相同 guppy版本低于basecall_model_version_id guppy的模型
    prefix_catch_li = [x for x in medaka_model_li if x.startswith(basecall_model_version_id_prefix)]
    # 若无结果 直接退出
    if not prefix_catch_li:
        exit("No medaka model found,check you data,exit!!!")
    # 若只有一个结果 直接采用
    if len(prefix_catch_li) == 1:
        medaka_model = prefix_catch_li[1]
        medaka_var_model = medaka_model.replace("_g","_variant_g")
        return [medaka_model, medaka_var_model]
    # 若多个结果 寻找前缀相同 guppy版本低于basecall_model_version_id guppy的模型
    min_medaka_model = ""
    min_medaka_var_model = ""
    min_medaka_guppy_version = 10000
    for medaka_model in prefix_catch_li:
        medaka_var_model = medaka_model.replace("_g","_variant_g")
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
    kmers = [seq[i:i+2] for i in range(len(seq) - 2 + 1)]
    return len(set(kmers)) < 3      


def mutation_classify(output_vcf: str,ref_seq: str) -> list:
    # 读取 VCF 文件中 SNP 位点（忽略 header）
    snps = []
    with VariantFile(output_vcf) as vcf:
        confidence = "High"
        for rec in vcf.fetch():
            chrom = rec.chrom
            pos = rec.pos
            ref = rec.ref
            alt = ','.join(str(a) for a in rec.alts)
            # qd = rec.info.get('QD', '.')
            # fs = rec.info.get('FS', '.')
            # mq = rec.info.get('MQ', '.')
            # af = rec.info.get('AF', '.')
            ref_dp, alt_dp = rec.info['DPS']
            # af = rec.info.get('AF', '.')
            total = ref_dp + alt_dp
            af = float_leave_1(alt_dp / total * 100) if total > 0 else 0
            # QUAL 过滤
            if rec.qual is not None and rec.qual < 2:
                continue

            # DP 过滤
            dp = rec.info.get("DP")
            if dp is None or dp < 1000:
                continue

            # GQ 过滤（从 FORMAT 字段的 SAMPLE 里取）
            sample = rec.samples[0]  # 假设只有一个样本
            gq = sample.get("GQ")
            if gq is None or gq < 3:
                continue
            # (QD < 2.0 || FS > 60.0 || MQ < 40.0)?
            # if qd < 2 or fs > 60 or mq < 40:
            #     confidence = "Low"
            # type:SNP/SV
            max_alt_len = max([len(a) for a in rec.alts])
            indel_num = abs(len(ref) - max_alt_len)
            type = "SNP" if indel_num <= 50 else "SV"
            if type == "SNP":
                snps.append({"confidence":confidence, "pos":pos, "ref":ref, "alt":alt, "af": af, "type":type})
            else:
                sv_type = "Duplication" if len(ref) < max_alt_len else "Deletion"
                snps.append({"confidence":confidence, "pos":pos, "ref":ref, "alt":alt,  "sv_type":sv_type, "sv_len":indel_num, "af": af, "type":type,})
                
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
    variant_li = []
    fq_path = f"{output_dir}/clean_reads.fq"
    model, var_model = get_medaka_model(fq_path)
    print(f"Use model:{model},{var_model}")
    # 使用medaka call SNP
    output_vcf = f"{output_dir}/medaka.annotated.vcf"
    concensus_cmd = f"medaka_consensus \
        -i {fq_path} \
        -d {output_dir}/ref.fa \
        -o {output_dir} -t 24 \
        -m {model}"
    
    variant_cmd = f"medaka_haploid_variant \
        -i {fq_path} \
        -r {output_dir}/ref.fa \
        -o {output_dir} -t 24 \
        -m {var_model}"
    if not Path(f"{output_dir}/consensus.fasta").exists():
        system(concensus_cmd)
    if not Path(f"{output_dir}/medaka.vcf").exists():
        system(variant_cmd)
    # 将mutation分类至high confidence 和 low confidence并存入字典
    variant_li = mutation_classify(output_vcf,ref_seq)
    return variant_li
    

def annotate_mutation_to_gbk(gbk_file: str, snp_li: list, output_file: str) -> None:
    # 读取 gbk 文件
    record = SeqIO.read(gbk_file, "genbank")
    # 添加 feature 注释
    for snp in snp_li:
        pos = snp["pos"]
        ref = snp["ref"]
        alt = snp["alt"]
        type = snp["type"]
        af = snp["af"]
        confidence = snp["confidence"]
        feature = SeqFeature(
            location=FeatureLocation(pos, pos + 1),
            type="variation",
            qualifiers={
                "note": [f"{type}: {ref}>{alt},AF:{af}({confidence} Confidence)"],
                "ref": ref,
                "alt": alt
            }
        )
        record.features.append(feature)
    # 写出新的 gbk 文件
    with open(output_file, "w") as out_handle:
        SeqIO.write(record, out_handle, "genbank")


def image_insert(doc: DocxTemplate, data_dict: dict, output_dir: str) -> None:
    read_len_distribution_image_path = f"{output_dir}/read_length_distribution.png"
    depth_distribution_image_path = f"{output_dir}/depth_per_base.png"
    data_dict['image1'] = InlineImage(doc,read_len_distribution_image_path, width=Mm(170))
    data_dict['image2'] = InlineImage(doc,depth_distribution_image_path, width=Mm(170))


def convert_to_pdf(docx_path: str, output_dir: str = "."):
    subprocess.run([
        "libreoffice", "--headless", "--convert-to", "pdf", "--outdir",
        output_dir, docx_path
    ])



def report_generate(qc_info_dict: dict, map_info_dict: dict, variant_li: list, output_dir: str) -> None:
    data_dict = {}
    # user input message
    sample_id = "C1413CJPG0-1"
    clone_id = "-"
    order = "-"
    proposal = "-"
    proposal_name = "-"
    data_dict["sample_id"] = sample_id
    data_dict["clone_id"] = clone_id
    data_dict["order"] = order
    data_dict["proposal"] = proposal
    data_dict["proposal_name"] = proposal_name
    # 将qc map snp 信息存入同一字典
    data_dict.update(qc_info_dict)
    # 格式化qc整数数字
    data_dict["number_of_bases"] = f"{qc_info_dict["number_of_bases"]:,}"
    data_dict["number_of_reads"] = f"{qc_info_dict["number_of_reads"]:,}"
    data_dict["median_read_length"] = f"{int(qc_info_dict["median_read_length"]):,}"
    data_dict["mean_read_length"] = f"{int(qc_info_dict["mean_read_length"]):,}"
    data_dict["read_length_stdev"] = f"{int(qc_info_dict["read_length_stdev"]):,}"
    data_dict["n50"] = f"{int(qc_info_dict["n50"]):,}"
    data_dict["ref_len"] = f"{int(qc_info_dict["ref_len"]):,}"
    data_dict.update(map_info_dict)
    data_dict["QC"] = "PASS" if float(qc_info_dict["median_qual"]) > 15 else "FAILED"
    data_dict["Mapping"] = "PASS" if float(map_info_dict["map_ratio"].replace("%","")) > 95 and float(map_info_dict["coverage"].replace("%","")) > 95 else "FAILED"
    data_dict["snp_list"] = [] if not variant_li else [x for x in variant_li if x["type"]=="SNP"]
    data_dict["mutation_count"] = len(data_dict["snp_list"])
    data_dict["sv_list"] = [] if not variant_li else [x for x in variant_li if x["type"]=="SV"]
    data_dict["sv_count"] = len(data_dict["sv_list"])
    print(data_dict)
    # 读取模版
    doc = DocxTemplate(r"./nanopore sequencing report-final.docx")
    # 插入图片
    image_insert(doc,data_dict,output_dir)
    # 渲染报告并转PDF
    doc.render(data_dict)
    doc.save(f"{output_dir}/report.docx")
    # convert(f"{output_dir}/report.docx")
    output_pdf = f"{output_dir}/report.pdf"
    convert_to_pdf(f"{output_dir}/report.docx",output_dir)
    

def main() -> None:
    t1 = time.time()
    # input file
    gbk_file = "../painted_fq/C1413CJPG0-1_1220-01-A07-B08.gbk"
    input_fq_file = "../painted_fq/C1413CJPG0-1_1220-01-A07-B08.fastq.gz"
    # output dir
    output_dir = "./test_1/"
    # QC and record info and plot plot_read_length_distribution png
    qc_info_dict = quality_check(input_fq_file, output_dir)
    # extract info and seq from gbk
    ref_seq = reference_info_from_gbk(gbk_file, output_dir)
    ref_len = len(ref_seq)
    # map to referrence by minimap2
    map2reference(output_dir)
    # extract map info from bam file and plot depth per base png
    map_info_dict = obtain_map_result(ref_len, f"{output_dir}/aln.bam",f"{output_dir}/depth_per_base.png")
    variant_li = process_mutation(ref_seq,output_dir)
    annotated_gbk = output_dir + Path(gbk_file).name.replace(".gbk","_annotated.gbk")
    annotate_mutation_to_gbk(gbk_file,variant_li,annotated_gbk)
    report_generate(qc_info_dict, map_info_dict, variant_li, output_dir)
    print(f"pipeline time cost:{round(time.time() - t1, 2)} s")
    
    
if __name__ == "__main__":
    main()
    