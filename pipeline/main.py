from collections import defaultdict
from os import system
from pathlib import Path
from Bio import SeqIO
from pysam import AlignmentFile, AlignedSegment
import pandas as pd


def extract_qc_info(qc_res_path: str) -> None:
    ...


def quality_check(fq_path: str, output_data_path: str) -> None:
    qc_res_path = f"{output_data_path}/qc_summary/"
    if not Path(qc_res_path).exists():
        Path(qc_res_path).mkdir(exist_ok=1,parents=1)
    # 过滤接头
    porechop_cmd = f"porechop -i {fq_path} -o {output_data_path}/adapter_trimmed.fq --threads 24 \
         --adapter_threshold 85 --extra_end_trim 20 \
         --discard_middle"
    # 运行NanoFilt 过滤长度低于阈值以及低质量数据
    nanofilt_cmd = f"NanoFilt -q 10 -l 500 --headcrop 50 --tailcrop 50 {output_data_path}/adapter_trimmed.fq > {output_data_path}/clean_reads.fq"
    # 运行Nanoplot 生成qc报告
    nanoplot_cmd = f"NanoPlot -t 24 --fastq {output_data_path}/clean_reads.fq -o {qc_res_path}"
    system(porechop_cmd)
    system(nanofilt_cmd)
    system(nanoplot_cmd)
    # 提取qc 信息
    extract_qc_info(qc_res_path)
    
    
def sgRNA_detective() -> None:
    # 根据注释文件寻找label为gRNA的位置并进行检测
    ...
    
    
def coverage_check() -> None:
    ...


def stop_codon_check() -> None:
    ...
    

def mismatch_check() -> None:
    ...
    

def map2reference(fq_path: str, reference_path: str) -> None:
    # 运行minimap2 
    system()
    
    
def parse_alignment_result(bam_file_path: str, res_table_path: str) -> None:
    # 处理比对结果 统计覆盖率 SNP
    ...
    

def main() -> None:
    quality_check("../painted_fq/C2931XKUG0-1_c-ps232691-1.fastq","./test/")
    
    
if __name__ == "__main__":
    main()