from collections import defaultdict
from os import system
from pathlib import Path
from Bio import SeqIO
from pysam import AlignmentFile, AlignedSegment
import pandas as pd


def quality_check(fq_path: str, clean_data_path: str) -> None:
    # 运行NanoFilt 并处理报告
    system(f"gunzip -c {fq_path} | NanoFilt -q 12 -l 1000 --headcrop 50 --tailcrop 50 | gzip > clean_reads.fq.gz")
    
    
def sgRNA_detective() -> None:
    ...
    
    
def coverage_check() -> None:
    ...


def stop_codon_check() -> None:
    ...
    

def mismatch_check() -> None:
    ...
    

def map2reference(fq_path: str, reference_path: str) -> None:
    # 运行minimap2 
    ...
    
    
def parse_alignment_result(bam_file_path: str, res_table_path: str) -> None:
    # 处理比对结果 统计覆盖率 SNP
    ...
    

def main() -> None:
    quality_check("","")
    
    
if __name__ == "__main__":
    main()