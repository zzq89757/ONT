from itertools import count
from sys import argv
import pysam
from collections import Counter

def count_bases(bam_file: str, chrom: str, pos: int) -> Counter:
    """
    统计 BAM 文件中某个位点（1-based）的不同碱基频次。
    
    参数:
        bam_file: BAM 文件路径
        chrom: 染色体名称（如 'chr1'）
        pos: 1-based 的碱基位置

    返回:
        dict 格式，如 {'A': 10, 'C': 3, 'T': 5, 'G': 2, 'N': 1}
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    base_counts = Counter()
    # 设定过滤阈值（可调）
    MIN_MAPQ = 5
    MIN_BASEQ = 10
    # 注意 pos-1 是 0-based
    for pileup_column in bam.pileup(chrom, pos - 1, pos, truncate=True):
        if pileup_column.pos == pos - 1:
            for pileup_read in pileup_column.pileups:
                if not pileup_read.is_del and not pileup_read.is_refskip:
                    aln = pileup_read.alignment
                    if aln.mapping_quality < MIN_MAPQ:
                        continue
                    baseq = aln.query_qualities[pileup_read.query_position]
                    if baseq < MIN_BASEQ:
                        continue
                    base = pileup_read.alignment.query_sequence[pileup_read.query_position]
                    base_counts[base.upper()] += 1
    bam.close()
    return base_counts


if __name__ == "__main__":
    if len(argv) != 2:
        print("Usage: python count_base.py <position>")
        exit(1)

    # 读取 BAM 文件并统计碱基频次
    # 注意：这里的 pos 是 1-based 的
    # 例如：python count_base.py 123456
    # 将统计在染色体上位置 123456 的碱基频次
    counts = count_bases(r"/mnt/ntc_data/wayne/Project/NTC/ONT/pipeline/./C1880JSBG0-1_A01//aln.bam", "vector", int(argv[1]))

    print(counts)