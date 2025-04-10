from pysam import FastqFile

import gzip
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from collections import defaultdict

def fastq_quality_terminal(fastq_file, terminal_length=500):
    """
    计算 FASTQ 文件中每个碱基位置的平均质量值。

    参数：
    - fastq_file: FASTQ 文件路径（支持 .gz 压缩格式）
    - terminal_length: 只统计前 terminal_length 个碱基和最后 terminal_length 个碱基，防止数据过长影响分析

    返回：
    - avg_qualities: 每个碱基位置的平均质量值列表
    """
    qualities = defaultdict(list)

    # 使用 pysam.FastqFile 读取 FASTQ 文件（支持 .gz）
    with FastqFile(fastq_file, "r") as fastq:
        print("🚀 解析 FASTQ 数据...")
        for record in tqdm(fastq):
            qual_scores = record.get_quality_array() # 获取质量分（Phred Score）
            # 只统计 terminal_length 以内的碱基
            if len(qual_scores) < 1000:continue
            for i in range(terminal_length):
                qualities[i].append(qual_scores[i])
            for i in range(terminal_length):
                qualities[i + terminal_length].append(qual_scores[len(qual_scores) - terminal_length + i])

    print("📊 计算平均质量...")
    avg_qualities = [np.mean(q) if q else 0 for q in qualities.values()]
    return avg_qualities

def plot_quality_terminal(avg_qualities: list, output_image_path: str):
    """
    绘制质量分布折线图
    """
    plt.figure(figsize=(20, 10))
    plt.plot(range(1, len(avg_qualities) + 1), avg_qualities, marker='o', linestyle='-', color='b')
    plt.xticks(range(0, len(avg_qualities) + 1, 50))
    plt.xlabel("Position")
    plt.ylabel("Phred Score")
    plt.title("Phred per base")
    plt.grid(True, linestyle="--", alpha=0.5)

    # 添加 Q20 和 Q30 阈值参考线
    # plt.axhline(y=20, color='r', linestyle='--', label="Q20")
    # plt.axhline(y=30, color='g', linestyle='--', label="Q30")
    plt.axhline(y=10, color='r', linestyle='--', label="Q10")
    plt.axhline(y=15, color='g', linestyle='--', label="Q10")

    # plt.legend()
    # plt.show()
    plt.savefig(output_image_path)
    

def main() -> None:
    fastq_path = "C5920SXXG0-1_241211LR167.fastq"
    # 统计每个位置的平均质量
    avg_qualities = fastq_quality_terminal(fastq_path)
    # 分别绘制总的质量分布图以及首尾200bp的局部分布图
    all_pos_image_path = fastq_path.replace("fastq","png")
    head_image_path = fastq_path.replace(".fastq","_head.png")
    tail_image_path = fastq_path.replace(".fastq","_tail.png")
    # plot_quality_profile(avg_qualities,all_pos_image_path)
    plot_quality_terminal(avg_qualities[:500],head_image_path)
    plot_quality_terminal(avg_qualities[-500:],tail_image_path)

# 🚀 运行主程序
if __name__ == "__main__":
    main()
