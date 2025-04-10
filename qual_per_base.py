from pysam import FastqFile

import gzip
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool, RLock
from sys import path
path.append("/mnt/ntc_data/wayne/Project/NTC/ONT/")
from qual_terminal import fastq_quality_terminal



def async_in_iterable_structure(fun, iterable_structure, cpus):
    def init(l):
        global lock
        lock = l

    lock = RLock()
    p = Pool(int(cpus), initializer=init, initargs=(lock,))
    # apply async in iterable structure
    for i in iterable_structure:
        p.apply_async(func=fun, args=(i,))
    p.close()
    p.join()
    

def fastq_quality_stats(fastq_file, max_length=200000):
    """
    计算 FASTQ 文件中每个碱基位置的平均质量值。

    参数：
    - fastq_file: FASTQ 文件路径（支持 .gz 压缩格式）
    - max_length: 只统计前 max_length 个碱基，防止数据过长影响分析

    返回：
    - avg_qualities: 每个碱基位置的平均质量值列表
    """
    qualities = []

    # 使用 pysam.FastqFile 读取 FASTQ 文件（支持 .gz）
    with FastqFile(fastq_file, "r") as fastq:
        print("🚀 解析 FASTQ 数据...")
        for record in tqdm(fastq):
            qual_scores = record.get_quality_array()[:max_length] # 获取质量分（Phred Score）

            # 确保列表足够长
            if len(qualities) < len(qual_scores):
                qualities.extend([[] for _ in range(len(qual_scores) - len(qualities))])

            # 只统计 max_length 以内的碱基
            for i in range(min(len(qual_scores), max_length)):
                qualities[i].append(qual_scores[i])

    print("📊 计算平均质量...")
    avg_qualities = [np.mean(q) if q else 0 for q in qualities]
    return avg_qualities

def plot_quality_profile(sample: str, avg_qualities: list, output_image_path: str):
    """
    绘制质量分布折线图
    """
    width = 20 
    if len(avg_qualities) > 2000:width = 40
    plt.figure(figsize=(20, 10))
    plt.plot(range(1, len(avg_qualities) + 1), avg_qualities, marker='o', linestyle='-', color='b')
    if len(avg_qualities) < 2000:plt.xticks(range(0, len(avg_qualities) + 1, 50))
    plt.xlabel("Position")
    plt.ylabel("Phred Score")
    plt.title(f"Phred per base in sample {sample}")
    plt.grid(True, linestyle="--", alpha=0.5)

    # 添加 Q20 和 Q30 阈值参考线
    # plt.axhline(y=20, color='r', linestyle='--', label="Q20")
    # plt.axhline(y=30, color='g', linestyle='--', label="Q30")
    plt.axhline(y=10, color='r', linestyle='--', label="Q10")
    plt.axhline(y=15, color='g', linestyle='--', label="Q15")

    # plt.legend()
    # plt.show()
    plt.savefig(output_image_path)
    

def run_plot(fastq_path) -> None:
    # 统计每个位置的平均质量
    avg_qualities = fastq_quality_stats(fastq_path)
    # 分别绘制总的质量分布图以及首尾200bp的局部分布图
    all_pos_image_path = "png/" + fastq_path.replace("fastq","png")
    head_image_path = "png/" + fastq_path.replace(".fastq","_head.png")
    tail_image_path = "png/" + fastq_path.replace(".fastq","_tail.png")
    sample = fastq_path.replace(".fastq","")
    plot_quality_profile(sample, avg_qualities,all_pos_image_path)
    
    # 统计两端每个位置的平均质量
    terminal_length = 500
    avg_qualities = fastq_quality_terminal(fastq_path, terminal_length)
    plot_quality_profile(sample, avg_qualities[:terminal_length],head_image_path)
    plot_quality_profile(sample, avg_qualities[-terminal_length:],tail_image_path)


def main() -> None:
    fq_li = ["C5920SXXG0-1_241211LR167.fastq", "C2931XKUG0-2_HR912.fastq", "C9940703G0-1_ZX55.fastq", "C7027854G0-1_AU381.fastq"]
    async_in_iterable_structure(run_plot, fq_li, 24)
    
# 🚀 运行主程序
if __name__ == "__main__":
    main()
