import matplotlib.pyplot as plt
import seaborn as sns
from pysam import index, depth, flagstat, FastqFile, VariantFile

def plot_read_length_distribution(fq_path: str, out_png: str):
    read_lengths = []
    with FastqFile(fq_path) as fq:
        for seq in fq:
            read_lengths.append(len(seq.sequence))
            if len(seq.sequence) == 79483:
                print(seq.sequence)
    print(max(read_lengths))
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
    
    
def main() -> None:
    fq_path = "/mnt/ntc_data/wayne/Project/NTC/ONT/pipeline/test_2/clean_reads.fq"
    out_png = "/mnt/ntc_data/wayne/Project/NTC/ONT/pipeline/read_len.png"
    plot_read_length_distribution(fq_path, out_png)
    

if __name__ == "__main__":
    main()