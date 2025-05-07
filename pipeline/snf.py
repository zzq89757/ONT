from pysam import index, depth, flagstat, FastqFile, VariantFile
from main import float_leave_1

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
            # qd = rec.info.get('QD', '.')
            # fs = rec.info.get('FS', '.')
            # mq = rec.info.get('MQ', '.')
            af = rec.info.get('AF', 0.2)
            precise = rec.info.get('PRECISE', 0)
            sv_len = abs(rec.info.get('SVLEN', 0))
            sv_type = rec.info.get('SVTYPE', '')
            gt, gq, dr, dv = str(rec).split("\t")[-1].split(":")
            # GQ过滤
            if int(gq) == 0:
                continue
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
            dp = rec.info.get("SUPPORT",0)
            if dp is None or dp < 1000:
                continue

            # STDEV_LEN 过滤
            std_len = rec.info.get("STDEV_LEN",0)
            if not std_len:
                continue
            sample = rec.samples[0]  # 假设只有一个样本
            gq = sample.get("GQ")
            if gq is None or gq < 3:
                continue
            # type:SNP/SV
            type = "SNP" if sv_len <= 50 else "SV"
            if type == "SNP":
                snps.append({"confidence":confidence, "pos":pos, "ref":ref, "alt":alt, "af": af, "type":type})
            else:
                sv_type = "Duplication" if sv_type == "DUP" else "Deletion"
                snps.append({"confidence":confidence, "pos":pos,  "sv_type":sv_type, "sv_len":sv_len, "af": af, "type":type,})
                
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

a = sniffles_mutation_classify("/mnt/ntc_data/wayne/Project/NTC/ONT/pipeline/test_1/snf.vcf")
print(a)