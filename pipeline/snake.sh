snakemake -np $(for id in 1_A01 2_A02 3_A03; do echo -n "/mnt/ntc_data/wayne/Project/NTC/ONT/pipeline/C1880JSBG0-$id/aln.bam "; done)
