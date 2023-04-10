import pandas as pd

df = pd.DataFrame(columns=["SampleName", "mode", "fq1", "fq2"])

df["SampleName"] = ["P1", "P2", "R", "S"]
df["mode"] = ["se", "pe", "spet", "spet"]
df["fq1"] = ['ParMaroc_G04_R1_L001.fastq.gz', 'ParC438_F04_R1_L001.fastq.gz', 'P11_F2_M_C438_R_AB_A10_R1_L001.fastq.gz', 'P11_F2_M_C438_S_AB_E10_R1_L001.fastq.gz']
df["fq2"] = ['ParMaroc_G04_R2_L001.fastq.gz', 'ParC438_F04_R2_L001.fastq.gz', 'P11_F2_M_C438_R_AB_A10_R2_L001.fastq.gz ', 'P11_F2_M_C438_S_AB_E10_R2_L001.fastq.gz ']

df.to_csv("samples.csv", sep='\t', index=False)