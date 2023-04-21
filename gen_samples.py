import pandas as pd

df = pd.DataFrame(columns=["SampleName", "mode", "fq1", "fq2"])

df["SampleName"] = ["P1", "P2", "R", "S"]
df["mode"] = ["pe", "pe", "spet", "spet"] #mode must be pe(paired end), se(single end), spet ou spet no umi
df["fq1"] = ['VL-T7_S7_L001_R1_001.fastq.gz',
            'VL-T8_S8_L001_R1_001.fastq.gz',
            'P1_F2_B_H33_R_AB_A01_R1_L001.fastq.gz;P1_F2_B_H33_R_AB_A01_R1_L002.fastq.gz;P1_F2_B_H33_R_AB_A01_R1_L003.fastq.gz;P1_F2_B_H33_R_AB_A01_R1_L004.fastq.gz;P1_F2_B_H33_R_CD_B01_R1_L001.fastq.gz;P1_F2_B_H33_R_CD_B01_R1_L002.fastq.gz;P1_F2_B_H33_R_CD_B01_R1_L003.fastq.gz;P1_F2_B_H33_R_CD_B01_R1_L004.fastq.gz;P1_F2_B_H33_R_EF_C01_R1_L001.fastq.gz;P1_F2_B_H33_R_EF_C01_R1_L002.fastq.gz;P1_F2_B_H33_R_EF_C01_R1_L003.fastq.gz;P1_F2_B_H33_R_EF_C01_R1_L004.fastq.gz;P1_F2_B_H33_R_GH_D01_R1_L001.fastq.gz;P1_F2_B_H33_R_GH_D01_R1_L002.fastq.gz;P1_F2_B_H33_R_GH_D01_R1_L003.fastq.gz;P1_F2_B_H33_R_GH_D01_R1_L004.fastq.gz', 
            'P1_F2_B_H33_S_AB_E01_R1_L001.fastq.gz;P1_F2_B_H33_S_AB_E01_R1_L002.fastq.gz;P1_F2_B_H33_S_AB_E01_R1_L003.fastq.gz;P1_F2_B_H33_S_AB_E01_R1_L004.fastq.gz;P1_F2_B_H33_S_CD_F01_R1_L001.fastq.gz;P1_F2_B_H33_S_CD_F01_R1_L002.fastq.gz;P1_F2_B_H33_S_CD_F01_R1_L003.fastq.gz;P1_F2_B_H33_S_CD_F01_R1_L004.fastq.gz;P1_F2_B_H33_S_EF_G01_R1_L001.fastq.gz;P1_F2_B_H33_S_EF_G01_R1_L002.fastq.gz;P1_F2_B_H33_S_EF_G01_R1_L003.fastq.gz;P1_F2_B_H33_S_EF_G01_R1_L004.fastq.gz;P1_F2_B_H33_S_GH_H01_R1_L001.fastq.gz;P1_F2_B_H33_S_GH_H01_R1_L002.fastq.gz;P1_F2_B_H33_S_GH_H01_R1_L003.fastq.gz;P1_F2_B_H33_S_GH_H01_R1_L004.fastq.gz']

df["fq2"] = ['VL-T7_S7_L001_R2_001.fastq.gz',
            'VL-T8_S8_L001_R2_001.fastq.gz', 
            'P1_F2_B_H33_R_AB_A01_R2_L001.fastq.gz;P1_F2_B_H33_R_AB_A01_R2_L002.fastq.gz;P1_F2_B_H33_R_AB_A01_R2_L003.fastq.gz;P1_F2_B_H33_R_AB_A01_R2_L004.fastq.gz;P1_F2_B_H33_R_CD_B01_R2_L001.fastq.gz;P1_F2_B_H33_R_CD_B01_R2_L002.fastq.gz;P1_F2_B_H33_R_CD_B01_R2_L003.fastq.gz;P1_F2_B_H33_R_CD_B01_R2_L004.fastq.gz;P1_F2_B_H33_R_EF_C01_R2_L001.fastq.gz;P1_F2_B_H33_R_EF_C01_R2_L002.fastq.gz;P1_F2_B_H33_R_EF_C01_R2_L003.fastq.gz;P1_F2_B_H33_R_EF_C01_R2_L004.fastq.gz;P1_F2_B_H33_R_GH_D01_R2_L001.fastq.gz;P1_F2_B_H33_R_GH_D01_R2_L002.fastq.gz;P1_F2_B_H33_R_GH_D01_R2_L003.fastq.gz;P1_F2_B_H33_R_GH_D01_R2_L004.fastq.gz', 
            'P1_F2_B_H33_S_AB_E01_R2_L001.fastq.gz;P1_F2_B_H33_S_AB_E01_R2_L002.fastq.gz;P1_F2_B_H33_S_AB_E01_R2_L003.fastq.gz;P1_F2_B_H33_S_AB_E01_R2_L004.fastq.gz;P1_F2_B_H33_S_CD_F01_R2_L001.fastq.gz;P1_F2_B_H33_S_CD_F01_R2_L002.fastq.gz;P1_F2_B_H33_S_CD_F01_R2_L003.fastq.gz;P1_F2_B_H33_S_CD_F01_R2_L004.fastq.gz;P1_F2_B_H33_S_EF_G01_R2_L001.fastq.gz;P1_F2_B_H33_S_EF_G01_R2_L002.fastq.gz;P1_F2_B_H33_S_EF_G01_R2_L003.fastq.gz;P1_F2_B_H33_S_EF_G01_R2_L004.fastq.gz;P1_F2_B_H33_S_GH_H01_R2_L001.fastq.gz;P1_F2_B_H33_S_GH_H01_R2_L002.fastq.gz;P1_F2_B_H33_S_GH_H01_R2_L003.fastq.gz;P1_F2_B_H33_S_GH_H01_R2_L004.fastq.gz']

df.to_csv("samples.csv", sep='\t', index=False)