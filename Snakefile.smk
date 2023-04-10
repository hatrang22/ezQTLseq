__author__ = "INRAE GAFL"
__license__ = "MIT"
__copyright__ = "INRAE, GAFL 2023"

import os
import sys
import pandas as pd
from snakemake.utils import min_version
# snakemake built-in report function requires min version 5.1
min_version("5.1.0")


#############################################################################
#
#     WGS pipeline for SPET 
#     1) preprocessing and quality control
#     2) mapping using bwa and alignement stats
#     3) SNPs calling using GATK. paralisation based on ref splitted by chrs (#samples X #chrs)
#
#############################################################################


##########################
## PREPROCESSING FUNCTIONS
##########################

# Get fastq R1. Could have several fastq file for 1 sample
def get_fq1(wildcards):
    tt=samples.loc[(wildcards.sample), ["fq1"]].dropna()
    fs=re.split('[ ;,]+',tt.item())
    ml=list()
    for f in fs:
      if f.endswith(".gz"):
        f.strip()
        ml.append(config["fq_dir"]+"/"+f)
    return ml
# Get fastq R2. Could have several fastq file for 1 sample
def get_fq2(wildcards):
    tt=samples.loc[(wildcards.sample), ["fq2"]].dropna()
    fs=re.split('[ ;,]+',tt.item())
    ml=list()
    for f in fs:
      if f.endswith(".gz"):
        f.strip()
        ml.append(config["fq_dir"]+"/"+f)
    return ml

# Standardize inputs:
def standardize_bam_remove_duplicate(rmd_config: str, choices=["y", "yes", "n", "no"]):
  if isinstance(rmd_config, str):
    rmd = rmd_config.lower()
    if not rmd in choices:
      raise ValueError(f"Unknown choice ({rmd}). Choices: {choices}")

    return rmd
  else:
    raise TypeError(f"Choice of deduplication must be string among {choices}")

##########################
## GLOBAL VARIABLES
##########################

# Standardize inputs
brmd = standardize_bam_remove_duplicate(config["bam_remove_duplicate"])  

# Read the sample file using pandas lib (sample names+ fastq names) and create index using the sample name
samples = pd.read_csv(config["samplesfile"], sep='\t', dtype=str, comment='#').set_index(["SampleName"], drop=False) 

# Build a chr list based on GENOME_FAI
CHRS=list()
with open(config["REFPATH"]+"/"+config["GENOME_FAI"],'r') as fh:
    #line = fh.readline()
    for line in fh:
      line=line.strip()
      fs=re.split(r'\t',line)
      CHRS.append(fs[0])
      #CHRS.append(line)
fh.close()


###################################################################################
rule final_outs:
    input:
        # expand("test_{mode}.txt", mode=samples['mode'])
        # expand("{outdir}/fastp/{sample[0]}_{sample[1]}_trim.fastq.gz", outdir=config["outdir"], sample=zip(samples['SampleName'], samples['mode']))
        # expand("{outdir}/fastp/{sample}_trim.fastq.gz", outdir=config["outdir"], sample=samples['SampleName'])
        expand("{outdir}/fastqc/{sample[0]}_{sample[1]}.OK.done", outdir=config["outdir"],sample=zip(samples['SampleName'], samples['mode']))
        # "{outdir}/multiqc/multiqc_report_fastqc.html".format(outdir=config["outdir"]),
        #expand("{outdir}/mapped/{sample}_unsorted.bam", outdir=config["outdir"],sample=samples['SampleName']),
        #expand("{outdir}/mapped/{sample}_unsorted.dedup.bam", outdir=config["outdir"],sample=samples['SampleName']),
        #"{outdir}/multiqc/multiqc_report_bam.html".format(outdir=config["outdir"]),
        #expand("{outdir}/mapped/{sample}_clipped.bam", outdir=config["outdir"],sample=samples['SampleName']),
        #"{refp}/{ref}.dict".format(refp=config["REFPATH"],ref=config["GENOME"]),
        #expand("{outdir}/variant/gatk_gvcf/{sample}-{mychr}.g.vcf.gz",outdir=config["outdir"],sample=samples['SampleName'], mychr=CHRS),
        #expand("{outdir}/variant/gvcf_{mychr}_list.map", outdir=config["outdir"], mychr=CHRS),
        #expand("{outdir}/variant/gatk_genomicsdb_{mychr}.ok", outdir=config["outdir"], mychr=CHRS),
        #"{outdir}/variant/gatk_all.vcf.gz".format(outdir=config["outdir"]),
        # "{outdir}/variant/gatk_all.score_snps.pdf".format(outdir=config["outdir"]),
        #"{outdir}/variant/gatk_all.filtered_snps.vcf.gz".format(outdir=config["outdir"]),
        #"{outdir}/variant/gatk_all.keep_biallele.vcf.gz".format(outdir=config["outdir"]),
        #"{outdir}/variant/gatk_all.keep_filter_dp.vcf.gz".format(outdir=config["outdir"]),
        #"{outdir}/variant/gatk_all.keep_snps.vcf.gz".format(outdir=config["outdir"]),
        #"{outdir}/variant/gatk_all.filter_P_snps.vcf.gz".format(outdir=config["outdir"]),
        # "{outdir}/variant/gatk_all.keep_snps.stats.txt".format(outdir=config["outdir"]),
        # "{outdir}/Rqtl/gatk_all.filter_bulk_rs.table".format(outdir=config["outdir"]),
        # "{outdir}/Rqtl/6.Takagi.jpg".format(outdir=config["outdir"]),
        # "{outdir}/Rqtl/qtlsT.csv".format(outdir=config["outdir"]),
        # "{outdir}/Rqtl/qtlsG.csv".format(outdir=config["outdir"]),
        # "{outdir}/Rqtl/filtered.csv".format(outdir=config["outdir"]),"{outdir}/Rqtl/stats.txt".format(outdir=config["outdir"])
     

##################################################################################
########################  PREPROCESSING  #########################################
##################################################################################

# 1-1) reads preprocessing for different modes :
#PREPROC
# Adaptor removal, trim read with fastp
#

rule fastp:
    input:
        R1 = get_fq1,
        R2 = get_fq2
    output:
        R = "{outdir}/fastp/{{sample}}_{{mode}}_trim.fastq.gz".format(outdir=config["outdir"]),
        R1 = "{outdir}/fastp/{{sample}}_{{mode}}_1_trim.fastq.gz".format(outdir=config["outdir"]),
        R2 = "{outdir}/fastp/{{sample}}_{{mode}}_2_trim.fastq.gz".format(outdir=config["outdir"])
    threads: 4
    message: "Running fastp on files \n"
    params:
        mode       = "{mode}",
        spl        = config["outdir"]+"/fastp/{sample}",
        fqdir      = config["fq_dir"],
        outdir     = config["outdir"],
        modules    = config["MODULES"],
        fastp_bin  = config["fastp_bin"],
        bind       = config["BIND"],
        json       = config["outdir"]+"/fastp/{sample}_trim.json",
        html       = config["outdir"]+"/fastp/{sample}_trim.html"
    shell:
       """
       if [[ "{params.mode}" != "pe" ]]; then
        singularity exec {params.bind} {params.fastp_bin} fastp \
        --stdin \
        -i <(zcat {input.R1}) \
        -o {output.R} \
        -h {params.html} \
        -j {params.json} \
        --cut_mean_quality 20 \
        --cut_window_size 4 \
        --low_complexity_filter \
        --complexity_threshold 30 \
        -w {threads}
        touch {output.R1}  # create pseudo-R1-file for inputs of the new rule
        touch {output.R2}  # create pseudo-R2-file for inputs of the new rule
       else
        singularity exec {params.bind} {params.fastp_bin} fastp \
        --stdin \
        -i <(zcat {input.R1}) \
        -I <(zcat {input.R2}) \
        -o {output.R1} \
        -O {output.R2} \
        -h {params.html} \
        -j {params.json} \
        --correction \
        --detect_adapter_for_pe \
        --cut_mean_quality 20 \
        --cut_window_size 4 \
        --low_complexity_filter \
        --complexity_threshold 30 \
        -w {threads}
        touch {output.R}  # create pseudo-R-file for inputs of the new rule
       fi 

        exit 0

       """

# 1-2) Check quality control using FastQC: for each individual
#QC
#fastqc SPET mode
#function return fastq file and adds the fastq path

#fastqc 
#function return fastq1 and fastq2 files and adds the fastq path
rule fastqc:
    input:
        R = "{outdir}/fastp/{{sample}}_{{mode}}_trim.fastq.gz".format(outdir=config["outdir"]),
        R1 = "{outdir}/fastp/{{sample}}_{{mode}}_1_trim.fastq.gz".format(outdir=config["outdir"]),
        R2 = "{outdir}/fastp/{{sample}}_{{mode}}_2_trim.fastq.gz".format(outdir=config["outdir"])
    output:
        # to avoid the output name generated by fastqc (regex on the name) we use a flag file
        "{outdir}/fastqc/{{sample}}_{{mode}}.OK.done".format(outdir=config["outdir"])
    threads: 2
    params:
        mode       = "{mode}",
        fqdir      = config["fq_dir"],
        outdir     = config["outdir"],
        fastqc_bin = config["fastqc_bin"],
        bind       = config["BIND"]
    shell:
        """
        #mkdir -p {params.outdir}/fastqc #snakemake create automaticly the folders
        if [[ "{params.mode}" != "pe" ]]; then
            singularity exec {params.bind} {params.fastqc_bin} fastqc -o {params.outdir}/fastqc -t {threads} {input.R} && touch {output}
            rm -rf {input.R1} {input.R2}  # remove pseudo-files of the previous rule
        else
            singularity exec {params.bind} {params.fastqc_bin} fastqc -o {params.outdir}/fastqc -t {threads} {input.R1} {input.R2} && touch {output}
            rm -rf {input.R}  # remove pseudo-file of the previous rule
        fi
        """

# 1-3) check quality control using FastQC: a global multiQC for the analysis
# multiQC on fastqc outputs
rule multiqc_fastqc:
    input:
        expand("{outdir}/fastqc/{sample}.OK.done", outdir=config["outdir"], sample=samples['SampleName'])
    output:
        #report("{outdir}/multiqc/multiqc_report_fastqc.html".format(outdir=config["outdir"]), caption="report/multiqc.rst", category="Quality control")
        "{outdir}/multiqc/multiqc_report_fastqc.html".format(outdir=config["outdir"])
    threads:
        1
    params:
        outdir      = config["outdir"],
        multiqc_bin = config["multiqc_bin"],
        bind        = config["BIND"]
    shell:
        """
        #mkdir -p {params.outdir}/multiqc #snakemake create automaticly the folders
        singularity exec {params.bind} {params.multiqc_bin} multiqc --filename {output} {params.outdir}/fastqc
        """
# 1-4) merge forward and reverse reads for PE :
# output:
# merged reads: _merged.fastq.gz
# unmerged: _1_umerged.fastq.gz + _2_umerged.fastq.gz
rule mergePE:
    input:
        R1 = "{outdir}/fastp/{{sample}}_1_trim.fastq.gz".format(outdir=config["outdir"]),
        R2 = "{outdir}/fastp/{{sample}}_2_trim.fastq.gz".format(outdir=config["outdir"]),
    output:
        R = "{outdir}/fastp/{{sample}}_merged.fastq.gz".format(outdir=config["outdir"]),
        R1 = "{outdir}/fastp/{{sample}}_1_umerged.fastq.gz".format(outdir=config["outdir"]),
        R2 = "{outdir}/fastp/{{sample}}_2_umerged.fastq.gz".format(outdir=config["outdir"])
    params:
        modules       = config["MODULES"],
        fastp_bin     = config["fastp_bin"],
        bind          = config["BIND"],
        json          = config["outdir"]+"/fastp/{sample}_trim.json",
        html          = config["outdir"]+"/fastp/{sample}_trim.html",
    threads:
        10
    shell:
        """
        singularity exec {params.bind} {params.fastp_bin} fastp \
        -i {input.R1} \
        -I {input.R2} \
        --merge \
        --correction \
        --merged_out {output.R} \
        --out1 {output.R1} \
        --out2 {output.R2} \
        --json={params.json} \
        --html={params.html} \
        --length_required 50 \
        --thread {threads}
        rm -f {params.html} {params.json}
        """

###############################################################################
########################  DNA mapping   #######################################
###############################################################################
# 2-1) reference indexation
rule bwa_index:
    input:
        genome = config["REFPATH"] + "/" + config["GENOME"]
    output:
       config["REFPATH"] + "/" + config["GENOME"] + ".amb",
       config["REFPATH"] + "/" + config["GENOME"] + ".ann",
       config["REFPATH"] + "/" + config["GENOME"] + ".bwt",
       config["REFPATH"] + "/" + config["GENOME"] + ".pac",
       config["REFPATH"] + "/" + config["GENOME"] + ".sa",
       config["REFPATH"] + "/" + config["GENOME"] + ".fai"
    params:
        bind         = config["BIND"],
        bwa_bin      = config["bwa_bin"],
        samtools_bin = config["samtools_bin"]
    shell:
        """
        singularity exec {params.bind} {params.bwa_bin} bwa index -a bwtsw -b 500000000 {input.genome}
        singularity exe {params.bind} {params.samtools_bin} samtools faidx {input.genome} # create .fai file 
        """
# 2-2) mapping

# if smode == "pe":
#     ruleorder: bwa_pe > bwa_spet_se
# else:
#     ruleorder: bwa_spet_se > bwa_pe


# ####################SPET or SE mode##################################
#   Requirement:
#       - config["REFPATH"]/config["GENOME"] indexed
#       - trimmed reads in config["outdir"]/fastp/
#   Input:
#       - config["REFPATH"] and config["GENOME"]
#       - {sample}_trim.fastq.gz
#   Output:
#       - {sample}_sorted.bam
#       - {sample}_sorted.bam.bai
# ############################################################
rule bwa_spet_se :
    input:
        R = "{outdir}/fastp/{{sample}}_trim.fastq.gz".format(outdir=config["outdir"]),
        #fake input used to force index building before alignement if not present
        idx = config["REFPATH"] + "/" + config["GENOME"] + ".bwt"
    output:
        bam   = "{outdir}/mapped/raw/{{sample}}_sorted.raw.bam".format(outdir=config["outdir"]),
        #bai   = "{outdir}/mapped/raw/{{sample}}_sorted.raw.bam.bai".format(outdir=config["outdir"]),
    params:
        outtmp       = "{outdir}/mapped/{{sample}}".format(outdir=config["outdir"]),
        idxbase      = config["REFPATH"] + "/" + config["GENOME"],
        bind         = config["BIND"],
        bwa_bin      = config["bwa_bin"],
        samtools_bin = config["samtools_bin"],
        rg           = "@RG\\tID:{sample}\\tSM:{sample}"
    threads: 16
    shell:
        """
        singularity exec {params.bind} {params.bwa_bin} bwa mem -t {threads} -K 100000000 -R '{params.rg}' {params.idxbase} {input.R}|singularity exec {params.bind} {params.samtools_bin} samtools view -h -F 2048 -o {output.bam} 
        |singularity exec {params.bind} {params.samtools_bin} samtools sort -@2 -m 6G -o {params.outtmp}.um.tmp.bam
        #merge bam
        singularity exec {params.bind} {params.samtools_bin} samtools merge {output.bam} 
        singularity exec {params.bind} {params.samtools_bin} samtools flagstat -@ {threads} {output.bam}
        
        """
#2 stage : first with unmerged PE and second with merged reads
# ############################################################
#   Requirement:
#       - config["REFPATH"]/config["GENOME"] indexed
#       - trimmed reads in config["outdir"]/fastp/
#   Input:
#       - config["REFPATH"] and config["GENOME"]
#       - {sample}_merged.fastq.gz
#       -{sample}_1_umerged.fastq.gz + {sample}_2_umerged.fastq.gz
#   Output:
#       - {sample}_sorted.bam
#       - {sample}_sorted.bam.bai
# ############################################################
rule bwa_pe :
    input:
        R = "{outdir}/fastp/{{sample}}_merged.fastq.gz".format(outdir=config["outdir"]),
        R1 = "{outdir}/fastp/{{sample}}_1_umerged.fastq.gz".format(outdir=config["outdir"]),
        R2 = "{outdir}/fastp/{{sample}}_2_umerged.fastq.gz".format(outdir=config["outdir"]),
        #fake input used to force index building before alignement if not present
        idx = config["REFPATH"] + "/" + config["GENOME"] + ".bwt"
    output:
        bam   = "{outdir}/mapped/raw/{{sample}}_sorted.raw.bam".format(outdir=config["outdir"]),
        #bai   = "{outdir}/mapped/raw/{{sample}}_sorted.raw.bam.bai".format(outdir=config["outdir"]),
    params:
        outtmp       = "{outdir}/mapped/raw/{{sample}}".format(outdir=config["outdir"]),
        idxbase      = config["REFPATH"] + "/" + config["GENOME"],
        bind         = config["BIND"],
        bwa_bin      = config["bwa_bin"],
        samtools_bin = config["samtools_bin"],
        rg           = "@RG\\tID:{sample}\\tSM:{sample}"
    threads: 8
    shell:
        """
        singularity exec {params.bind} {params.bwa_bin} bwa mem \
        -t {threads} \
        -K 100000000 \
        -R '{params.rg}' \
        {params.idxbase} \
        {input.R1} {input.R2} \
        | singularity exec {params.bind} {params.samtools_bin} samtools view -h -F 2048 \
        | singularity exec {params.bind} {params.samtools_bin} samtools sort -@2 -m 6G -o {params.outtmp}.um.tmp.bam -

        singularity exec {params.bind} {params.bwa_bin} bwa mem \
        -t {threads} \
        -K 100000000 \
        -R '{params.rg}' \
        {params.idxbase} \
        {input.R} \
        | singularity exec {params.bind} {params.samtools_bin} samtools view -h -F 2048 \
        |singularity exec {params.bind} {params.samtools_bin} samtools sort -@2 -m 6G -o {params.outtmp}.m.tmp.bam -
        # merge bams
        singularity exec {params.bind} {params.samtools_bin} samtools merge {output.bam} {params.outtmp}.m.tmp.bam {params.outtmp}.um.tmp.bam && rm -f {params.outtmp}.m.tmp.bam {params.outtmp}.um.tmp.bam
        singularity exec {params.bind} {params.samtools_bin} samtools flagstat -@ {threads} {output.bam}
        """

# 2-3) mark end remove duplicate PCR 

# if smode == "spet":
#     ruleorder: remove_duplicates_spet > keep_duplicates > remove_duplicates_pe
# elif smode == "pe" and brmd.startswith("y"):
#     ruleorder: remove_duplicates_pe > keep_duplicates > remove_duplicates_spet
# # elif smode == "pe" and brmd.startswith("n"):
# #     ruleorder: keep_duplicates > remove_duplicates_pe > remove_duplicates_spet
# else: 
#     ruleorder: keep_duplicates > remove_duplicates_pe > remove_duplicates_spet



##https://github.com/tecangenomics/nudup
rule remove_duplicates_spet:
    input:
        bam   = "{outdir}/mapped/raw/{{sample}}_sorted.raw.bam".format(outdir=config["outdir"]),
        R     = get_fq2
    output:
        dedup = "{outdir}/mapped/raw/{{sample}}.dedup.bam".format(outdir=config["outdir"]),
        bam   = "{outdir}/mapped/{{sample}}_sorted.bam".format(outdir=config["outdir"]),
        bai   = "{outdir}/mapped/{{sample}}_sorted.bam.bai".format(outdir=config["outdir"])
    params:
        samtools_bin = config["samtools_bin"],
        bind     = config["BIND"],
        outdir   = config["outdir"],
        #nudup_bin = config["nudup_bin"]
    threads: 4
    message:"Remove duplicate.\n"
    shell:
        """
        singularity exec nudup_v2.3.3.sif nudup.py -f {input.R} --out {output.dedup} {input.bam} \
        |singularity exec {params.bind} {params.samtools_bin} samtools sort -@2 -m 6G -o {output.bam}

        #samtools index
        singularity exec {params.bind} {params.samtools_bin}  samtools index -@4 {output.bam}   
        """
rule remove_duplicates_pe:
    input:
        bam   = "{outdir}/mapped/raw/{{sample}}_sorted.raw.bam".format(outdir=config["outdir"]),
        #bai   = "{outdir}/mapped/raw/{{sample}}_sorted.raw.bam.bai".format(outdir=config["outdir"]),
    output:
        bam   = "{outdir}/mapped/{{sample}}_sorted.bam".format(outdir=config["outdir"]),
        bai   = "{outdir}/mapped/{{sample}}_sorted.bam.bai".format(outdir=config["outdir"]),
    params:
        tempdir=config["outdir"] + "/mapped/tmpdir",
        bind         = config["BIND"],
        samtools_bin = config["samtools_bin"],
        sambamba_bin= config["sambamba_bin"], # sambamba for markdup
    threads: 4
    message: "Mark/remove PCR duplicates.\n"
    shell:
        """
        ulimit -n 4048
        mkdir -p {params.tempdir}
        singularity exec {params.bind} {params.samtools_bin} samtools rmdup -s {input.bam} {output.bam}
        #singularity exec {params.bind} {params.sambamba_bin} sambamba markdup --tmpdir={params.tempdir} -t {threads} {input.bam} {output.bam} && \
        rm -f {input.bam} && \
        singularity exec {params.bind} {params.samtools_bin} samtools index -@ {threads} {output.bam}
        #--remove-duplicates #remove duplicates instead of just marking them
        """

rule keep_duplicates:
    input:
        bam   = "{outdir}/mapped/raw/{{sample}}_sorted.raw.bam".format(outdir=config["outdir"]),
        #bai   = "{outdir}/mapped/raw/{{sample}}_sorted.raw.bam.bai".format(outdir=config["outdir"]),
    output:
        bam   = "{outdir}/mapped/{{sample}}_sorted.bam".format(outdir=config["outdir"]),
        bai   = "{outdir}/mapped/{{sample}}_sorted.bam.bai".format(outdir=config["outdir"]),
    params:
        bind         = config["BIND"],
        samtools_bin = config["samtools_bin"],
    threads: 4
    message: "Keep PCR duplicates.\n"
    shell:
        """
        mv {input.bam} {output.bam}
        singularity exec {params.bind} {params.samtools_bin} samtools index -@ {threads} {output.bam}
        """

# 2_3) bam stats 
rule bam_stats:
    input:
        bam   = "{outdir}/mapped/{{sample}}_sorted.bam".format(outdir=config["outdir"])
    output:
        stats = "{outdir}/mapped/{{sample}}.stats.txt".format(outdir=config["outdir"])
    params:
        outdir       = config["outdir"],
        bind         = config["BIND"],
        samtools_bin = config["samtools_bin"]
    threads: 1
    shell:
        """
        singularity exec {params.bind} {params.samtools_bin} samtools stats -@ {threads} {input.bam} > {output.stats}
        """
# 2-4) multiqc for bams 
rule multiqc_bam:
    input:
        expand("{outdir}/mapped/{sample}.stats.txt", outdir=config["outdir"], sample=samples['SampleName'])
    output:
        "{outdir}/multiqc/multiqc_report_bam.html".format(outdir=config["outdir"])
    threads:
        1
    params:
        outdir      = config["outdir"]+"/mapped/*.stats.txt",
        multiqc_bin = config["multiqc_bin"],
        bind        = config["BIND"]
    shell:
        """
        singularity exec {params.bind} {params.multiqc_bin} multiqc --filename {output} {input}
        """ 
# 2-6) clip probe
## sort clip the first 40pb to remove probederived sequence
##by using BamUtil
rule rm_probe:
    input:
       bam  = "{outdir}/mapped/{{sample}}_sorted.bam".format(outdir=config["outdir"])
    output:
       clip = "{outdir}/mapped/{{sample}}_clipped.bam".format(outdir=config["outdir"])
    params:
       bind        = config["BIND"],
       outdir      = config["outdir"],
       bamutil_bin = config["bamUtil_bin"],
       probe_length = config["probe_length"]
    threads: 4
    shell:
       """
       singularity exec {params.bind} {params.bamutil_bin} bam TrimBam {input.bam} {output.clip} --clip -L {params.probe_length}
       """
################################################################################
###########################  SNP calling using  GATK ###########################
################################################################################
# 3) SNP caller using GATK4
#gatk:4 rules
# 3-1) create a reference dictionnary the ref must not be a sl and a ref.fa => ref.dict (not ref.fa.dict)
rule gatk4_ref_dict:
    input:
        ref  = "{refp}/{ref}".format(refp=config["REFPATH"],ref=config["GENOME"]),
        #fai  = "{refp}/{ref}".format(refp=config["REFPATH"],ref=config["GENOME_FAI"]),
    output:
        dic  = "{refp}/{ref}.dict".format(refp=config["REFPATH"],ref=config["GENOME"]) #config["REFPATH"] + "/" + config["GENOME"] +".dict"
    params:
        bind = config["BIND"],
        samtools_bin = config["samtools_bin"],
        gatk4_bin = config["gatk4_bin"]
    shell:
        """
        singularity exec {params.bind} {params.gatk4_bin} gatk CreateSequenceDictionary -R {input.ref}
        #singularity exec {params.bind} {params.samtools_bin} samtools faidx {input.ref} #create fasta index file
        #singularity exec {params.bind} {params.samtools_bin} samtools dict -o {output.dic} {input.ref}
        """
# ------------------  parrallel by chr -----------------
# 3-2) HaplotypeCaller for each sample and each chr

# if smode == "spet":
#     ruleorder: gatk4_hc_spet > gatk4_hc_pe_se
# else:
#     ruleorder: gatk4_hc_pe_se > gatk4_hc_spet


rule gatk4_hc_spet:
    input:
        #"{outdir}/variant/chrslist.OK".format(outdir=config["outdir"]),
        bam     = "{outdir}/mapped/{{sample}}_clipped.bam".format(outdir=config["outdir"]),
        #refdict = "{refp}/{ref}.dict".format(refp=config["REFPATH"],ref=config["GENOME"]),
    output:
        gvcf="{outdir}/variant/gatk_gvcf/{{sample}}-{{mychr}}.g.vcf.gz".format(outdir=config["outdir"])
    params:
        ch="{mychr}",
        ref = config["REFPATH"] + "/" + config["GENOME"],
        bind = config["BIND"],
        samtools_bin = config["samtools_bin"],
        gatk4_bin = config["gatk4_bin"]
    threads: 10
    shell:
        """
        singularity exec {params.bind} {params.gatk4_bin} \
        gatk --java-options '-Xmx24G' HaplotypeCaller \
        --reference {params.ref} \
        --input {input.bam} \
        -ERC GVCF -L {params.ch} \
        --output {output.gvcf}
        exit 0
        """
rule gatk4_hc_pe_se:
    input:
        #"{outdir}/variant/chrslist.OK".format(outdir=config["outdir"]),
        bam     = "{outdir}/mapped/{{sample}}_sorted.bam".format(outdir=config["outdir"]),
        refdict = "{refp}/{ref}.dict".format(refp=config["REFPATH"],ref=config["GENOME"]),
    output:
        gvcf="{outdir}/variant/gatk_gvcf/{{sample}}-{{mychr}}.g.vcf.gz".format(outdir=config["outdir"])
    params:
        ch="{mychr}",
        ref = config["REFPATH"] + "/" + config["GENOME"],
        bind = config["BIND"],
        samtools_bin = config["samtools_bin"],
        gatk4_bin = config["gatk4_bin"]
    threads: 1
    shell:
        """
        singularity exec {params.bind} {params.gatk4_bin} \
        gatk --java-options '-Xmx24G -XX:ParallelGCThreads={threads}' HaplotypeCaller \
        --reference {params.ref} \
        --input {input.bam} \
        -ERC GVCF -L {params.ch} \
        --output {output.gvcf}
        exit 0
        """

# 3-3) Generate a map (text file) of gvcf files for each chr
rule gatk4_gvcf_map_file:
    input:
        gvcfs = expand("{outdir}/variant/gatk_gvcf/{sample}-{{mychr}}.g.vcf.gz", outdir=config["outdir"], sample=samples['SampleName']),
        
    output:
        gvcfmap      = "{outdir}/variant/gvcf_{{mychr}}_list.map".format(outdir=config["outdir"]),
    params:
        ch="{mychr}",
        #csv         = expand("{sample}\t{outdir}/variant/gatk_gvcf/{sample}-{{mychr}}.g.vcf.gz", outdir=config["outdir"], sample=samples['SampleName']),
        csv         = expand("{sample}\t{outdir}/variant/gatk_gvcf/{sample}", outdir=config["outdir"], sample=samples['SampleName']),
        bamsp       = "{outdir}/mapped".format(outdir=config["outdir"]),
        outlist     = config["outdir"]
    threads: 1
    run:
        ch=params.ch
        x=params.csv
        f = open(output.gvcfmap, "w")
        for i in x:
            f.write("{l}-{mchr}.g.vcf.gz\n".format(l=i,mchr=ch))
        f.close()



# 3-4) build a datastore of gvcf for each chr

# see: https://gatk.broadinstitute.org/hc/en-us/articles/360040096732-GenomicsDBImport
rule gatk4_gdb:
    input:
        gvcf = expand("{outdir}/variant/gatk_gvcf/{sample}-{{mychr}}.g.vcf.gz",outdir=config["outdir"],sample=samples['SampleName']),
        gvcfmap = "{outdir}/variant/gvcf_{{mychr}}_list.map".format(outdir=config["outdir"])
    output:
        "{outdir}/variant/gatk_genomicsdb_{{mychr}}.ok".format(outdir=config["outdir"])
    params:
        ch="{mychr}",
        tempdir=config["outdir"] + "/variant/tmpgatkdir",
        gdb  = config["outdir"] + "/variant/" + "GenomicsDB_{mychr}",
        ref = config["REFPATH"] + "/" + config["GENOME"],
        bind = config["BIND"],
        samtools_bin = config["samtools_bin"],
        gatk4_bin = config["gatk4_bin"]
    threads: 12
    shell:
        """
        mkdir -p {params.tempdir}
        singularity exec {params.bind} {params.gatk4_bin} \
        gatk --java-options '-Xmx8G' GenomicsDBImport \
        --genomicsdb-workspace-path {params.gdb} \
        --batch-size 200 \
        -L {params.ch} \
        --sample-name-map {input.gvcfmap} \
        --reader-threads {threads} \
        --tmp-dir={params.tempdir} && touch {output}
        """
# 3-5) genotype (GenotypeGVCFs) for all samples for each chr
# https://gatk.broadinstitute.org/hc/en-us/articles/360047218551-GenotypeGVCFs
rule gatk4_gc:
    input:
        gdb="{outdir}/variant/gatk_genomicsdb_{{mychr}}.ok".format(outdir=config["outdir"])
    output:
        vcf="{outdir}/variant/gatk_{{mychr}}_genotyped.vcf.gz".format(outdir=config["outdir"])
        # WARNING if {outdir}/variant/gatk_{{mychr}}.vcf.gz is NOT working we need {{mychr}}_something
    params:
        tempdir=config["outdir"] + "/tmpgatkdir",
        gdb  = config["outdir"] + "/variant/" + "GenomicsDB_{mychr}",
        ref = config["REFPATH"] + "/" + config["GENOME"],
        bind = config["BIND"],
        samtools_bin = config["samtools_bin"],
        gatk4_bin = config["gatk4_bin"]
    threads: 1
    message: "GATK4 genotype_variants vcf\n"
    shell:
        """
        mkdir -p {params.tempdir}
        singularity exec {params.bind} {params.gatk4_bin} \
        gatk --java-options '-Xmx8G' GenotypeGVCFs \
        --variant gendb://{params.gdb} \
        --reference {params.ref} \
        --output {output.vcf} \
        --tmp-dir={params.tempdir}
        """
# 3-6) merge all the vcf produced by chr
rule combinevcf:
    input:
        #"{outdir}/variant/chrslist.OK".format(outdir=config["outdir"]),
        vcf=expand("{outdir}/variant/gatk_{mychr}_genotyped.vcf.gz",outdir=config["outdir"],mychr=CHRS),
        refdict = config["REFPATH"] + "/" + config["GENOME"] +".dict",
    output:
        gvcf = "{outdir}/variant/gatk_all.vcf.gz".format(outdir=config["outdir"])
    params:
        mlist = expand(" -I {outdir}/variant/gatk_{mychr}_genotyped.vcf.gz",outdir=config["outdir"],mychr=CHRS),
        ref = config["REFPATH"] + "/" + config["GENOME"],
        bind = config["BIND"],
        samtools_bin = config["samtools_bin"],
        gatk4_bin = config["gatk4_bin"]
    shell:
        """
        singularity exec {params.bind} {params.gatk4_bin} \
        gatk --java-options '-Xmx16G' GatherVcfs  {params.mlist} -O {output.gvcf}
        singularity exec {params.bind} {params.gatk4_bin} \
        gatk --java-options '-Xmx16G' IndexFeatureFile -I {output.gvcf}
        #-I vcf_list_to_merge.txt
        """
# -------------------- SNPs selection and filtering --------------------------------
# 3-7) select SNPs only
rule gatk4_select_snps_variants:
    input:
        vcf = "{outdir}/variant/gatk_all.vcf.gz".format(outdir=config["outdir"])
    output:
        vcf="{outdir}/variant/gatk_all.raw_snps.vcf.gz".format(outdir=config["outdir"])
    params:
        ref = config["REFPATH"] + "/" + config["GENOME"],
        bind = config["BIND"],
        samtools_bin = config["samtools_bin"],
        gatk4_bin = config["gatk4_bin"]
    threads: 10
    message: "GATK4 select only SNPs variants vcf\n"
    shell:
        """
        singularity exec {params.bind} {params.gatk4_bin} \
        gatk --java-options '-Xmx8G -XX:ParallelGCThreads={threads}' SelectVariants \
        --select-type-to-include SNP \
        --variant {input.vcf} \
        --reference {params.ref} \
        --output {output.vcf}
        """
# 3-8) hard filtering as outlined in GATK docs
# (https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set)
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
# QualByDepth (QD)
# FisherStrand (FS)
# StrandOddsRatio (SOR)
# RMSMappingQuality (MQ)
# MappingQualityRankSumTest (MQRankSum)
# ReadPosRankSumTest (ReadPosRankSum)
#
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112?id=2806
# For reference, here are some basic filtering thresholds to improve upon.
# QD < 2.0
# QUAL < 30.0
# SOR > 3.0
# FS > 60.0
# MQ < 40.0
# MQRankSum < -12.5
# ReadPosRankSum < -8.0
rule gatk4_filterSnps:
    input:
        vcf="{outdir}/variant/gatk_all.raw_snps.vcf.gz".format(outdir=config["outdir"])
    output:
        vcf="{outdir}/variant/gatk_all.filtered_snps.vcf.gz".format(outdir=config["outdir"])
    params:
        vcftmp ="{outdir}/variant/vcf.snp.tmp.vcf.gz".format(outdir=config["outdir"]),
        ref = config["REFPATH"] + "/" + config["GENOME"],
        bind = config["BIND"],
        samtools_bin = config["samtools_bin"],
        gatk4_bin = config["gatk4_bin"],
        QUAL_filter = config["snp_QUAL_filter"],
        QD_filter = config["snp_QD_filter"],
        FS_filter = config["snp_FS_filter"],
        MQ_filter = config["snp_MQ_filter"],
        SOR_filter = config["snp_SOR_filter"],
        MQRankSum_filter = config["snp_MQRankSum_filter"],
        ReadPosRankSum_filter = config["snp_ReadPosRankSum_filter"]
    threads: 2
    shell:
        """
        singularity exec {params.bind} {params.gatk4_bin} \
        gatk --java-options '-Xmx8G' VariantFiltration \
        --variant {input.vcf} \
        --reference {params.ref} \
        --output {params.vcftmp} \
        -filter-name "QUAL_filter" -filter "QUAL {params.QUAL_filter}" \
        -filter-name "QD_filter" -filter "QD {params.QD_filter}" \
        -filter-name "FS_filter" -filter "FS {params.FS_filter}" \
        -filter-name "MQ_filter" -filter "MQ {params.MQ_filter}" \
        -filter-name "SOR_filter" -filter "SOR {params.SOR_filter}" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum {params.MQRankSum_filter}" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum {params.ReadPosRankSum_filter}" \
        --genotype-filter-expression "DP<3" \
        --genotype-filter-name "LowDP3" \

        singularity exec {params.bind} {params.gatk4_bin} \
        gatk --java-options '-Xmx8G' SelectVariants \
        --variant {params.vcftmp} \
        --set-filtered-gt-to-nocall \
        --exclude-filtered \
        --reference {params.ref} \
        --output {output.vcf} \


        rm -f {params.vcftmp} {params.vcftmp}.tbi
        """

# 3-9) rules for creating quality score graphs for raw SNPs.
# First we do a table with the quality score using the VariantsToTable option
rule gatk4_snps_raw_score_table:
    input:
        vcf="{outdir}/variant/gatk_all.raw_snps.vcf.gz".format(outdir=config["outdir"])
    output:
        csv="{outdir}/variant/gatk_all.score_raw_snps.csv".format(outdir=config["outdir"])
    params:
        ref = config["REFPATH"] + "/" + config["GENOME"],
        bind = config["BIND"],
        samtools_bin = config["samtools_bin"],
        gatk4_bin = config["gatk4_bin"]
    threads: 2
    shell:
        """
        singularity exec {params.bind} {params.gatk4_bin} \
        gatk --java-options '-Xmx8G' VariantsToTable \
        -V {input.vcf} \
        -R {params.ref} \
        -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
        -O {output.csv}
        """
# 3-10) rules for creating quality score graphs for filtered SNPs.
# First we do a table with the quality score using the VariantsToTable option
rule gatk4_snps_filtered_score_table:
    input:
        vcf="{outdir}/variant/gatk_all.filtered_snps.vcf.gz".format(outdir=config["outdir"])
    output:
        csv="{outdir}/variant/gatk_all.score_filtered_snps.csv".format(outdir=config["outdir"])
    params:
        ref = config["REFPATH"] + "/" + config["GENOME"],
        bind = config["BIND"],
        samtools_bin = config["samtools_bin"],
        gatk4_bin = config["gatk4_bin"]
    threads: 2
    shell:
        """
        singularity exec {params.bind} {params.gatk4_bin} \
        gatk --java-options '-Xmx8G -XX:ParallelGCThreads={threads}' VariantsToTable \
        -V {input.vcf} \
        -R {params.ref} \
        -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
        -O {output.csv}
        """

## 3-11) plots run Rscript to create quality score graphs to determine filters for snps vcf
## Rscript "gatk_scores_qual_raw_vs_filtered.R" must be download
rule gatk4_snps_quality_graphe_pdf:
    input:
        rawsnps="{outdir}/variant/gatk_all.score_raw_snps.csv".format(outdir=config["outdir"]),
        filteredsnps = "{outdir}/variant/gatk_all.score_filtered_snps.csv".format(outdir=config["outdir"])
    output:
        pdf="{outdir}/variant/gatk_all.score_snps.pdf".format(outdir=config["outdir"])
    params:
        bind = config["BIND"],
        R_bin = config["R_bin"]
    shell:
        """
        singularity exec {params.bind} {params.R_bin} \
        Rscript --vanilla gatk_scores_qual_raw_vs_filtered.R {input.rawsnps} {input.filteredsnps} "SNPs" {output.pdf}
        """

################################END GATK #################################################


# 1) --------- we select variant from the previous filter and keep only biallelic samples ---------------

## 4-1 ---------------keep biallelic ------------------
rule keep_biallelic:
    input:
        snp = "{outdir}/variant/gatk_all.filtered_snps.vcf.gz".format(outdir=config["outdir"])
    output:
        biallel = "{outdir}/variant/gatk_all.keep_biallele.vcf.gz".format(outdir=config["outdir"])
    params:
        ref = config["REFPATH"] + "/" + config["GENOME"],
        bind = config["BIND"],
        gatk4_bin = config["gatk4_bin"]
    threads: 1
    shell:
        """
        singularity exec {params.bind} {params.gatk4_bin} \
        gatk --java-options '-Xmx8G -XX:ParallelGCThreads={threads}' SelectVariants \
        --variant {input.snp} \
        --exclude-filtered \
        --restrict-alleles-to BIALLELIC \
        --output {output.biallel}
        """
## 4-2 ---------- filter in depth -----------------
rule filter_by_dp:
    input:
        snp = "{outdir}/variant/gatk_all.keep_biallele.vcf.gz".format(outdir=config["outdir"])
    output:
        allele = "{outdir}/variant/gatk_all.keep_filter_dp.vcf.gz".format(outdir=config["outdir"])
    params:
        ref = config["REFPATH"] + "/" + config["GENOME"],
        bind = config["BIND"],
        gatk4_bin = config["gatk4_bin"],
        min_dp = config["min_dp"],
        max_dp = config["max_dp"],
        dp_p = config["dp_p"]
    threads: 1
    shell:
        """
        singularity exec {params.bind} {params.gatk4_bin} \
        gatk --java-options '-Xmx8G -XX:ParallelGCThreads={threads}' SelectVariants \
        --variant {input.snp} \
        -select "vc.getGenotype('P1').getDP()>{params.dp_p} && vc.getGenotype('P2').getDP() > {params.dp_p} && (vc.getGenotype('R').getDP() > {params.min_dp} && vc.getGenotype('R').getDP() < {params.max_dp}) && (vc.getGenotype('S').getDP() > {params.min_dp} && vc.getGenotype('S').getDP() < {params.max_dp})" \
        --output {output.allele}
        """

## 4-3) ------------- keep only SNP without missing genotypes --------------
rule gatk4_no_missing:
    input:
        vcf = "{outdir}/variant/gatk_all.keep_filter_dp.vcf.gz".format(outdir=config["outdir"])
    output:
        keep_genotype="{outdir}/variant/gatk_all.keep_snps.vcf.gz".format(outdir=config["outdir"])
    params:
        ref = config["REFPATH"] + "/" + config["GENOME"],
        bind = config["BIND"],
        gatk4_bin = config["gatk4_bin"]
    threads: 1
    shell:
        """
        singularity exec {params.bind} {params.gatk4_bin} \
        gatk --java-options '-Xmx8G -XX:ParallelGCThreads={threads}' SelectVariants \
        --variant {input.vcf} \
        --exclude-filtered \
        --max-nocall-fraction 0.00 \
        --output {output.keep_genotype}
        """

## 4-4) ---------keep polymorphism with the parents -----------
rule keep_variant_Parents:
    input:
        snp="{outdir}/variant/gatk_all.keep_snps.vcf.gz".format(outdir=config["outdir"])
    output:
        polymorph="{outdir}/variant/gatk_all.filter_P_snps.vcf.gz".format(outdir=config["outdir"])
    params:
        ref = config["REFPATH"] + "/" + config["GENOME"],
        bind = config["BIND"],
        gatk_bin = config["gatk_bin"]
    shell:
        """
        singularity exec {params.bind} {params.gatk_bin} \
        java -jar /opt/jvarkit_github/jvarkit/dist/vcffilterjdk.jar \
        -e 'return !variant.getGenotype("P1").sameGenotype(variant.getGenotype("P2"));' -o {output.polymorph} {input.snp}#keep only the polymorph different with the parents 
        """
## 4-5) ---------keep and extract R and S bulks ----------------
rule keep_extract_bulks:
     input:
        parents="{outdir}/variant/gatk_all.filter_P_snps.vcf.gz".format(outdir=config["outdir"])
     output:
        bulks="{outdir}/variant/gatk_all.filter_bulk.vcf.gz".format(outdir=config["outdir"])
     params:
        bcftools_bin = config["bcftools_bin"],
        bind = config["BIND"]
     shell:
        """
        singularity exec {params.bind} {params.bcftools_bin} bcftools view -s R,S {input.parents} -Oz -o {output.bulks}
        singularity exec {params.bind} {params.bcftools_bin} bcftools index {output.bulks} -t
        """

## 4-6) -------------- statistics -------------------------------------
rule gatk_stats:
    input:
        snp = "{outdir}/variant/gatk_all.filter_bulk.vcf.gz".format(outdir=config["outdir"])
    output:
        stat = "{outdir}/variant/gatk_all.keep_snps.stats.txt".format(outdir=config["outdir"])
    params:
        outdir       = config["outdir"],
        idxbase      = config["REFPATH"] + "/" + config["GENOME"],
        bind         = config["BIND"],
        bwa_bin      = config["bwa_bin"],
        bcftools_bin = config["bcftools_bin"]
    threads: 2
    shell:
        """
        singularity exec {params.bind} {params.bcftools_bin} bcftools stats -s P1,P2,R,S {input.snp} > {output.stat}
        """

## 4-7 ----------create a table file for Rscript ---------------
rule vcf2table:
    input:
        bulks= "{outdir}/variant/gatk_all.filter_bulk.vcf.gz".format(outdir=config["outdir"])
    output:
        rs="{outdir}/Rqtl/gatk_all.filter_bulk_rs.table".format(outdir=config["outdir"])
    params:
        ref = config["REFPATH"] + "/" + config["GENOME"],
        bind = config["BIND"],
        samtools_bin = config["samtools_bin"],
        gatk4_bin = config["gatk4_bin"]
    threads: 1
    shell:
        """
        singularity exec {params.bind} {params.gatk4_bin} \
        gatk VariantsToTable \
        -V {input.bulks} \
        -F CHROM -F POS -F REF -F ALT -GF AD -GF DP -GF GQ -GF PL \
        -O {output.rs}
        """
## 4-8 ---------run the script Rqtl -----------------
rule do_qtlseqR:
    input:
        snpfile="{outdir}/Rqtl/gatk_all.filter_bulk_rs.table".format(outdir=config["outdir"])
    output:
        jpg6="{outdir}/Rqtl/6.Takagi.jpg".format(outdir=config["outdir"]),
        qtlsT="{outdir}/Rqtl/qtlsT.csv".format(outdir=config["outdir"]),
        qtlsG="{outdir}/Rqtl/qtlsG.csv".format(outdir=config["outdir"]),
        filtered="{outdir}/Rqtl/filtered.csv".format(outdir=config["outdir"]),
        stats="{outdir}/Rqtl/stats.txt".format(outdir=config["outdir"]),
        plots="{outdir}/Rqtl/7.bulk_plots.pdf".format(outdir=config["outdir"])
    params:
        R_bulk_size = config["R_bulk_size"],
        S_bulk_size = config["S_bulk_size"],
        nb_takagi_reps = config["nb_takagi_reps"],
        filter_threshold = config["filter_threshold"],
        window_size = config["window_size"],
        false_discovery_rate_G = config["false_discovery_rate_G"],
        min_depth_in_bulk = config["min_depth_in_bulk"],
        max_depth_in_bulk = config["max_depth_in_bulk"],
        bind = config["BIND"],
        QTL_bin = config["QTL_bin"]
    shell:
        """
        singularity exec {params.bind} {params.QTL_bin} Rscript do_qtl.R \
        {input.snpfile} \
        {output.stats} \
        {output.filtered} \
        {output.qtlsT} \
        {output.qtlsG} \
        {output.plots} \
        {params.R_bulk_size} \
        {params.S_bulk_size} \
        {params.nb_takagi_reps} \
        {params.filter_threshold} \
        {params.window_size} \
        {params.false_discovery_rate_G} \
        {output.jpg6} \
        {params.min_depth_in_bulk} \
        {params.max_depth_in_bulk}
        """