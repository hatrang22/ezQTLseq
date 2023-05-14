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
#     QTLSeq pipeline for different sequencing mode: Paired-End(PE), Single-End(SE), Single Primer Enrichment Technology(SPET)
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
    ml=list()
    if list(tt):  # if len(tt) > 0
        fs=re.split('[ ;,]+', tt.item())
        for f in fs:
            if f.endswith(".gz"):
                f.strip()
                ml.append(config["fq_dir"]+"/"+f)
    return ml

# Standardize samples
def standardize_samples(df, valid_mode=["pe", "se", "spet", "spetnoumi"], require_col=["SampleName", "mode", "fq1", "fq2"]):
    # check column names (all values of require_col must be presented in df)
    if not all([col_name in df.columns for col_name in require_col]):
        raise ValueError(f"Missing column(s) of the sample file. Required columns: {require_col}")

    df['mode'] = df['mode'].apply(lambda x: x.lower())  # convert mode to lowercase
    # check valid mode (all sample modes of df must be a valid mode)
    if not all([mode in valid_mode for mode in df['mode'].to_list()]):
        raise ValueError(f"Unvalid mode found in the sample file. Valid mode: {valid_mode}")

    # check if missing forward read file(s)
    if any(df.isna()['fq1']):
        raise ValueError("Missing forward read file(s) in the sample file")

    # check if missing reverse read (or umi) file(s) for pe (or spet)
    if any(df[df["mode"]=="pe"].isna()['fq2']):
        raise ValueError("Missing reverse read file(s) for paired end in the sample file")
    if any(df[df["mode"]=="spet"].isna()['fq2']):
        raise ValueError("Missing UMI file(s) for SPET in the sample file")

    return df

# Generate .fai file if not exist
def gen_fai():	
    path_fai = config['REFPATH'] + "/" + config['GENOME_FAI'] 	
    bind = config["BIND"]	
    samtools_bin = config["samtools_bin"]	
    genome = config["REFPATH"] + "/" + config["GENOME"] 	
    if not os.path.exists(path_fai):	
        shell("singularity exec {bind} {samtools_bin} samtools faidx {genome}")  # create .fai file
    

##########################
## GLOBAL VARIABLES
##########################

# Read the sample file using pandas lib (sample names+ fastq names) and create index using the sample name
samples = pd.read_csv(config["samplesfile"], sep='\t', dtype=str, comment='#').set_index(["SampleName"], drop=False)
# Standardize samples
samples = standardize_samples(samples)

# Generate .fai file if not exist
gen_fai()

# Build a chr list based on GENOME_FAI
CHRS=list()
with open(config["REFPATH"]+"/"+config["GENOME_FAI"],'r') as fh:
    for line in fh:
      line=line.strip()
      fs=re.split(r'\t',line)
      CHRS.append(fs[0])
fh.close()

###################################################################################
rule final_outs:
    input:
        # expand("{outdir}/fastp/{sample[0]}_{sample[1]}_trim.fastq.gz", outdir=config["outdir"], sample=zip(samples['SampleName'], samples['mode'])),
        # expand("{outdir}/fastqc/{sample[0]}_{sample[1]}.OK.done", outdir=config["outdir"],sample=zip(samples['SampleName'], samples['mode'])),
        "{outdir}/multiqc/multiqc_report_fastqc.html".format(outdir=config["outdir"]),
        # expand("{outdir}/mapped/raw/{sample[0]}_{sample[1]}_sorted.raw.bam", outdir=config["outdir"],sample=zip(samples['SampleName'], samples['mode'])),
        # expand("{outdir}/mapped/raw/{sample[0]}_{sample[1]}_sorted.raw.bam.bai", outdir=config["outdir"],sample=zip(samples['SampleName'], samples['mode'])),
        # expand("{outdir}/mapped/{sample[0]}_{sample[1]}_sorted.bam", outdir=config["outdir"],sample=zip(samples['SampleName'], samples['mode'])),
        # expand("{outdir}/mapped/{sample[0]}_{sample[1]}.stats.txt", outdir=config["outdir"],sample=zip(samples['SampleName'], samples['mode'])),
        "{outdir}/multiqc/multiqc_report_bam.html".format(outdir=config["outdir"]),
        #"{refp}/{ref}.dict".format(refp=config["REFPATH"],ref=os.path.splitext(config["GENOME"])[0])
        # expand("{outdir}/variant/gatk_gvcf/{sample[0]}_{sample[1]}-{mychr}.g.vcf.gz",outdir=config["outdir"],sample=zip(samples['SampleName'], samples['mode']), mychr=CHRS),
        # expand("{outdir}/variant/gvcf_{mychr}_list.map", outdir=config["outdir"], mychr=CHRS),
        # expand("{outdir}/variant/gatk_genomicsdb_{mychr}.ok", outdir=config["outdir"], mychr=CHRS),
        # "{outdir}/variant/gatk_all.vcf.gz".format(outdir=config["outdir"]),
        # "{outdir}/variant/gatk_all.score_snps.pdf".format(outdir=config["outdir"]),
        # "{outdir}/variant/gatk_all.filtered_snps.vcf.gz".format(outdir=config["outdir"]),
        # "{outdir}/variant/gatk_all.keep_biallele.vcf.gz".format(outdir=config["outdir"]),
        # "{outdir}/variant/gatk_all.keep_filter_dp.vcf.gz".format(outdir=config["outdir"]),
        # "{outdir}/variant/gatk_all.keep_snps.vcf.gz".format(outdir=config["outdir"]),
        # "{outdir}/variant/gatk_all.filter_P_snps.vcf.gz".format(outdir=config["outdir"]),
        # "{outdir}/variant/gatk_all.keep_snps.stats.txt".format(outdir=config["outdir"]),
        # "{outdir}/Rqtl/gatk_all.filter_bulk_rs.table".format(outdir=config["outdir"]),
        "{outdir}/Rqtl/6.Takagi.jpg".format(outdir=config["outdir"]),
        "{outdir}/Rqtl/qtlsT.csv".format(outdir=config["outdir"]),
        "{outdir}/Rqtl/qtlsG.csv".format(outdir=config["outdir"]),
        "{outdir}/Rqtl/filtered.csv".format(outdir=config["outdir"]),"{outdir}/Rqtl/stats.txt".format(outdir=config["outdir"])
     

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
    run:
        if params.mode != "pe":  # if not pe, namely, se, spet, spetnoumi
            shell(
                    "singularity exec {params.bind} {params.fastp_bin} fastp \
                    --stdin \
                    -i <(zcat {input.R1}) \
                    -o {output.R} \
                    -h {params.html} \
                    -j {params.json} \
                    --cut_mean_quality 20 \
                    --cut_window_size 4 \
                    --low_complexity_filter \
                    --complexity_threshold 30 \
                    -w {threads}"
                )
            shell("touch {output.R1}")  # create pseudo-R1-file for inputs of the new rule
            shell("touch {output.R2}")  # create pseudo-R2-file for inputs of the new rule
        else:
            shell(
                    "singularity exec {params.bind} {params.fastp_bin} fastp \
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
                    -w {threads}"
                )
            shell("touch {output.R}")  # create pseudo-R-file for inputs of the new rule


# 1-2) Check quality control using FastQC: for each individual
#QC
#fastqc for differents modes
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
    run:
        if params.mode != "pe":  # if not pe, namely, se, spet, spetnoumi
            shell("singularity exec {params.bind} {params.fastqc_bin} fastqc -o {params.outdir}/fastqc -t {threads} {input.R} && touch {output}")
        else:
            shell("singularity exec {params.bind} {params.fastqc_bin} fastqc -o {params.outdir}/fastqc -t {threads} {input.R1} {input.R2} && touch {output}")

# 1-3) check quality control using FastQC: a global multiQC for the analysis
# multiQC on fastqc outputs
rule multiqc_fastqc:
    input:
        expand("{outdir}/fastqc/{sample[0]}_{sample[1]}.OK.done", outdir=config["outdir"], sample=zip(samples['SampleName'], samples['mode']))
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
        singularity exec {params.bind} {params.multiqc_bin} multiqc --filename {output} {params.outdir}/fastqc
        """


###############################################################################
########################  DNA mapping   #######################################
###############################################################################
# 2-1) reference indexation
rule bwa_index:
    input:
        genome = config["REFPATH"] + "/" + config["GENOME"] 
    output:
       config["REFPATH"] + "/" + config["GENOME"]  + ".amb",
       config["REFPATH"] + "/" + config["GENOME"]  + ".ann",
       config["REFPATH"] + "/" + config["GENOME"]  + ".bwt",
       config["REFPATH"] + "/" + config["GENOME"]  + ".pac",
       config["REFPATH"] + "/" + config["GENOME"]  + ".sa"
    params:
        bind         = config["BIND"],
        bwa_bin      = config["bwa_bin"],
        samtools_bin = config["samtools_bin"]
    shell:
        """
        singularity exec {params.bind} {params.bwa_bin} bwa index -a bwtsw -b 500000000 {input.genome}
        """

# 2-2) mapping
# ############################################################
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

rule bwa_mapping_wow_merge:
    input:
        R = "{outdir}/fastp/{{sample}}_{{mode}}_trim.fastq.gz".format(outdir=config["outdir"]),
        R1 = "{outdir}/fastp/{{sample}}_{{mode}}_1_trim.fastq.gz".format(outdir=config["outdir"]),
        R2 = "{outdir}/fastp/{{sample}}_{{mode}}_2_trim.fastq.gz".format(outdir=config["outdir"]),
        #fake input used to force index building before alignement if not present
        idx = config["REFPATH"] + "/" + config["GENOME"]  + ".bwt"
    output:
        bam   = "{outdir}/mapped/raw/{{sample}}_{{mode}}_sorted.raw.bam".format(outdir=config["outdir"]),
        bai   = "{outdir}/mapped/raw/{{sample}}_{{mode}}_sorted.raw.bam.bai".format(outdir=config["outdir"]),
    params:
        modules           = config["MODULES"],
        fastp_bin         = config["fastp_bin"],
        json              = config["outdir"]+"/fastp/{sample}_{mode}_trim.json",
        html              = config["outdir"]+"/fastp/{sample}_{mode}_trim.html",
        mode              = "{mode}",
        outtmp            = "{outdir}/mapped/{{sample}}_{{mode}}".format(outdir=config["outdir"]),
        idxbase           = config["REFPATH"] + "/" + config["GENOME"] ,
        bind              = config["BIND"],
        bwa_bin           = config["bwa_bin"],
        samtools_bin      = config["samtools_bin"],
        rg                = "@RG\\tID:{sample}\\tSM:{sample}"
    threads: 8
    run:
        if params.mode != "pe":  # if not pe, namely, se, spet, spetnoumi
            shell(
                    "singularity exec {params.bind} {params.bwa_bin} bwa mem \
                    -t {threads} \
                    -K 100000000 \
                    -R '{params.rg}' \
                    {params.idxbase} \
                    {input.R} \
                    |singularity exec {params.bind} {params.samtools_bin} samtools view -h -F 2048 \
                    |singularity exec {params.bind} {params.samtools_bin} samtools sort -@2 -m 6G -o {output.bam}"
                )
            shell("singularity exec {params.bind} {params.samtools_bin} samtools flagstat -@ {threads} {output.bam}")
            shell("singularity exec {params.bind} {params.samtools_bin} samtools index -@ {threads} -b {output.bam}")
            shell(" rm -rf {input.R1} {input.R2}")  # remove pseudo-files of the previous rule

        else:
            #merge PE
            shell(
                    "singularity exec {params.bind} {params.fastp_bin} fastp \
                    -i {input.R1} \
                    -I {input.R2} \
                    --merge \
                    --correction \
                    --merged_out {params.outtmp}.merged.tmp.fastq.gz \
                    --out1 {params.outtmp}_1.unmerged.tmp.fastq.gz \
                    --out2 {params.outtmp}_2.unmerged.tmp.fastq.gz \
                    --json={params.json} \
                    --html={params.html} \
                    --length_required 50 \
                    --thread {threads}"
                )
            shell("rm -f {params.html} {params.json}")
            
            #mapping
            #2 stage : first with unmerged PE and second with merged reads
            shell(
                    "singularity exec {params.bind} {params.bwa_bin} bwa mem \
                    -t {threads} \
                    -K 100000000 \
                    -R '{params.rg}' \
                    {params.idxbase} \
                    {params.outtmp}_1.unmerged.tmp.fastq.gz {params.outtmp}_2.unmerged.tmp.fastq.gz \
                    | singularity exec {params.bind} {params.samtools_bin} samtools view -h -F 2048 \
                    | singularity exec {params.bind} {params.samtools_bin} samtools sort -@2 -m 6G -o {params.outtmp}.um.tmp.bam -"
                )
            shell(
                    "singularity exec {params.bind} {params.bwa_bin} bwa mem \
                    -t {threads} \
                    -K 100000000 \
                    -R '{params.rg}' \
                    {params.idxbase} \
                    {params.outtmp}.merged.tmp.fastq.gz \
                    |singularity exec {params.bind} {params.samtools_bin} samtools view -h -F 2048 \
                    |singularity exec {params.bind} {params.samtools_bin} samtools sort -@2 -m 6G -o {params.outtmp}.m.tmp.bam -"
                )
            shell("rm -rf {input.R}")  # remove pseudo-file of the previous rule

            # merge bams
            shell("singularity exec {params.bind} {params.samtools_bin} samtools merge {output.bam} {params.outtmp}.m.tmp.bam {params.outtmp}.um.tmp.bam && rm -f {params.outtmp}.m.tmp.bam {params.outtmp}.um.tmp.bam")
            shell("singularity exec {params.bind} {params.samtools_bin} samtools flagstat -@ {threads} {output.bam}")
            shell("singularity exec {params.bind} {params.samtools_bin} samtools index -@ {threads} -b {output.bam}")
            shell("rm -rf {params.outtmp}.merged.tmp.fastq.gz")
            shell("rm -rf {params.outtmp}_1.unmerged.tmp.fastq.gz  {params.outtmp}_2.unmerged.tmp.fastq.gz")


# 2-3) mark end remove duplicate PCR 
## For SPET with UMI: https://github.com/tecangenomics/nudup

rule remove_or_keep_duplicate: 
    input:
        bam   = "{outdir}/mapped/raw/{{sample}}_{{mode}}_sorted.raw.bam".format(outdir=config["outdir"]),
        R     = get_fq2
    output:
        bam   = "{outdir}/mapped/{{sample}}_{{mode}}_sorted.bam".format(outdir=config["outdir"]),
        bai   = "{outdir}/mapped/{{sample}}_{{mode}}_sorted.bam.bai".format(outdir=config["outdir"])
    params:
        tmpdir      = config["outdir"] + "/mapped/tmpdir",
        bamutil_bin  = config["bamUtil_bin"],
        probe_length = config["probe_length"],
        samtools_bin = config["samtools_bin"],
        bind         = config["BIND"],
        outdir       = config["outdir"],
        sambamba_bin = config["sambamba_bin"], # sambamba for markdup
        bam_rm_dup   = config["bam_remove_duplicate"],
        R2_concat    = "{outdir}/mapped/{{sample}}_R2.fastq.gz".format(outdir=config["outdir"]),
        nudup_bin    = config["nudup_bin"],
        mode         = "{mode}",
        prefix       = "{outdir}/mapped/{{sample}}_{{mode}}".format(outdir=config["outdir"]),
        outtmp       = "{outdir}/mapped/{{sample}}_{{mode}}".format(outdir=config["outdir"]),
        dedup        = "{outdir}/mapped/{{sample}}_{{mode}}.sorted.dedup.bam".format(outdir=config["outdir"]),
    threads: 10
    message: "Mark/remove or keep PCR duplicates.\n"
    run:
        if params.mode == "spet":
            shell("cat {input.R} > {params.R2_concat}")
            shell("singularity exec {params.bind} {params.nudup_bin} nudup.py -f {params.R2_concat} -o {params.prefix} {input.bam} && rm {params.R2_concat}")
            shell("singularity exec {params.bind} {params.bamutil_bin} bam TrimBam {params.dedup} {params.outtmp}.clipped.tmp.bam --clip -L {params.probe_length}")
            shell("singularity exec {params.bind} {params.samtools_bin} samtools sort -o {output.bam} {params.outtmp}.clipped.tmp.bam")
            shell("singularity exec {params.bind} {params.samtools_bin} samtools index -@ {threads} -b {output.bam}")
        elif params.mode == "pe" and params.bam_rm_dup == "y":
            shell("ulimit -n 4048")
            shell("mkdir -p {params.tmpdir}")
            #singularity exec {params.bind} {params.samtools_bin} samtools rmdup -s {input.bam} {output.bam}
            shell(
                    "singularity exec {params.bind} {params.sambamba_bin} sambamba markdup --remove-duplicates --tmpdir={params.tmpdir} -t {threads} {input.bam} {output.bam} && \
                    rm -f {input.bam} && \
                    singularity exec {params.bind} {params.samtools_bin} samtools index -@ {threads} -b {output.bam}"
                )
            #--remove-duplicates #remove duplicates instead of just marking them
            shell("rm -rf {params.tmpdir}")
        elif params.mode == "spetnoumi":
            shell("mv {input.bam} {output.bam}")
            shell("singularity exec {params.bind} {params.bamutil_bin} bam TrimBam {output.bam} {params.outtmp}.clipped.tmp.bam --clip -L {params.probe_length}")
            shell("singularity exec {params.bind} {params.samtools_bin} samtools sort -o {output.bam} {params.outtmp}.clipped.tmp.bam")
            shell("singularity exec {params.bind} {params.samtools_bin} samtools index -@ {threads} -b {output.bam}")
        else:
            shell("mv {input.bam} {output.bam}")
            shell("singularity exec {params.bind} {params.samtools_bin} samtools index -@ {threads} -b {output.bam}")
        

# 2_3) bam stats 
rule bam_stats:
    input:
        bam   = "{outdir}/mapped/{{sample}}_{{mode}}_sorted.bam".format(outdir=config["outdir"]),
    output:
        stats = "{outdir}/mapped/{{sample}}_{{mode}}.stats.txt".format(outdir=config["outdir"])
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
        expand("{outdir}/mapped/{sample[0]}_{sample[1]}.stats.txt", outdir=config["outdir"], sample=zip(samples['SampleName'], samples['mode']))
    output:
        "{outdir}/multiqc/multiqc_report_bam.html".format(outdir=config["outdir"])
    threads:
        2
    params:
        outdir      = config["outdir"]+"/mapped/*.stats.txt",
        multiqc_bin = config["multiqc_bin"],
        bind        = config["BIND"]
    shell:
        """
        singularity exec {params.bind} {params.multiqc_bin} multiqc --filename {output} {input}
        """ 
################################################################################
###########################  SNP calling using  GATK ###########################
################################################################################
# 3) SNP caller using GATK4
#gatk:4 rules
# 3-1) create a reference dictionnary the ref must not be a sl and a ref.fa => ref.dict (not ref.fa.dict)
rule gatk4_ref_dict:
    input:
        ref  = "{refp}/{ref}".format(refp=config["REFPATH"],ref=config["GENOME"])
    output:
        dic  = "{refp}/{ref}.dict".format(refp=config["REFPATH"],ref=os.path.splitext(config["GENOME"])[0])
    params:
        bind = config["BIND"],
        samtools_bin = config["samtools_bin"],
        gatk4_bin = config["gatk4_bin"]
    threads: 1
    shell:
        """
        singularity exec {params.bind} {params.gatk4_bin} gatk CreateSequenceDictionary -R {input.ref} -O {output.dic}
        """
# ------------------  parrallel by chr -----------------
# 3-2) HaplotypeCaller for each sample and each chr
rule gatk4_hc:
    input:
        bam     = "{outdir}/mapped/{{sample}}_{{mode}}_sorted.bam".format(outdir=config["outdir"]),
        #refdict = "{refp}/{ref}.dict".format(refp=config["REFPATH"],ref=config["GENOME"])
    output:
        gvcf="{outdir}/variant/gatk_gvcf/{{sample}}_{{mode}}-{{mychr}}.g.vcf.gz".format(outdir=config["outdir"])
    params:
        ch="{mychr}",
        ref = config["REFPATH"] + "/" + config["GENOME"] ,
        bind = config["BIND"],
        samtools_bin = config["samtools_bin"],
        gatk4_bin = config["gatk4_bin"]
    threads: 16
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

# 3-3) Generate a map (text file) of gvcf files for each chr
rule gatk4_gvcf_map_file:
    input:
        gvcfs = expand("{outdir}/variant/gatk_gvcf/{sample[0]}_{sample[1]}-{{mychr}}.g.vcf.gz", outdir=config["outdir"], sample=zip(samples['SampleName'], samples['mode'])),   
    output:
        gvcfmap      = "{outdir}/variant/gvcf_{{mychr}}_list.map".format(outdir=config["outdir"]),
    params:
        ch          ="{mychr}",
        #csv         = expand("{sample}\t{outdir}/variant/gatk_gvcf/{sample}-{{mychr}}.g.vcf.gz", outdir=config["outdir"], sample=samples['SampleName']),
        csv         = expand("{sample[0]}\t{outdir}/variant/gatk_gvcf/{sample[0]}_{sample[1]}", outdir=config["outdir"], sample=zip(samples['SampleName'],samples['mode'])),
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
        #gvcfs = expand("{outdir}/variant/gatk_gvcf/{sample[0]}-{{mychr}}_{sample[1]}.g.vcf.gz", outdir=config["outdir"], sample=zip(samples['SampleName'], samples['mode'])),
        gvcfmap = "{outdir}/variant/gvcf_{{mychr}}_list.map".format(outdir=config["outdir"])
    output:
        "{outdir}/variant/gatk_genomicsdb_{{mychr}}.ok".format(outdir=config["outdir"])
    params:
        ch="{mychr}",
        tmpdir=config["outdir"] + "/variant/tmpgatkdir",
        gdb  = config["outdir"] + "/variant/" + "GenomicsDB_{mychr}",
        ref = config["REFPATH"] + "/" + config["GENOME"] ,
        bind = config["BIND"],
        samtools_bin = config["samtools_bin"],
        gatk4_bin = config["gatk4_bin"]
    threads: 12
    shell:
        """
        mkdir -p {params.tmpdir}
        singularity exec {params.bind} {params.gatk4_bin} \
        gatk --java-options '-Xmx8g' GenomicsDBImport \
        --genomicsdb-workspace-path {params.gdb} \
        --batch-size 200 \
        -L {params.ch} \
        --sample-name-map {input.gvcfmap} \
        --reader-threads {threads} \
        --tmp-dir={params.tmpdir} && touch {output}
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
        tmpdir=config["outdir"] + "/tmpgatkdir",
        gdb  = config["outdir"] + "/variant/" + "GenomicsDB_{mychr}",
        ref = config["REFPATH"] + "/" + config["GENOME"] ,
        bind = config["BIND"],
        samtools_bin = config["samtools_bin"],
        gatk4_bin = config["gatk4_bin"]
    threads: 2
    message: "GATK4 genotype_variants vcf\n"
    shell:
        """
        mkdir -p {params.tmpdir}
        singularity exec {params.bind} {params.gatk4_bin} \
        gatk --java-options '-Xmx8G' GenotypeGVCFs \
        --variant gendb://{params.gdb} \
        --reference {params.ref} \
        --output {output.vcf} \
        --tmp-dir={params.tmpdir}
        """
# 3-6) merge all the vcf produced by chr
rule combinevcf:
    input:
        #"{outdir}/variant/chrslist.OK".format(outdir=config["outdir"]),
        vcf=expand("{outdir}/variant/gatk_{mychr}_genotyped.vcf.gz",outdir=config["outdir"],mychr=CHRS),
        refdict = config["REFPATH"] + "/" + os.path.splitext(config["GENOME"])[0] +".dict",
    output:
        gvcf = "{outdir}/variant/gatk_all.vcf.gz".format(outdir=config["outdir"])
    params:
        mlist = expand(" -I {outdir}/variant/gatk_{mychr}_genotyped.vcf.gz",outdir=config["outdir"],mychr=CHRS),
        ref = config["REFPATH"] + "/" + config["GENOME"],
        bind = config["BIND"],
        samtools_bin = config["samtools_bin"],
        gatk4_bin = config["gatk4_bin"]
    threads: 2
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
        ref = config["REFPATH"] + "/" + config["GENOME"] ,
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
        ref = config["REFPATH"] + "/" + config["GENOME"] ,
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
        ref = config["REFPATH"] + "/" + config["GENOME"] ,
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
    threads: 2
    shell:
        """
        singularity exec {params.bind} {params.R_bin} \
        Rscript --vanilla gatk_scores_qual_raw_vs_filtered.R {input.rawsnps} {input.filteredsnps} "SNPs" {output.pdf}
        """

################################ END GATK #################################################


# 1) --------- we select variant from the previous filter and keep only biallelic samples ---------------

## 4-1 ---------------keep biallelic ------------------
rule keep_biallelic:
    input:
        snp = "{outdir}/variant/gatk_all.filtered_snps.vcf.gz".format(outdir=config["outdir"])
    output:
        biallel = "{outdir}/variant/gatk_all.keep_biallele.vcf.gz".format(outdir=config["outdir"])
    params:
        ref = config["REFPATH"] + "/" + config["GENOME"] ,
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
        ref = config["REFPATH"] + "/" + config["GENOME"] ,
        bind = config["BIND"],
        gatk4_bin = config["gatk4_bin"],
        min_dp = config["min_dp"],
        max_dp = config["max_dp"],
        dp_p = config["dp_p"]
    threads: 2
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
        ref = config["REFPATH"] + "/" + config["GENOME"] ,
        bind = config["BIND"],
        gatk4_bin = config["gatk4_bin"]
    threads: 2
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
        ref = config["REFPATH"] + "/" + config["GENOME"] ,
        bind = config["BIND"],
        gatk_bin = config["gatk_bin"]
    threads: 2
    shell:
        """
        singularity exec {params.bind} {params.gatk_bin} \
        java -jar /opt/jvarkit_github/jvarkit/dist/vcffilterjdk.jar \
        -e 'return !variant.getGenotype("P1").sameGenotype(variant.getGenotype("P2"));' -o {output.polymorph} {input.snp}
        #keep only the polymorph different with the parents 
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
    threads: 2
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
        idxbase      = config["REFPATH"] + "/" + config["GENOME"] ,
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
        ref = config["REFPATH"] + "/" + config["GENOME"] ,
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
    threads: 1
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
