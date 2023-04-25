# QTLseq Analysis Pipeline User Guide 
Welcome to the QTLseq Analysis Pipeline User Guide! This guide will provide you with step-by-step instructions on how to use the QTLseq analysis pipeline with different sequencing modes. The QTLseq analysis pipeline is designed to identify quantitative trait loci (QTLs) associated with complex traits in genetic data.


## Introduction to QTLseq Analysis Pipeline
The QTLseq analysis pipeline is a computational workflow for identifying QTLs associated with complex traits using next-generation sequencing data. It involves several steps, including read mapping, SNP calling, QTL mapping, and statistical analysis. The pipeline is designed to be flexible and customizable to accommodate different sequencing modes, such as paired-end, single-end, Single Primer Enrichment Technology (SPET) with or without UMI sequencing.


## Setting Up Your Pipeline
To set up the QTLseq analysis pipeline in Git for different sequencing modes, follow these general steps:

### Step 1: Clone the Pipeline Repository
1. Open a terminal window and navigate to the directory where you want to clone the QTLseq analysis pipeline Git repository.
2. Run the following command to clone the repository:

```bash
git clone https://forgemia.inra.fr/gafl/users/ha_trang_phung/pipeline_qtlSeq.git
```
### Step 2: Generate the Sample file
1. Look for a sample file template provided by the pipeline, which contains the necessary information for configuring the pipeline  for your specific sequencing data. The sample file must be CSV format



 
2. Open the copied sample file using a text editor or Excel, and fill in the required information for your specific sequencing data. 
    - In the "SampleName" column,  "P1", "P2", "R", and "S" are used as the unique identifiers for each of the four samples
    - In the "mode" column, you can specify one of the four available modes for your sequencing data.These modes may include paired-end(pe), single-end(se), SPET(spet), SPET without UMI(spet_no_umi). 
    - The "fq1" column stands for read1(R1/ForwardRead) and the "fq2" column stands for read2(R2/ReverseRead) in the case of Paired End or UMI in the case of SPET. 
Note that the QTLseq analysis pipeline is capable of processing multiple input files for each sample simultaneously, if needed.To specify multiple input files for a sample in the sample file, simply separate them with a semicolon (;) in the "fq1" and "fq2" column

        **IMPORTANT**:Please make sure the order of input files specified in the "fq1" column of the sample file need to corresponds correctly with the order of input files in the "fq2", if you are working with paired-end sequencing data or SPET with UMI sequencing data. For single-end reads or SPET without UMI, you can simply put the file names in the 'fq1' column and leave the 'fq2' column empty.
        
    - You can also generate automatically your sample file by using the _gen_samples.py_. The _gen_sample.py_ file is a Python script that takes input parameters such as SampleName, mode and fq1 et fq2 inputs files 
...........................
```bash
python gen_samples.py
```
### Step 3: Configuring the Config file

The config file contains the paths or URLs to the sample file, data directory, and reference genome, as well as various parameters for SNP calling and QTLseq R. It may also include container paths and parameters for binding Singularity containers for software dependencies.

In the config file, you can edit the sample file path, reference genome path, and other parameters to customize the analysis pipeline according to your needs.



### Step 4: Running the analysis pipeline
.

## Clean outputs
```bash
rm -rf logs/ out/ .snakemake/
```

## Quicktest with dryrun

```bash
snakemake --configfile config.yaml -s Snakefile.smk -np
```

## Run pipeline

Direct run:

```bash
./run_ezqtlseq_pipeline.slurm
```

Submit job:

```bash
sbatch run_ezqtlseq_pipeline.slurm
```

The pipeline workflow is available [here](dag.pdf).
