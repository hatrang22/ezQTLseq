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

(TODO...)

## Generate samples

```bash
python gen_samples.py
```

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
