#!/bin/bash
# run ezQTLseq snakemake workflow
# Pipeline features:
# 1) read preprocessing, reads QC,
# 2) mapping with bwa,mapping QC
# 3) SNP calling and filtering using GATK
# 4) do QTLseqR analysis
#
# Snakemake features: fastq from csv file, config, modules, SLURM
# use of cluster: SLURM, sample file
# using a config file in yaml
# using ressources CPU and memory
# Using cluster with SLURM cluster parametrer in a cluster.json
# singularity images
# Samples (samples names and fastq files) in csv file

#SBATCH --job-name=ezqtlseq # job name (-J)
#SBATCH --time="96:00:00" #max run time "hh:mm:ss" or "dd-hh:mm:ss" (-t)
#SBATCH --cpus-per-task=1 # max nb of cores (-c)
#SBATCH --ntasks=1 #nb of tasks
#SBATCH --mem=6G # max memory (-m)
#SBATCH --output=ezqtlseq.out #stdout (-o)

########################## On genotoul ###############################
# Uncomment the module load for genotoul
## snakemake 5.8
#module load bioinfo/snakemake-5.8.1
## for singularity
#module load system/singularity-3.5.3
######################################################################

########################## On ifb-core ###############################
## Uncomment the next line
#module load snakemake
######################################################################


# All rules included in one snakemake file:
#RULES=Snakefile_para_gatk_PEmerged.smk
RULES=Snakefile.smk

# config file to be edited
CONFIG=config.yaml

# slurm directive by rule (can be edited if needed)
CLUSTER_CONFIG=cluster.json

# sbatch directive to pass to snakemake
CLUSTER='sbatch --mem={cluster.mem} -t {cluster.time} -c {cluster.cores} -J {cluster.jobname} -o logs/{cluster.out}'
# Maximum number of jobs to be submitted at a time (see cluster limitation)
MAX_JOBS=500

# Full clean up:
#rm -fr .snakemake logs *.out *.html *.log *.pdf
# Clean up only the .snakemake
rm -fr .snakemake

# log directory is mandatory (see $CLUSTER) else slurm jobs failed but not the master job
mkdir -p logs

# Generate the dag files
# With samples
#snakemake --configfile $CONFIG -s $RULES --dag | dot -Tpdf > dag.pdf
# Only rules
snakemake --configfile $CONFIG -s $RULES --rulegraph | dot -Tpdf > dag.pdf

# Dry run (simulation)
#snakemake --configfile $CONFIG -s $RULES -np -j $MAX_JOBS --cluster-config $CLUSTER_CONFIG --cluster "$CLUSTER" > snakemake_dryrun.out

# Full run (if everething is ok: uncomment it)
snakemake --configfile $CONFIG -s $RULES -p -j $MAX_JOBS --cluster-config $CLUSTER_CONFIG --cluster "$CLUSTER" 


# If latency problem add to the run:
# --latency-wait 60

exit 0
