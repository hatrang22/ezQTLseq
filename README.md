# Pipeline QTLseq

In development...

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
