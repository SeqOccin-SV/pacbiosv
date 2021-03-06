# PacBio SV detection pipeline

## Running PacBio SV detection pipeline on your dataset

To run this snakemake pipeline, you need to have access to an installation of python3 with snakemake 5 module installed.

#### Setting python env to have access to PacBio helper script

```bash
# check for module load on genotoul
conda env create -p ./pbsv -f envs/pbsv_env.yaml
conda activate ./pbsv/
```

## Running the pipeline

### Snakemake

First, modify samples.tsv/bamSamples.tsv and config.yaml

Example samples.tsv
|sample|path|
|-------|----|
|name_CLR|1.subreads.bam,2.subreads.bam|
|name_HiFi1|1.ccs.bam,1.ccs.bam|
|name_HiFi2|1.fq.gz,2.fq.gz|

Example bamSamples.tsv

|sample|bam_path|
|-------|----|
|name|path/mapped.bam|

then use the launch script

```bash
sbatch -j 4 ./launch.pbsv.smkj
```

##### Summary of the bash commands executed by the pipeline

###### Using pbmm2 to run minimap2 with preset for CCS
```bash
pbmm2 align {ref} {reads} {output.bam} --sort [--preset CCS [--sample {sample} --rg '@RG\tID:movie{sample}']]
samtools index {output.bam}
```

###### Running Pacbio Detection
```bash
pbsv discover {output.bam} {output.svsig.gz}
pbsv call {ref} {output.svsig.gz} {output.vcf}
```

