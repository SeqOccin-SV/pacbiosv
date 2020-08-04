# PacBio SV detection pipeline

### Running PacBio SV detection pipeline on your dataset

#### Setting python env to have access to PacBio helper script

# check for module load on genotoul
conda env create -p ./pbsv -f envs/pbsv_env.yaml
conda activate ./pbsv/

#### Running the pipeline

##### Snakemake

```bash
sbatch -j 4 ./launch.pbsv.smkj
```

##### Bash command summary

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

