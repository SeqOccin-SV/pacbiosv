#!/usr/bin/env python

# ref = '/home/adf/seqoccin/data/reference/hs37d5_hsa10.fa'
import pandas as pd
import pathlib
from pprint import pprint

# VAR
samples = pd.read_table(config['samples'], comment='#').set_index('sample', drop=False)
# ~ seqtypes = ['subreads','bamccs','fastqccs']

###############################################################
### Methods
###############################################################

def get_type(wildcards):
	print('Inside get_type\n')
	pprint(wildcards)
	type = 0
	if (wildcards.input.file.endswith('subreads.fofn') ):
		type = 1
	elif (wildcards.input.file.endswith('bamccs.fofn') ):
		type = 2
	elif (wildcards.input.file.endswith('fastqccs.fofn') ):
		type = 3
	return type

def get_preset(wildcards):
	return '--preset CCS' if (get_type(wildcards) == 2 or get_type(wildcards) == 3) else ''

def get_info(wildcards):
	return '--sample {sample} --rg {params.rg}' if (get_type(wildcards) == 3) else ''

def get_files(wildcards):
	''' Parse the samples.txt from config file
	'''
	files = samples.loc[wildcards.sample, 'path'].split(',')
	return(files)

def get_suffix(string):
	''' Return path file suffixes
	'''
	return(''.join(pathlib.Path(string).suffixes))

def get_threads(rule, default):
	cluster_config = snakemake.workflow.cluster_config
	if rule in cluster_config and "threads" in cluster_config[rule]:
		return cluster_config[rule]["threads"]
	elif "default" in cluster_config and "threads" in cluster_config["default"]:
		return cluster_config["default"]["threads"]
	else:
		return default
	
def write_fofn(wildcards):
	''' write a fofn file based on samples.txt information
	the file is named depending on submitted datatypes for CLR and CCS differentiation
	'''
	# ~ pprint('in Write_fofn')
	
	files = get_files(wildcards)
	fofnName = ''
	type = 0
	suffixes = map(get_suffix, files)
	suffixSet = list(set(suffixes))
	
	pathlib.Path("fofn").mkdir(parents=True, exist_ok=True)
	
	if (len(suffixSet) > 1):
		raise Exception('''You got several datatypes for sample {}. 
	You have to analyze CLR and HiFi data separatedly.'''.format(wildcards.sample))
	else:
		if (suffixSet[0].endswith('subreads.bam')):
			fofnName = 'fofn/'+wildcards.sample+'_subreads.fofn'
			type = 1
		elif (suffixSet[0].endswith('ccs.bam')):
			fofnName = 'fofn/'+wildcards.sample+'_bamccs.fofn'
			type = 2
		elif (suffixSet[0].endswith('fastq.gz') or suffixSet[0].endswith('fq.gz') or suffixSet[0].endswith('fastq') or suffixSet[0].endswith('fq')):
			fofnName = 'fofn/'+wildcards.sample+'_fastqccs.fofn'
			type = 3
	
	with open(fofnName, 'w') as fh:
		for f in files:
			fh.write(f+"\n")
			
	return(fofnName)

###############################################################
### Wildcards
###############################################################
# ~ wildcard_constraints:
	# ~ seqtype="[subreads|bamccs|fastqccs]"

###############################################################
### Rules
###############################################################

rule all:
	''' Generic all rule to launch the full pipeline
	'''
	input:
	# do not get the default rule thingy
		expand("calling/{sample}-pbmm2.vcf", sample=samples.index)

# Could add a rule to produce an index for reference 
# Test if it improves computation time

# ~ rule write_fofn:
	# ~ ''' rule to produce a fofn file for pbmm2
	# ~ '''
	# ~ input:
		# ~ write_fofn
	# ~ output:
		# ~ fofn="fofn/{sample}.fofn"
	# ~ log:
		# ~ stdout="logs/fofn/{sample}.out",
		# ~ stderr="logs/fofn/{sample}.log"
	# ~ shell:
		# ~ ' '

rule pbmm2:
	''' rule to align reads with minimap2 wrapper from fastq file
	generally, fastq file input is used for HiFi reads coming out of CCS
	'''
	input:
		file=write_fofn
	output:
		# ~ bam=temp("mapping/{sample}-pbmm2.bam")
		bam=temp("mapping/{sample}-pbmm2.bam")
	log:
		stdout="logs/pbmm2/{sample}.out",
		stderr="logs/pbmm2/{sample}.log"
	benchmark:
		"bench/{sample}.pbmm2.benchmark.txt"
	params:
		rg = '@RG\\tID:movie{sample}\\tSM:{sample}',
	conda:
		'envs/pbsv_env.yaml'
	threads: get_threads('pbmm2',10)
	script:
		'scripts/pbmm2.py'

# rule bam_sort:
# 	''' Rule to sort bam file from pbmm2
# 	This is done in a separate rule instead of using pbmm2 --sort
# 	for memory usage optimisation
# 	'''
# 	input:
# 		"mapping/{sample}-pbmm2.bam"
# 	output:
# 		"mapping/{sample}-sorted-pbmm2.bam"
# 	log:
# 		"logs/samtools/{sample}-sort.log"
# 	conda:
# 		'envs/pbsv_env.yaml'
# 	threads: get_threads('bam_sort',10)
# 	shell:
# 		"export TMPDIR=./ ;"
# 		"samtools sort -@ {threads} -o {output} {input} 2> {log}"
 
rule bam_index:
	input:
		"mapping/{sample}-pbmm2.bam"
	output:
		"mapping/{sample}-pbmm2.bam.bai"
	log:
		"logs/samtools/{sample}-index.log"
	conda:
		'envs/pbsv_env.yaml'
	threads: get_threads('bam_index',4)
	shell:
		"samtools index -@ {threads} {input} 2> {log}"

rule pbsv_discover:
	''' first rule for sv detection, use bam to look for regions with possible variants
	'''
	input:
		bam="mapping/{sample}-pbmm2.bam",
		bai="mapping/{sample}-pbmm2.bam.bai"
	output:
		"calling/{sample}-pbmm2.svsig.gz"
	log:
		"logs/pbsv/{sample}_discover.log"
	conda:
		'envs/pbsv_env.yaml'
	threads: get_threads('pbsv_discover',4)
	shell:
		"pbsv discover {input.bam} {output}"
		" 2> {log}"

rule pbsv_call:
	''' second rule for sv detection, use .svsig.gz information to call variants
	'''
	input:
		"calling/{sample}-pbmm2.svsig.gz"
	output:
		"calling/{sample}-pbmm2.vcf"
	log:
		"logs/pbsv/{sample}_call.log"
	benchmark:
		"bench/{sample}.pbsvCall.benchmark.txt"
	conda:
		'envs/pbsv_env.yaml'
	threads: get_threads('pbsv_call',10)
	shell:
		"pbsv call -j {threads} "+config['ref']+" {input} {output}"
                # specific parameter for CCS, easier to consider a SV" --ccs"
		" 2> {log}"
