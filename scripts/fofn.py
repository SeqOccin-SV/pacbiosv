#/bin/env python3

import os
import pandas as pd
import pathlib
from pprint import pprint

###################################################################
### Method
###################################################################

def get_type(fofn):
	print('Inside get_type\n')
	pprint(fofn)
	type = 0
	if (fofn.endswith('subreads.fofn') ):
		type = 1
	elif (fofn.endswith('ccs.fofn') ):
        	type = 2
	elif (fofn.endswith('fastqccs.fofn') ):
        	type = 3
	return type

def get_files(pdObject):
	''' Parse the samples.txt from config file
	'''
	files = pdObject.loc[pdObject.sample, 'path'].split(',')
	return(files)

def get_suffix(string):
	''' Return path file suffixes
	'''
	return(''.join(pathlib.Path(string).suffixes))
	
###################################################################
### Script
###################################################################

# ~ def write_fofn(wildcards):
# ~ ''' write a fofn file based on samples.txt information
# ~ the file is named depending on submitted datatypes for CLR and CCS differentiation
# ~ '''

# VAR
# Present in Snakefile. Dont know how to pass it to script
samples = pd.read_table(snakemake.config['samples']).set_index('sample', drop=False)

files = get_files(samples)
fofnName = ''
type = 0
suffixes = map(get_suffix, files)
suffixSet = list(set(suffixes))

pathlib.Path("fofn").mkdir(parents=True, exist_ok=True)

if (len(suffixSet) > 1):
	raise Exception('''You got several datatypes for sample {}. 
You have to analyze CLR and HiFi data separatedly.'''.format(snakemake.sample))
else:
	if (suffixSet[0].endswith('subreads.bam')):
		fofnName = 'fofn/'+snakemake.sample+'_subreads.fofn'
		type = 1
	elif (suffixSet[0].endswith('ccs.bam')):
		fofnName = 'fofn/'+snakemake.sample+'_bamccs.fofn'
		type = 2
	elif (suffixSet[0].endswith('fastq.gz') or suffixSet[0].endswith('fq.gz') or suffixSet[0].endswith('fastq') or suffixSet[0].endswith('fq')):
		fofnName = 'fofn/'+snakemake.sample+'_fastqccs.fofn'
		type = 3

with open(fofnName, 'w') as fh:
	for f in files:
		fh.write(f+"\n")
		
return(fofnName)

