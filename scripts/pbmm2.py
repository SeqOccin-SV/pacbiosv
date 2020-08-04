#/bin/env python3

import os
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


###################################################################
### Script
###################################################################

command = ' '.join(['pbmm2 align', snakemake.config['ref'], snakemake.input.file, snakemake.output.bam, '--sort'])

if (get_type(snakemake.input.file)==2): # only if CCS
	command += ' --preset CCS'

if (get_type(snakemake.input.file)==3): # only if from fastq, not ok for subreads.bam
	command += ' '.join(' --sample', snakemake.sample,'--rg', snakemake.params.rg)

command += ' '.join([' >', snakemake.log.stdout, '2>', snakemake.log.stderr])

print(command)

os.system(command)

