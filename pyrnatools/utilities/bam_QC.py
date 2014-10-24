#!/usr/bin/python

########################################################################
# 23 September 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################


##Requires RSeQC which I have uninstalled due to issues with pysam and bx multiple installations. Fix these problems and this will work!
import subprocess
import sys, re, os
import ConfigParser
import itertools
import argparse
from collections import defaultdict
import rnaseq_toolkit
from multiprocessing import Pool

#All these come from the RSeqQC module

def total_qc_of_bam(bam_file, genome, paired, conditions):
	outdir = conditions[bam_file]
	command = "bam_stat.py -i {} > {}".format(bam_file, outdir+"/read_mapping_statistics.txt")
	#subprocess.call(command, shell=True)
	#This program is used to calculate reads mapping statistics from provided BAM file. This script determines "uniquely mapped reads" from mapping quality, 
	#which quality the probability that a read is misplaced (Do NOT confused with sequence quality, sequence quality measures the probability that a base-calling was wrong) .
	command2 = "clipping_profile.py -i {} -o {}".format(bam_file, outdir+"/clipping_profile")
	#subprocess.call(command2, shell=True)
	#This program is used to estimate clipping profile of RNA-seq reads from BAM or SAM file. 
	#Note that to use this funciton, CIGAR strings within SAM/BAM file should have "S" operation (This means your reads aligner should support clipped mapping).
	command3 = "infer_experiment.py -r {} -i {} > {}".format(genome, bam_file, outdir+"/experiment_type.txt")
	subprocess.call(command3, shell=True)
	#Consult http://rseqc.sourceforge.net/#download-rseqc to understand results
	if paired:
		command4 = "inner_distance.py -i {} -o {} -r {}".format(bam_file, outdir+"/inner_distance", genome)
		subprocess.call(command4.split(), shell=False)
	command5 = "read_distribution.py  -i {} -r {} > {}".format(bam_file, genome, outdir+"/read_distribution.txt")
	subprocess.call(command5, shell=True)
	command6 = "read_duplication.py -i {} -o {}".format(bam_file, outdir+"/read_duplication")
	subprocess.call(command6.split(), shell=False)
	command7 = "read_GC.py -i {} -o {}".format(bam_file, outdir+"/GC_content")
	subprocess.call(command7.split(), shell=False)
	command8 = "read_quality.py -i {} -o {}".format(bam_file, outdir+"/read_quality")
	subprocess.call(command8.split(), shell=False)
	command9 = "read_NVC.py -i {} -o {}".format(bam_file, outdir+"/read_NVC")
	subprocess.call(command9.split(), shell=False)

def genebody():
	command3 = "geneBody_coverage.py -r {}.housekeeping.bed -i {} -o {}".format(genome, bam_list, output)

	#Read coverage over gene body. This module is used to check if reads coverage is uniform and if there is any 5'/3' bias. 
	#This module scales all transcripts to 100 nt and calculates the number of reads covering each nucleotide position. 
	#Finally, it generates plots illustrating the coverage profile along the gene body.

def ConfigSectionMap(section):
	dict1 = {}
	options = Config.options(section)
	for option in options:
		try:
			dict1[option] = Config.get(section, option)
			if dict1[option] == -1:
				DebugPrint("skip: %s" % option)
		except:
			print("exception on %s!" % option)
			dict1[option] = None
	return dict1

def total_qc_fc(args):
	return total_qc_of_bam(*args)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='QC for bam files based on RSeqQC Package\n')
	parser.add_argument('-c','--CONFIG', help='Config file containing [[Output_directory]\nKey is bam file and value is Output_directory', required=False)
	parser.add_argument('-a', help='Perform Genebody coverage plot of bam files')
	parser.add_argument('-paired', help='Are samples Paired end', action="store_true")
	parser.add_argument('-gen', help='Genome, option are mm10 or hg19')
	parser.add_argument('-t', help='Threads, default=4. Keep it low', default=4)
	args = vars(parser.parse_args())
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["CONFIG"])
	#Read design matrix and create list of conditions and directories
	conditions = ConfigSectionMap("Output_directory")
	path = os.path.dirname(rnaseq_toolkit.__file__)
	reference = path + "/data/{}_Enembl_ens.bed".format(args["gen"])
	for key in conditions:
		subprocess.call(["mkdir", conditions[key]])
	pool = Pool(int(args["t"]))
	pool.map(total_qc_fc, itertools.izip(list(conditions.keys()), itertools.repeat(reference), itertools.repeat(args["paired"]), itertools.repeat(conditions))) 
	pool.close()
	pool.join()