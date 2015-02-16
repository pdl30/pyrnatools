#!/usr/bin/python

########################################################################
# 12 Jan 2015
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import sys, re, os
import argparse
import ConfigParser

def ConfigSectionMap(section, Config):
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

def miso_index_gtf(gtf):
	name = re.sub(".gtf$", "", gtf)
	output_dir = "indexed_{}".format(name)
	#GTF needs to be in GFF format
	#Use gtf2gff3.pl and gtf2gff3.cfg in pyrnatools/tools directory for conversion
	#Not neat for ensembl, still issues with source problems
	command = "index_gff --index {} {}".format(gtf, output_dir)
	subprocess.call(command.split())
	return output_dir

def miso_run(conditions, index, insert, sd, read_len):
	for sample in sorted(conditions):
		name = os.path.basename(sample)
		output_dir = re.sub(".bam$", "_miso", name)
	#	if os.path.isdir(output_dir):
	#		print "Miso directory already exists, will not overwrite!"
		if insert:
			command = "miso --run {} {} --output-dir {} --read-len {} --paired-end {} {}".format(index, sample, output_dir, read_len, insert, sd)
			#print command
			command2 = "summarize_miso --summarize-samples {} {}".format(output_dir, output_dir)
			print command2
		else:
			command = "miso --run {} {} --output-dir {} --read-len {}".format(index, sample, output_dir, read_len)
			#print command
			command2 = "summarize_miso --summarize-samples {} {}".format(output_dir, output_dir)
			print command2
			#subprocess.call(command.split())

def dexseq_prep(dexseq_dir, gtf, conditions, paired, orientation):
	#This prepares annotation and counts bam files
	count_program = os.path.join(dexseq_dir, "python_scripts/dexseq_count.py")
	anno_program = os.path.join(dexseq_dir, "python_scripts/deseq_prepare_annotation.py")
	#gff = re.sub(".gtf$", ".gff", gtf)
	#command1 = "python {} {} {}".format(anno_program, gtf, gff)

	for sample in sorted(conditions):
		bam_name = os.path.basename(sample)
		output = re.sub(".bam$", "_dexseq.count", bam_name)
		if paired:
			p = "yes"
		else:
			p = "no"
		command = "python {} -f bam -p {} -s {} {} {} {}".format(count_program, p, orientation, gtf, sample, output)
		subprocess.call(command.split())


def spliceR():
	rscript = "library(spliceR)"

def mats():


def main():
	parser = argparse.ArgumentParser(description='Differential expression for RNA-seq experiments. Runs DESEQ2 by default\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")

	miso_parser = subparsers.add_parser('miso', help="Runs MISO")
	miso_parser.add_argument('-c','--config', help='Config file containing parameters, please see documentation for usage!', required=True)
	miso_parser.add_argument('-g','--gtf', help='GTF file for indexing', required=False)
	miso_parser.add_argument('-i','--index', help='Already indexed GTF directory', required=False)
	miso_parser.add_argument('-n','--insert', help='Insert size', default=None, required=False)
	miso_parser.add_argument('-s','--sd', help='SD', required=False)
	miso_parser.add_argument('-r','--len', help='read length', required=False)

	dexseq_parser = subparsers.add_parser('dexseq', help="Runs DEXSEQ")
	dexseq_parser.add_argument('-c','--config', help='Config file containing bam files, please see documentation for usage!', required=True)
	dexseq_parser.add_argument('-g','--gtf', help='GTF file formatted by dexseq', required=True)
	dexseq_parser.add_argument('-p', action='store_true', help='Use if samples are paired end. Will find sd and insert size for bam files', required=False)
	dexseq_parser.add_argument('-o','--orientation', help='Options are yes, no or reverse. Test First!!!', required=True)

	msat_parser = subparsers.add_parser('mats', help="Runs MATS")
	dexseq_parser = subparsers.add_parser('spliceR', help="Runs spliceR")
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	if args["subparser_name"] == "miso":
		Config = ConfigParser.ConfigParser()
                Config.optionxform = str
                Config.read(args["config"])

        #Read design matrix and create list of conditions and directories
                conditions = ConfigSectionMap("Conditions", Config)

		if args["gtf"]:
			index = miso_index_gtf(args["gtf"])
		elif args["index"]:
			index = args["index"]
		else:
			raise Exception("Please provide either GTF or indexed GTF!")
		miso_run(conditions, args["index"], args["insert"], args["sd"], args["len"])
	elif args["subparser_name"] == "dexseq":
		Config = ConfigParser.ConfigParser()
		Config.optionxform = str
		Config.read(args["config"])

	#Read design matrix and create list of conditions and directories
		conditions = ConfigSectionMap("Conditions", Config)
		dexseq_dir = "/raid/home/patrick/R/x86_64-pc-linux-gnu-library/3.1/DEXSeq"
		dexseq_prep(dexseq_dir, args["gtf"], conditions, args["p"], args["orientation"])
	elif args["subparser_name"] == "mats":
		mats = "/home/patrick/Programs/MATS.3.0.8/RNASeq-MATS.py"
	elif args["subparser_name"] == "spliceR":
		

main()

#MATS OPTIONS
mats_option = """-s1 rep1_1[:rep1_2][,rep2_1[:rep2_2]]*	FASTQ file(s) for the sample_1. For the paired-end data, two files must be colon separated and replicates must be in a comma separated list (Only if using fastq)
-s2 rep1_1[:rep1_2][,rep2_1[:rep2_2]]*	FASTQ file(s) for the sample_2. For the paired-end data, two files must be colon separated and replicates must be in a comma separated list (Only if using fastq)
-b1 s1_rep1.bam[,s1_rep2.bam]	Mapping results for the sample_1 in bam format. Replicates must be in a comma separated list (Only if using bam)
-b2 s2.rep1.bam[,s2.rep2.bam]	Mapping results for the sample_2 in bam format. Replicates must be in a comma separated list (Only if using bam)
-t readType	Type of read used in the analysis. readType is either 'paired' or 'single'. 'paired' is for paired-end data and 'single' is for single-end data
-len <int>	The length of each read
-gtf gtfFile	An annotation of genes and transcripts in GTF format
-bi bowtieIndexBase	The basename of the bowtie indexes (ebwt files). The base name does not include the first period. For example, use hg19 for hg19.1.ebwt. (Only if using fastq)
-o outDir	The output directory
Optional:
-a <int>	The "anchor length" used in TopHat. At least “anchor length” NT must be mapped to each end of a given junction. The default is 8
-r1 <float>[,<float>]*	The insert sizes of sample_1 data. This applies only for the paired-end data. Replicates are separated by comma. The default is 15 for each replicate
-r2 <float>[,<float>]*	The insert sizes of sample_2 data. This applies only for the paired-end data. Replicates are separated by comma. The default is 15 for each replicate
-sd1 <float>[,<float>]*	The standard deviations for the r1. Replicates are separated by comma. The default is 70 for each replicate
-sd2 <float>[,<float>]*	The standard deviations for the r2. Replicates are separated by comma. The default is 70 for each replicate
-c <float>	The cutoff splicing difference. The cutoff used in the null hypothesis test for differential splicing. The default is 0.0001 for 0.01% difference. Valid: 0 ≤ cutoff < 1
-analysis analysisType	Type of analysis to perform. analysisType is either 'P' or 'U'. 'P' is for paired analysis and 'U' is for unpaired analysis. The default is 'U'
-expressionChange <float>	Filters out AS events with whose gene expression levels differs more than the given cutoff fold change between the two samples. Valid: fold change > 1.0. The default is 10000.0
-keepTemp	Enables MATS to keep its temporary files generated during the run. Temporary files are generally for debugging. The default is to delete all temporary files."""