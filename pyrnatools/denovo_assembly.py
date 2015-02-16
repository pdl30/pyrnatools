#!/usr/bin/python

########################################################################
# 30 July 2014
# Patrick Lombard, Centre for Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import os, re, sys
import argparse
import ConfigParser
import subprocess

def run_cufflinks(idict, threads, libtype, gtf=None):
	for key in idict:
		bam_name = os.path.basename(key)
		output = re.sub(".bam$", "_clinks", bam_name)
		#outdir = key+"_clinks"
		if gtf:	
			command = "cufflinks -p {} -g {} --library-type {} -o {} {}".format(threads, gtf, libtype, output, key)
		else:
			command = "cufflinks -p {} --library-type {} -o {} {}".format(threads, output, key)
		subprocess.call(command.split())

def run_cuffmerge(idict, threads, fasta, gtf):
	tmp = open("tmp_cuff_input.txt", "w")
	for key in idict:
		bam_name = os.path.basename(key)
		output = re.sub(".bam$", "_clinks/transcripts.gtf", bam_name)
		tmp.write("{}\n".format(output)),
	tmp.close()
	command = "cuffmerge -p {} -g {} -o cuffmerge/ -s {} tmp_cuff_input.txt".format(threads, gtf, fasta)
	subprocess.call(command.split())

def run_cuffdiff(idict, comp1, comp2, threads, mergedgtf, outdir):
	list1 = 1


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

def main():
	parser = argparse.ArgumentParser(description='Assembly for RNA-seq experiments.\n')
	parser.add_argument('-c','--config', help='Config file containing Conditions which are bam/sam files', required=False)
	parser.add_argument('-g','--gtf', help='GTF file,  If not provided, will do denovo alignment', required=False)
	parser.add_argument('-l','--library', help='Cufflinks option for library type, program default is unstranded', default='fr-unstranded', required=False)
	parser.add_argument('-f','--fasta', help='Reference Fasta', required=False)
	parser.add_argument('-t','--threads', help='Number of threads, default=4', default=4, required=False)
	parser.add_argument('-o','--outdir', help='Output directory for merged transcripts', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	conditions = ConfigSectionMap("Conditions", Config)
	#Run cufflinks
	#run_cufflinks(conditions, args["threads"], args["library"], args["gtf"])
	run_cuffmerge(conditions, args["threads"], args["fasta"], args["gtf"])

	#comparisons = ConfigSectionMap("Comparisons", Config)
	#for comp in comparisons:
#		c = comparisons[comp].split(",")
		#comps = [x.strip(' ') for x in c]

main()