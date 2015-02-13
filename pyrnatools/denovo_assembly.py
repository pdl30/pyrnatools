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
		outdir = key+"_clinks"
		if gtf:
			command = "cufflinks -p {} -g {} --library-type {} -o {} {}".format(threads, libtype, outdir, key)
		else:
			command = "cufflinks -p {} --library-type {} -o {} {}".format(threads, outdir, key)
	subprocess.call(command.split())

def run_cuffmerge(idict, threads, fasta, outdir=None):
	tmp = open("tmp_cuff_input.txt", "w")
	for key in idict:
		outdir = key+"_clinks"
		tmp.write("{}\n".format(outdir)),
	tmp.close()
	if outdir:
		command = "cuffmerge -p {} -o {} -s {} tmp_cuff_input.txt".format(threads, outdir, fasta)
	else:
		command = "cuffmerge -p {} -s {} tmp_cuff_input.txt".format(threads, fasta)
	subprocess.call(command.split())

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
	conditions = ConfigSectionMap("Conditions")

	#Run cufflinks
	if args["library"]:
		run_cufflinks(conditions, args["threads"], args["library"], args["gtf"])
#	else:
#		run_cufflinks(conditions, args["threads"])

#	run_cuffmerge(conditions, args["threads"], args["fasta"], args["outdir"])
	#else:
	#	run_cuffmerge(conditions, args["threads"], args["fasta"])