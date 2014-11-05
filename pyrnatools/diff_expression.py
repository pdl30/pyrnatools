#!/usr/bin/python

########################################################################
# 19 May 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import sys, re, os
import ConfigParser
import itertools
import argparse
from collections import defaultdict
from pyrnatools.tools import gfold, deseq2
from multiprocessing import Pool
import pysam

def create_design_for_R(idict):
	output = open("tmp_design.txt", "w")
	output.write("sampleName\tfileName\tcondition\n"),
	for key in sorted(idict.keys()):
		name = re.sub(".bam$", "", key)
		count = re.sub(".bam$", ".count", key)
		output.write("{}\t{}\t{}\n".format(name, count, idict[key]))
	output.close()

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

def run_rcode(rscript, name):
	rcode = open(name, "w")
	rcode.write(rscript)
	rcode.close()
	try:
		subprocess.call(['Rscript', name])
	except:
		error("Error in running {}\n".format(name))
		error("Error: %s\n" % str(sys.exc_info()[1]))
		error( "[Exception type: %s, raised in %s:%d]\n" % ( sys.exc_info()[1].__class__.__name__, 
		os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
		traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
		sys.exit(1)

def main():
	parser = argparse.ArgumentParser(description='Differential expression for RNA-seq experiments. Runs DESEQ2 by default\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	deseq2_parser = subparsers.add_parser('deseq2', help="Runs DESEQ2")
	gfold_parser = subparsers.add_parser('gfold', help="Runs GFOLD")
	count_parser = subparsers.add_parser('count', help="Prints Counts from DESEQ")
	count_parser.add_argument('-c','--config', help='Config file containing parameters, please see documentation for usage!', required=False)
	deseq2_parser.add_argument('-c','--config', help='Config file containing parameters, please see documentation for usage!', required=False)
	gfold_parser.add_argument('-c','--config', help='Config file containing parameters, please see documentation for usage!', required=False)
	deseq2_parser.add_argument('-i','--input', help='Combined counts file from HTSeq or pyrna_count.py',required=True)
	deseq2_parser.add_argument('-p','--padj', help='Option for DESEQ2', default=0.05, required=False)
	gfold_parser.add_argument('-a','--alt', help='Use HTseq counts faked files. Assumes normalisation is already in place. Requires a standard file extension.', required=False)
	count_parser.add_argument('-i','--input', help='Combined counts file',required=True)
	args = vars(parser.parse_args())

	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])

	#Read design matrix and create list of conditions and directories
	conditions = ConfigSectionMap("Conditions", Config)
	

	if args["subparser_name"] == "count":
		create_design_for_R(conditions)
		rscript = deseq2.print_norm_counts(args["input"])
		run_rcode(rscript, "get_counts.R")
	if args["subparser_name"] == "gfold":
		comparisons = ConfigSectionMap("Comparisons", Config)
		if args["alt"]:
			for comp in comparisons:
				c = comparisons[comp].split(",")
				comps = [x.strip(' ') for x in c]
				gfold.run_gfold_diff(args["input"], conditions, comps[0], comps[1], alt=args["alt"])
		else:
			for comp in comparisons:
				c = comparisons[comp].split(",")
				comps = [x.strip(' ') for x in c]
				if args["n"]:
					gfold.gfold_housekeeper(comps[0], comps[1])
					gfold.run_gfold_diff(conditions, comps[0], comps[1],True)
				else:
					gfold.run_gfold_diff(conditions, comps[0], comps[1])

	if args["subparser_name"] == "deseq2":
		create_design_for_R(conditions)
		comparisons = ConfigSectionMap("Comparisons", Config)
		for comp in comparisons:
			c = comparisons[comp].split(",")
			comps = [x.strip(' ') for x in c]
			rscript = deseq2.write_deseq(conditions, comps[0], comps[1], args["padj"]) ##Needs changing!!!
			run_rcode(rscript, "deseq2_rcode.R")
