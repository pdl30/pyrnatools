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

def annotate_sam(bam_file, gtf_file, stranded, outfile=None):
	print "==> Counting sam file...\n"
	if outfile:
		htout = open(outfile, "w")
	else:
		count_file = re.sub(".bam$", ".count", bam_file)
		htout = open(count_file,"w")
	command = ["htseq-count", "--mode=union", "--stranded={}".format(stranded), "--quiet", "-f", "bam", bam_file,  gtf_file]
	subprocess.call(command, stdout=htout)

def join_counts(idict, outfile):
	data = defaultdict(list)
	output = open(outfile, "w")
	output.write("ID"),
	for bam in sorted(idict):
		name = re.sub(".bam$", "", bam) 
		output.write("\t{}".format(name)),
		with open(name+".count") as f:
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				if word[0].startswith("__"):
					pass
				else:
					data[word[0]].append(word[1])
	output.write("\n"),
	for key2 in sorted(data):
		data2 = data[key2]
		output.write(key2+"\t" + "\t".join(data2) + "\n"),
	output.close()
	
def run_gfold_count(args):
	return gfold.run_gfold_c(*args)
def anno_function(args):
	return annotate_sam(*args)

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
	parser = argparse.ArgumentParser(description='Counts features from BAM files\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	htseq_parser = subparsers.add_parser('htseq', help="Runs HTseq count")
	gfold_parser = subparsers.add_parser('gfold', help="Runs GFOLD Count")
	htseq_parser.add_argument('-c','--config', help='Config file containing [Conditions], please see documentation for usage!\nPlease use a unique name for every input bam file!', required=False)
	gfold_parser.add_argument('-c','--config', help='Config file containing [Conditions], please see documentation for usage!\nPlease use a unique name for every input bam file!', required=False)
	htseq_parser.add_argument('-i','--input', help='Input bam file', required=False)
	gfold_parser.add_argument('-i','--input', help='Input bam file', required=False)
	htseq_parser.add_argument('-g','--gtf', help='GTF file', required=False)
	gfold_parser.add_argument('-g','--gtf', help='GTF file', required=False)
	gfold_parser.add_argument('-n',action='store_true', help='Gapdh Normlisation', required=False)
	htseq_parser.add_argument('-s','--stranded', help='Option for HTSeq', default="yes", required=False)
	htseq_parser.add_argument('-t','--threads', help='Number of threads, default=8', default=8, required=False)
	gfold_parser.add_argument('-t','--threads', help='Number of threads, default=8', default=8, required=False)
	htseq_parser.add_argument('-o','--outfile', help='Output counts file. If using config, will be matrix, else will be single output', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	if args["subparser_name"] == "gfold":
		if args["config"]:
			Config = ConfigParser.ConfigParser()
			Config.optionxform = str
			Config.read(args["config"])

			#Read design matrix and create list of conditions and directories
			conditions = ConfigSectionMap("Conditions", Config)
			pool = Pool(int(args["threads"]))
			pool.map(run_gfold_count, itertools.izip(list(conditions.keys()), itertools.repeat(args["gtf"]))) ##Running annotation in parallel
			pool.close()
			pool.join()
		elif args["input"]:
			gfold.run_gfold_c(args["input"], args["gtf"])

	elif args["subparser_name"] == "htseq":
		if args["config"]:
			Config = ConfigParser.ConfigParser()
			Config.optionxform = str
			Config.read(args["config"])

			#Read design matrix and create list of conditions and directories
			conditions = ConfigSectionMap("Conditions", Config)

			pool = Pool(int(args["threads"]))
			pool.map(anno_function, itertools.izip(list(conditions.keys()), itertools.repeat(args["gtf"]), itertools.repeat(args["stranded"]))) ##Running annotation in parallel
			pool.close()
			pool.join()	
			join_counts(conditions, args["outfile"])
		elif args["input"]:
			annotate_sam(args["input"], args["gtf"], args["stranded"], args["outfile"])