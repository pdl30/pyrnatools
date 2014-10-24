#!/usr/bin/python

########################################################################
# 28 July 2014
# Patrick Lombard, Centre for Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import os, re, sys
import argparse

def read_dir(idir):
	ifiles = [f for f in os.listdir(idir) if f.endswith("_analysis.tsv")]
	data = {}
	for ifile in ifiles:
		name = re.sub("_deseq2_analysis.tsv", "", ifile)
		data[name] = {}
		with open(ifile) as f:
			next(f)
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				data[name][word[0]] = (word[2], word[5], word[6])
	return data

def read_gfold_dir(idir):
	ifiles = [f for f in os.listdir(idir) if f.endswith(".diff")]
	data = {}
	for ifile in ifiles:
		name = re.sub(".diff", "", ifile)
		data[name] = {}
		with open(ifile) as f:
			for i in range(0,10):
				next(f)
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				data[name][word[0]] = word[2]
#	print data
	return data

def read_counts(count_file, idict):
	with open(count_file) as f:
		header = next(f)
		header  = header.rstrip()
		print header,
		for key in sorted(idict):
			print "\t{} LFC".format(key),
			print "\t{} Pvalue".format(key),
			print "\t{} AdjPval".format(key),
		print "\n",
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if word[0].startswith("_"):
				pass
			else:
				result = []
				for key in sorted(idict):
					res = idict[key].get(word[0], "NA")
					result.append(res[0])
					result.append(res[1])
					result.append(res[2])
				if len(result) > 0:
					print "{}".format(line),
					for k in result:
						print "\t{}".format(k),
					print "\n",

def read_gfold_counts(count_file, idict):
	with open(count_file) as f:
		header = next(f)
		header  = header.rstrip()
		print header,
		for key in sorted(idict):
			print "\t{} Gfold LFC".format(key),
		print "\n",
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if word[0].startswith("_"):
				pass
			else:
				result = []
				for key in sorted(idict):
					res = idict[key][word[0]]
					result.append(res)
				if len(result) > 0:
					print "{}".format(line),
					for k in result:
						print "\t{}".format(k),
					print "\n",

def main():
	parser = argparse.ArgumentParser(description='Creates RNA-seq table from DESEQ/GFOLD results\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	deseq2_parser = subparsers.add_parser('deseq2', help="Runs DESEQ2")
	gfold_parser = subparsers.add_parser('gfold', help="Runs GFOLD")

	deseq2_parser.add_argument('-m','--matrix', help='Counts Matrix', required=True)
	deseq2_parser.add_argument('-i','--idir', help='Directory containing DESEQ2 results', required=True)

	gfold_parser.add_argument('-m','--matrix', help='Counts Matrix', required=True)
	gfold_parser.add_argument('-i','--idir', help='Directory containing Gfold results', required=True)

	args = vars(parser.parse_args())
	if args["subparser_name"] == "deseq2":
		deseq_data = read_dir(args["idir"])
		read_counts(args["matrix"], deseq_data)
	elif args["subparser_name"] == "gfold":
		gfold_data = read_gfold_dir(args["idir"])
		read_gfold_counts(args["matrix"], gfold_data)