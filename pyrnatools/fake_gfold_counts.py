#!/usr/bin/python

########################################################################
# 09 Oct 2014
# Patrick Lombard, Centre for Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import os, re, sys
import argparse

def htseq_to_gfold(ilist, gfold_file):
	gfold_data = {}
	with open(gfold_file) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			gfold_data[word[0]] = (word[1], word[3])
	for ifile in ilist:
		name = ifile.strip("_htseq.tsv")
		output = open(name+"_htgold.tsv", "w")
		with open(ifile) as f:
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				m = re.search("__", word[0])
				if m:
					pass
				else:
					output.write("{}\t{}\t{}\t{}\t{}\n".format(word[0], gfold_data[word[0]][0], word[1], gfold_data[word[0]][1], 10)),
		output.close()

def main():
	parser = argparse.ArgumentParser(description='Takes HTseq output and creates a fake gfold file.\n ')
	parser.add_argument('-i', '--input', help='Input directory containing files with _htseq.tsv ending', required=True)
	parser.add_argument('-e', '--ex', help='Example gfold count file processed using same GTF', required=True)
	args = vars(parser.parse_args())
	ifiles = [f for f in os.listdir(args["input"]) if f.endswith("_htseq.tsv")]
	htseq_to_gfold(ifiles, args["ex"])