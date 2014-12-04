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

def run_gfold_c(bam_file, sample_dict, gtf_file):
	count_file = re.sub(".bam$", "_gfold.count", bam_file)
	command = "samtools view {0} | gfold count -ann {1} -tag stdin -o {2}".format(bam_file, gtf_file, count)
	subprocess.call(command, shell=True)

def gfold_housekeeper(cond1, cond2):
	#Not used anymore!
	#Here it is Gapdh but can be anything!
	count1 = cond1+"_gfold.count"
	count2 = cond2 +"_gfold.count"
	with open(count1) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if word[1] == "Gapdh":
				gapdh1 = float(word[2])
				print gapdh1
	with open(count2) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if word[1] == "Gapdh":
				gapdh2 = float(word[2])
	norm1 = 100000/float(gapdh1)
	norm2 = 100000/float(gapdh2)
	output1 = open(cond1+"_gapdh.count", "w")
	output2 = open(cond2+"_gapdh.count", "w")
	with open(count1) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			norm = float(word[2])*norm1
			output1.write("{}\t{}\t{}\t{}\t{}\n".format(word[0], word[1], norm, word[3], word[4])),
	with open(count2) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			norm = float(word[2])*norm2
			output2.write("{}\t{}\t{}\t{}\t{}\n".format(word[0], word[1], norm, word[3], word[4])),
	output1.close()
	output2.close()

def run_gfold_diff(sample_dict, cond1, cond2, norm=False, alt=None):
	if norm:
		command2 = "gfold diff -s1 {0} -s2 {1} -suf _gapdh.count -o {0}_vs_{1}.diff -norm NO".format(cond1, cond2)
	elif alt:
		command2 = "gfold diff -s1 {0} -s2 {1} -suf {2} -o {0}_vs_{1}.diff -norm NO".format(cond1, cond2, alt)
	else:
		command2 = "gfold diff -s1 {0} -s2 {1} -suf _gfold.count -o {0}_vs_{1}.diff".format(cond1, cond2)
	print command2
	subprocess.call(command2.split())