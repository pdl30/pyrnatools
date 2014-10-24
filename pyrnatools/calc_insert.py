#!/usr/bin/python

########################################################################
# 15 May 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import re, os, sys
import argparse

picard = "/home/patrick/Programs/picard-tools-1.70"

def head(ifile, ofile):
	output = open(ofile, "w")
	with open(ifile) as myfile:
		head = [next(myfile) for x in xrange(100000)]
	size = 0
	read1 = head[1]
	bases = read1.rstrip()
	for h in head:
		output.write(h),
	output.close()
	return len(list(bases))
		
def calc_insert(fq1, fq2, index):
	f = open('/dev/null', 'w')
	name = re.sub("_1.fq", "", fq1)
	read_size = head("{}_1.fq".format(name), "{0}_trun1.fq".format(name))
	head("{}_2.fq".format(name), "{0}_trun2.fq".format(name))
	c3 = "bowtie --sam -I 0 -X 500 {0} -1 {1}_trun1.fq -2 {1}_trun2.fq > {1}_trun.sam".format(index, name)
	subprocess.call(c3, shell=True, stderr=f)
	os.remove("{0}_trun1.fq".format(name))
	os.remove("{0}_trun2.fq".format(name))
	c4 = "java -Xmx2g -jar {0}/SortSam.jar INPUT={1}_trun.sam OUTPUT={1}_trunsorted.sam SORT_ORDER=coordinate".format(picard, name)
	subprocess.call(c4, shell=True, stderr=f)
	os.remove("{0}_trun.sam".format(name))
	c5 = "java -Xmx2g -jar {0}/CollectInsertSizeMetrics.jar INPUT={1}_trunsorted.sam HISTOGRAM_FILE={1}_hist.pdf OUTPUT={1}_result.txt".format(picard, name)
	subprocess.call(c5, shell=True, stderr=f)
	os.remove("{0}_trunsorted.sam".format(name))

	lines = open("{}_result.txt".format(name)).readlines()
	for i, line in enumerate(lines):
		if i == 7:
			info= line.rstrip().split("\t")
			r_size = float(info[4]) - 2*int(read_size)
			print "mate inner distance = {}\nmate standard deviation = {}\n".format(r_size, info[5]),
	os.remove("{}_hist.pdf".format(name))
	os.remove("{}_result.txt".format(name))

def main():
	parser = argparse.ArgumentParser(description='Calculates insert size for paired end data\n ')
	parser.add_argument('-p', '--pair', help='Paired end fastqs. Please put them in order!', required=False, nargs='+')
	parser.add_argument('-i', '--index', help='Bowtie1 Index', required=True)
	args = vars(parser.parse_args())
	calc_insert(args["pair"][0], args["pair"][1], args["index"])