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
import subprocess
import ConfigParser
from multiprocessing import Pool, Manager
import itertools

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

def downloader(gsm):
	old_path = os.getcwd()
#	download2 = "wget -c -nv -q 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='{0}'&targ=self&view=quick&form=text&token=_blank' -O tmp/{0}.soft".format(gsm)
	#subprocess.call(download2, shell=True)
	if not os.path.isdir(gsm):
		os.mkdir(gsm)
	f = open("tmp/{}.soft".format(gsm), "r")
	lines = f.readlines()
	sra = None
	gsm_sra = []
	os.chdir(gsm)
	new_path = os.getcwd()
	for line in lines:
		line = line.rstrip()
		if line.startswith("!Sample_supplementary_file"):
			sra_path = line.split("ByExp/sra/")
			if len(sra_path) > 1:
				sra = sra_path[1]
				download3 = "wget -r -c --no-verbose -N -nd ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/{}".format(sra)
				subprocess.call(download3, shell=True)
	sras = [f for f in os.listdir(new_path) if f.endswith(".sra")]
	for sra in sras:
		command0 = "fastq-dump --split-3 {}".format(sra)
		subprocess.call(command0, shell=True)
		os.remove(sra)
	fqs = [f for f in os.listdir(new_path) if f.endswith(".fastq")]
	command = "cat"
	for fq in fqs:
		command += " {}".format(fq)
	command += " > {}.fastq".format(gsm)
	print command
	subprocess.call(command, shell=True)
	os.chdir(old_path)

def download_function(args):
	return downloader(*args)

def main():
	parser = argparse.ArgumentParser(description='Various clustering for RNA-seq experiments using DESEQ2 counts\n')
	parser.add_argument('-c','--config', help='Config file containing parameters, please see documentation for usage!', required=False)
	parser.add_argument('-t','--threads', help='Number of threads', default=20, required=False)
	args = vars(parser.parse_args())
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	if not os.path.isdir("tmp"):
		os.mkdir("tmp/")
	conditions = ConfigSectionMap("Conditions", Config)
	pool = Pool(int(args["threads"]))
	
	pool.map(download_function, itertools.izip(list(conditions.keys())))
	pool.close()
	pool.join()	
	
main()