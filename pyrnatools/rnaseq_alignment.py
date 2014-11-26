#!/usr/bin/python

########################################################################
# 15 May 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import re
import argparse
import os, sys

def run_fastqc(fq1):
	command = ["fastqc", "{}".format(fq1)] #outdir must exist!
	subprocess.call(command)

def find_adapters(fq):
	adapters = []
	name = re.sub(".fastq", "", fq)
	report = name+"_fastqc/fastqc_data.txt"
	flist = open(report).readlines()
	parsing = False
	for line in flist:
		if line.startswith(">>Overrepresented sequences\tfail"):
			parsing = True
		elif line.startswith(">>END_MODULE"):
			parsing = False
		if parsing:
			if line.startswith(">>"):
				continue
			if line.startswith("#"):
				continue
			else:
				line = line.rstrip()
				word = line.split("\t")
				if word[3] != "No Hit":
					adapters.append(word[0])
	return adapters

def cut_adapters(adapters, fq1, outdir, rev_adapters=None, fq2=None):

	adapt1 = ""
	for i in adapters:
		adapters = "-a {} ".format(i)
		adapt1 = adapters+adapt1
	if rev_adapters:
		adapt2 = ""
		for i in rev_adapters:
			adapters = "-a {} ".format(i)
			adapt2 = adapters+adapt2

	if rev_adapters:
		command1 = "cutadapt -q 20 {0} --minimum-length=10 --paired-output {1}/tmp.2.fastq -o {1}/tmp.1.fastq {2} {3}".format(adapt1, outdir, fq1, fq2)
		p = subprocess.Popen(command1.split())
		p.communicate()
		command2 = "cutadapt -q 20 {0} --minimum-length=10 --paired-output {1}/trimmed_1.fastq -o {1}/trimmed_2.fastq {1}/tmp.2.fastq {1}/tmp.1.fastq".format(adapt2, outdir)
		p = subprocess.Popen(command2.split())
		p.communicate()
		cleanup = ["rm", "{0}/tmp_2.fastq".format(outdir), "{0}/tmp_1.fastq".format(outdir)]
		subprocess.call(cleanup)
	else:
		command1 = "cutadapt -q 20 {0} --minimum-length=10 -o {1}/trimmed.fastq {2}".format(adapt1, outdir, fq1)
		p = subprocess.Popen(command1.split())
		p.communicate()

def paired_tophat(fastq1, fastq2, index, gtf, outdir, insert, sd, threads):
	report = outdir+'/'+'tophat_report.txt'
	report1_o = open(report, "wb")
	uniq = "tophat2 -r {5} --mate-std-dev {6} -p {7} --GTF {0} -o {1} {2} {3} {4}".format(gtf, outdir, index, fastq1, fastq2, insert, sd, threads)
	p = subprocess.Popen(uniq.split(), stderr=report1_o)
	p.communicate()

def single_tophat(fastq, index, gtf, outdir, threads):
	report = outdir+'/'+'tophat_report.txt'
	report1_o = open(report, "wb")
	uniq = "tophat2 --GTF {0} -p {4} -o {1} {2} {3}".format(gtf, outdir, index, fastq, threads)
	p = subprocess.Popen(uniq.split(), stderr=report1_o)
	p.communicate()


def main():
	parser = argparse.ArgumentParser(description='ChIP-seq bowtie wrapper\n ')
	parser.add_argument('-f', '--fastq', help='Single end fastq', required=False)
	parser.add_argument('-p', '--pair', help='Paired end fastqs. Please put them in order!', required=False, nargs='+')
	parser.add_argument('-i', '--index', help='Genome Index', required=True)
	parser.add_argument('-g', '--gtf', help='GTF file', required=True)
	parser.add_argument('-a', '--insert', help='Insert size for paried end, default=50', default=50, required=False)
	parser.add_argument('-b', '--sd', help='Insert size SD for paried end, default=20', default=20, required=False)
	parser.add_argument('-t', '--threads', help='Number of threads', default=1, required=False)
	#parser.add_argument('-aligner',, help='Which aligner to use, options are Tophat/Star, default is Tophat [STAR not implemented yet]', required=False, default="Tophat")
	parser.add_argument('-o', '--out', help='Name of results directory', required=True)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	path = os.getcwd()
	
	if os.path.isdir(args["out"]):
		print("Results directory already exists!")
	else:
		subprocess.call(["mkdir", args["out"]])

	#if args["ALIGNER"] == "Tophat":
	if args["pair"]:
		fq1 = args["pair"][0]
		fq2 = args["pair"][1]
		print "==> Running FastQC...\n"
		run_fastqc(fq1)
		run_fastqc(fq2)
		fwd_adapt = find_adapters(fq1)
		rev_adapt = find_adapters(fq2)
		if fwd_adapt or rev_adapt:
			print "==> Removing adapters...\n"
			cut_adapters(fwd_adapt, fq1, args["out"], fq2=fq2, rev_adapters=rev_adapt)
			fq1 = args["out"]+"/trimmed_1.fastq" 
			fq2 = args["out"]+"/trimmed_2.fastq"
		print "==> Aligning fastq's...\n"
		paired_tophat(fq1, fq2, args["index"], args["gtf"], args["out"], args["insert"], args["sd"], args["threads"])
	elif args["fastq"]:
		fq1 = args["fastq"]
		print "==> Running FastQC...\n"
		run_fastqc(fq1)
		adapt = find_adapters(fq1)
		if adapt:
			print "==> Removing adapters...\n"
			cut_adapters(adapt, fq1, args["out"])
			fq1 = args["out"]+"/trimmed.fastq" 
		print "==> Aligning fastq's...\n"
		single_tophat(args["fastq"], args["index"], args["gtf"], args["out"], args["threads"])
