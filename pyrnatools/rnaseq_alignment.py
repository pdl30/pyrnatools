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
from pyrnatools import calc_insert
import math

def run_fastqc(fq1):
	print "Running FastQC on {}\n".format(fq1)
	command = ["fastqc", "-q", "{}".format(fq1)] #outdir must exist!
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

def cut_adapters(paired, adapters, fq1, outdir, rev_adapters=None, fq2=None):
	f = open('{}/trim_report.txt'.format(outdir), 'w')
	adapt1 = ""
	for i in adapters:
		adapters = "-a {} ".format(i)
		adapt1 = adapters+adapt1
	
	if paired:
		command1 = "cutadapt -q 20 {0} --minimum-length=20 --paired-output {1}/tmp_2.fastq -o {1}/tmp_1.fastq {2} {3}".format(adapt1, outdir, fq1, fq2)
		p = subprocess.Popen(command1.split(), stdout=f)
		p.communicate()
		if rev_adapters:
			adapt2 = ""
			for i in rev_adapters:
				adapters = "-a {} ".format(i)
				adapt2 = adapters+adapt2			
			command2 = "cutadapt -q 20 {0} --minimum-length=20 --paired-output {1}/trimmed_1.fastq -o {1}/trimmed_2.fastq {1}/tmp_2.fastq {1}/tmp_1.fastq".format(adapt2, outdir)
			p = subprocess.Popen(command2.split(), stdout=f)
			p.communicate()
		else:
			command2 = "cutadapt -q 20 --minimum-length=20 --paired-output {0}/trimmed_1.fastq -o {0}/trimmed_2.fastq {0}/tmp_2.fastq {0}/tmp_1.fastq".format(outdir)
			p = subprocess.Popen(command2.split(), stdout=f)
			p.communicate()
		cleanup = ["rm", "{0}/tmp_2.fastq".format(outdir), "{0}/tmp_1.fastq".format(outdir)]
		subprocess.call(cleanup)
	else:
		command1 = "cutadapt -q 20 {0} --minimum-length=20 -o {1}/trimmed.fastq {2}".format(adapt1, outdir, fq1)
		p = subprocess.Popen(command1.split(), stdout=f)
		p.communicate()

def paired_tophat(fastq1, fastq2, index, gtf, outdir, insert, sd, threads):
	report = outdir+'/'+'tophat_report.txt'
	report1_o = open(report, "wb")
	uniq = "tophat2 -r {5} --mate-std-dev {6} -p {7} --GTF {0} -o {1} {2} {3} {4} --no-coverage-search".format(gtf, outdir, index, fastq1, fastq2, insert, sd, threads)
	p = subprocess.Popen(uniq.split(), stderr=report1_o)
	p.communicate()

def single_tophat(fastq, index, gtf, outdir, threads):
	report = outdir+'/'+'tophat_report.txt'
	report1_o = open(report, "wb")
	uniq = "tophat2  --no-coverage-search --GTF {0} -p {4} -o {1} {2} {3}".format(gtf, outdir, index, fastq, threads)
	p = subprocess.Popen(uniq.split(), stderr=report1_o)
	p.communicate()

def single_star(fastq, index, threads):
	command1 = "STAR --genomeDir {} --readFilesIn {} --runThreadN {} --outFilterMismatchNmax 2".format(index, fastq, threads)
	star_unique("Aligned.out.sam", "uniquely_mappable.sam")
	command3 = "samtools view -bS uniquely_mappable.sam > uniquely_mappable.bam"
	command2 = "samtools sort uniquely_mappable.bam uniquely_mappable_sorted"
	subprocess.call(command1, shell=True)
	subprocess.call(command2, shell=True)
	subprocess.call(command3, shell=True)
	os.remove("uniquely_mappable.sam")
	os.remove("uniquely_mappable.bam")

def paired_star(fastq1, fastq2, index, threads):
	command1 = "STAR --genomeDir {} --readFilesIn {} {} --runThreadN {} --outFilterIntronMotifs None --outFilterMismatchNmax 2".format(index, fastq1, fastq2, threads)
	star_unique("Aligned.out.sam", "uniquely_mappable.sam")
	command3 = "samtools view -bS uniquely_mappable.sam > uniquely_mappable.bam"
	command2 = "samtools sort uniquely_mappable.bam uniquely_mappable_sorted"
	subprocess.call(command1, shell=True)
	subprocess.call(command2, shell=True)
	subprocess.call(command3, shell=True)
	os.remove("uniquely_mappable.sam")
	os.remove("uniquely_mappable.bam")

def star_unique(samfile, outfile):
	output=  open(outfile, "w")
	with open(samfile) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if line.startswith("@"):
				output.write("{}\n".format(line)),
				continue
			if len(word) > 10:
				if word[4] == "255":
					if word[11] == "NH:i:1":
						output.write("{}\n".format(line)),
						
def main():
	parser = argparse.ArgumentParser(description='RNA-seq alignment\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	tophat_parser = subparsers.add_parser('tophat', help="Tophat alignment")
	star_parser = subparsers.add_parser('star', help="STAR alignment")

	tophat_parser.add_argument('-f', '--fastq', help='Single end fastq', required=False)
	tophat_parser.add_argument('-p', '--pair', help='Paired end fastqs. Please put them in order!', required=False, nargs='+')
	tophat_parser.add_argument('-i', '--index', help='Genome Index', required=True)
	tophat_parser.add_argument('-g', '--gtf', help='GTF file', required=True)
	tophat_parser.add_argument('-a', '--insert', help='Insert size for paried end, default=50', default=50, required=False)
	tophat_parser.add_argument('-b', '--sd', help='Insert size SD for paried end, default=20', default=20, required=False)
	tophat_parser.add_argument('-t', '--threads', help='Number of threads', default=1, required=False)
	tophat_parser.add_argument('-c', help='Will find sd and insert automatically', action='store_true', required=False)
	tophat_parser.add_argument('-o', '--out', help='Name of results directory', required=True)

	star_parser.add_argument('-f', '--fastq', help='Single end fastq', required=False)
	star_parser.add_argument('-p', '--pair', help='Paired end fastqs. Please put them in order!', required=False, nargs='+')
	star_parser.add_argument('-i', '--index', help='Genome Index Folder', required=True)
	star_parser.add_argument('-t', '--threads', help='Number of threads', default=1, required=False)
	star_parser.add_argument('-o', '--out', help='Name of results directory', required=True)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	path = os.getcwd()
	
	if os.path.isdir(args["out"]):
		print("Results directory already exists!")
	else:
		subprocess.call(["mkdir", args["out"]])

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
			cut_adapters(True, fwd_adapt, fq1, args["out"], fq2=fq2, rev_adapters=rev_adapt)
			fq1 = args["out"]+"/trimmed_1.fastq" 
			fq2 = args["out"]+"/trimmed_2.fastq"
			fq1 = os.path.abspath(fq1)
			fq2 = os.path.abspath(fq2)
		else:
			fq1 = os.path.abspath(fq1)
			fq2 = os.path.abspath(fq2)

	else:
		fq1 = args["fastq"]
		print "==> Running FastQC...\n"
		run_fastqc(fq1)
		adapt = find_adapters(fq1)
		if adapt:
			print "==> Removing adapters...\n"
			cut_adapters(False, adapt, fq1, args["out"])
			fq1 = args["out"]+"/trimmed.fastq" 
			fq1 = os.path.abspath(fq1)
		else:
		#	subprocess.call("mv {} {}".format(fq1, args["out"]), shell=True) #This shouldn't be done. 
			fq1 = os.path.abspath(fq1)

	if args["subparser_name"] == "tophat":
		print "==> Running Tophat...\n"
		if args["pair"]:
			if args["c"]:
				fq = False
				if args["pair"][0].endswith(".fq"):
					fq = True
				name = calc_insert.run_bowtie(fq1, fq2, args["index"], fq)
				calc_insert.sort_index(name)
				insert = calc_insert.pysam_insert_size(name)
				paired_tophat(fq1, fq2, args["index"], args["gtf"], args["out"], int(insert[0]), int(insert[1]), args["threads"])
			else:
				paired_tophat(fq1, fq2, args["index"], args["gtf"], args["out"], args["insert"], args["sd"], args["threads"])
		elif args["fastq"]:
			single_tophat(fq1, args["index"], args["gtf"], args["out"], args["threads"])
	elif args["subparser_name"] == "star":
		print "==> Running STAR...\n"
		os.chdir(args["out"])
		if args["pair"]:
			paired_star(fq1, fq2, args["index"], args["threads"])
		elif args["fastq"]:
			single_star(fq1, args["index"], args["threads"])