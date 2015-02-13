#!/usr/bin/python

########################################################################
# 12 Jan 2015
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import sys, re, os
import argparse

def miso_index_gtf(gtf):
	name = re.sub(".gtf$", "", gtf)
	output_dir = "indexed_{}".format(name)
	#GTF needs to be in GFF format
	#Use gtf2gff3.pl and gtf2gff3.cfg in pyrnatools/tools directory for conversion
	#Not neat for ensembl, still issues with source problems
	command = "index_gff --index {} {}".format(gtf, output_dir)
	subprocess.call(command.split())
	return output_dir

def miso_run(conditions, index):
	for sample in sorted(conditions):
		output_dir = re.sub(".bam$", "_miso", sample)
		if os.path.isdir(output_dir):
			print "Miso directory already exists, will not overwrite!"
		else:
			command = "miso --run {} {} --output-dir {} --read-len {}".format(index, sample, output_dir, read_len)
			subprocess.call(command.split())

def dexseq_prep(dexseq_dir, gtf, conditions, paired, orientation):
	#This prepares annotation and counts bam files
	count_program = os.path.join(dexseq_dir, "python_scripts/dexseq_count.py")
	anno_program = os.path.join(dexseq_dir, "python_scripts/deseq_prepare_annotation.py")
	gff = re.sub(".gtf$", ".gff", gtf)
	command1 = "python {} {} {}".format(anno_program, gtf, gff)

	for sample in sorted(conditions):
		output = re.sub(".bam$", "_dexseq.count", sample)
		if paired:
			p = "yes"
		else:
			p = "no"
		command = "python {} -p {} -s {} {} {} {}".format(count_program, p, orientation, gff, sample, output)
		subprocess.call(command.split())

def main():
	parser = argparse.ArgumentParser(description='Differential expression for RNA-seq experiments. Runs DESEQ2 by default\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")

	miso_parser = subparsers.add_parser('miso', help="Runs MISO")
	miso_parser.add_argument('-c','--config', help='Config file containing parameters, please see documentation for usage!', required=True)
	miso_parser.add_argument('-g','--gtf', help='GTF file for indexing', required=False)
	miso_parser.add_argument('-i','--index', help='Already indexed GTF directory', required=False)
	miso_parser.add_argument('-p', action='store_true', help='Use if samples are paired end. Will find sd and insert size for bam files', required=False)

	dexseq_parser = subparsers.add_parser('dexseq', help="Runs DEXSEQ")
	dexseq_parser.add_argument('-c','--config', help='Config file containing parameters, please see documentation for usage!', required=True)
	dexseq_parser.add_argument('-g','--gtf', help='GTF file', required=True)
	dexseq_parser.add_argument('-p', action='store_true', help='Use if samples are paired end. Will find sd and insert size for bam files', required=False)
	dexseq_parser.add_argument('-o','--orientation', help='Options are yes, no or reverse. Test First!!!', required=True)

	msat_parser = subparsers.add_parser('mmsat', help="Runs MMSAT")
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	if args["subparser_name"] == "miso":
		if args["gtf"]:
			index = miso_index_gtf(args["gtf"])
		elif args["index"]:
			index = args["index"]
		else:
			raise Exception("Please provide either GTF or indexed GTF!")
		F
	elif args["subparser_name"] == "dexseq":
		dexseq_dir = "/raid/home/patrick/R/x86_64-pc-linux-gnu-library/3.1/DEXSeq"
		dexseq_prep(dexseq_dir, args["gtf"], conditions, args["p"], args["orientation"])