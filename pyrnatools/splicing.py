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
import ConfigParser

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

def miso_index_gtf(gtf):
	name = re.sub(".gtf$", "", gtf)
	output_dir = "indexed_{}".format(name)
	#GTF needs to be in GFF format
	#Use gtf2gff3.pl and gtf2gff3.cfg in pyrnatools/tools directory for conversion
	#Not neat for ensembl, still issues with source problems
	command = "index_gff --index {} {}".format(gtf, output_dir)
	subprocess.call(command.split())
	return output_dir

def miso_run(conditions, index, insert, sd, read_len):
	for sample in sorted(conditions):
		name = os.path.basename(sample)
		output_dir = re.sub(".bam$", "_miso", name)
	#	if os.path.isdir(output_dir):
	#		print "Miso directory already exists, will not overwrite!"
		if insert:
			command = "miso --run {} {} --output-dir {} --read-len {} --paired-end {} {}".format(index, sample, output_dir, read_len, insert, sd)
			#print command
			command2 = "summarize_miso --summarize-samples {} {}".format(output_dir, output_dir)
			print command2
		else:
			command = "miso --run {} {} --output-dir {} --read-len {}".format(index, sample, output_dir, read_len)
			#print command
			command2 = "summarize_miso --summarize-samples {} {}".format(output_dir, output_dir)
			print command2
			#subprocess.call(command.split())

def dexseq_prep(dexseq_dir, gtf, conditions, paired, orientation):
	#This prepares annotation and counts bam files
	count_program = os.path.join(dexseq_dir, "python_scripts/dexseq_count.py")
	anno_program = os.path.join(dexseq_dir, "python_scripts/deseq_prepare_annotation.py")
	#gff = re.sub(".gtf$", ".gff", gtf)
	#command1 = "python {} {} {}".format(anno_program, gtf, gff)

	for sample in sorted(conditions):
		bam_name = os.path.basename(sample)
		output = re.sub(".bam$", "_dexseq.count", bam_name)
		if paired:
			p = "yes"
		else:
			p = "no"
		command = "python {} -f bam -p {} -s {} {} {} {}".format(count_program, p, orientation, gtf, sample, output)
		subprocess.call(command.split())


def spliceR():
	rscript = "library(spliceR)"

def mats(conditions, comp1, comp2, gtf, insert, sd, rlen, mats_program):
	rev_cond = reverse_dict(conditions)
	list1 = ",".join(list(rev_cond[comp1]))
	list2 = ",".join(list(rev_cond[comp2]))
	outdir = "{}_vs_{}".format(comp1, comp2)
	a = ",{}".format(insert) * (len(list(rev_cond[comp1])) - 1)
	b = ",{}".format(insert) * (len(list(rev_cond[comp2])) - 1)
	c = ",{}".format(sd) * (len(list(rev_cond[comp1])) - 1)
	d = ",{}".format(sd) * (len(list(rev_cond[comp2])) - 1)
	a = "{}".format(insert) + a
	b = "{}".format(insert) + b
	c = "{}".format(sd) + c
	d = "{}".format(sd) + d
	if insert:
		command = "python {} -b1 {} -b2 {} -gtf {} -o {} -t paired -len {} -r1 {} -sd1 {} -r2 {} -sd2 {} -analysis P".format(mats_program, list1, list2, 
			gtf, outdir, rlen, a, c, b, d)
		subprocess.call(command.split())


def reverse_dict(idict):
	inv_map = {}
	for k, v in idict.iteritems():
		inv_map[v] = inv_map.get(v, [])
		inv_map[v].append(k)
	return inv_map

def main():
	parser = argparse.ArgumentParser(description='Differential expression for RNA-seq experiments. Runs DESEQ2 by default\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")

	miso_parser = subparsers.add_parser('miso', help="Runs MISO")
	miso_parser.add_argument('-c','--config', help='Config file containing parameters, please see documentation for usage!', required=True)
	miso_parser.add_argument('-g','--gtf', help='GTF file for indexing', required=False)
	miso_parser.add_argument('-i','--index', help='Already indexed GTF directory', required=False)
	miso_parser.add_argument('-n','--insert', help='Insert size', default=None, required=False)
	miso_parser.add_argument('-s','--sd', help='SD', required=False)
	miso_parser.add_argument('-r','--len', help='read length', required=False)

	dexseq_parser = subparsers.add_parser('dexseq', help="Runs DEXSEQ")
	dexseq_parser.add_argument('-c','--config', help='Config file containing bam files, please see documentation for usage!', required=True)
	dexseq_parser.add_argument('-g','--gtf', help='GTF file formatted by dexseq', required=True)
	dexseq_parser.add_argument('-p', action='store_true', help='Use if samples are paired end. Will find sd and insert size for bam files', required=False)
	dexseq_parser.add_argument('-o','--orientation', help='Options are yes, no or reverse. Test First!!!', required=True)

	mats_parser = subparsers.add_parser('mats', help="Runs MATS")
	mats_parser.add_argument('-c','--config', help='Config file containing bam files, please see documentation for usage!', required=True)
	mats_parser.add_argument('-g','--gtf', help='GTF file formatted by dexseq', required=True)
	mats_parser.add_argument('-n','--insert', help='Insert size, use for paired end data. In order to make things simple, will assume all bams are the same', default=None, required=False)
	mats_parser.add_argument('-s','--sd', help='SD', required=False)
	mats_parser.add_argument('-r','--len', help='read length', required=False)
	
	splice_parser = subparsers.add_parser('spliceR', help="Runs spliceR")
	splice_parser.add_argument('-c','--config', help='Config file containing bam files, please see documentation for usage!', required=True)
	splice_parser.add_argument('-g','--gtf', help='GTF file formatted by dexseq', required=True)
	splice_parser.add_argument('-o','--output', help='Output directory', required=True)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	conditions = ConfigSectionMap("Conditions", Config)
	
	if args["subparser_name"] == "miso":	
		if args["gtf"]:
			index = miso_index_gtf(args["gtf"])
		elif args["index"]:
			index = args["index"]
		else:
			raise Exception("Please provide either GTF or indexed GTF!")
		miso_run(conditions, args["index"], args["insert"], args["sd"], args["len"])

	elif args["subparser_name"] == "dexseq":
		conditions = ConfigSectionMap("Conditions", Config)
		dexseq_dir = "/raid/home/patrick/R/x86_64-pc-linux-gnu-library/3.1/DEXSeq"
		dexseq_prep(dexseq_dir, args["gtf"], conditions, args["p"], args["orientation"])

	elif args["subparser_name"] == "mats":
		mats_program = "/home/patrick/Programs/MATS.3.0.8/RNASeq-MATS.py"
		comparisons = ConfigSectionMap("Comparisons", Config)
		for comp in comparisons:
			c = comparisons[comp].split(",")
			comps = [x.strip(' ') for x in c]
			mats(conditions, comps[0], comps[1], args["gtf"], args["insert"], args["sd"], args["len"], mats_program)

	elif args["subparser_name"] == "spliceR":
		rcode = spliceR()

main()

