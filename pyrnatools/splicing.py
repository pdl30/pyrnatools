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
from multiprocessing import Pool
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

def miso_run(conditions, index, event, insert, sd, read_len):
	for sample in sorted(conditions):
		name = os.path.basename(sample)
		output_dir = re.sub(".bam$", "_miso", name)
		if insert:
			command = "miso --run {} {} --output-dir {}_{} --read-len {} --paired-end {} {}".format(index, sample, output_dir, event, read_len, insert, sd)
			command2 = "summarize_miso --summarize-samples {0}_{1} {0}_{1}".format(output_dir, event)
			subprocess.call(command.split())
			subprocess.call(command2.split())
		else:
			command = "miso --run {} {} --output-dir {}_{} --read-len {}".format(index, sample, output_dir,event, read_len)
			command2 = "summarize_miso --summarize-samples {0}_{1} {0}_{1}".format(output_dir, event)
			subprocess.call(command.split())
			subprocess.call(command2.split())

def dexseq_prep(sample, dexseq_dir, gtf, paired, orientation):
	#This prepares annotation and counts bam files. Not sure if working correctly!!!
	count_program = os.path.join(dexseq_dir, "python_scripts/dexseq_count.py")
	bam_name = os.path.basename(sample)
	output = re.sub(".bam$", "_dexseq.count", bam_name)
	dex_out = open("count_output.txt", "a")
	if paired:
		p = "yes"
	else:
		p = "no"
	command = "python {} -f bam -p {} -s {} {} {} {}".format(count_program, p, orientation, gtf, sample, output)
	#python dexseq_count.py -s reverse Homo_sapiens.GRCh37.74.ucsc_dexseq.gff  /accepted_hits.bam 75K.DIF.1_1.fq.mm10.count -f bam -p yes
	subprocess.call(command.split(), stdout=dex_out)

def create_design_for_R(idict):
	output = open("tmp_design.txt", "w")
	output.write("sampleName\tfileName\tcondition\n"),
	for key in sorted(idict.keys()):
		bam_name = os.path.basename(key)
		name = re.sub(".bam$", "", bam_name)
		count = re.sub(".bam$", "_dexseq.count", bam_name)
		output.write("{}\t{}\t{}\n".format(name, count, idict[key]))
	output.close()

def dexseq_run(conditions, comp1, comp2, gtf):
	rscript = "suppressPackageStartupMessages(library(DEXSeq))\n"
	rscript += "pdata <- read.table('tmp_design.txt', header=T)\n"
	rscript += "counts1 <- pdata[which(pdata[,3] == '{}'),]\n".format(comp1)
	rscript += "counts2 <- pdata[which(pdata[,3] == '{}'),]\n".format(comp2)
	rscript += "new_pdata <- rbind(counts1, counts2)\n"
	rscript += "dxd <- DEXSeqDataSetFromHTSeq(as.character(new_pdata[,2]), sampleData=new_pdata, design= ~ sample + exon + condition:exon, flattenedfile='{}')\n".format(gtf)
	rscript += "dxd = estimateSizeFactors( dxd )\n"
	rscript += "dxd = estimateDispersions( dxd )\n"
	rscript += "dxd = testForDEU( dxd )\n"
	rscript += "dxr1 = DEXSeqResults( dxd )\n"
	rscript += "DEXSeqHTML( dxr1, FDR=0.1, color=c('#FF000080', '#0000FF80') )\n"
	rscript += "write.table(dxr1, file='{}_vs_{}.tsv', sep='\\t', quote=F)\n".format(comp1, comp2)
	return rscript

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
		print command
		#subprocess.call(command.split())

def reverse_dict(idict):
	inv_map = {}
	for k, v in idict.iteritems():
		inv_map[v] = inv_map.get(v, [])
		inv_map[v].append(k)
	return inv_map

def run_rcode(rscript, name):
	rcode = open(name, "w")
	rcode.write(rscript)
	rcode.close()
	try:
		subprocess.call(['Rscript', name])
	except:
		error("Error in running {}\n".format(name))
		error("Error: %s\n" % str(sys.exc_info()[1]))
		error( "[Exception type: %s, raised in %s:%d]\n" % ( sys.exc_info()[1].__class__.__name__, 
		os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
		traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
		sys.exit(1)

def splice_graph_prep(bam_file):
	orig_path = os.getcwd()
	bam_name = os.path.basename(bam_file)
	ofolder = bam_name.strip(".bam$")
	if not os.path.isdir(ofolder):
		os.mkdir(ofolder)
	os.chdir(ofolder)
	gtf = "/home/patrick/72_roberto_splicing/annotation/Homo_sapiens.GRCh37.74.ucsc.gtf"
	fa = "/home/patrick/Reference_Genomes/hg19/UCSC/Chrom_fa/ucsc_hg19.fa"
	command = "sam_filter.py {} /home/patrick/Programs/SpliceGrapher-0.2.4/classifiers/Homo_sapiens.zip -f {} -m {} -o {}_filt.sam\n".format(bam_file, fa, gtf, ofolder)
	print command
	subprocess.call(command.split())
	if not os.path.isdir("predictions"):
		os.mkdir("predictions")
	if not os.path.isdir("forms"):
		os.mkdir("forms")
	command2 = "predict_graphs.py {}_filt.sam -m {} -d predictions".format(ofolder, gtf)
	subprocess.call(command2.split())
	command3 = "find_splice_forms.py {}_filt.sam -m {} -d forms".format(ofolder, gtf) #Then use spliceg_to_gtf.py to get usable gene models
	subprocess.call(command2.split())
	os.chdir(orig_path)

def dexseq_prep_fun(args):
	return dexseq_prep(*args)

def spliceg_fun(args):
	return splice_graph_prep(*args)

def main():
	parser = argparse.ArgumentParser(description='Overview of a few programs for Splicing analysis\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")

	miso_parser = subparsers.add_parser('miso', help="Runs MISO. Does not run compare_samples. Please provide already processed GTF files using index_gff")
	miso_parser.add_argument('-c','--config', help='Config file containing Conditions', required=True)
	miso_parser.add_argument('-i','--index', help='Already indexed GTF directory', required=False)
	miso_parser.add_argument('-n','--insert', help='Insert size', default=None, required=False)
	miso_parser.add_argument('-s','--sd', help='SD', required=False)
	miso_parser.add_argument('-r','--len', help='read length', required=False)
	miso_parser.add_argument('-a','--name', help='AS event name', required=False)

	dexseq_parser = subparsers.add_parser('dexseq', help="Runs DEXSEQ, run each comparison in a different directory")
	dexseq_parser.add_argument('-c','--config', help='Config file containing Conditions and Comparisons, please see documentation for usage!', required=True)
	dexseq_parser.add_argument('-r', action='store_true', help='Will run counting instead of differential analysis', required=False)
	dexseq_parser.add_argument('-g','--gtf', help='GTF file formatted by dexseq using deseq_prepare_annotation.py', required=True)
	dexseq_parser.add_argument('-t','--threads', help='threads, default=1', default=1, required=False)
	dexseq_parser.add_argument('-p', action='store_true', help='Use if samples are paired end.', required=False)
	dexseq_parser.add_argument('-o','--orientation', help='Options are yes, no or reverse. Test First!!!', required=False)

	mats_parser = subparsers.add_parser('mats', help="Runs MATS")
	mats_parser.add_argument('-c','--config', help='Config file containing bam files, please see documentation for usage!', required=True)
	mats_parser.add_argument('-g','--gtf', help='GTF file formatted by dexseq', required=True)
	mats_parser.add_argument('-n','--insert', help='Insert size, use for paired end data. In order to make things simple, will assume all bams are the same', default=None, required=False)
	mats_parser.add_argument('-s','--sd', help='SD', required=False)
	mats_parser.add_argument('-r','--len', help='read length', required=False)

	splice_parser = subparsers.add_parser('spliceG', help="Runs Splicegrapher")
	splice_parser.add_argument('-c','--config', help='Config file containing bam files, please see documentation for usage!', required=True)
	splice_parser.add_argument('-t','--threads', help='threads, default=1', default=1, required=False)
	#splice_parser.add_argument('-g','--gtf', help='GTF file formatted by dexseq', required=True)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	conditions = ConfigSectionMap("Conditions", Config)
	
	if args["subparser_name"] == "miso":	
		miso_run(conditions, args["index"], args["name"], args["insert"], args["sd"], args["len"])

	elif args["subparser_name"] == "spliceG":
		pool = Pool(int(args["threads"]))
		pool.map(spliceg_fun, itertools.izip(list(conditions.keys()))) ##Running annotation in parallel
		pool.close()
		pool.join()
		#for bam in sorted(conditions):
		#	splice_graph_prep(bam)

	elif args["subparser_name"] == "dexseq":
		dexseq_dir = "/raid/home/patrick/R/x86_64-pc-linux-gnu-library/3.1/DEXSeq"
		if args["r"]:
			pool = Pool(int(args["threads"]))
			pool.map(dexseq_prep_fun, itertools.izip(list(conditions.keys()), itertools.repeat(dexseq_dir), itertools.repeat(args["gtf"]), itertools.repeat(args["p"]),
				itertools.repeat(args["orientation"]))) ##Running annotation in parallel
			pool.close()
			pool.join()
		else:
			create_design_for_R(conditions)
			comparisons = ConfigSectionMap("Comparisons", Config)
			for comp in comparisons:
					c = comparisons[comp].split(",")
					comps = [x.strip(' ') for x in c]
					rscript = dexseq_run(conditions, comps[0], comps[1], args["gtf"])
					run_rcode(rscript, "dexseq.R")

	elif args["subparser_name"] == "mats":
		mats_program = "/home/patrick/Programs/rMATS.3.0.9/RNASeq-MATS.py"
		comparisons = ConfigSectionMap("Comparisons", Config)
		for comp in comparisons:
			c = comparisons[comp].split(",")
			comps = [x.strip(' ') for x in c]
			mats(conditions, comps[0], comps[1], args["gtf"], args["insert"], args["sd"], args["len"], mats_program)
