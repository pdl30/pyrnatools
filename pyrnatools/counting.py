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
from pyrnatools.tools import gfold, deseq2
from multiprocessing import Pool
import tempfile

def annotate_sam(bam_file, gtf_file, stranded, outdir):
	print "==> Counting sam file...\n"
	if bam_file.endswith(".bam"):
		name = os.path.basename(bam_file)
		count_file = re.sub(".bam$", ".count", name)
		command = "htseq-count --stranded={} --quiet -f bam {} {} > {}/{}".format(stranded, bam_file,  gtf_file, outdir, count_file)
	else:
		name = os.path.basename(bam_file)
		count_file = re.sub(".sam$", ".count", name) #Just in case supplied file is sam!
		command = "htseq-count --stranded={} --quiet -f sam {} {} > {}/{}".format(stranded, bam_file,  gtf_file, outdir, count_file)
	subprocess.call(command, shell=True)

def join_counts(idict, outdir):
	data = defaultdict(list)
	output = open("{}/combined_counts.tsv".format(outdir), "w")
	output.write("ID"),
	for bam in sorted(idict):
		if bam.endswith(".bam"):
			name = os.path.basename(bam)
			count_file = re.sub(".bam$", ".count", name)
		else:
			name = os.path.basename(bam)
			count_file = re.sub(".sam$", ".count", name)
		output.write("\t{}".format(bam)),
		with open("{}/{}".format(outdir, count_file)) as f:
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				if word[0].startswith("__"):
					pass
				else:
					data[word[0]].append(word[1])
	output.write("\n"),
	for key2 in sorted(data):
		data2 = data[key2]
		output.write(key2+"\t" + "\t".join(data2) + "\n"),
	output.close()

def featurecounts(conditions, threads, gtf_file, stranded, paired, outfile, bam=None):
	# -s  0 (unstranded), 1 (stranded) and 2 (reversely stranded). 0 by default
	command = "featureCounts -a {} -T {} -o {}".format(gtf_file, threads, outfile)
	command = command.split()
	if paired:
		command.append(" -p")
	if stranded == "yes":
		command.append(" -s 1")
	elif stranded == "no":
		command.append("-s 0")
	elif stranded == "reverse":
		command.append("-s 2")
	if bam:
		command.append(bam)
	else:
		command.extend(sorted(list(conditions.keys())))
	print command
	subprocess.call(command)
	#eatureCounts -p -a /home/patrick/Reference_Genomes/mm10/Ensembl/76/Mus_musculus.GRCm38.76_ucsc.gtf -o tmp.count 
def run_gfold_count(args):
	return gfold.run_gfold_c(*args)
def anno_function(args):
	return annotate_sam(*args)

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

def convert_gtf_to_ucsc(igtf):
	a = tempfile.NamedTemporaryFile(delete=False)
	with open(igtf) as f:
		for line in f:
			new_chr = None
			line = line.rstrip()
			word = line.split("\t")
			if line.startswith("#"):
				pass
			if re.match(r"^\d", word[0]):
				new_chr = "chr" + word[0]
			elif re.match(r"^X", word[0]):
				new_chr = "chrX"
			elif re.match(r"^Y", word[0]):
				new_chr = "chrY"
			elif word[0] == "MT":
				new_chr = "chrM"
			else:
				pass
			if new_chr:
				a.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(new_chr, word[1], word[2], word[3], word[4], word[5], word[6], word[7], word[8])),
	a.close()
	return a.name

def infer(bam, genome):
	path1 = "/home/patrick/Reference_Genomes/pyngspipe_references/"
	if genome == "hg19":
		refbed = path1 + "hg19/hg19_Ensembl.bed"
	elif genome == "mm10":
		refbed = path1 + "mm10/mm10_Ensembl.bed"
	infercommand = "infer_experiment.py -i {} -r {} > infer_res.txt".format(bam, refbed)
	subprocess.call(infercommand, shell=True)

def main():
	parser = argparse.ArgumentParser(description='Counts features from BAM files\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")

	htseq_parser = subparsers.add_parser('htseq', help="Runs HTseq count")
	htseq_parser.add_argument('-c','--config', help='Config file containing [Conditions] with bam files as keys and colnames as values', required=False)
	htseq_parser.add_argument('-i','--input', help='Input bam file', required=False)
	htseq_parser.add_argument('-g','--gtf', help='GTF file', required=False)
	htseq_parser.add_argument('-s','--stranded', help='Option for HTSeq, default=yes', default="yes", required=False)
	htseq_parser.add_argument('-t','--threads', help='Number of threads, default=8', default=8, required=False)
	htseq_parser.add_argument('-o','--output', help='Output counts file directory, default is current directory', required=False) #Will output all counts files and combined file if specified
	htseq_parser.add_argument('-e', action='store_true', help='Converts the GTF to UCSC format', required=False)

	gfold_parser = subparsers.add_parser('gfold', help="Runs GFOLD Count")
	gfold_parser.add_argument('-c','--config', help='Config file containing [Conditions], please see documentation for usage!\nPlease use a unique name for every input bam file!', required=False)
	gfold_parser.add_argument('-i','--input', help='Input bam file', required=False)
	gfold_parser.add_argument('-n',action='store_true', help='Gapdh Normlisation', required=False)
	gfold_parser.add_argument('-g','--gtf', help='GTF file', required=False)
	gfold_parser.add_argument('-t','--threads', help='Number of threads, default=8', default=8, required=False)
	gfold_parser.add_argument('-o','--output', help='Output counts file directory, default is current directory', required=False)

	feat_parser = subparsers.add_parser('feat', help="Runs featureCount")
	feat_parser.add_argument('-c','--config', help='Config file containing [Conditions] with bam files as keys and colnames as values', required=False)
	feat_parser.add_argument('-i','--input', help='Input bam file', required=False)
	feat_parser.add_argument('-g','--gtf', help='GTF file', required=False)
	feat_parser.add_argument('-s','--stranded', help='Option for featueCount, default=yes', default="yes", required=False)
	feat_parser.add_argument('-t','--threads', help='Number of threads, default=8', default=8, required=False)
	feat_parser.add_argument('-o','--output', help='Output counts file', required=False) #Will output all counts files and combined file if specified
	feat_parser.add_argument('-p', action='store_true', help='Samples are paired end', required=False)
	#feat_parser.add_argument('-e', action='store_true', help='Converts the GTF to UCSC format', required=False)

	infer_parser = subparsers.add_parser('infer', help="Runs infer_experiment.py")
	infer_parser.add_argument('-i','--input', help='Input bam file', required=True)
	infer_parser.add_argument('-g','--genome', help='Options are hg19/mm10', required=True)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	
	if args["output"]:
		output = args["output"]
	else:
		output = os.getcwd()
	
	if args["subparser_name"] == "gfold":
		if args["config"]:
			Config = ConfigParser.ConfigParser()
			Config.optionxform = str
			Config.read(args["config"])

			#Read design matrix and create list of conditions and directories
			conditions = ConfigSectionMap("Conditions", Config)
			pool = Pool(int(args["threads"]))
			pool.map(run_gfold_count, itertools.izip(list(conditions.keys()), itertools.repeat(args["gtf"]), itertools.repeat(output))) ##Running annotation in parallel
			pool.close()
			pool.join()
		elif args["input"]:
			gfold.run_gfold_c(args["input"], args["gtf"], output)
	elif args["subparser_name"] == "feat":
		if args["config"]:
			Config = ConfigParser.ConfigParser()
			Config.optionxform = str
			Config.read(args["config"])
			conditions = ConfigSectionMap("Conditions", Config)
			featurecounts(conditions, int(args["threads"]), args["gtf"], args["stranded"], args["p"], args["output"])
		elif args["input"]:
			featurecounts(None, int(args["threads"]), args["gtf"], args["stranded"], args["p"], args["output"], bam=args["input"])

	elif args["subparser_name"] == "htseq":
		if args["e"]:
			gtf = convert_gtf_to_ucsc(args["gtf"])
		else:
			gtf = args["gtf"]

		if args["config"]:
			Config = ConfigParser.ConfigParser()
			Config.optionxform = str
			Config.read(args["config"])

			#Read design matrix and create list of conditions and directories
			conditions = ConfigSectionMap("Conditions", Config)

			pool = Pool(int(args["threads"]))
			pool.map(anno_function, itertools.izip(list(conditions.keys()), itertools.repeat(gtf), itertools.repeat(args["stranded"]), itertools.repeat(output))) ##Running annotation in parallel
			pool.close()
			pool.join()	
			join_counts(conditions, output)
		elif args["input"]:
			annotate_sam(args["input"], gtf, args["stranded"], output)
	elif args["subparser_name"] == "infer":
		infer(args["bam"], args["genome"])