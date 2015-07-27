#!/usr/bin/python

########################################################################
# 15 May 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import sys, os, re
import argparse
import pyrnatools
import pysam
import pybedtools
import pkg_resources
import ConfigParser
from multiprocessing import Pool, Manager

def convert_bam_bed(bam, name, paired, outdir):
	count = 0
	print "==> Converting bam to bed...\n"
#	if aligner=="T":
	outbam = open("{}/{}.unique.bam".format(outdir, name), "wb")
	filtered_bam = pysam.view( "-bq 50", bam) ##Filters for uniquely aligned reads!
	for read in filtered_bam:
		count += 1 
		outbam.write(read)

	inbam = pybedtools.BedTool("{}/{}.unique.bam".format(outdir, name))
	bed = inbam.bam_to_bed(split=True)
	bed.saveas("{}/{}.BED".format(outdir, name))
	#STAR conversion
	#elif aligner=="S":
#		samfile = pysam.Samfile(name+".bam", "rb")
#		for alignedread in samfile.fetch():
#			count += 1
#		samfile.close()
#		inbam = pybedtools.BedTool(name+".bam")
#		bed = inbam.bam_to_bed(split=True)
#		bed.saveas(name+".BED")
	if paired:
		count /= 2
	return count

def change_ens_ucsc_for_bed(name, outdir):
	outbed2 = open("{}/{}_ucsc.BED".format(outdir, name), "w")
	print "==> Converting Ensembl to UCSC chromosomes...\n"
	with open("{}/{}.BED".format(outdir, name)) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
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
			outbed2.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(new_chr, word[1], word[2], word[3], word[4], word[5])),
	outbed2.close()

##Must include scaling!
def genomeCoverage(name, genome, outdir, rpm=None, split=False):
	print "==> Converting bed to bedGraph...\n"
	if split:
		if rpm:
			command = "genomeCoverageBed -bg -strand + -scale {} -i {}/{}.BED -g {} > {}/{}_pos_rpm.bedGraph".format(rpm, outdir, name, genome, outdir, name)
			subprocess.call(command, shell=True)
			command = "genomeCoverageBed -bg -strand - -scale {} -i {}/{}.BED -g {} > {}/{}_neg_rpm.bedGraph".format(rpm, outdir, name, genome, outdir, name)
			subprocess.call(command, shell=True)
			output1 = "{}/{}_pos_rpm.bedGraph".format(outdir, name)
			output2 = "{}/{}_neg_rpm.bedGraph".format(outdir, name)
		else:
			command = "genomeCoverageBed -bg -strand + -i {}/{}.BED -g {} > {}/{}_pos.bedGraph".format(outdir, name, genome, outdir, name)
			subprocess.call(command, shell=True)
			command = "genomeCoverageBed -bg -strand - -i {}/{}.BED -g {} > {}/{}_neg.bedGraph".format(outdir, name, genome, outdir, name)
			subprocess.call(command, shell=True)
			output1 = "{}/{}_pos.bedGraph".format(outdir, name)
			output2 = "{}/{}_neg.bedGraph".format(outdir, name)
		output = [output1, output2]
	else:
		if rpm:
			command = "genomeCoverageBed -bg -scale {} -i {}/{}.BED -g {} > {}/{}_rpm.bedGraph".format(rpm, outdir, name, genome, outdir, name)
			subprocess.call(command, shell=True)
			output = "{}/{}_rpm.bedGraph".format(outdir, name)
		else:
			command = "genomeCoverageBed -bg -i {} -g {}/{}.BED > {}/{}_rpm.bedGraph".format(outdir, name, genome, outdir, name)
			output = "{}/{}.bedGraph".format(outdir, name)
	return output

def bedgraphtobigwig(bedgraph, chrom, split):
	print "==> Converting bedGraph to bigWig...\n"
	if split:
		for bedg in bedgraph:
			bw = re.sub(".bedGraph$", ".bw", bedg)
			command = ["bedGraphToBigWig", bedg, chrom, bw]
			subprocess.call(command)
	else:	
		bw = re.sub(".bedGraph$", ".bw", bedgraph)
		command = ["bedGraphToBigWig", bedgraph, chrom, bw]
		subprocess.call(command)

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

def main():
	parser = argparse.ArgumentParser(description="Processes RNA-seq samples to bigWig tracks.\nIf Tophat2 is specified, this will pull out the uniquely mapped reads\nOtherwiseit is assumed that the bam file is already uniquely mapped!")
	parser.add_argument('-c', '--config', help='Contains [Conditions] with bam files as keys. Use either this for multiple files or input for one at a time', required=False)
	parser.add_argument('-i', '--input', help='Bam file from aligner etc.', required=False)
	parser.add_argument('-p', action='store_true', help='Use if samples are paired end.', required=False)
	parser.add_argument('-g', '--genome', help='Genome the samples are aligned to, options include mm10/mm9/hg19', required=True)
	parser.add_argument('-s', action='store_true', help='Split tracks by strand', required=False) 
	parser.add_argument('-rpm', action='store_true', help='Scale to RPM', required=False) 
	parser.add_argument('-ens', action='store_true', help='If samples are aligned to ensembl genome, convert to UCSC coordinates', required=False) 
	parser.add_argument('-o', '--outdir', help='Output directory', required=True)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	chrom = pkg_resources.resource_filename('pyrnatools', 'data/{}.chrom.sizes'.format(args["genome"]))
	if not os.path.isfile(chrom):
		raise Exception("Unsupported Genome!")
	
	if args["input"]:#
		name = os.path.basename(args["input"])
		name = re.sub(".bam$", "", name)
		unique_reads = convert_bam_bed(args["input"], name, args["p"], args["outdir"])

		if args["ens"]:
			change_ens_ucsc_for_bed(name, args["outdir"])
			name = name+"_ucsc"
		if args["rpm"]:
			scale = float(1000000)/int(unique_reads)
			bedgraph = genomeCoverage(name, chrom, args["outdir"], rpm=scale, split=args["s"])	
		else:
			bedgraph = genomeCoverage(name, chrom, args["outdir"], split=args["s"])	
		bedgraphtobigwig(bedgraph, chrom, args["s"])
	
	elif args["config"]:
		Config = ConfigParser.ConfigParser()
		Config.optionxform = str
		Config.read(args["config"])

		conditions = ConfigSectionMap("Conditions", Config)
		
		for key in conditions:
			name = os.path.basename(args["input"])
			name = re.sub(".bam$", "", name)
			unique_reads = convert_bam_bed(bam, name, args["p"], args["outdir"])
			
			if args["ens"]:
				change_ens_ucsc_for_bed(name, args["outdir"])
				name = name+"_ucsc"
			if args["rpm"]:
				scale = float(1000000)/int(unique_reads)
				bedgraph = genomeCoverage(name, args["genome"], args["outdir"], rpm=scale, split=args["s"])	
			else:
				bedgraph = genomeCoverage(name, args["genome"], args["outdir"], split=args["s"])	
			bedgraphtobigwig(bedgraph, chrom, args["s"])
	
