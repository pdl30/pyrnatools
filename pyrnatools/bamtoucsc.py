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

def convert_bam_bed(name, paired):
	count = 0
	print "==> Converting bam to bed...\n"
#	if aligner=="T":
	outbam = open(name+".unique.bam", "wb")
	filtered_bam = pysam.view( "-bq 50", name+".bam") ##Filters for uniquely aligned reads!
	for read in filtered_bam:
		count += 1 
		outbam.write(read)

	inbam = pybedtools.BedTool(name+".unique.bam")
	bed = inbam.bam_to_bed(split=True)
	bed.saveas(name+".BED")
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

def change_ens_ucsc_for_bed(name):
	outbed2 = open(name+"_ucsc.BED", "w")
	print "==> Converting Ensembl to UCSC chromosomes...\n"
	with open(name+".BED") as f:
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

##Must include scaling!
def genomeCoverage(name, house=None, rpm=None):
	print "==> Converting bed to bedGraph...\n"
	inbed = pybedtools.BedTool(name+".BED")
	if house:
		outcov = inbed.genome_coverage(bg=True, genome='mm10', scale=house)
		output = name+"_house.bedGraph"
	elif rpm:
		outcov = inbed.genome_coverage(bg=True, genome='mm10', scale=rpm)
		output = name+"_rpm.bedGraph"
	else:
		outcov = inbed.genome_coverage(bg=True, genome='mm10')
		output = name+".bedGraph"
	outcov.saveas(output)
	return output

def bedgraphtobigwig(bedgraph, chrom):
	bw = re.sub(".bedGraph$", ".bw", bedgraph)
	print "==> Converting bedGraph to bigWig...\n"
	command = ["bedGraphToBigWig", bedgraph, chrom, bw]
	subprocess.call(command)

def normalise_to_housekeeper(count_file):
	print "==> Normalising to Housekeeper...\n"
	with open(count_file) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if word[0] == "ENSMUSG00000057666":
				housekeeper = int(word[1])	
		print housekeeper
	return housekeeper

def main():
	parser = argparse.ArgumentParser(description="Processes RNA-seq samples to bigWig tracks.\nIf Tophat2 is specified, this will pull out the uniquely mapped reads\nOtherwiseit is assumed that the bam file is already uniquely mapped!")
	parser.add_argument('-i', '--input', help='Bam file from aligner etc.', required=True)
	parser.add_argument('-p', action='store_true', help='Use if samples are paired end.', required=False)
	parser.add_argument('-g', '--genome', help='Genome the samples are aligned to, options include mm10/mm9/hg19', required=True)
	parser.add_argument('-rpm', action='store_true', help='Scale to RPM', required=False) 
	parser.add_argument('-a', '--house', help='Housekeeper normalisation. Input file is HTSEQ-count file containing gene for normalisation on first line', required=False)
	parser.add_argument('-ens', action='store_true', help='If samples are aligned to ensembl genome, convert to UCSC coordinates', required=False) 
	args = vars(parser.parse_args())

	chrom = pkg_resources.resource_filename('pyrnatools', 'data/{}.chrom.sizes'.format(args["genome"]))
	if not os.path.isfile(chrom):
		raise Exception("Unsupported Genome!")
	
	name = re.sub(".bam$", "", args["input"])

	unique_reads = convert_bam_bed(name, args["p"])
	#print unique_reads
	if args["ens"]:
		change_ens_ucsc_for_bed(name)
		name = name+"_ucsc"

	if args["rpm"]:
		scale = float(1000000)/int(unique_reads)
		bedgraph = genomeCoverage(name, rpm=scale)	
	elif args["house"]:
		house = normalise_to_housekeeper(args["house"])
		scale = float(1000)/int(house) #Works and checked
		bedgraph = genomeCoverage(name, house=scale)
	else:
		bedgraph = genomeCoverage(name)
	
	bedgraphtobigwig(bedgraph, chrom)
