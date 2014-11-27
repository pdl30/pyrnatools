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
import shutil

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
	download2 = "wget -c -nv -q 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='{0}'&targ=self&view=quick&form=text&token=_blank' -O tmp/{0}.soft".format(gsm)
	subprocess.call(download2, shell=True)
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
	subprocess.call(command, shell=True)
	#Remove old SRR fastqs
	for fq in fqs:
		os.remove(fq)
	os.chdir(old_path)

def create_staging_table(accession, out):
	output = open(out, "w")
	header = """id\tCT_General\tCT_Subtype\tbto\tfactor\tfactorGeneId\tdetails\tGSE\tGSM\tgroupname\trepository\texp_type\tfilename\tcontrol\tSRX\tpvaluemerged\t
		paired\tmapping_genome\ttotal_raw_reads\ttrimmed\treads_after_trimming\tuniquely_mappable_reads\tFinal_peak_number\taligner_used\traw_data_code\tIn_house_code\t
		Completion_date\tAdded_to_compendium\tstatus\tprocessing_comments\talignment_report\torganism\tsubmitter\n"""
	output.write(header),
	command = "wget -c -nv -q -O tmp.txt http://www.ebi.ac.uk/arrayexpress/files/{0}/{0}.sdrf.txt".format(accession)
	command = "wget -c -nv -q -O tmp2.txt http://www.ebi.ac.uk/arrayexpress/files/{0}/{0}.idf.txt".format(accession)
	subprocess.call(command, shell=True)

	f = open("tmp.txt", "r")
	lines = f.readlines()
	header = lines[0].rstrip()
	head = header.split("\t")
	cell_type =  [i for i, x in enumerate(head) if x.startswith("Characteristics") and "cell type" in x]
	cell_line =  [i for i, x in enumerate(head) if x.startswith("Characteristics") and "cell line" in x]
	organism =  [i for i, x in enumerate(head) if x.startswith("Characteristics") and "organism" in x]
	library_layout = [i for i, x in enumerate(head) if x.startswith("Comment") and "LIBRARY_LAYOUT" in x]
	library_type = [i for i, x in enumerate(head) if x.startswith("Comment") and "LIBRARY_STRATEGY" in x]
	names = [i for i, x in enumerate(head) if "Assay Name" in x]

	for i in xrange(1, len(lines)):
		line = lines[i]
		line = line.rstrip()
		word = line.split("\t")
		#print cell_type, cell_line, organism, library_layout, library_type
		print word[cell_type[0]], word[cell_line[0]], word[organism[0]], word[library_layout[0]], word[library_type[0]]

	with open("tmp2.txt") as f:
		for line in f:
			if line.startswith("Protocol Description"):
				description = re.sub("Protocol Description\t", "", line)


def download_ebi(accession, staging):
	#Example EBI Format
	#Source Name	Comment[ENA_SAMPLE]	Material Type	Provider	Characteristics[organism]	Characteristics[specimen with known storage state]	
	#Characteristics[cell line]	Characteristics[strain]	Characteristics[cell type]	
	#Protocol REF	Protocol REF	Protocol REF	Protocol REF	Extract Name	Material Type	
	#Comment[LIBRARY_LAYOUT]	Comment[LIBRARY_SOURCE]	Comment[LIBRARY_STRATEGY]
	command = "wget -b -o tmp.txt http://www.ebi.ac.uk/arrayexpress/files/{0}/{0}.sdrf.txt".format(accession)
	with open("tmp.txt") as f:
		header = next(f)
		header = header.rstrip()
		head = header.split("\t")
		cell_type = filter(lambda x:re.search(r'cell_type', x), names)
		for line in f:
			line = line.s

def download_function(args):
	return downloader(*args)

def main():
	parser = argparse.ArgumentParser(description='Downloads data from GEO and Arrayexpress\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	geo_parser = subparsers.add_parser('geo', help="Use if downloading data from GEO")
	ebi_parser = subparsers.add_parser('ebi', help="Use if downloading data from Arrayexpress")

	geo_parser.add_argument('-c','--config', help='Input ConfigParser file containing [Conditions] with GSMs as keys', required=True)
	geo_parser.add_argument('-t','--threads', help='Number of threads, default=20', default=20, required=False)

	ebi_parser.add_argument('-a', '--accession', help='Arrayexpress accession number', required=True )
	ebi_parser.add_argument('-s', action='store_true', help='This will create a staging table and won\'t process the samples', required=False )
	ebi_parser.add_argument('-o', '--output', help='Output staging table', required=False )
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	if args["subparser_name"] == "geo":
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
		shutil.rmtree('tmp/')
	elif args["subparser_name"]	 == "ebi":
		if args["s"]:
			create_staging_table(args["accession"], args["output"])
		else:
			download_ebi(args["accession"])

main()