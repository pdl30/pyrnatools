#!/usr/bin/python

########################################################################
# 15 May 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import re, os, sys
import argparse
import math
import pkg_resources

def head(ifile, ofile):
	output = open(ofile, "w")
	with open(ifile) as myfile:
		head = [next(myfile) for x in xrange(100000)]
	size = 0
	read1 = head[1]
	bases = read1.rstrip()
	for h in head:
		output.write(h),
	output.close()
	return len(list(bases))
		
def run_bowtie(fq1, fq2, index, fq):
	f = open('/dev/null', 'w')
	if fq:
		name = re.sub("_1.fq", "", fq1)
	else:
		name = re.sub("_1.fastq", "", fq1)
	read_size = head(fq1, "{0}_trun1.fq".format(name))
	head(fq2, "{0}_trun2.fq".format(name))
	c3 = "bowtie --sam -I 0 -X 500 {0} -1 {1}_trun1.fq -2 {1}_trun2.fq > {1}_trun.sam".format(index, name)
	subprocess.call(c3, shell=True, stderr=f)
	return name

def getmeanval(dic,maxbound=-1):
	nsum=0;  n=0;	
	for (k,v) in dic.items():
		if maxbound!=-1 and k>maxbound:
			continue
		nsum=nsum+k*v;
		n=n+v;
	meanv=nsum*1.0/n;
	nsum=0; n=0;
	for (k,v) in dic.items():
		if maxbound!=-1 and k>maxbound:
			continue;
		nsum=nsum+(k-meanv)*(k-meanv)*v;
		n=n+v;
		varv=math.sqrt(nsum*1.0/(n-1));
	return [meanv,varv];
 
def get_insert(name):
	samfile = name+"_trun.sam"
	plrdlen={};
	plrdspan={};
	objmrl=re.compile('([0-9]+)M$');
	objmtj=re.compile('NH:i:(\d+)');
	nline=0
	sam = open(samfile, "r")
	for lines in sam:
		field=lines.strip().split();
		nline=nline+1;
		if len(field)<12:
			continue;
		try:
			mrl=objmrl.match(field[5]);
			if mrl==None: # ignore non-perfect reads
				continue;
			readlen=int(mrl.group(1));
			if readlen in plrdlen.keys():
				plrdlen[readlen]=plrdlen[readlen]+1;
			else:
				plrdlen[readlen]=1;
			if field[6]!='=':
				continue;
			dist=int(field[8]);
			if dist<=0: # ignore neg dist
				continue;
			mtj=objmtj.search(lines);

			if dist in plrdspan.keys():
				plrdspan[dist]=plrdspan[dist]+1;
			else:
				plrdspan[dist]=1;
		except ValueError:
			continue;
	if len(plrdspan)==0:
		print('No qualified paired-end reads found. Are they single-end reads?');
	else:
		maxv=max(plrdspan,key=plrdspan.get);
		spanval=getmeanval(plrdspan,maxbound=maxv*3);
		return spanval

def infer_command(name, refbed):
	samfile = name+"_trun.sam"
	infercommand = "infer_experiment.py -i {} -r {} > infer_res.txt".format(samfile, refbed)
	print infercommand
	subprocess.call(infercommand, shell=True)
	reverse = read_infer()
	return reverse

def read_infer():
	with open("infer_res.txt") as f:
		for line in f:
			line = line.rstrip()
			if line.startswith("Fraction of reads explained by \"1++,1--,2+-,2-+\": "):
				per1 = line.lstrip("Fraction of reads explained by \"1++,1--,2+-,2-+\": ")
			elif line.startswith("Fraction of reads explained by \"1+-,1-+,2++,2--\": "):
				per2 = line.lstrip("Fraction of reads explained by \"1+-,1-+,2++,2--\": ")
	if float(per1) > float(per2):
		reverse = False
	else:
		reverse = True
	return reverse

#def main():
if __name__ == "__main__":	
	parser = argparse.ArgumentParser(description='Collection of tools for paired end RNA-seq data\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	insert_parser = subparsers.add_parser('insert', help="Calculates insert size for paired end data")
	type_parser = subparsers.add_parser('type', help="Will infer experiment orientation")
	type_parser.add_argument('-p', '--pair', help='Paired end fastqs. Please put them in order!', required=True, nargs='+')
	type_parser.add_argument('-i', '--index', help='Bowtie1 Index, must be in UCSC format!', required=True)
	type_parser.add_argument('-g', '--genome', choices=["mm10", "hg19"], help="Genome of sample, options are hg19 and mm10")
	insert_parser.add_argument('-p', '--pair', help='Paired end fastqs. Please put them in order!', required=True, nargs='+')
	insert_parser.add_argument('-i', '--index', help='Bowtie1 Index', required=True)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	
	fq = False
	if args["pair"][0].endswith(".fq"):
		fq = True
	name = run_bowtie(args["pair"][0], args["pair"][1], args["index"], fq)
	print name
	if args["subparser_name"] == "insert":
		insert = get_insert(name)
		print "mate inner distance = {}\nmate standard deviation = {}".format(insert[0], insert[1])
	elif args["subparser_name"] == "type":
		path1 = pkg_resources.resource_filename('pyrnapipe', 'data/')
		refbed = path1 + "{}_Ensembl.bed".format(args["genome"])
		reverse = infer_command(name, refbed)
		if reverse:
			print "Orientation = reverse"
		else:
			print "Orientation = yes"