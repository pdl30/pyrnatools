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

def write_deseq(ifile, sample_dict, cond1, cond2, padj):
	print "==> Running differental expression analysis...\n"
	rscript =  "library(DESeq2)\n"
	rscript += "pdata <- read.table('tmp_design.txt', header=T)\n"
	#Make sure they match!
	rscript += "counts <- read.table('{}', sep='\\t', header=T, row.names=1)\n".format(ifile)
	rscript += "rnaseq_dds <- DESeqDataSetFromMatrix(countData = counts, colData = data.frame(pdata), design = ~ condition)\n"
	rscript += "rnaseq_dds$condition <- factor(rnaseq_dds$condition, levels=unique(pdata[,3]))\n"
	rscript += "rnaseq_dds <- DESeq(rnaseq_dds)\n"	
	rscript += "rnaseq_res <- results(rnaseq_dds, contrast=c('condition','{0}','{1}'))\n".format(cond1, cond2)
	rscript += "rnaseq_sig <- rnaseq_res[which(rnaseq_res$padj <= {}),]\n".format(padj)
	rscript += "write.table(rnaseq_sig, file='{0}_vs_{1}_deseq2_significant.tsv', sep='\\t', quote=F)\n".format(cond1, cond2)
	rscript += "write.table(rnaseq_res, file='{0}_vs_{1}_deseq2_analysis.tsv', sep='\\t', quote=F)\n".format(cond1, cond2)
	return rscript



#def html_report():
#	rscript =  "library(ReportingTools)\n"
#	rscript += "library(DESeq2)\n"
#	rscript += "conditions <- c(rep("case",3), rep("control", 3))"
#	mockRna.dse <- DESeqDataSetFromMatrix(countData = mockRnaSeqData, colData = as.data.frame(conditions), design = ~ conditions)
#	colData(mockRna.dse)$conditions <- factor(colData(mockRna.dse)$conditions, levels=c("control", "case"))
#	 ## Get a DESeqDataSet object
#	mockRna.dse <- DESeq(mockRna.dse)Now the results can be written to a report using theDESeqDataSetobject.
#	des2Report <- HTMLReport(shortName ='RNAseq_analysis_with_DESeq2', title ='RNA-seq analysis of differential expression using DESeq2',reportDirectory = "./reports")
#	publish(mockRna.dse,des2Report, pvalueCutoff=0.05, annotation.db="org.Mm.eg.db", factor = colData(mockRna.dse)$conditions, reportDir="./reports")
#	finish(des2Report)
