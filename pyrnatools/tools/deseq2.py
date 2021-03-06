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

def write_deseq(ifile, sample_dict, cond1, cond2, padj, f, outdir, tmpdesign):
	print "==> Running differental expression analysis...\n"
	rscript =  "suppressMessages(library(DESeq2))\n"
	rscript +=  "suppressMessages(library(RColorBrewer)); suppressMessages(library(ggplot2)); suppressMessages(library(gplots))\n"
	rscript += "pdata <- read.table('{}', header=T)\n".format(tmpdesign)
	#Make sure they match!
	rscript += "counts <- read.table('{}', sep='\\t', header=T, row.names=1)\n".format(ifile)
	if f:
		rscript += "counts <- counts[,6:dim(counts)[2]]\n"
	else:
		pass
	rscript += "rnaseq_dds <- DESeqDataSetFromMatrix(countData = counts, colData = data.frame(pdata), design = ~ condition)\n"
	rscript += "rnaseq_dds$condition <- factor(rnaseq_dds$condition, levels=unique(pdata[,3]))\n"
	rscript += "rnaseq_dds <- DESeq(rnaseq_dds)\n"	
	rscript += "rnaseq_res <- results(rnaseq_dds, contrast=c('condition','{0}','{1}'))\n".format(cond1, cond2)
	rscript += "rnaseq_sig <- rnaseq_res[which(rnaseq_res$padj <= {}),]\n".format(padj)
	rscript += "norm_counts <- counts(rnaseq_dds, normalized=TRUE)\n"
	rscript += "rownames(norm_counts) <- rownames(counts); colnames(norm_counts) <- colnames(counts)\n"
	rscript += "write.table(norm_counts, file='{2}/deseq2_counts.tsv', sep='\\t', quote=F)\n".format(cond1, cond2, outdir)
	rscript += "write.table(rnaseq_sig, file='{2}/{0}_vs_{1}_deseq2_significant.tsv', sep='\\t', quote=F)\n".format(cond1, cond2, outdir)
	rscript += "write.table(rnaseq_res, file='{2}/{0}_vs_{1}_deseq2_analysis.tsv', sep='\\t', quote=F)\n".format(cond1, cond2, outdir)
	rscript += "rld<-rlog(rnaseq_dds); colnames(rld) <- pdata$sampleName; rlogMat<-assay(rld); distsRL<-dist(t(assay(rld)))\n"
	rscript += "pdf('{0}/Sample-RLD-plots.pdf'); hmcol<-colorRampPalette(brewer.pal(9,'GnBu'))(100); mat<-as.matrix(distsRL)\n".format(outdir)
	rscript += "hc<-hclust(distsRL); par(cex.main=1);\n";
	rscript += "heatmap.2(mat,Rowv=as.dendrogram(hc),symm=TRUE,trace='none',col=rev(hmcol),margin=c(13,13),main='Sample-to-sample distances',cexRow=1,cexCol=1)\n"
	rscript += "data<-plotPCA(rld,intgroup=c('condition'),returnData=TRUE, ntop = nrow(rld)); percentVar<-round(100*attr(data,'percentVar'))\n"
	rscript += "ggplot(data,aes(PC1, PC2,color=condition,label=name))+geom_point(size=3)+xlab(paste0('PC1: ',percentVar[1],'% variance'))+ylab(paste0('PC2: ',percentVar[2],'% variance')) +geom_text(colour = 'black', alpha = 0.8, size = 2)\n"
	rscript += "dev.off()\n"

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
