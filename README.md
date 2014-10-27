#pyrnatools 


### Installation

Clone this repository and then:

```bash
$ cd pyrnatools/
$ python setup.py install --user
```

This will install the scripts in the pyrnatools/scripts directory. For more information on the individual scripts, use the --help command after each script. 

##Core Pipeline

 - pyrna_align.py -> Wrapper for FASTQC, cutadapt and tophat2 sequence aligner. 
 - pyrna_diff_exp.py -> Wrapper for Htseq-count and deseq2/gfold differential expression analysis
 - pyrna_viz.py -> Converts bam file to UCSC formatted bigWig

##Additional Tools
 - pyrna_insert_size.py -> Calculates insert size for paired end RNA-seq data
 - pyrna_diff_summary.py -> Simple script that takes a counts matrix and a directory of DESEQ2/GFOLD results and create a summary matrix
 - pyrna_select_sigs.py -> Takes pyrna_diff_summary.py and prints out significant genes
 - pyrna_denovo.py -> Cufflinks wrapper for denovo assembly
 - pyrna_htseq_to_gfold.py -> Take a HTseq-count output and converts it to GFOLD formatted counts file
 
