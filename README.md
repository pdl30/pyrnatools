#pyrnatools 


### Installation

Clone this repository and then:

```bash
$ cd pyrnatools/
$ python setup.py install --user
```

This will install the scripts in the pyrnatools/scripts directory. For more information on the individual scripts, use the --help command after each script. 

##Core Pipeline

 - pyrna_align.py -> Wrapper for FASTQC, cutadapt and tophat2 sequence aligner
 - pyrna_count.py -> Wrapper for Htseq-count and gfold-count
 - pyrna_diff.py -> Wrapper for deseq2/gfold differential expression analysis
 - pyrna_ucsc.py -> Converts bam file to UCSC formatted bigWig

##Additional Tools
 - pyrna_pair_tools.py -> Infers sample orientation and calculates insert size for paired end RNA-seq data
 - pyrna_denovo.py -> Cufflinks wrapper for denovo assembly
 - pyrna_download.py -> Tool for downloading samples from GEO

