#!/bin/bash
# This script defines environment variables.
# It's sourced by all shell script 
# and is also included by the Makefile
# The syntax needs to be both bash and GNU Make compatible.

#################
#   Folders	#
#################
export BASEDIR=/path/to/repo/analysis
export BIN=/path/to/repo/bin
export DATA=/path/to/repo/data
export RESULTS=/path/to/repo/results

#################
#   Binaries	#
#################
# Bedtools
export INTERSECT_BED=/path/to/bedtools/bin/intersectBed
export OVERLAP=/path/to/bedtools/bin/getOverlap
export GROUP_BY=/path/to/bedtools/bin/groupBy
export CLUSTER=/path/to/bedtools/bin/clusterBed
export MERGE=/path/to/bedtools/bin/mergeBed
export CLOSEST=/path/to/bedtools/bin/closestBed
export SAMPLE=/path/to/bedtools/bin/sampleFile
export BEDTOBAM=/path/to/bedtools/bin/bedToBam
export B12ToB6=/path/to/bedtools/bin/bed12ToBed6
export GETFASTA=/path/to/bedtools/bin/fastaFromBed
export MAP=/path/to/bedtools/bin/mapBed
export GENOME_COVERAGE_BED=/path/to/bedtools/bin/genomeCoverageBed

# CPAT
CPAT=/path/to/cpat.py

# Pinstripe
export PSTR="mono /path/to/pinstripe.exe"

# Tuxedo suite
export CUFFLINKS=/path/to/cufflinks
export CUFFNORM=/path/to/cuffnorm
export TOPHAT=/path/to/tophat2
export BOWTIE2=/path/to/bowtie2-2.1.0
export CUFFMERGE=/path/to/cuffmerge
export GFFREAD=/path/to/gffread
export CUFFQUANT=/path/to/cuffquant

# Samtools
export SAMTOOLS=/path/to/samtools

# UCSC tools
export LIFTOVER=/path/to/liftOver
export BEDGRAPH_TO_BIGWIG=/path/to/bedGraphToBigWig
export BIGWIG_TO_BEDGRAPH=/path/to/bigWigToBedGraph
export OVERLAPSELECT=/path/to/overlapSelect
export gtfToGenePred=/path/to/gtfToGenePred
export genePredToBed=/path/to/genePredToBed
export twoBitToFa=/path/to/twoBitToFa

# NCBI BLAST
export BLASTN=/path/to/blastn
export MAKEBLASTDB=/path/to/makeblastdb
export CONVERT2BLASTMASK=/path/to/convert2blastmask

# PFAM SCAN
export pfam_scan=/path/to/pfam_scan.pl
export hmmpress=/path/to/HMMER
export transeq=/path/to/EMBOSS

# WATER AND NEEDLE
export NEEDLE=/path/to/EMBOSS
export WATER=/path/to/EMBOSS

