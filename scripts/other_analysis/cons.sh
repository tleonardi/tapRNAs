#!/bin/bash
set -xe -o pipefail
shopt -s extglob
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh

# table browser of the UCSC genome browser
mkdir -p $BASEDIR/downstream/conservation
if [[ ! -f $DATA/hg38.phastCons100way.bedg ]]; then
	$BIGWIG_TO_BEDGRAPH $DATA/hg38.phastCons100way.bw $DATA/hg38.phastCons100way.bedg
fi


# Merge all the exons 
$B12ToB6 -i $BASEDIR/nc/posConNCwithCodingPartner.bed | sort -k1,1 -k2,2n -k3,3n | $MERGE -c 4 -o distinct -i - > $BASEDIR/downstream/conservation/merged_pcRNA_exons.bed
$B12ToB6 -i $BASEDIR/data/Gencode.v21.lincRNAs.bed | sort -k1,1 -k2,2n -k3,3n | $MERGE -c 4 -o distinct -i - > $BASEDIR/downstream/conservation/merged_lncRNA_exons.bed
$B12ToB6 -i $BASEDIR/coding/Gencode-current_hsa_coding.bed | sort -k1,1 -k2,2n -k3,3n | $MERGE -c 4 -o distinct -i - > $BASEDIR/downstream/conservation/merged_coding_exons.bed


$MAP -a $BASEDIR/downstream/conservation/merged_pcRNA_exons.bed -b $DATA/hg38.phastCons100way.bedg -c 4 -o mean -null 0 > $BASEDIR/downstream/conservation/merged_pcRNA_exons.bed.map
$MAP -a $BASEDIR/downstream/conservation/merged_lncRNA_exons.bed -b $DATA/hg38.phastCons100way.bedg -c 4 -o mean -null 0 > $BASEDIR/downstream/conservation/merged_lncRNA_exons.bed.map
$MAP -a $BASEDIR/downstream/conservation/merged_coding_exons.bed -b $DATA/hg38.phastCons100way.bedg -c 4 -o mean -null 0 > $BASEDIR/downstream/conservation/merged_coding_exons.bed.map


Rscript <(cat <<EOF
library(ggplot2)
a <- read.delim("$BASEDIR/downstream/conservation/merged_pcRNA_exons.bed.map", header=F)
b <- read.delim("$BASEDIR/downstream/conservation/merged_lncRNA_exons.bed.map", header=F)
c <- read.delim("$BASEDIR/downstream/conservation/merged_coding_exons.bed.map", header=F)
a\$type <- "pcRNAs"
b\$type <- "lncRNAs"
c\$type <- "coding"
comb <-rbind(a,b,c)
pdf("$BASEDIR/downstream/conservation/plot.pdf")
ggplot(comb, aes(x=V5, fill=type)) + geom_density(alpha=0.2) + xlab("Mean exon conservation score")
dev.off()
EOF)

