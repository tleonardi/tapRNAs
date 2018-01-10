#!/bin/bash
set -xe -o pipefail
export LC_ALL=C
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../include.sh

PC=$BASEDIR/../results/posConNCwithCodingPartner.bed
# Needs to be exported for R
export DIR=$BASEDIR/downstream/cageSupport

mkdir -p $DIR


wget -O - "http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_coord.bed.gz" | gunzip > $DIR/hg19.cage_peak_phase1and2combined_coord.bed
wget -O - "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" | gunzip > $DIR/hg19ToHg38.over.chain
$LIFTOVER $DIR/hg19.cage_peak_phase1and2combined_coord.bed $DIR/hg19ToHg38.over.chain $DIR/hg38.cage_peak_coord.bed.unsorted $DIR/hg38.cage_peak_coord.bed_unmapped

sort -k1,1 -k2,2n -k3,3n $DIR/hg38.cage_peak_coord.bed.unsorted > $DIR/hg38.cage_peak_coord.bed

# Print 5p exon
awk 'BEGIN{OFS=FS="\t"}{
        split($11,SIZES,",");
        if($6=="+"){
                print $1,$2,$2,".",$5,$6
        }
        else if($6=="-"){
                print $1,$3,$3,".",$5,$6
        }
}' $PC | sort -k1,1 -k2,2n -k3,3n -k6,6 | uniq > $DIR/pcRNAs_TSS.bed

$CLOSEST -D a -t first -s -a $DIR/pcRNAs_TSS.bed -b $DIR/hg38.cage_peak_coord.bed > $DIR/pcRNAs_TSS_closestCage.bed


Rscript <(cat <<'EOF'
setwd(Sys.getenv("DIR"))
library(ggplot2)
library(reshape2)
library(dplyr)
a <- read.delim("pcRNAs_TSS_closestCage.bed", header=F)
d <- data.frame(Dist=seq(0,1000,10))
for(i in d$Dist){
	d[d$Dist==i, "N"] <- nrow(a[abs(a$V16)<=i,])
}
d <- mutate(d, Freq=N/nrow(a))
pdf("plot.pdf")
print(ggplot(d, aes(x=Dist, y=Freq)) + geom_line() + geom_point(size=1) + ylim(c(0,1)) +xlab("Distance of TSS from closest CAGE tag\n(upstream or downstream)") + ylab("Fraction of pcRNAs") +scale_x_continuous(breaks=seq(0,1000,100), labels=paste(seq(0,1000,100),"bp", sep="")) + theme_bw() +annotate("text", x = 100, y = 1, label = paste("Median: ", median(abs(a$V16)), "bp", sep="")) + annotate("text", x = 100, y = 0.95, label = paste("3rd quartile: ", quantile(abs(a$V16),0.75), "bp", sep="")))
dev.off()
EOF
)


