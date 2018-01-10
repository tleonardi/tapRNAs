#!/bin/bash
set -xe -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh


PC=$BASEDIR/../results/pc_annotation.txt

DIR=$BASEDIR/downstream/pairwiseCons
mkdir -p $DIR

n=0
(while read -r line
do
	((n++))
	A=$(echo $line | awk '{print $1}')
	B=$(echo $line | awk '{print $2}')
	echo -ne "$n\t$A\t$B\t"
	$NEEDLE -stdout -auto <(grep -A 1 $A $BASEDIR/../results/posConNCwithCodingPartner.fasta) <(grep -A 1 ${B}_${A} $BASEDIR/../results/posConnedPromoterMouseTS_nr.fasta) | grep '# Identity' |  perl -pe 's/# Identity:\s+([0-9]+)\/([0-9]+)\s\(.+/$1\t$2/'
done < <(cut -f 1,12 $PC | awk '{FS=OFS="\t"}{for(i=1;i<=split($2,A,",");i++){print $1,A[i]}}' | tail -n +2))> $DIR/pairWiseCons.txt


awk 'BEGIN{OFS=FS="\t"}{print $2,$3,$4,$5,$4/$5}' $DIR/pairWiseCons.txt | sort -k1,1 -k5,5nr | sort -k1,1 -u > $DIR/pairWiseCons_highest.txt



# Merge human pc-associated transcripts to their 
# correponsing mouse transcripts
join -1 1 -2 2 -t $'\t' <(cut -f2 $PC | sort -k1,1 | uniq) <(cut -f1,2 GencodeV21_info_tr.txt | sort -k2,2) | sort -k2,2 | \
	join -1 2 -2 1 -t $'\t' - <(cut -f2,3 $BASEDIR/data/hg38ToMm10_orthologs_mart.txt | sort -k1,1) | sort -k3,3 | \
	join -1 3 -2 2 -t $'\t' - <(cut -f1,2 GencodeM4_info_tr.txt | sort -k2,2) > $DIR/human2mouseCoding.txt

# Extract Fasta for human coding
cut -f2 $DIR/human2mouseCoding.txt | sort -u | $PSTR join --idCol 1 --columns 1 -x stdin $BASEDIR/data/GencodeV21.bed | cut -f1-12 > $DIR/humanCoding.bed
bedtools getfasta -split -name -s -fi $BASEDIR/data/hg38.fa -bed $DIR/humanCoding.bed -fo $DIR/humanCoding.fasta

# Extract Fasta for mouse coding
cut -f4 $DIR/human2mouseCoding.txt | sort -u | $PSTR join --idCol 1 --columns 1 -x stdin $BASEDIR/data/GencodeM4.bed | cut -f1-12 > $DIR/mouseCoding.bed
bedtools getfasta -split -name -s -fi $BASEDIR/data/mm10.fa -bed $DIR/mouseCoding.bed -fo $DIR/mouseCoding.fasta


n=0
(while read -r line
do
        ((n++))
        A=$(echo $line | awk '{print $2}')
        B=$(echo $line | awk '{print $4}')
	C=$(echo $line | awk '{print $3}')
        echo -ne "$n\t$C\t$A\t$B\t"
        $NEEDLE -stdout -auto <(grep -A 1 $A $DIR/humanCoding.fasta) <(grep -A 1 ${B} $DIR/mouseCoding.fasta) | grep '# Identity' |  perl -pe 's/# Identity:\s+([0-9]+)\/([0-9]+)\s\(.+/$1\t$2/'
done < $DIR/human2mouseCoding.txt )> $DIR/pairWiseConsCoding.txt


awk 'BEGIN{OFS=FS="\t"}{print $2,$5/$6}' $DIR/pairWiseConsCoding.txt | sort -k1,1 -k2,2nr | sort -k1,1 -u > $DIR/pairWiseConsCoding_highest.txt


Rscript <(cat <<EOF
 library(ggplot2)

 pc <- read.delim("$BASEDIR/../pc_annotation.txt", header=T)
 a <- read.delim("$BASEDIR/downstream/pairwiseCons/pairWiseCons_highest.txt", header=F)
 a <- a[,c(1,5)]
 colnames(a) <- c("V1", "V2")
 b <- read.delim("$BASEDIR/downstream/pairwiseCons/pairWiseConsCoding_highest.txt", header=F)
 a$Type <- "pcRNAs"
 b$Type <- "pcCoding"
 aa <- rbind(a,b)
 ggplot(aa, aes(V2, fill=Type)) + geom_density(colour="black", alpha=0.5)
 ggplot(aa, aes(y=V2, x=Type)) + geom_boxplot()


 pc <- merge(pc, a, by.x="Name", by.y="V1")
 pdf("$BASEDIR/downstream/pairwiseCons/plot.pdf")
 ggplot(pc, aes(V5)) + geom_density() + xlab("Pairwise global conservation")
 dev.off()
EOF)

