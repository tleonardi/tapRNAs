#!/bin/bash

set -xe -o pipefail
export LC_ALL=C
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh


PC=$BASEDIR/../results/pc_annotation.txt

DIR=$BASEDIR/downstream/promoterCons
mkdir -p $DIR





$PSTR join -idCol 1 --columns 2 $BASEDIR/data/GencodeV21_info_tr.txt $BASEDIR/../results/posConCodingTranscripts.bed | awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$13,$5,$6,$7,$8,$9,$10,$11,$12}' > $DIR/posConCodingTranscripts_GeneID.bed


# Non coding and Coding promoters bins:
# Prom: +/- 500
# Up: from -1500 to -500  
# Down: from 500 to 1500
awk 'BEGIN{OFS="\t"}{if ($6=="+") {print $1,$2-500,$2+500,$4,$5,$6} else {print $1,$3-500,$3+500,$4,$5,$6}}'   $DIR/posConCodingTranscripts_GeneID.bed | uniq > $DIR/posConCodingTranscripts_GeneID_promoters.bed
awk 'BEGIN{OFS="\t"}{if ($6=="+") {print $1,$2-1500,$2-500,$4,$5,$6} else {print $1,$3+500,$3+1500,$4,$5,$6}}' $DIR/posConCodingTranscripts_GeneID.bed | uniq > $DIR/posConCodingTranscripts_GeneID_promoters_up.bed
awk 'BEGIN{OFS="\t"}{if ($6=="+") {print $1,$2+500,$2+1500,$4,$5,$6} else {print $1,$3-1500,$3-500,$4,$5,$6}}' $DIR/posConCodingTranscripts_GeneID.bed | uniq > $DIR/posConCodingTranscripts_GeneID_promoters_down.bed

awk 'BEGIN{OFS="\t"}{if ($6=="+") {print $1,$2-500,$2+500,$4,$5,$6} else {print $1,$3-500,$3+500,$4,$5,$6}}'   $BASEDIR/../results/posConNCwithCodingPartner.bed > $DIR/posConNCwithCodingPartner_promoters.bed
awk 'BEGIN{OFS="\t"}{if ($6=="+") {print $1,$2-1500,$2-500,$4,$5,$6} else {print $1,$3+500,$3+1500,$4,$5,$6}}' $BASEDIR/../results/posConNCwithCodingPartner.bed > $DIR/posConNCwithCodingPartner_promoters_up.bed
awk 'BEGIN{OFS="\t"}{if ($6=="+") {print $1,$2+500,$2+1500,$4,$5,$6} else {print $1,$3-1500,$3-500,$4,$5,$6}}' $BASEDIR/../results/posConNCwithCodingPartner.bed > $DIR/posConNCwithCodingPartner_promoters_down.bed


# Get FASTA
bedtools getfasta -name -s -fi $BASEDIR/data/hg38.fa -bed $DIR/posConCodingTranscripts_GeneID_promoters.bed -fo $DIR/posConCodingTranscripts_GeneID_promoters.fa
bedtools getfasta -name -s -fi $BASEDIR/data/hg38.fa -bed $DIR/posConCodingTranscripts_GeneID_promoters_up.bed -fo $DIR/posConCodingTranscripts_GeneID_promoters_up.fa
bedtools getfasta -name -s -fi $BASEDIR/data/hg38.fa -bed $DIR/posConCodingTranscripts_GeneID_promoters_down.bed -fo $DIR/posConCodingTranscripts_GeneID_promoters_down.fa

bedtools getfasta -name -s -fi $BASEDIR/data/hg38.fa -bed  $DIR/posConNCwithCodingPartner_promoters.bed -fo $DIR/posConNCwithCodingPartner_promoters.fa
bedtools getfasta -name -s -fi $BASEDIR/data/hg38.fa -bed  $DIR/posConNCwithCodingPartner_promoters_up.bed -fo $DIR/posConNCwithCodingPartner_promoters_up.fa
bedtools getfasta -name -s -fi $BASEDIR/data/hg38.fa -bed  $DIR/posConNCwithCodingPartner_promoters_down.bed -fo $DIR/posConNCwithCodingPartner_promoters_down.fa


# For each pair of PC->coding gene
n=0
(while read -r line
do
        ((n++))
	# A: PC ID
        A=$(echo $line | awk '{print $1}')
	# B: Coding gene ID
        B=$(echo $line | awk '{print $2}')

	# Grep the sequence of the pc
	SEQA=$(grep -A1 $A $DIR/posConNCwithCodingPartner_promoters.fa | tail -1)

	# Loop through each of the sequences (i.e. transcripts) of the coding gene
	# nb stores the index of the coding transcript so that we can pull out the same
	# from _up and _down
	nb=0
	while read -r SEQB
	do
		((nb++))
		# Needle promoter
		NEED_PROM=$($NEEDLE -stdout -auto <(echo -e ">SEQA\n$SEQA") <(echo -e ">SEQB\n$SEQB") | grep '# Identity' |  perl -pe 's/# Identity:\s+([0-9]+)\/([0-9]+)\s\((.+)%\)\n/\t$3/')

		SEQB_RAND=$(grep -v ">ENSG" $DIR/posConCodingTranscripts_GeneID_promoters.fa | shuf | awk 'NR<11{print ">RAND"NR"\n"$0}')
		NEED_RAND=($($NEEDLE -stdout -auto <(echo -e ">SEQA\n$SEQA") <(echo -e "$SEQB_RAND") | grep '# Identity' |  perl -pe 's/# Identity:\s+([0-9]+)\/([0-9]+)\s\((.+)%\)\n/\t$3/'))
		
		# Extract _up and _down in position nb
		SEQA_up=$(grep -A 1 $A $DIR/posConNCwithCodingPartner_promoters_up.fa | grep -v ">")
		SEQB_up=$(grep -A 1 $B $DIR/posConCodingTranscripts_GeneID_promoters_up.fa | grep -v ">" | awk -v nb=$nb 'NR==nb')
		NEED_UP=$($NEEDLE -stdout -auto <(echo -e ">SEQA\n$SEQA_up") <(echo -e ">SEQB\n$SEQB_up") | grep '# Identity' |  perl -pe 's/# Identity:\s+([0-9]+)\/([0-9]+)\s\((.+)%\)\n/\t$3/')

		SEQB_UP_RAND=$(grep -v ">ENSG" $DIR/posConCodingTranscripts_GeneID_promoters_up.fa | shuf | awk 'NR<11{print ">RAND"NR"\n"$0}')
	        NEED_UP_RAND=($($NEEDLE -stdout -auto <(echo -e ">SEQA\n$SEQA") <(echo -e "$SEQB_UP_RAND") | grep '# Identity' |  perl -pe 's/# Identity:\s+([0-9]+)\/([0-9]+)\s\((.+)%\)\n/\t$3/'))

		SEQA_down=$(grep -A 1 $A $DIR/posConNCwithCodingPartner_promoters_down.fa | grep -v ">")
		SEQB_down=$(grep -A 1 $B $DIR/posConCodingTranscripts_GeneID_promoters_down.fa | grep -v ">" | awk -v nb=$nb 'NR==nb')
                NEED_DOWN=$($NEEDLE -stdout -auto <(echo -e ">SEQA\n$SEQA_down") <(echo -e ">SEQB\n$SEQB_down") | grep '# Identity' |  perl -pe 's/# Identity:\s+([0-9]+)\/([0-9]+)\s\((.+)%\)\n/\t$3/')

		SEQB_DOWN_RAND=$(grep -v ">ENSG" $DIR/posConCodingTranscripts_GeneID_promoters_down.fa | shuf | awk 'NR<11{print ">RAND"NR"\n"$0}')
                NEED_DOWN_RAND=($($NEEDLE -stdout -auto <(echo -e ">SEQA\n$SEQA") <(echo -e "$SEQB_DOWN_RAND") | grep '# Identity' |  perl -pe 's/# Identity:\s+([0-9]+)\/([0-9]+)\s\((.+)%\)\n/\t$3/'))


		echo -e "$A\t$B$NEED_PROM$NEED_UP$NEED_DOWN\tReal"
		for i in $(seq 0 9); do echo -e "$A\t$B\t${NEED_RAND[i]}\t${NEED_UP_RAND[i]}\t${NEED_DOWN_RAND[i]}\tRandom"; done
	done < <(grep -A1 $B $DIR/posConCodingTranscripts_GeneID_promoters.fa | grep -v ">")
echo $n >&2
done < <(cut -f 1,2 $PC | tail -n +2)
)> $DIR/pairWiseCons.txt

awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,($3+$4+$5)/3}' $DIR/pairWiseCons.txt | sort -k1,1 -k2,2 -k7,7nr | sort -k1,1 -k6,6 -u > $DIR/pairWiseCons_best.txt


Rscript <(cat <<EOF
library(ggplot2)
library(reshape2)
library(dplyr)
a <- read.delim("$DIR/pairWiseCons_best.txt", header=F)
colnames(a) <- c("Tx", "Gene", "Prom", "Up", "Down", "Type", "Tot")
a <- a[,-7]
a<- melt(a) 

b<-read.delim("analysis/nc/posConNCwithCodingPartner_andContext_expr.bedx", header=F) %>% select(V4, V15)

a <- merge(a, b, by.x="Tx", by.y="V4")

pdf("$DIR/plot.pdf")
ggplot(a, aes(x=variable, y=value)) + geom_boxplot()
dev.off()
EOF)

