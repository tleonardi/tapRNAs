#!/bin/bash
set -e -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh



P=8

OUT=$BASEDIR/transcriptomes/cufflinks
mkdir -p $OUT/hsa
mkdir -p $OUT/logs/



# Make a GTF with only CDS and exon features
	if [[ ! -f "$BASEDIR/data/GencodeV21_exons-cds.gtf" ]]; then
		awk 'BEGIN{OFS=FS="\t"} !/^#/{if($3 == "exon" || $3 == "CDS"){print}}' $BASEDIR/data/GencodeV21.gtf > $BASEDIR/data/GencodeV21_exons-cds.gtf
	fi
	if [[ ! -f "$BASEDIR/data/GencodeM4_exons_cds.gtf" ]]; then
		awk 'BEGIN{OFS=FS="\t"} !/^#/{if($3 == "exon" || $3 == "CDS"){print}}' $BASEDIR/data/GencodeM4.gtf > $BASEDIR/data/GencodeM4_exons_cds.gtf
	fi

# Create annotation of rRNA and mtRNA
	if [[ ! -f "$BASEDIR/data/GencodeV21_mask.gtf" ]]; then
		awk 'NR>5{if($14 ~ /Mt_/ || $14 ~ "rRNA"){print}}' $BASEDIR/data/GencodeV21_exons-cds.gtf > $BASEDIR/data/GencodeV21_mask.gtf
	fi
	if [[ ! -f "$BASEDIR/data/GencodeM4_mask.gtf" ]]; then
		awk 'NR>5{if($14 ~ /Mt_/ || $14 ~ "rRNA"){print}}' $BASEDIR/data/GencodeM4_exons_cds.gtf > $BASEDIR/data/GencodeM4_mask.gtf
	fi

# Create fasta indexes
	if [[ ! -f "$BASEDIR/data/hg38.fa.fai" ]]; then
		$SAMTOOLS faidx $BASEDIR/data/hg38.fa

	fi
	if [[ ! -f "$BASEDIR/data/mm10.fa.fai" ]]; then
		$SAMTOOLS faidx $BASEDIR/data/mm10.fa
	fi

# Define stranded datasets and datasets the require high amount of memory
STRANDED_MOUSE="SRR496249$|SRR496250$|SRR496221$|SRR496222$|SRR530639$|SRR530640$"
STRANDED_HUMAN="SRR307919$|SRR317040$|SRR307925$|SRR317039$|SRR307911$|SRR307912$|SRR307923$|SRR307924$"
STRANDED=$(echo "$STRANDED_MOUSE|$STRANDED_HUMAN")

# Human
for i in $BASEDIR/transcriptomes/tophat/hsa/*; do
	SAMP_NAME=$(basename $i)
	if [[ ! -d $OUT/hsa/$SAMP_NAME ]] && [[ -f $i/accepted_hits.bam ]]; then
        	if [[ $SAMP_NAME =~ $STRANDED ]]; then
                	LIBTYPE="fr-firststrand";
        	else
                	LIBTYPE="fr-unstranded";
        	fi
		BAM_SIZE=$(du -b "$i/accepted_hits.bam" | cut -f 1)
		if [[ $BAM_SIZE -ge 2147483648 ]]; then
			MEM=50000;
		else
			MEM=30000;
		fi
		bsub -q research-rh6 -J c$SAMP_NAME -oo $OUT/logs/$SAMP_NAME.log -M $MEM -n $P -R "'"rusage[mem=$MEM]"'"  $CUFFLINKS -p $P --library-type $LIBTYPE -F 0.05 --multi-read-correct --frag-bias-correct $BASEDIR/data/hg38.fa -M $BASEDIR/data/GencodeV21_mask.gtf -g $BASEDIR/data/GencodeV21_exons-cds.gtf -o $OUT/hsa/$SAMP_NAME $i/accepted_hits.bam
	fi
done;


mkdir -p $OUT/mmu
mkdir -p $OUT/logs/

# Mouse
for i in $BASEDIR/transcriptomes/tophat/mmu/*; do
	SAMP_NAME=$(basename $i)
	if [[ ! -d $OUT/mmu/$SAMP_NAME ]] && [[ -f $i/accepted_hits.bam ]]; then
       		if [[ $SAMP_NAME =~ $STRANDED ]]; then
       			LIBTYPE="fr-firststrand";
       		else
       			LIBTYPE="fr-unstranded";
       		fi		
		BAM_SIZE=$(du -b "$i/accepted_hits.bam" | cut -f 1)
		if [[ $BAM_SIZE -ge 2147483648 ]]; then
			MEM=50000;
		else
			MEM=30000;
		fi
		bsub -q research-rh6 -J c$SAMP_NAME -oo $OUT/logs/$SAMP_NAME.log -M $MEM -n $P -R "'"rusage[mem=$MEM]"'"  $CUFFLINKS -p $P --library-type $LIBTYPE -F 0.05 --multi-read-correct --frag-bias-correct $BASEDIR/data/mm10.fa -M $BASEDIR/data/GencodeM4_mask.gtf -g $BASEDIR/data/GencodeM4_exons_cds.gtf -o $OUT/mmu/$SAMP_NAME $i/accepted_hits.bam
	fi
done;

touch $BASEDIR/.assemble
