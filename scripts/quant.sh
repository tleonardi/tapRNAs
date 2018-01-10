#!/bin/bash
set -e -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh


P=4

OUT=$BASEDIR/transcriptomes/cuffquant
mkdir -p $OUT/hsa
mkdir -p $OUT/logs/


# Define stranded datasets 
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
		bsub -q research-rh6 -J $SAMP_NAME -oo $OUT/logs/$SAMP_NAME.log -M 20000 -n $P -R "'"rusage[mem=20000]"'"  $CUFFQUANT -p $P --library-type $LIBTYPE --multi-read-correct --frag-bias-correct $BASEDIR/data/hg38.fa -M $BASEDIR/data/GencodeV21_mask.gtf -o $OUT/hsa/$SAMP_NAME $BASEDIR/transcriptomes/cuffmerge/hsa/merged.gtf $i/accepted_hits.bam
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
		bsub -q research-rh6 -J $SAMP_NAME -oo $OUT/logs/$SAMP_NAME.log -M 20000 -n $P -R "'"rusage[mem=20000]"'"  $CUFFQUANT -p $P --library-type $LIBTYPE --multi-read-correct --frag-bias-correct $BASEDIR/data/mm10.fa -M $BASEDIR/data/GencodeM4_mask.gtf -o $OUT/mmu/$SAMP_NAME $BASEDIR/transcriptomes/cuffmerge/mmu/merged.gtf $i/accepted_hits.bam
	fi
done;

touch $BASEDIR/.quant
