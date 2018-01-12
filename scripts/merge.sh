#!/bin/bash
set -e -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh


CUFF_OUT=$BASEDIR/transcriptomes/cufflinks
# Number of CPUs for LSF
P=4

# Create output folders
mkdir -p $BASEDIR/transcriptomes/cuffmerge/logs
mkdir -p $BASEDIR/transcriptomes/cuffmerge/hsa
mkdir -p $BASEDIR/transcriptomes/cuffmerge/mmu

# Create a list of cufflinks-assembled human transcriptomes
for i in $CUFF_OUT/hsa/*; do 
	SAMP=$(basename $i)
	if [[ -f $i/transcripts.gtf ]]; then
		echo $i/transcripts.gtf;
	fi
done > $CUFF_OUT/hsa.list

# Run cuffmerge for human
bsub -q research-rh6 -J hsa_merge -oo $BASEDIR/transcriptomes/cuffmerge/logs/hsa_cuffmerge.log -K -n $P -M 40000 -R 'rusage[mem=40000]'  $CUFFMERGE -p $P -g $BASEDIR/data/GencodeV21_exons-cds.gtf -o $BASEDIR/transcriptomes/cuffmerge/hsa $CUFF_OUT/hsa.list &


# Create a list of cufflinks-assembled mouse transcriptomes
for i in $CUFF_OUT/mmu/*; do 
	if [[ -f $i/transcripts.gtf ]]; then
		echo $i/transcripts.gtf;
	fi
done > $CUFF_OUT/mmu.list

# Run cuffmerge for mouse
bsub -q research-rh6 -J mmu_merge -oo $BASEDIR/transcriptomes/cuffmerge/logs/mmu_cuffmerge.log -K -n $P -M 20000 -R 'rusage[mem=20000]'  $CUFFMERGE -p $P -g $BASEDIR/data/GencodeM4_exons_cds.gtf -o $BASEDIR/transcriptomes/cuffmerge/mmu $CUFF_OUT/mmu.list & 

wait
# Convert cuffmerge GTF outputs to BED
bsub -q research-rh6 -J hsa_gtf2bed -oo $BASEDIR/transcriptomes/cuffmerge/hsa/gtf2bed.log -K -M 10000 -n $P -R 'rusage[mem=10000]' "$PSTR gtfToBed -p $P -a gene_id,oId,class_code $BASEDIR/transcriptomes/cuffmerge/hsa/merged.gtf | sort -k1,1 -k2,2n > $BASEDIR/transcriptomes/cuffmerge/hsa/merged.bedx" &
bsub -q research-rh6 -J mmu_gtf2bed -oo $BASEDIR/transcriptomes/cuffmerge/mmu/gtf2bed.log -K -M 5000 -n $P -R 'rusage[mem=5000]' "$PSTR gtfToBed -p $P -a gene_id,oId,class_code $BASEDIR/transcriptomes/cuffmerge/mmu/merged.gtf | sort -k1,1 -k2,2n > $BASEDIR/transcriptomes/cuffmerge/mmu/merged.bedx" &

wait
cut -f1-12 $BASEDIR/transcriptomes/cuffmerge/hsa/merged.bedx > $BASEDIR/transcriptomes/cuffmerge/hsa/merged.bed
cut -f1-12 $BASEDIR/transcriptomes/cuffmerge/mmu/merged.bedx > $BASEDIR/transcriptomes/cuffmerge/mmu/merged.bed
touch $BASEDIR/.merge
