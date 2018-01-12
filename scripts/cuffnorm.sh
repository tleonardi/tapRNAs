#!/bin/bash
set -e -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh

# Number of CPUs for LSF
P=16

OUT=$BASEDIR/transcriptomes/cuffnorm
mkdir -p $OUT/hsa
mkdir -p $OUT/mmu
mkdir -p $OUT/logs/

# Create sample file for Human
echo -e "sample_id\tgroup_label" > $OUT/hsa/samples.txt
(for i in $BASEDIR/transcriptomes/cuffquant/hsa/*; do
	SAMP_NAME=$(basename $i)
	if [[ -f $i/abundances.cxb ]]; then
		SAMP_LABEL=$(grep  $SAMP_NAME $BASEDIR/../samples_matrix.txt | cut -f2)
		echo -e "$i/abundances.cxb\t$SAMP_LABEL"
	fi
done)>>$OUT/hsa/samples.txt;


# Create sample file for Mouse
echo -e "sample_id\tgroup_label" > $OUT/mmu/samples.txt;
(for i in $BASEDIR/transcriptomes/cuffquant/mmu/*; do
        SAMP_NAME=$(basename $i)
        if [[ -f $i/abundances.cxb ]]; then
                SAMP_LABEL=$(grep  $SAMP_NAME $BASEDIR/../samples_matrix.txt | cut -f2)
                echo -e "$i/abundances.cxb\t$SAMP_LABEL"
        fi
done)>>$OUT/mmu/samples.txt;


bsub -q research-rh6 -J hsa_norm -oo $OUT/logs/hsa_cuffnorm.log -M 80000 -n $P -R "'"rusage[mem=80000]"'"  $CUFFNORM --output-format cuffdiff --use-sample-sheet -p $P --library-type fr-unstranded -o $OUT/hsa $BASEDIR/transcriptomes/cuffmerge/hsa/merged.gtf $OUT/hsa/samples.txt
bsub -q research-rh6 -J mmu_norm -oo $OUT/logs/mmu_cuffnorm.log -M 80000 -n $P -R "'"rusage[mem=80000]"'"  $CUFFNORM --output-format cuffdiff --use-sample-sheet -p $P --library-type fr-unstranded -o $OUT/mmu $BASEDIR/transcriptomes/cuffmerge/mmu/merged.gtf $OUT/mmu/samples.txt


touch $BASEDIR/.cuffnorm
