#!/bin/bash
set -xe -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh


# Note on memory requirments from UCSC:
# Please note that the wigToBigWig utility uses a lot of memory; somewhere on the order of 1.5 times more memory than the uncompressed wiggle input file. 
# The bedGraphToBigWig utility uses somewhere on the order of 1/4 times more RAM than the uncompressed bedGraph input file

# Human
mkdir -p $BASEDIR/transcriptomes/coverage/hsa
mkdir -p $BASEDIR/transcriptomes/coverage/logs

for i in $BASEDIR/transcriptomes/tophat/hsa/*; do
	SAMP_NAME=$(basename $i)
	if [[ ! -f $BASEDIR/transcriptomes/coverage/hsa/$SAMP_NAME.bedg ]]; then
		bsub -q research-rh6 -J $SAMP_NAME -K -oo $BASEDIR/transcriptomes/coverage/logs/$SAMP_NAME.log -M 5000 -R 'rusage[mem=5000]' "$GENOME_COVERAGE_BED -split -bg -ibam $i/accepted_hits.bam -g $BASEDIR/data/chromInfo_hg38.txt > $BASEDIR/transcriptomes/coverage/hsa/$SAMP_NAME.bedg && $BEDGRAPH_TO_BIGWIG $BASEDIR/transcriptomes/coverage/hsa/$SAMP_NAME.bedg $BASEDIR/data/chromInfo_hg38.txt $BASEDIR/transcriptomes/coverage/hsa/$SAMP_NAME.bw"&
	fi
done;

wait

for i in $BASEDIR/transcriptomes/coverage/hsa/*bw; do
	NAME=$(basename $i .bw);
	TISSUE=$(grep $NAME $BASEDIR/../samples_matrix.txt | cut -f2)
	LONG=$(grep $NAME $BASEDIR/../samples_matrix.txt | cut -f2-7|sed 's/\t/  /g')
	echo track $NAME > $BASEDIR/transcriptomes/coverage/hsa/$NAME.txt
	echo bigDataUrl $NAME.bw >> $BASEDIR/transcriptomes/coverage/hsa/$NAME.txt
	echo shortLabel RNA-Seq ${TISSUE}-${NAME} >> $BASEDIR/transcriptomes/coverage/hsa/$NAME.txt
	echo longLabel $LONG >> $BASEDIR/transcriptomes/coverage/hsa/$NAME.txt
	echo parent RNA-SeqComposite on >> $BASEDIR/transcriptomes/coverage/hsa/$NAME.txt
	echo type bigWig >> $BASEDIR/transcriptomes/coverage/hsa/$NAME.txt
	echo >> $BASEDIR/transcriptomes/coverage/hsa/$NAME.txt
done	


#Mouse 
mkdir -p $BASEDIR/transcriptomes/coverage/mmu
for i in $BASEDIR/transcriptomes/tophat/mmu/*; do
	SAMP_NAME=$(basename $i)
	if [[ ! -f $BASEDIR/transcriptomes/coverage/mmu/$SAMP_NAME.bedg ]]; then
		bsub -q research-rh6 -J $SAMP_NAME -K -oo $BASEDIR/transcriptomes/coverage/logs/$SAMP_NAME.log -M 5000 -R 'rusage[mem=5000]' "$GENOME_COVERAGE_BED -split -bg -ibam $i/accepted_hits.bam -g $BASEDIR/data/chromInfo_mm10.txt > $BASEDIR/transcriptomes/coverage/mmu/$SAMP_NAME.bedg && $BEDGRAPH_TO_BIGWIG $BASEDIR/transcriptomes/coverage/mmu/$SAMP_NAME.bedg $BASEDIR/data/chromInfo_mm10.txt $BASEDIR/transcriptomes/coverage/mmu/$SAMP_NAME.bw"&
	fi
done;

wait

for i in $BASEDIR/transcriptomes/coverage/mmu/*bw; do
	NAME=$(basename $i .bw);
	TISSUE=$(grep $NAME $BASEDIR/../samples_matrix.txt | cut -f2)
	LONG=$(grep $NAME $BASEDIR/../samples_matrix.txt | cut -f2-7|sed 's/\t/  /g')
	echo track $NAME > $BASEDIR/transcriptomes/coverage/mmu/$NAME.txt
	echo bigDataUrl $NAME.bw >> $BASEDIR/transcriptomes/coverage/mmu/$NAME.txt
	echo shortLabel RNA-Seq ${TISSUE}-${NAME} >> $BASEDIR/transcriptomes/coverage/mmu/$NAME.txt
	echo longLabel $LONG >> $BASEDIR/transcriptomes/coverage/mmu/$NAME.txt
	echo parent RNA-SeqComposite on >> $BASEDIR/transcriptomes/coverage/mmu/$NAME.txt
	echo type bigWig >> $BASEDIR/transcriptomes/coverage/mmu/$NAME.txt
	echo >> $BASEDIR/transcriptomes/coverage/mmu/$NAME.txt
	
done;

touch $BASEDIR/.coverage
