#!/bin/bash
set -xe -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh


mkdir -p $BASEDIR/hub
mkdir -p $BASEDIR/hub/hg38
mkdir -p $BASEDIR/hub/mm10

{ cat << EOF
hub pcRNAs
shortLabel pcRNAs project $(date +"%m-%d-%y")
longLabel Several datasets relative to the pcRNA project
genomesFile genomes.txt
email tl344@ebi.ac.uk
EOF
}>$BASEDIR/hub/hub.txt

{ cat << EOF
genome hg38
trackDb hg38/trackDb.txt

genome mm10
trackDb mm10/trackDb.txt
EOF
}>$BASEDIR/hub/genomes.txt


fetchChromSizes hg38 > $BASEDIR/hub/hg38.chrom.sizes
fetchChromSizes mm10 > $BASEDIR/hub/mm10.chrom.sizes

bedToBigBed $BASEDIR/transcriptomes/cuffmerge/hsa/merged.bed $BASEDIR/hub/hg38.chrom.sizes $BASEDIR/hub/hg38/hsa_merged_transcripts.bb
bedToBigBed $BASEDIR/transcriptomes/cuffmerge/mmu/merged.bed $BASEDIR/hub/mm10.chrom.sizes $BASEDIR/hub/mm10/mmu_merged_transcripts.bb


#bedToBigBed $BASEDIR/nc/posConNCwithCodingPartner_nr.bed $BASEDIR/hub/hg38.chrom.sizes $BASEDIR/hub/hg38/hsa_pcRNAs.bb
bedToBigBed $BASEDIR/nc/posConNCwithCodingPartner.bed $BASEDIR/hub/hg38.chrom.sizes $BASEDIR/hub/hg38/hsa_pcRNAs_redundant.bb


cp $BASEDIR/transcriptomes/coverage/hsa/*.bw $BASEDIR/hub/hg38
cp $BASEDIR/transcriptomes/coverage/mmu/*.bw $BASEDIR/hub/mm10

bedToBigBed $BASEDIR/nc/posConnedPromoterMouseTS.bed $BASEDIR/hub/mm10.chrom.sizes $BASEDIR/hub/mm10/mmu_pcRNAs.bb
bedToBigBed $BASEDIR/downstream/hic/data/GM12878_looplist.bed $BASEDIR/hub/hg38.chrom.sizes $BASEDIR/hub/hg38/GM12878_looplist.bb

{ cat <<EOF
track mergedTranscripts
bigDataUrl hsa_merged_transcripts.bb
shortLabel Merged transcripts
longLabel All de novo and reference transcript after cuffmerge
type bigBed 12 .
visibility hide
color 0,0,0

track pcRNAsALL
bigDataUrl hsa_pcRNAs_redundant.bb
shortLabel pcRNAs (all)
longLabel All pcRNAs (all isoforms)
type bigBed 12 .
visibility hide
color 255,140,0

track RNA-SeqComposite
compositeTrack on
shortLabel pcRNA RNA-Seq
longLabel pcRNA RNA-Seq data composite track
type bigWig

track GM12878-loops
bigDataUrl GM12878_looplist.bb
shortLabel GM12878-loops
longLabel HiC loops in GM12878
type bigBed 12 .
visibility hide
color 0,0,0

EOF
} > $BASEDIR/hub/hg38/trackDb.txt 

for i in $BASEDIR/transcriptomes/coverage/hsa/*.txt; do cat $i; done >> $BASEDIR/hub/hg38/trackDb.txt


{ cat <<EOF
track mergedTranscripts
bigDataUrl mmu_merged_transcripts.bb
shortLabel Merged transcripts
longLabel All de novo and reference transcript after cuffmerge
type bigBed 12 .
visibility hide
color 0,0,0

track pcRNAs
bigDataUrl mmu_pcRNAs.bb
shortLabel pcRNAs
longLabel All pcRNAs
type bigBed 12 .
visibility hide
color 255,140,0

track RNA-SeqComposite
compositeTrack on
shortLabel pcRNA RNA-Seq
longLabel pcRNA RNA-Seq data composite track
type bigWig

EOF
}> $BASEDIR/hub/mm10/trackDb.txt
for i in $BASEDIR/transcriptomes/coverage/mmu/*.txt; do cat $i; done >> $BASEDIR/hub/mm10/trackDb.txt
