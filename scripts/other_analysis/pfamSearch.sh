#!/bin/bash
set -xe -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../include.sh


OUT=$BASEDIR/downstream/blastx
mkdir -p $OUT

export PERL5LIB=/path/to/pfam_scanDir:$PERL5LIB

wget -O - ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz | gunzip > $OUT/Pfam-A.hmm
wget -O - ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz | gunzip > $OUT/Pfam-A.hmm.dat
wget -O - ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz | gunzip > $OUT/active_site.dat

$hmmpress $OUT/Pfam-A.hmm

$GETFASTA -name -split -s -fi $BASEDIR/data/hg38.fa -bed $BASEDIR/nc/posConNCwithCodingPartner.bed | perl -pe 's/::chr.+//' > $OUT/pcRNAs.fa

$transeq -sequence $OUT/pcRNAs.fa -outseq $OUT/pcRNAs_tx.fa -frame F 

$pfam_scan -fasta $OUT/pcRNAs_tx.fa -dir $OUT -cpu 8 > $OUT/pfmaScan.out

awk '!/#|^$/' $OUT/pfmaScan.out | perl -pe 's/ +/\t/g' > $OUT/pfmaScan.clean.out

cut -f1 $OUT/pfmaScan.clean.out | perl -pe 's/_[0-9]+$//' | sort -u > $OUT/pfmaScan_hits.txt

join -1 1 -2 1 -t $'\t' <(cut -f4,13,26 $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx | sort -k1,1) <(sort -k1,1 $OUT/pfmaScan_hits_manualAnnot.txt) 
