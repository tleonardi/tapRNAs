#!/bin/bash
set -xe -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh



if [[ -z $RESULTS ]]; then
	echo "RESULTS folder is not set. Are you running the script through the Makefile?";
	exit 1;
fi

if [ -d "$RESULTS" ]; then
	echo "The RESULTS directory already exists. Delete it or rename it.";
	exit 1;
fi


mkdir $RESULTS

# All mm10 nc transcripts
#	Novel TS with no ORF, more than 1 exon and no overlap to protein coding exons
#	plus gencode RNAs with no annoted ORF and more than one exon
cp $BASEDIR/nc/mm10CombinedNCTranscripts.bed $RESULTS

# All hg38 nc transcripts
#	Gencode spliced, no annotated ORF, no overlap with coding exons
#	plus all nover transcripts with no ORF, spliced and no overlap with gencode coding exons
cp $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap.bed $RESULTS/GencodeV21PlusNovel_ncSplice_noCodeOlap.bed


# BED of human pcRNAs
cp $BASEDIR/nc/posConNCwithCodingPartner.bed $RESULTS

# BED of mouse pcRNAs
cp $BASEDIR/nc/posConnedPromoterMouseTS.bed $RESULTS

# BED of human pcRNA-associated coding transcripts
cp $BASEDIR/downstream/hic/coding.bed $RESULTS/posConCodingTranscripts.bed

# Final list of pcRNAs (red and nr)
cp $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx $RESULTS
cp $BASEDIR/../R/report/pc_annotation.txt $RESULTS


# Human and mouse de novo transcriptomes
cp $BASEDIR/transcriptomes/cuffmerge/hsa/merged.bed $RESULTS/hsa_merged.bed
cp $BASEDIR/transcriptomes/cuffmerge/mmu/merged.bed $RESULTS/mmu_merged.bed

# Make FASTA file of human and mouse (redundant) pcRNAs
$GETFASTA -name -split -s -fi $BASEDIR/data/hg38.fa -bed $BASEDIR/nc/posConNCwithCodingPartner.bed -fo $RESULTS/posConNCwithCodingPartner.fasta
$GETFASTA -name -split -s -fi $BASEDIR/data/mm10.fa -bed <(awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$13"_"$4,$5,$6,$7,$8,$9,$10,$11,$12}' $BASEDIR/nc/posConnedPromoterMouseTS_nr.bedx) -fo $RESULTS/posConnedPromoterMouseTS_nr.fasta

# Make FASTA file of human codings
# $GETFASTA -name -split -s -fi $BASEDIR/data/hg38.fa -bed $BASEDIR/coding/PConCodingPartners.bed  -fo $RESULTS/PConCodingPartners.fasta

head -1 $BASEDIR/../R/report/pc_annotation.txt | awk 'BEGIN{OFS=FS="\t"}{print $0,"Nanostring"}'> $RESULTS/pc_annotation.csv
join -1 1 -2 4 -t $'\t' <(cat $BASEDIR/../R/report/pc_annotation.txt | tail -n +2 | sort -k1,1) <(cat $BASEDIR/nc/posConNCwithCodingPartner.bed | sort -k4,4) | awk 'BEGIN{OFS=FS="\t"}{print "<a href=http://genome-euro.ucsc.edu/cgi-bin/hgTracks?org=human&db=hg38&position="$21":"$22-2500"-"$23+2500" target=\"_blank\">"$1"</a>",$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20} ' | sort -k3,3 -t $'\t' | join -1 3 -2 1 -t $'\t' -a 1 - <(cut -f4,5 $DATA/nanostring/nano_annot.txt | sort -u | sort -k1,1) | awk 'BEGIN{OFS=FS="\t"}{if($21==""){$21="NotAvail"}else{$21="<a href=nanostring.html#"tolower($21)" target=\"_blank\">Exp</a>/<a href=nanostring.html#"tolower($21)"-1 target=\"_blank\">Correlation</a>"} if($9=="_"){$9="0"}else{$9="1"}print $2,$3,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21}' >> $RESULTS/pc_annotation.csv

# Save some HiC data files
mkdir -p $RESULTS/HiC
cp $BASEDIR/downstream/hic/endPoints/pcPeak_endPoint.bedx $RESULTS/HiC
touch $BASEDIR/.save

