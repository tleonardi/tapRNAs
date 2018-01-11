#!/bin/bash
set -xe -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh


P=16
TMP=/tmp/$LSB_JOBID

# Prepare tmp folder and download data
mkdir -p $TMP


# Extended blasted promoters 500bp each way
awk 'BEGIN{OFS="\t"}{print $1,$2-500,$3+500,$4,0,$6}' $BASEDIR/nc/hsaNcVsMm10_2_gt100.bed > $BASEDIR/nc/hsaNcVsMm10_2_gt100_ext500.bed

# Select only 5' exon of mm10CombinedNCTranscripts.bed (nc transcripts in combined gencode and novel datasets)
awk 'BEGIN{OFS="\t"}{split($11,SZ,",");if ($6=="+") {print $1,$2,$2+SZ[1],$4,$5,$6} else {print $1,$3-SZ[$10],$3,$4,$5,$6}}' $BASEDIR/nc/mm10CombinedNCTranscripts.bed > $BASEDIR/nc/mm10CombinedNCTranscripts_5p.bed

# Select only 5' exon of all mouse nc, spliced RNAs
awk 'BEGIN{OFS="\t"}{split($11,SZ,",");if ($6=="+") {print $1,$2,$2+SZ[1],$4,$5,$6} else {print $1,$3-SZ[$10],$3,$4,$5,$6}}' $BASEDIR/nc/mm10_allmrna_spliced_NC_uniqID.bed > $BASEDIR/mouse/mm10_allmrna_5p.bed


# Select only extended promoters that overlap with at least one nc transcript
$INTERSECT_BED -wo -a $BASEDIR/nc/hsaNcVsMm10_2_gt100_ext500.bed  -b $BASEDIR/nc/mm10CombinedNCTranscripts_5p.bed | cut -f1-12 > $BASEDIR/nc/hsaNcVsMm10_2_gt100_ext500_NCTranscripts_5p-Overlap.bed

$INTERSECT_BED -v -a $BASEDIR/nc/hsaNcVsMm10_2_gt100_ext500.bed  -b $BASEDIR/nc/mm10CombinedNCTranscripts_5p.bed > $BASEDIR/nc/hsaNcVsMm10_2_gt100_ext500_NCTranscripts_5p-NOOverlap.bed

$INTERSECT_BED -wo -a $BASEDIR/nc/hsaNcVsMm10_2_gt100_ext500_NCTranscripts_5p-NOOverlap.bed -b $BASEDIR/mouse/mm10_allmrna_5p.bed > $BASEDIR/nc/hsaNcVsMm10_2_gt100_ext500_allmrna_5p-Overlap.bed

$PSTR join --rowPerMatch -x -idCol 10 --columns 4 $BASEDIR/nc/hsaNcVsMm10_2_gt100_ext500_NCTranscripts_5p-Overlap.bed $BASEDIR/nc/mm10CombinedNCTranscripts.bed > $BASEDIR/nc/mm10ConsPromWithNCTranscript.bed
$PSTR join --rowPerMatch -x -idCol 10 --columns 4 $BASEDIR/nc/hsaNcVsMm10_2_gt100_ext500_allmrna_5p-Overlap.bed $BASEDIR/nc/mm10_allmrna_spliced_NC_uniqID.bed > $BASEDIR/nc/mm10ConsPromWithAllmrnaNC.bed


# Fields 1:12 mouse ncRNA (ID field 4) with conserved promoter corresponding to human ncRNA with ID in field 13
# In same cases the same ncRNA has overlaps with multiple conserved promoters of the same human ncRNA
# To avoid counting these cases twice we sort -u
cat $BASEDIR/nc/mm10ConsPromWithNCTranscript.bed $BASEDIR/nc/mm10ConsPromWithAllmrnaNC.bed | sort -u | sort -k1,1 -k2,2n -k3,3n > $BASEDIR/nc/mm10ConsPromWithNCTranscript+AllmrnaNC.bed


# Add information on the closest protein coding ID (with human syntenic ID) for the mouse RNA
join -1 4 -2 1 -t $'\t' <(sort -k4,4 $BASEDIR/nc/mm10ConsPromWithNCTranscript+AllmrnaNC.bed) <(sort -k1,1 $BASEDIR/nc/mm10ncRNAToClosestCodingGene_humanSynt.txt) > $BASEDIR/nc/mm10ConsPromWithNCTranscript+AllmrnaNC_mm10Cod2Nc.bedx

# Add information on the human protein coding
# OUTPUT:
# 1: Human nc ID
# 2: Mouse nc ID
# 3-13: Mouse ncRNA BED line
# 14-17: Mouse positional comparison line
# 18-21: Human positional comparison line
join -1 13 -2 1 -t $'\t' <(sort -k13,13 $BASEDIR/nc/mm10ConsPromWithNCTranscript+AllmrnaNC_mm10Cod2Nc.bedx) <(sort -k1,1 $BASEDIR/nc/hg38ncRNAToClosestCodingGene.txt) > $BASEDIR/nc/mm10ConsPromWithNCTranscript+AllmrnaNC_mm10Cod2Nc_hg38Cod2Nc.bedx


# Select only transcripts for which the human and mouse coding match
# and reformat the file:
# 1 Human pcRNA ID
# 2 Mouse pcRNA ID
# 3 Human coding ID
# 4 Human orientation
# 5 Mouse orientation
# 6 Same orientation flag
# 7 Sense or antisense overlap or BIDIR flag
# 8 Distance in human
awk 'function orientation(x){if(x=="DS-S" || x=="US-S" || x=="OLAP"){return("S")}else if(x=="Asense" || x=="BIDIR" || x=="DS-AS" || x=="US-AS"){return("AS")}}
BEGIN{OFS=FS="\t"}$17==$21 && orientation($14)==orientation($18){if($14==$18){SAME=1}else{SAME=0}if($18 == "Asense" || $18 == "OLAP" || $18 == "BIDIR"){OLAP=1}else{OLAP=0}print $1,$2,$17,$18,$14,SAME,OLAP,$20}' $BASEDIR/nc/mm10ConsPromWithNCTranscript+AllmrnaNC_mm10Cod2Nc_hg38Cod2Nc.bedx > $BASEDIR/nc/mm10posConRNAs_humanID_codingID.bedx 

# Annotate Human pcRNA - Human coding gene - Mouse pcRNAs
# For each pcRNA keep only the closest coding
# Choice criteria (in order of precedence)
# Overlap (i.e. Asense and OLAP) -> Distance
cat $BASEDIR/nc/mm10posConRNAs_humanID_codingID.bedx | sort -k1,1 -k3,3 -k4,4 -k5,5 -k6,6n -k7,7n -k8,8n | $GROUP_BY -g 1,3,4,7,8 -c 2,5 -o distinct,distinct | sort -k1,1 -k4,4nr -k5,5n | sort -k1,1 -u  > $BASEDIR/nc/posCon_human2Coding2Mouse.txt



###Get set of transcripts in human coords
# Join (reporting only matching rows) the annotation of posCons with the nc transcriptome
$PSTR join -idCol 1 -x --columns 2 $BASEDIR/nc/posCon_human2Coding2Mouse.txt $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap.bed  > $BASEDIR/nc/posConNCwithCodingPartner.bedx
cut -f1-12 $BASEDIR/nc/posConNCwithCodingPartner.bedx > $BASEDIR/nc/posConNCwithCodingPartner.bed





# Compile annotation of mouse pcRNAs
cat $BASEDIR/nc/posCon_human2Coding2Mouse.txt | awk 'BEGIN{OFS=FS="\t"}{N=split($6,A,","); for(i=1; i<=N; i++){print $1,A[i]}}' | sort -u > $BASEDIR/nc/posConnedPromoterMouseTS_list.txt
$PSTR join -x --rowPerMatch -idCol 2 --columns 1 $BASEDIR/nc/posConnedPromoterMouseTS_list.txt $BASEDIR/nc/mm10CombinedNCTranscripts+allmrna_spliced_NC_uniqID.bed | awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$13,$5,$6,$7,$8,$9,$10,$11,$12,$4}' > $BASEDIR/nc/posConnedPromoterMouseTS.bedx

$PSTR findBest --ncFraction 0.5 -d $BASEDIR/nc/posConnedPromoterMouseTS.bedx | sort -k4,4 > $BASEDIR/nc/posConnedPromoterMouseTS_nr.bedx
# Make BED file
cut -f1-12 $BASEDIR/nc/posConnedPromoterMouseTS_nr.bedx | sort -k1,1 -k2,2n > $BASEDIR/nc/posConnedPromoterMouseTS.bed
$GROUP_BY -g 4 -c 13 -o distinct -i $BASEDIR/nc/posConnedPromoterMouseTS_nr.bedx > $BASEDIR/nc/posConnedPromoterMouseTS_nr.txt


# Remove redundancy from human posCons and conver to to standard bed
$PSTR  findBest --ncFraction 0.5 $BASEDIR/nc/posConNCwithCodingPartner.bedx | sort -k1,1 -k2,2n > $BASEDIR/nc/posConNCwithCodingPartner_nr.bedx
cut -f1-12 $BASEDIR/nc/posConNCwithCodingPartner_nr.bedx > $BASEDIR/nc/posConNCwithCodingPartner_nr.bed

rm -rf $TMP
md5sum $BASEDIR/nc/posConNCwithCodingPartner.bedx $BASEDIR/nc/posConNCwithCodingPartner_nr.bedx > $BASEDIR/.identify_find

