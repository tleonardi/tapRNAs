#!/bin/bash
set -xe -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh

TMP=/tmp/$LSB_JOBID

# Prepare tmp folder
mkdir -p $TMP


#how many are lincRNA
$INTERSECT_BED -u -s -split -a $BASEDIR/nc/posConNCwithCodingPartner.bed -b $BASEDIR/data/Gencode.v21.lincRNAs.bed | awk 'BEGIN{OFS="\t"}{print $4,1}' > $BASEDIR/nc/posConLinc.txt


#Get posCon pair
cut -f4,13 $BASEDIR/nc/posConNCwithCodingPartner.bedx > $BASEDIR/nc/posConIDPairs.txt


#Get set of coding genes with posCon RNA and their common names
join -1 13 -2 1 -t $'\t' <(sort -k13,13 $BASEDIR/nc/posConNCwithCodingPartner.bedx) <(cut -f2,9 $BASEDIR/data/GencodeV21_info_tr.txt | sort -k1,1 -k2,2 -u) | awk 'BEGIN{OFS=FS="\t"}{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$1,$14}' > $BASEDIR/nc/posConNCwithCodingPartnerUnigene.bedx

#get miRNA overlapped by nc for later annotation of nc/coding set
echo "select transcriptId,geneName from wgEncodeGencodeAttrsV20 where geneName like '%miR%'" | mysql hg38 --user=genome --host=genome-mysql.cse.ucsc.edu | awk 'BEGIN{OFS=FS="\t"}{sub("\\.[0-9]+$", "", $1); print $1,$2}' > $BASEDIR/nc/miRID.txt 
$PSTR join -x $BASEDIR/nc/miRID.txt $BASEDIR/data/Gencode-current_hsa.bed | awk ' BEGIN{OFS="\t"}{print $1,$2,$3,$13,$5,$6,$7,$8,$9,$10,$11,$12}' > $BASEDIR/nc/miRID.bed
$INTERSECT_BED -s -split -wb -a $BASEDIR/nc/posConNCwithCodingPartner.bed -b $BASEDIR/nc/miRID.bed | cut -f4,16 | sort -k1,1 | $GROUP_BY -g 1 -o distinct -c 2 > $BASEDIR/nc/miRInPosCon.txt

# Save pcRNA context
mkdir -p $BASEDIR/nc/context
cut -f 1,3,5 $BASEDIR/nc/posCon_human2Coding2Mouse.txt > $BASEDIR/nc/context/pcRNAContexts.txt



### Get list of transcripts in Gencode annotations (20% similarity)
$INTERSECT_BED -split -s -u -f 0.2 -a $BASEDIR/nc/posConNCwithCodingPartner.bed -b $BASEDIR/data/Gencode-current_hsa.bed | awk 'BEGIN{OFS="\t"}{print $4,0}' > $BASEDIR/nc/pcRNA_NovTSGenc_20pc.txt



#########################################################
# PREPARE PC PROMOTERS
#########################################################
# Select promoters of each pcRNA
awk 'BEGIN{OFS="\t"}{if ($6=="+") {print $1,$2-2000,$2+2000,$4,$5,$6} else {print $1,$3-2000,$3+2000,$4,$5,$6}}' $BASEDIR/nc/posConNCwithCodingPartner.bedx | sort -k1,1 -k2,2n > $BASEDIR/nc/posConNC_promoters.bed

# Merge the overlapping (or book ended) promoters
$MERGE -c 4,6 -o distinct,distinct -i $BASEDIR/nc/posConNC_promoters.bed | awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,"Prom_"NR,$4,$5}' > $BASEDIR/nc/posConNC_promoters_clusters.bedx
awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$4,"0",$6}' $BASEDIR/nc/posConNC_promoters_clusters.bedx > $BASEDIR/nc/posConNC_promoters_clusters.bed
awk 'BEGIN{OFS=FS="\t"}{N=split($5,PC,","); for(i=1;i<=N;i++){print PC[i],$4}}' $BASEDIR/nc/posConNC_promoters_clusters.bedx > $BASEDIR/nc/posConNC_promoters_clusters.txt

awk 'BEGIN{OFS=FS="\t"}{print $1,$2,"-", "-"}' $BASEDIR/nc/posConNC_promoters_clusters.txt > $BASEDIR/nc/posConProm2CodProm.txt


#########################################################
### Compile into final annotation
#########################################################
$PSTR join -empty 0 -columns 2,3 $BASEDIR/nc/context/pcRNAContexts.txt $BASEDIR/nc/posConNCwithCodingPartnerUnigene.bedx | $PSTR join -empty 0 $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodingGeneOlap.txt stdin | $PSTR join -empty 1 $BASEDIR/nc/pcRNA_NovTSGenc_20pc.txt stdin | $PSTR join -empty 0 $BASEDIR/nc/posConLinc.txt stdin | $PSTR join -empty "_" $BASEDIR/nc/miRInPosCon.txt stdin | $PSTR join -idCol 4 -columns 13,14,15,16 $BASEDIR/nc/ncRNAcmb_mouseCons.bedx stdin | $PSTR join -columns 2 $BASEDIR/nc/posConnedPromoterMouseTS_nr.txt stdin | $PSTR join --idCol 1 --columns 2,4 $BASEDIR/nc/posConProm2CodProm.txt stdin > $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx


# Make a table that associate human genes to transcripts (used by prepare_data.R)
join -1 1 -2 1 -a 1 -t $'\t' <(cut -f4 $BASEDIR/coding/Gencode-current_hsa_coding.bed | sort -k1,1) <(cut -f1-2 $BASEDIR/data/GencodeV21_info_tr.txt| sort -k1,1) > $BASEDIR/coding/GencodeV21_coding_trID2geneID.txt


md5sum $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx > $BASEDIR/.annotate
