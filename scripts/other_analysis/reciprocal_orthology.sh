#!/bin/bash
set -xe -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh

# Number of CPUs
P=16

TMP=/tmp/$LSB_JOBID

# Prepare tmp folder and download data
mkdir -p $TMP


mkdir -p $BASEDIR/downstream/reciprocal_orthology

cut -f13 $BASEDIR/nc/posConnedPromoterMouseTS_nr.bedx | sort -u > $BASEDIR/downstream/reciprocal_orthology/mouse_pcRNAs.txt



(
$PSTR join --columns 1 --idCol 1 -x $BASEDIR/downstream/reciprocal_orthology/mouse_pcRNAs.txt $BASEDIR/nc/mm10CombinedNCTranscripts.bed | awk 'BEGIN{OFS="\t"}{if ($6=="+") {print $1,$2-500,$2,$4,$5,$6} else {print $1,$3,$3+500,$4,$5,$6}}'
$PSTR join --columns 1 --idCol 1 -x $BASEDIR/downstream/reciprocal_orthology/mouse_pcRNAs.txt $BASEDIR/nc/mm10_allmrna_spliced_NC_uniqID.bed | awk 'BEGIN{OFS="\t"}{if ($6=="+") {print $1,$2-500,$2,$4,$5,$6} else {print $1,$3,$3+500,$4,$5,$6}}'
)>$BASEDIR/downstream/reciprocal_orthology/mouse_pcRNA_promoters.bed
$PSTR getDna $BASEDIR/data/mm10.fa $BASEDIR/downstream/reciprocal_orthology/mouse_pcRNA_promoters.bed > $BASEDIR/downstream/reciprocal_orthology/mouse_pcRNA_promoters.fa


# Make a mask file based on the soft-masking of the genome
# and then create a BLAST DB
$CONVERT2BLASTMASK -in $BASEDIR/data/hg38.fa -masking_algorithm repeat -masking_options "repeatmasker and tandem repeats from UCSC" -outfmt maskinfo_asn1_bin -out $BASEDIR/data/hg38_mask.asnb 
$MAKEBLASTDB -in $BASEDIR/data/hg38.fa -mask_data $BASEDIR/data/hg38_mask.asnb -dbtype nucl

### Blast mouse promoter regions against hg38
# outfmt 6 indicates the output format. The columns are as follows:
# 	qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
$BLASTN -task blastn -db $BASEDIR/data/hg38.fa -out $BASEDIR/downstream/reciprocal_orthology/MmuNcVsHg38.blast  -query $BASEDIR/downstream/reciprocal_orthology/mouse_pcRNA_promoters.fa -outfmt 6 -evalue 0.001 -num_threads $P -db_soft_mask 40 -lcase_masking

### Clean and filter the results of blastn of ncRNA promoter from mouse to human
# we keep only promoter with match >100nt and Eval<1e-10. We then print in BED format the human coords. Convert blast (1 based) results back into bed format (0 based)
# We keep multimapping promoters
awk 'BEGIN{OFS="\t"}$4>=100 && $11<1E-10 {if ($10>$9) {print $2,$9-1,$10,$1,0,"."} else {print $2,$10-1,$9,$1,0,"."}}'  $BASEDIR/downstream/reciprocal_orthology/MmuNcVsHg38.blast > $BASEDIR/downstream/reciprocal_orthology/MmuNcVsHg38_gt100.bed

# Extended blasted promoters 500bp each way
awk 'BEGIN{OFS="\t"}{print $1,$2-500,$3+500,$4,0,$6}' $BASEDIR/downstream/reciprocal_orthology/MmuNcVsHg38_gt100.bed > $BASEDIR/downstream/reciprocal_orthology/MmuNcVsHg38_gt100_ext500.bed

# Extract first exon of human ncRNAs
awk 'BEGIN{OFS="\t"}{split($11,SZ,",");if ($6=="+") {print $1,$2,$2+SZ[1],$4,$5,$6} else {print $1,$3-SZ[$10],$3,$4,$5,$6}}' $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap.bed > $BASEDIR/downstream/reciprocal_orthology/Gencode-currentPlusNovel_ncSplice_noCodeOlap_5p.bed

$INTERSECT_BED -wo -a $BASEDIR/downstream/reciprocal_orthology/MmuNcVsHg38_gt100_ext500.bed -b $BASEDIR/downstream/reciprocal_orthology/Gencode-currentPlusNovel_ncSplice_noCodeOlap_5p.bed | cut -f4,10 | sort -u > $BASEDIR/downstream/reciprocal_orthology/mouse2human_blasted.txt


#join -1 1 -2 2 -a 1 -t $'\t' <(cut -f4,13,25,26 $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx | sort -k1,1) <(sort -k2,2 $BASEDIR/downstream/reciprocal_orthology/mouse2human_blasted.txt) | head





rm -rf $TMP

