#!/bin/bash
set -xe -o pipefail
shopt -s extglob

PC=$BASEDIR/nc/posConNCwithCodingPartner.bedx
HIC=$BASEDIR/downstream/hic
DIR=$BASEDIR/downstream/ctcf

mkdir -p $DIR/data

if [[ ! -f $DIR/data/tfbs.bed ]]; then
	wget -O - http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz | gunzip > $DIR/data/tfbs.bed
fi

awk '$4=="CTCF"' $DIR/data/tfbs.bed > $DIR/data/ctcf_hg19.bed

$LIFTOVER -bedPlus=5 $DIR/data/ctcf_hg19.bed $BASEDIR/data/hg19ToHg38.over.chain $DIR/data/ctcf.unsorted.bed $DIR/data/ctcf.unmapped

sort -k1,1 -k2,2n -k3,3n $DIR/data/ctcf.unsorted.bed > $DIR/data/ctcf.bed




#################################################################
#  How many pcRNAs are associated with CTCF sites in their pomoter?
#################################################################

mkdir -p $DIR/CtcfOverlapPromoter

$INTERSECT_BED -wo -a $HIC/posConAnnot_PROMOTER+-2000.bed -b $DIR/data/ctcf.bed | bedtools groupby -i - -g 4 -c 4 -o count > $DIR/CtcfOverlapPromoter/pc_with_ctcf_in_prom.txt
$INTERSECT_BED -wo -a $HIC/Gencode.v21.linc.spliced_PROMOTER+-2000.bed -b $DIR/data/ctcf.bed | bedtools groupby -i - -g 4 -c 4 -o count > $DIR/CtcfOverlapPromoter/lincs_with_ctcf_in_prom.txt
$INTERSECT_BED -wo -a $HIC/coding_PROMOTER+-2000.bed -b $DIR/data/ctcf.bed | bedtools groupby -i - -g 4 -c 4 -o count > $DIR/CtcfOverlapPromoter/pcCoding_with_ctcf_in_prom.txt
$INTERSECT_BED -wo -a $BASEDIR/coding/Gencode-current_hsa_coding_PROMOTER+-2000.bed -b $DIR/data/ctcf.bed | bedtools groupby -i - -g 4 -c 4 -o count > $DIR/CtcfOverlapPromoter/gencode_coding_with_ctcf_in_prom.txt

cat <<EOT >$DIR/CtcfOverlapPromoter/ctcf_contact_summary.txt
Category WithLoop Total
pcRNA $(wc -l $DIR/CtcfOverlapPromoter/pc_with_ctcf_in_prom.txt | awk '{print $1}') $(wc -l $HIC/posConAnnot_PROMOTER+-2000.bed | awk '{print $1}')
lincRNAs $(wc -l $DIR/CtcfOverlapPromoter/lincs_with_ctcf_in_prom.txt | awk '{print $1}') $(wc -l $HIC/Gencode.v21.linc.spliced_PROMOTER+-2000.bed | awk '{print $1}')
pcCoding $(wc -l $DIR/CtcfOverlapPromoter/pcCoding_with_ctcf_in_prom.txt | awk '{print $1}') $(wc -l $HIC/coding_PROMOTER+-2000.bed | awk '{print $1}')
GencodeCoding $(wc -l $DIR/CtcfOverlapPromoter/gencode_coding_with_ctcf_in_prom.txt | awk '{print $1}') $(wc -l $BASEDIR/coding/Gencode-current_hsa_coding_PROMOTER+-2000.bed | awk '{print $1}')
EOT



#################################################################
#       ALL: Make heatmaps of CTCF overlap
#################################################################
# bamCoverage and computeMatrix are part of deepTools
mkdir -p $DIR/pcCoverageByCtcfAllLines


# Convert the loops to big wig
awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,"CTCF_"NR}' $DIR/data/ctcf.bed > $DIR/data/ctcf.bed4
$BEDTOBAM -i $DIR/data/ctcf.bed4 -g $BASEDIR/data/hg38.chromSizes > $DIR/data/ctcf.bam.unsorted

$SAMTOOLS sort $DIR/data/ctcf.bam.unsorted $DIR/data/ctcf
$SAMTOOLS index $DIR/data/ctcf.bam $DIR/data/ctcf.bam.bai
bamCoverage -b $DIR/data/ctcf.bam -o $DIR/data/ctcf.bw


# Compute matrices
mkdir -p $DIR/pcCoverageByCtcfAllLines/matrix/
computeMatrix scale-regions -b 20000 -a 20000 -m 10000 --regionsFileName $HIC/pcCoverageByLoops/annotationByPosition/posConAnnot.bed --sortRegions "no" --outFileNameMatrix $DIR/pcCoverageByCtcfAllLines/matrix/peaks_pc_matrix.txt --scoreFileName $DIR/data/ctcf.bw --missingDataAsZero --outFileName /dev/null -p 1

Rscript $BIN/dtCompare.R --files $DIR/pcCoverageByCtcfAllLines/matrix/peaks_pc_matrix.txt --labels Ctcf --profileStdErr --noHeat --profW 10 --profH 10 --outFile $DIR/pcCoverageByCtcfAllLines/CTCF_coverage

convert -density 300 $DIR/pcCoverageByCtcfAllLines/CTCF_coverage_profile.pdf -quality 100 $DIR/pcCoverageByCtcfAllLines/CTCF_coverage_profile.png

