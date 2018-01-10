#!/bin/bash
set -xe -o pipefail
shopt -s extglob
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh




PC=$BASEDIR/nc/posConNCwithCodingPartner.bedx
HIC=$BASEDIR/downstream/hic
P=8
DIR=$BASEDIR/downstream/hic-tads
mkdir -p $DIR/data


CELLLINES="HMEC HUVEC NHEK K562 HeLa KBM7 IMR90 GM12878"

for i in $CELLLINES; do
        if [[ ! -f $DIR/data/${i}_domainlist_HG19.txt ]]; then
                if [[ $i == "GM12878" ]]; then
                        URL="GM12878_primary+replicate"
                else
                        URL=$i
                fi
                wget --limit-rate=100k -O - "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${URL}_Arrowhead_domainlist.txt.gz" | gunzip > $DIR/data/${i}_domainlist_HG19.txt
        fi
        # Convert them to bed and lift over
        awk 'BEGIN{OFS=FS="\t"}{print "chr"$1,$2,$6,".",0,".",$2,$2,"0,0,0","2",$3-$2","$6-$5,0","$5-$2}' $DIR/data/${i}_domainlist_HG19.txt | tail -n +2 > $DIR/data/${i}_domainlist_HG19.bed
        $LIFT $DIR/data/${i}_domainlist_HG19.bed $BASEDIR/data/hg19ToHg38.over.chain $DIR/data/${i}_domainlist_unsorted.bed $DIR/data/${i}_domainlist_unmapped
        sort -k1,1 -k2,2n -k3,3n $DIR/data/${i}_domainlist_unsorted.bed > $DIR/data/${i}_domainlist.bed
done

for i in $CELLLINES; do
	cat $DIR/data/${i}_domainlist.bed | sort -k1,1 -k2,2n -k3,3n | cut -f1-3 > $DIR/data/allLines_domainlist.bed
done


#################################################################
#  How many pcRNAs are associated with TADs in the promoter?
#################################################################

mkdir -p $DIR/TADOverlapAllLinesPromoter
for n in $(seq 1000 1000 50000); do
	
	# Annotate TAD ends
	awk -v n=${n} 'BEGIN{OFS=FS="\t"}{print $1,$2-n,$2+n,"TAD"NR"-5p\n"$1,$3-n,$3+n,"TAD"NR"-3p"}' $DIR/data/allLines_domainlist.bed | sort -k1,1 -k2,2n | $MERGE -c 4 -o collapse  -i - > $DIR/TADOverlapAllLinesPromoter/allLines_TADends_${n}kb.bed
	
	$INTERSECT_BED -wo -a $HIC/posConAnnot_PROMOTER+-2000.bed -b $DIR/TADOverlapAllLinesPromoter/allLines_TADends_${n}kb.bed | $GROUP_BY -i - -g 4 -c 4 -o count > $DIR/TADOverlapAllLinesPromoter/pc_with_TADs_in_prom_${n}kb.txt
	$INTERSECT_BED -wo -a $HIC/Gencode.v21.linc.spliced_PROMOTER+-2000.bed -b $DIR/TADOverlapAllLinesPromoter/allLines_TADends_${n}kb.bed | $GROUP_BY -i - -g 4 -c 4 -o count > $DIR/TADOverlapAllLinesPromoter/lincs_with_TADs_in_prom_${n}kb.txt
	$INTERSECT_BED -wo -a $HIC/coding_PROMOTER+-2000.bed -b $DIR/TADOverlapAllLinesPromoter/allLines_TADends_${n}kb.bed | $GROUP_BY -i - -g 4 -c 4 -o count > $DIR/TADOverlapAllLinesPromoter/pcCoding_with_TADs_in_prom_${n}kb.txt
	$INTERSECT_BED -wo -a $BASEDIR/coding/Gencode-current_hsa_coding_PROMOTER+-2000.bed -b $DIR/TADOverlapAllLinesPromoter/allLines_TADends_${n}kb.bed | $GROUP_BY -i - -g 4 -c 4 -o count > $DIR/TADOverlapAllLinesPromoter/gencode_coding_with_TADs_in_prom_${n}kb.txt
	
cat <<EOT >$DIR/TADOverlapAllLinesPromoter/TADs_contact_summary_${n}kb.txt
Category WithLoop Total
pcRNA $(wc -l $DIR/TADOverlapAllLinesPromoter/pc_with_TADs_in_prom_${n}kb.txt | awk '{print $1}') $(wc -l $HIC/posConAnnot_PROMOTER+-2000.bed | awk '{print $1}')
lincRNAs $(wc -l $DIR/TADOverlapAllLinesPromoter/lincs_with_TADs_in_prom_${n}kb.txt | awk '{print $1}') $(wc -l $HIC/Gencode.v21.linc.spliced_PROMOTER+-2000.bed | awk '{print $1}')
pcCoding $(wc -l $DIR/TADOverlapAllLinesPromoter/pcCoding_with_TADs_in_prom_${n}kb.txt | awk '{print $1}') $(wc -l $HIC/coding_PROMOTER+-2000.bed | awk '{print $1}')
GencodeCoding $(wc -l $DIR/TADOverlapAllLinesPromoter/gencode_coding_with_TADs_in_prom_${n}kb.txt | awk '{print $1}') $(wc -l $BASEDIR/coding/Gencode-current_hsa_coding_PROMOTER+-2000.bed | awk '{print $1}')
EOT
done

(for n in $(seq 1000 1000 50000); do
	awk -v n=${n} 'BEGIN{OFS="\t"; FS=" "}!/^Category/{print n,$1,$2,$3}' $DIR/TADOverlapAllLinesPromoter/TADs_contact_summary_${n}kb.txt
done)>$DIR/TADOverlapAllLinesPromoter/TADs_contact_summary_all.txt


### Domain coverage by pcRNAs
# Convert the pcRNAs to big wig
$BEDTOBAM -i  $DIR/data/allLines_TADends.bed -g $BASEDIR/data/hg38.chromSizes >  $DIR/data/allLines_TADends.bam.unsorted
$SAMTOOLS sort $DIR/data/allLines_TADends.bam.unsorted $DIR/data/allLines_TADends
$SAMTOOLS index $DIR/data/allLines_TADends.bam $DIR/data/allLines_TADends.bam.bai
bamCoverage -b $DIR/data/allLines_TADends.bam -o $DIR/data/allLines_TADends.bw

mkdir -p $DIR/pcCoveragebyTADends/matrix/
computeMatrix scale-regions -b 20000 -a 20000 -m 10000 --regionsFileName $HIC/pcCoverageByLoops/annotationByPosition/posConAnnot.bed --sortRegions "no" --outFileNameMatrix $DIR/pcCoveragebyTADends/matrix/peaks_pc_matrix.txt --scoreFileName $DIR/data/allLines_TADends.bw --missingDataAsZero --outFileName /dev/null -p $P
/homes/tl344/tom/bin/R312/bin/Rscript ~/tom/my_scripts/dtCompare.R --files $DIR/pcCoveragebyTADends/matrix/peaks_pc_matrix.txt --labels Loops --profileStdErr --noHeat --profW 10 --profH 10 --outFile $DIR/pcCoveragebyTADends/pcCoverageByTADs
convert -density 300 $DIR/pcCoveragebyTADends/pcCoverageByTADs_profile.pdf -quality 100 $DIR/pcCoveragebyTADends/pcCoverageByTADs_profile.png



touch touch $BASEDIR/.tad
