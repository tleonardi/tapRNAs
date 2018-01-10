#!/bin/bash
set -xe -o pipefail
shopt -s extglob
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh


PC=$BASEDIR/nc/posConNCwithCodingPartner.bedx
P=16
DIR=$BASEDIR/downstream/hic
mkdir -p $DIR/data


#################################################################
#               PREPARE DATA
#################################################################

# Find transcript IDs for the PC-associated coding genes
join -1 2 -2 1 -t $'\t' <(cut -f1-2 $BASEDIR/data/GencodeV21_info_tr.txt | sort -k2,2) <(cut -f13 $PC | sort -u) > $DIR/pc_coding_txIDs.txt

# Annotate pcRNA-associated coding genes 
$PSTR join --idCol 2 -x --columns 2 $DIR/pc_coding_txIDs.txt $BASEDIR/data/Gencode-current_hsa.bed | cut -f 1-12 > $DIR/coding.bed

# Annotate Gencode spliced lincs with no coding overlap
awk '$10>1' $BASEDIR/data/Gencode.v21.lincRNAs.bed | $OVERLAPSELECT -strand -overlapBases=20 -nonOverlapping -inFmt=bed $BASEDIR/coding/Gencode-current_hsa_codingExons.bed stdin stdout > $DIR/Gencode.v21.linc.spliced.bed


# Annotate pcRNA promoters
awk 'BEGIN{OFS=FS="\t"}{if($6=="+"){print $1,$2-2000,$2+2000,$4,$5,$6}if($6=="-"){print $1,$3-2000,$3+2000,$4,$5,$6}}' $PC | sort -k1,1 -k2,2n -k3,3n | bedtools merge -c 4 -o distinct -i - > $DIR/posConAnnot_PROMOTER+-2000.bed
# Annotate Gencode linc promoters
awk 'BEGIN{OFS=FS="\t"}{if($6=="+"){print $1,$2-2000,$2+2000,$4,$5,$6}if($6=="-"){print $1,$3-2000,$3+2000,$4,$5,$6}}' $DIR/Gencode.v21.linc.spliced.bed | sort -k1,1 -k2,2n -k3,3n | bedtools merge -c 4 -o distinct -i - > $DIR/Gencode.v21.linc.spliced_PROMOTER+-2000.bed
# Annotate Gencode coding promoters
awk 'BEGIN{OFS=FS="\t"}$1!="chrM"{if($6=="+"){print $1,$2-2000,$2+2000,$4,$5,$6}if($6=="-"){print $1,$3-2000,$3+2000,$4,$5,$6}}' $BASEDIR/coding/Gencode-current_hsa_coding.bed | sort -k1,1 -k2,2n -k3,3n | bedtools merge -c 4 -o distinct -i - > $BASEDIR/coding/Gencode-current_hsa_coding_PROMOTER+-2000.bed
# Annotate pcRNA-associated coding promoters
awk 'BEGIN{OFS=FS="\t"}{if($6=="+"){print $1,$2-2000,$2+2000,$4,$5,$6}if($6=="-"){print $1,$3-2000,$3+2000,$4,$5,$6}}' $DIR/coding.bed | sort -k1,1 -k2,2n -k3,3n | bedtools merge -c 4 -o distinct -i - > $DIR/coding_PROMOTER+-2000.bed

# Make BED file of all pcRNAs for Deeptools
mkdir -p $DIR/pcCoverageByLoops/annotationByPosition/
cut -f1-6 $PC > $DIR/pcCoverageByLoops/annotationByPosition/posConAnnot_allPC.bed
echo "#pcRNAs" >> $DIR/pcCoverageByLoops/annotationByPosition/posConAnnot_allPC.bed

# Make BED file of all pcRNA-coding for Deeptools
cut -f1-6 $DIR/coding.bed > $DIR/pcCoverageByLoops/annotationByPosition/coding.bed
echo "#pcCoding" >> $DIR/pcCoverageByLoops/annotationByPosition/coding.bed

# Make BED file of randome Gencode lincRNAs for Deeptools
bedtools sample -n 5000 -seed 383847 -i  $DIR/Gencode.v21.linc.spliced.bed | cut -f1-6 > $DIR/pcCoverageByLoops/annotationByPosition/posConAnnot_OnlyGencodeLinc.bed
echo "#GencodeLinc" >> $DIR/pcCoverageByLoops/annotationByPosition/posConAnnot_OnlyGencodeLinc.bed

# Make BED file of randome Gencode coding for Deeptools
bedtools sample -n 5000 -seed 383847 -i $BASEDIR/coding/Gencode-current_hsa_coding.bed  | cut -f1-6 > $DIR/pcCoverageByLoops/annotationByPosition/posConAnnot_OnlyGencodeCoding.bed
echo "#GencodeCoding" >> $DIR/pcCoverageByLoops/annotationByPosition/posConAnnot_OnlyGencodeCoding.bed



# Annotate the different pcRNA categories
mkdir -p  $DIR/pcAnnotByOrientation
for i in $(cut -f15 $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx | sort -u); do
        awk -v I=${i} 'BEGIN{OFS=FS="\t"} $15==I{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}'  $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx > $DIR/pcAnnotByOrientation/posConAnnot_Only${i}.bed
	echo "#$i" >> $DIR/pcAnnotByOrientation/posConAnnot_Only${i}.bed
done

# And their promoters
for i in $DIR/pcAnnotByOrientation/posConAnnot_!(*promoter*).bed; do
	NAME=$(basename $i .bed)
	awk 'BEGIN{OFS=FS="\t"}{if($6=="+"){print $1,$2-5000,$2+5000,$4,$5,$6}if($6=="-"){print $1,$3-5000,$3+5000,$4,$5,$6}}' $i | sort -k1,1 -k2,2n -k3,3n | bedtools merge -c 4 -o distinct -i - > $DIR/pcAnnotByOrientation/${NAME}_promoter.bed;
done



#########################################################
#			Analysis			#	
#########################################################


CELLLINES="HMEC HUVEC NHEK K562 HeLa KBM7 IMR90 GM12878"

for i in $CELLLINES; do
	if [[ ! -f $DIR/data/${i}_looplist_HG19.txt ]]; then
		if [[ $i == "GM12878" ]]; then
			URL="GM12878_primary+replicate"
		else
			URL=$i
		fi
		wget --limit-rate=100k -O - "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${URL}_HiCCUPS_looplist.txt.gz" | gunzip > $DIR/data/${i}_looplist_HG19.txt
	fi
	# Convert them to bed and lift over
	awk 'BEGIN{OFS=FS="\t"}{print "chr"$1,$2,$6,".",0,".",$2,$2,"0,0,0","2",$3-$2","$6-$5,0","$5-$2}' $DIR/data/${i}_looplist_HG19.txt | tail -n +2 > $DIR/data/${i}_looplist_HG19.bed
	$LIFTOVER $DIR/data/${i}_looplist_HG19.bed $BASEDIR/data/hg19ToHg38.over.chain $DIR/data/${i}_looplist_unsorted.bed $DIR/data/${i}_looplist_unmapped
	sort -k1,1 -k2,2n -k3,3n $DIR/data/${i}_looplist_unsorted.bed > $DIR/data/${i}_looplist.bed
done

rm -f $DIR/data/allLines_looplist.bed
cat $DIR/data/*_looplist.bed | sort -k1,1 -k2,2n -k3,3n > $DIR/data/allLines_looplist.bed

# Annotate the loop starting points (precisely, both starts and ends)
awk 'BEGIN{OFS=FS="\t"}{split($11,LEN,",");print $1,$2,$2+LEN[1],$1"-"$3-LEN[2]"-"$3"\n"$1,$3-LEN[2],$3,$1"-"$2"-"$2+LEN[1]}' $DIR/data/allLines_looplist.bed | sort -u | sort -k1,1 -k2,2n -k3,3n >  $DIR/allLines_loop_start_end_named.bedx

# Annotate the 12.903 distinct peaks (they are only 12851 because some did not lift over)
$MERGE -i $DIR/allLines_loop_start_end_named.bedx > $DIR/allLinespeaks.bed

# Annotate peaks distinctly
for i in $CELLLINES; do
	FILE=$DIR/data/${i}_looplist.bed
	awk 'BEGIN{OFS=FS="\t"}{split($11,LEN,",");print $1,$2,$2+LEN[1],$1"-"$3-LEN[2]"-"$3"\n"$1,$3-LEN[2],$3,$1"-"$2"-"$2+LEN[1]}' $FILE | sort -k1,1 -k2,2n -k3,3n > $DIR/${i}_loop_start_end_named.bedx
	$MERGE -i $DIR/${i}_loop_start_end_named.bedx > $DIR/${i}_peaks.bed
done
 

#################################################################
#  How many pcRNAs are associated with loops in the promoter?
#################################################################

mkdir -p $DIR/loopOverlapAllLinesPromoter

$INTERSECT_BED -wo -a $DIR/posConAnnot_PROMOTER+-2000.bed -b $DIR/allLinespeaks.bed | bedtools groupby -i - -g 4 -c 4 -o count > $DIR/loopOverlapAllLinesPromoter/pc_with_loops_in_prom.txt
$INTERSECT_BED -wo -a $DIR/Gencode.v21.linc.spliced_PROMOTER+-2000.bed -b $DIR/allLinespeaks.bed | bedtools groupby -i - -g 4 -c 4 -o count > $DIR/loopOverlapAllLinesPromoter/lincs_with_loops_in_prom.txt
$INTERSECT_BED -wo -a $DIR/coding_PROMOTER+-2000.bed -b $DIR/allLinespeaks.bed | bedtools groupby -i - -g 4 -c 4 -o count > $DIR/loopOverlapAllLinesPromoter/pcCoding_with_loops_in_prom.txt
$INTERSECT_BED -wo -a $BASEDIR/coding/Gencode-current_hsa_coding_PROMOTER+-2000.bed -b $DIR/allLinespeaks.bed | bedtools groupby -i - -g 4 -c 4 -o count > $DIR/loopOverlapAllLinesPromoter/gencode_coding_with_loops_in_prom.txt

cat <<EOT >$DIR/loopOverlapAllLinesPromoter/loops_contact_summary.txt
Category WithLoop Total
pcRNA $(wc -l $DIR/loopOverlapAllLinesPromoter/pc_with_loops_in_prom.txt | awk '{print $1}') $(wc -l $DIR/posConAnnot_PROMOTER+-2000.bed | awk '{print $1}')
lincRNAs $(wc -l $DIR/loopOverlapAllLinesPromoter/lincs_with_loops_in_prom.txt | awk '{print $1}') $(wc -l $DIR/Gencode.v21.linc.spliced_PROMOTER+-2000.bed | awk '{print $1}')
pcCoding $(wc -l $DIR/loopOverlapAllLinesPromoter/pcCoding_with_loops_in_prom.txt | awk '{print $1}') $(wc -l $DIR/coding_PROMOTER+-2000.bed | awk '{print $1}')
GencodeCoding $(wc -l $DIR/loopOverlapAllLinesPromoter/gencode_coding_with_loops_in_prom.txt | awk '{print $1}') $(wc -l $BASEDIR/coding/Gencode-current_hsa_coding_PROMOTER+-2000.bed | awk '{print $1}')
EOT


#################################################################
#        ALL: How many pcRNAs are associated with loops?
#################################################################

mkdir -p $DIR/loopOverlapAllLines
# Intersect loops with pcRNAs
$INTERSECT_BED -wo -a $PC -b $DIR/allLinespeaks.bed | bedtools groupby -i - -g 4 -c 4 -o count > $DIR/loopOverlapAllLines/pc_with_loops.txt
$INTERSECT_BED -wo -split -a $PC -b $DIR/allLinespeaks.bed | bedtools groupby -i - -g 4 -c 4 -o count > $DIR/loopOverlapAllLines/pc_with_exonic_loops.txt

# Do the same for gencode lincRNAs as control
$INTERSECT_BED -wo -a $DIR/Gencode.v21.linc.spliced.bed -b $DIR/allLinespeaks.bed | bedtools groupby -i - -g 4 -c 4 -o count > $DIR/loopOverlapAllLines/lincs_with_loops.txt
$INTERSECT_BED -wo -split -a $DIR/Gencode.v21.linc.spliced.bed -b $DIR/allLinespeaks.bed | bedtools groupby -i - -g 4 -c 4 -o count > $DIR/loopOverlapAllLines/lincs_with_exonic_loops.txt

# Do the same for pcRNA protein coding genes
$INTERSECT_BED -wo -a $DIR/coding.bed -b $DIR/allLinespeaks.bed | bedtools groupby -i - -g 4 -c 4 -o count > $DIR/loopOverlapAllLines/pcCoding_with_loops.txt
$INTERSECT_BED -wo -split -a $DIR/coding.bed -b $DIR/allLinespeaks.bed | bedtools groupby -i - -g 4 -c 4 -o count >$DIR/loopOverlapAllLines/pcCoding_with_exonic_loops.txt

# Do the same for all Gencode coding
$INTERSECT_BED -wo -a $BASEDIR/coding/Gencode-current_hsa_coding.bed -b $DIR/allLinespeaks.bed | bedtools groupby -i - -g 4 -c 4 -o count > $DIR/loopOverlapAllLines/gencode_coding_with_loops.txt
$INTERSECT_BED -wo -split -a $BASEDIR/coding/Gencode-current_hsa_coding.bed -b $DIR/allLinespeaks.bed | bedtools groupby -i - -g 4 -c 4 -o count > $DIR/loopOverlapAllLines/gencode_coding_with_exonic_loops.txt

# Do it for pcRNAs, by orientation
mkdir -p $DIR/loopOverlapAllLines/byOrientation
for i in $DIR/pcAnnotByOrientation/posConAnnot_!(*promoter*).bed; do
        NAME=$(basename $i .bed)
        $INTERSECT_BED -wo -a $i -b $DIR/allLinespeaks.bed | bedtools groupby -i - -g 4 -c 4 -o count > $DIR/loopOverlapAllLines/byOrientation/${NAME}_with_loops.txt
        $INTERSECT_BED -wo -split -a $i -b $DIR/allLinespeaks.bed | bedtools groupby -i - -g 4 -c 4 -o count > $DIR/loopOverlapAllLines/byOrientation/${NAME}_with_exonic_loops.txt
done

# Make a summary table
cat <<EOT >$DIR/loopOverlapAllLines/loops_contact_summary.txt
Category WithLoop WithExonicLoop Total
pcRNA $(wc -l $DIR/loopOverlapAllLines/pc_with_loops.txt| awk '{print $1}') $(wc -l $DIR/loopOverlapAllLines/pc_with_exonic_loops.txt| awk '{print $1}') $(wc -l $PC | awk '{print $1}')
lincRNAs $(wc -l $DIR/loopOverlapAllLines/lincs_with_loops.txt| awk '{print $1}') $(wc -l $DIR/loopOverlapAllLines/lincs_with_exonic_loops.txt| awk '{print $1}') $(wc -l $DIR/Gencode.v21.linc.spliced.bed | awk '{print $1}')
pcCoding $(wc -l $DIR/loopOverlapAllLines/pcCoding_with_loops.txt| awk '{print $1}') $(wc -l $DIR/loopOverlapAllLines/pcCoding_with_exonic_loops.txt| awk '{print $1}') $(wc -l $DIR/coding.bed | awk '{print $1}')
GencodeCoding $(wc -l $DIR/loopOverlapAllLines/gencode_coding_with_loops.txt| awk '{print $1}') $(wc -l $DIR/loopOverlapAllLines/gencode_coding_with_exonic_loops.txt| awk '{print $1}') $(wc -l $BASEDIR/coding/Gencode-current_hsa_coding.bed | awk '{print $1}')
Asense $(wc -l $DIR/loopOverlapAllLines/byOrientation/posConAnnot_OnlyAsense_with_loops.txt | awk '{print $1}') $(wc -l $DIR/loopOverlapAllLines/byOrientation/posConAnnot_OnlyAsense_with_exonic_loops.txt | awk '{print $1}') $(wc -l $DIR/pcAnnotByOrientation/posConAnnot_OnlyAsense.bed | awk '{print $1}')
Bidir $(wc -l $DIR/loopOverlapAllLines/byOrientation/posConAnnot_OnlyBIDIR_with_loops.txt | awk '{print $1}') $(wc -l $DIR/loopOverlapAllLines/byOrientation/posConAnnot_OnlyBIDIR_with_exonic_loops.txt | awk '{print $1}') $(wc -l $DIR/pcAnnotByOrientation/posConAnnot_OnlyBIDIR.bed | awk '{print $1}')
DS-AS $(wc -l $DIR/loopOverlapAllLines/byOrientation/posConAnnot_OnlyDS-AS_with_loops.txt | awk '{print $1}') $(wc -l $DIR/loopOverlapAllLines/byOrientation/posConAnnot_OnlyDS-AS_with_exonic_loops.txt | awk '{print $1}') $(wc -l $DIR/pcAnnotByOrientation/posConAnnot_OnlyDS-AS.bed | awk '{print $1}')
DS-S $(wc -l $DIR/loopOverlapAllLines/byOrientation/posConAnnot_OnlyDS-S_with_loops.txt | awk '{print $1}') $(wc -l $DIR/loopOverlapAllLines/byOrientation/posConAnnot_OnlyDS-S_with_exonic_loops.txt | awk '{print $1}') $(wc -l $DIR/pcAnnotByOrientation/posConAnnot_OnlyDS-S.bed | awk '{print $1}')
Olap $(wc -l $DIR/loopOverlapAllLines/byOrientation/posConAnnot_OnlyOLAP_with_loops.txt | awk '{print $1}') $(wc -l $DIR/loopOverlapAllLines/byOrientation/posConAnnot_OnlyOLAP_with_exonic_loops.txt | awk '{print $1}') $(wc -l $DIR/pcAnnotByOrientation/posConAnnot_OnlyOLAP.bed | awk '{print $1}')
US-AS $(wc -l $DIR/loopOverlapAllLines/byOrientation/posConAnnot_OnlyUS-AS_with_loops.txt | awk '{print $1}') $(wc -l $DIR/loopOverlapAllLines/byOrientation/posConAnnot_OnlyUS-AS_with_exonic_loops.txt | awk '{print $1}') $(wc -l $DIR/pcAnnotByOrientation/posConAnnot_OnlyUS-AS.bed | awk '{print $1}')
US-S $(wc -l $DIR/loopOverlapAllLines/byOrientation/posConAnnot_OnlyUS-S_with_loops.txt | awk '{print $1}') $(wc -l $DIR/loopOverlapAllLines/byOrientation/posConAnnot_OnlyUS-S_with_exonic_loops.txt | awk '{print $1}') $(wc -l $DIR/pcAnnotByOrientation/posConAnnot_OnlyUS-S.bed | awk '{print $1}')
EOT





#################################################################
#        Intersect with HMM
#################################################################

mkdir -p $DIR/endPoints/
$INTERSECT_BED -wa -wb -a $DIR/posConAnnot_PROMOTER+-2000.bed -b $DIR/allLines_loop_start_end_named.bedx | awk 'BEGIN{OFS=FS="\t"}{split($8, a, "-"); print a[1],a[2],a[3],$4}' | sort -u | sort -k1,1 -k2,2n -k3,3n | $MERGE -c 4 -o distinct -i - > $DIR/endPoints/pcPeak_endPoint.bedx
$INTERSECT_BED -wa -wb -a $DIR/Gencode.v21.linc.spliced_PROMOTER+-2000.bed -b $DIR/allLines_loop_start_end_named.bedx | awk 'BEGIN{OFS=FS="\t"}{split($8, a, "-"); print a[1],a[2],a[3],$4}' | sort -u | sort -k1,1 -k2,2n -k3,3n | $MERGE -c 4 -o distinct -i - > $DIR/endPoints/lincRNAPeak_endPoint.bedx
# Download and liftover HMM data
HMMLINES="Gm12878 H1hesc Hepg2 Hmec Hsmm Huvec K562 Nhek Nhlf"

for i in $HMMLINES; do
	if [[ ! -f $DIR/data/${i}_HMM.bed ]]; then
 	       wget --limit-rate=100k -O - "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmm${i}HMM.bed.gz" | gunzip > $DIR/data/${i}_HMM_HG19.bed
	       $LIFTOVER $DIR/data/${i}_HMM_HG19.bed $BASEDIR/data/hg19ToHg38.over.chain $DIR/data/${i}_HMM_unsorted.bed $DIR/data/${i}_HMM_unmapped
	       sort -k1,1 -k2,2n -k3,3n $DIR/data/${i}_HMM_unsorted.bed > $DIR/data/${i}_HMM.bed
	fi
	# Intersect the end point with chromHMM
	$INTERSECT_BED -wo -a $DIR/endPoints/pcPeak_endPoint.bedx -b $DIR/data/${i}_HMM.bed > $DIR/pcRNA_loop_endPoints_hmm${i}.bed
	awk 'BEGIN{OFS=FS="\t"}{print $1$2$3$4,$3-$2,$14,$8}' $DIR/pcRNA_loop_endPoints_hmm${i}.bed >  $DIR/pcRNA_loop_endPoints_hmm${i}_table.txt
	# Intersect all end-points with chrom HMM for background
	$INTERSECT_BED -wo -a $DIR/endPoints/lincRNAPeak_endPoint.bedx -b $DIR/data/${i}_HMM.bed > $DIR/genome_loop_endPoints_hmm${i}.bed
	awk 'BEGIN{OFS=FS="\t"}{print $1$2$3$4,$3-$2,$14,$8}'  $DIR/genome_loop_endPoints_hmm${i}.bed > $DIR/genome_loop_endPoints_hmm${i}_table.txt
done





#################################################################
#       ALL: Make heatmaps of loop overlap
#################################################################

mkdir -p $DIR/pcCoverageByLoopsAllLines


# Convert the loops to big wig
awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,"Peak_"NR}' $DIR/allLinespeaks.bed > $DIR/allLinespeaks.bed4
$BEDTOBAM -i $DIR/allLinespeaks.bed4 -g $BASEDIR/data/hg38.chromSizes > $DIR/allLinespeaks.bam.unsorted
$SAMTOOLS sort $DIR/allLinespeaks.bam.unsorted $DIR/allLinespeaks
$SAMTOOLS index $DIR/allLinespeaks.bam $DIR/allLinespeaks.bam.bai
bamCoverage -b $DIR/allLinespeaks.bam -o $DIR/AllLinespeaks.bw

# Combine pcRNA annotation by category
cat $DIR/pcCoverageByLoops/annotationByPosition/posConAnnot_allPC.bed $DIR/pcCoverageByLoops/annotationByPosition/coding.bed $DIR/pcCoverageByLoops/annotationByPosition/posConAnnot_OnlyGencodeLinc.bed $DIR/pcCoverageByLoops/annotationByPosition/posConAnnot_OnlyGencodeCoding.bed > $DIR/pcCoverageByLoops/annotationByPosition/posConAnnot.bed
cat $DIR/pcAnnotByOrientation/posConAnnot_Only!(*promoter.bed) | cut -f1-6 > $DIR/pcCoverageByLoops/annotationByPosition/posConAnnot_pcByPos.bed

# Compute matrices
computeMatrix scale-regions -b 20000 -a 20000 -m 10000 --regionsFileName $DIR/pcCoverageByLoops/annotationByPosition/posConAnnot.bed --sortRegions "no" --outFileNameMatrix $DIR/pcCoverageByLoopsAllLines/matrix/peaks_pc_matrix.txt --scoreFileName $DIR/AllLinespeaks.bw --missingDataAsZero --outFileName /dev/null -p $P
/homes/tl344/tom/bin/R312/bin/Rscript ~/tom/my_scripts/dtCompare.R --files $DIR/pcCoverageByLoopsAllLines/matrix/peaks_pc_matrix.txt --labels Loops --profileStdErr --noHeat --profW 10 --profH 10 --outFile $DIR/pcCoverageByLoopsAllLines/pcCoverageByLoops
convert -density 300 $DIR/pcCoverageByLoopsAllLines/pcCoverageByLoops_profile.pdf -quality 100 $DIR/pcCoverageByLoopsAllLines/pcCoverageByLoops_profile.png


computeMatrix scale-regions -b 20000 -a 20000 -m 10000 --regionsFileName $DIR/pcCoverageByLoops/annotationByPosition/posConAnnot_pcByPos.bed --sortRegions "no" --outFileNameMatrix $DIR/pcCoverageByLoopsAllLines/matrix/peaks_pc_matrix_pcByPos.txt --scoreFileName $DIR/AllLinespeaks.bw --missingDataAsZero --outFileName /dev/null -p $P
/homes/tl344/tom/bin/R312/bin/Rscript ~/tom/my_scripts/dtCompare.R --files $DIR/pcCoverageByLoopsAllLines/matrix/peaks_pc_matrix_pcByPos.txt --labels Loops --profileStdErr --noHeat --profW 10 --profH 10 --outFile $DIR/pcCoverageByLoopsAllLines/pcCoverageByLoops_pcByPos
convert -density 300 $DIR/pcCoverageByLoopsAllLines/pcCoverageByLoops_pcByPos_profile.pdf -quality 100 $DIR/pcCoverageByLoopsAllLines/pcCoverageByLoops_pcByPos_profile.png


# Repeat for promoter
computeMatrix reference-point --referencePoint TSS -b 10000 -a 10000 --regionsFileName $DIR/pcCoverageByLoops/annotationByPosition/posConAnnot.bed --sortRegions "no" --outFileNameMatrix $DIR/pcCoverageByLoopsAllLines/matrix/peaks_pc_matrix_promoter.txt --scoreFileName $DIR/AllLinespeaks.bw --missingDataAsZero --outFileName /dev/null -p $P
/homes/tl344/tom/bin/R312/bin/Rscript ~/tom/my_scripts/dtCompare.R --files $DIR/pcCoverageByLoopsAllLines/matrix/peaks_pc_matrix_promoter.txt --labels Loops --profileStdErr --noHeat --profW 10 --profH 10 --outFile $DIR/pcCoverageByLoopsAllLines/pcCoverageByLoop_promoters
convert -density 300 $DIR/pcCoverageByLoopsAllLines/pcCoverageByLoop_promoters_profile.pdf -quality 100 $DIR/pcCoverageByLoopsAllLines/pcCoverageByLoop_promoters_profile.png


touch $BASEDIR/.hic
