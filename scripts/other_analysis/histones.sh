#!/bin/bash
set -xe -o pipefail
shopt -s extglob
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh
export LC_ALL=C

P=4
TMP=/tmp/$LSB_JOBID
mkdir -p $TMP

# Download histone mod data
mkdir -p $BASEDIR/downstream/histones/data/hg19
CD=$BASEDIR/downstream/histones/
fetchChromSizes hg38 > $BASEDIR/data/hg38.chromSizes
# Download data
for i in $(awk '!/#/{print $1}' $DATA/histone_chipseq_data.txt); do
	NAME=$(basename $i .bw | sed 's/wgEncodeBroadHistone//' | sed -E 's/StdAln.+//')
	if [[ ! -f $CD/data/hg19/${NAME}.bw ]]; then
		wget -O - $i > $CD/data/hg19/${NAME}.bw
	fi;
done;


function convert () { 
	TMP=$CD/tmp/$LSB_JOBID
	mkdir -p $TMP
	NAME=$(basename $1 .bw);
	$BIGWIG_TO_BEDGRAPH $CD/data/hg19/${NAME}.bw $TMP/${NAME}.bedGraph.hg19
	$LIFTOVER -bedPlus=3 $TMP/${NAME}.bedGraph.hg19 $BASEDIR/data/hg19ToHg38.over.chain $TMP/${NAME}.bedGraph.hg38.unsorted $TMP/${NAME}.bedGraph.hg19_unmapped
	rm $TMP/${NAME}.bedGraph.hg19_unmapped $TMP/${NAME}.bedGraph.hg19
	sort -k1,1 -k2,2n -k3,3n -S 20G -T $TMP $TMP/${NAME}.bedGraph.hg38.unsorted > $TMP/${NAME}.bedGraph.hg38
	rm $TMP/${NAME}.bedGraph.hg38.unsorted
	$MAP -d -1 -c 4 -o median -i $TMP/${NAME}.bedGraph.hg38 > $TMP/${NAME}.bedGraph.hg38.merged
	rm $TMP/${NAME}.bedGraph.hg38
	$BEDGRAPH_TO_BIGWIG $TMP/${NAME}.bedGraph.hg38.merged $BASEDIR/data/hg38.chromSizes $CD/data/${NAME}.bw
	rm $TMP/${NAME}.bedGraph.hg38.merged
}
export CD
export LIFTOVER
export BASEDIR
export -f convert

mkdir -p $CD/data/logs
for i in $CD/data/hg19/*.bw; do
	NAME=$(basename $i .bw);
	if [[ ! -f $CD/data/${NAME}.bw ]]; then 
		rm -f $CD/data/logs/$(basename $i).log
		bsub -K -q research-rh6 -M 20000 -oo $CD/data/logs/$(basename $i).log -n 1 -R 'rusage[mem=20000]' convert $i &
	fi	
done;

wait



# Make matrices and heatmaps by position

mkdir -p $CD/annotationByPosition
# Prepare annot by class
for i in $(cut -f 15 $BASEDIR/../results/posConNCwithCodingPartner_andContext_expr.bedx | sort -u); do
	awk -v I=$i 'BEGIN{OFS=FS="\t"} $15==I{print $1,$2,$3,$4,$5,$6}'  $BASEDIR/../results/posConNCwithCodingPartner_andContext_expr.bedx > $CD/annotationByPosition/posConAnnot_Only$i.bed
	echo "#$i" >> $CD/annotationByPosition/posConAnnot_Only$i.bed
done

# Make annotation of 100 random Gencode lincs
$SAMPLE -n 1000 -seed 383847 -i $BASEDIR/data/Gencode.v21.lincRNAs.bed > $CD/annotationByPosition/posConAnnot_OnlyGencodeLinc.bed
echo "#GencodeLinc" >> $CD/annotationByPosition/posConAnnot_OnlyGencodeLinc.bed

# Make annotation of 100 random Gencode coding genes
$SAMPLE -n 1000 -seed 383847 -i $BASEDIR/coding/Gencode-current_hsa_coding.bed > $CD/annotationByPosition/posConAnnot_OnlyGencodeCoding.bed
echo "#GencodeCoding" >> $CD/annotationByPosition/posConAnnot_OnlyGencodeCoding.bed

# Merge all files together
cat $CD/annotationByPosition/posConAnnot_Only*.bed > $CD/annotationByPosition/posConAnnot.bed


# Compute matrices
mkdir -p $CD/annotationByPosition/matrix/logs
for i in $CD/data/*.bw; do
	NAME=$(basename $i .bw);
	if [[ ! -f $CD/annotationByPosition/matrix/${NAME}.txt ]]; then
		bsub -q research-rh6 -M 2000 -K -o $CD/annotationByPosition/matrix/logs/${NAME}.log -n 1 -R 'rusage[mem=2000]' computeMatrix reference-point --referencePoint TSS --regionsFileName $CD/annotationByPosition/posConAnnot.bed --sortRegions "no" --outFileNameMatrix $CD/annotationByPosition/matrix/${NAME}.txt --scoreFileName $i --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --binSize 10 --missingDataAsZero --outFileName /dev/null -p 1 &
	fi
done

wait

# Make heatmaps
TIS="Gm12878 H1hesc K562 Hsmm"
for CELL in $TIS; do
	LABELS="H3k4me3"
	SECLABELS="$CELL"
	FILES="$CD/annotationByPosition/matrix/${CELL}H3k4me3.txt"
	for FILE in $(ls $CD/annotationByPosition/matrix/$CELL!(*H3k4me3).txt); do
		MOD=$(basename $FILE .txt | sed "s/$CELL//" | sed "s/StdSig.*//" | sed "s/UcdSig.*//");
			LABELS="$LABELS $MOD"
			SECLABELS="$SECLABELS $CELL"
			FILES="$FILES $FILE"
	done
	FILEf=""
		Rscript $BIN/dtCompare.R --noHeat --files $FILES --labels $LABELS  --profileStdErr --heatH 14 --heatW 14 --profW 10 --profH 10 --sortByFirst --title $CELL --profileSplitLabels --outFile $CD/annotationByPosition/$CELL
done;


# As above, but keep al pcRNAs together
mkdir -p $CD/annotationAllPc
(cut -f1-6 $BASEDIR/../results/posConNCwithCodingPartner_andContext_expr.bedx 
echo "#pcRNAs"
# posConAnnot_OnlyGencodeCoding.bed and posConAnnot_OnlyGencodeLinc.bed already contain the "#" annotation at the end
cat $CD/annotationByPosition/posConAnnot_OnlyGencodeCoding.bed
cat $CD/annotationByPosition/posConAnnot_OnlyGencodeLinc.bed
)>$CD/annotationAllPc/posConAnnot.bed


mkdir -p $CD/annotationAllPc/matrix/logs
for i in $CD/data/*.bw; do
        NAME=$(basename $i .bw);
        if [[ ! -f $CD/annotationAllPc/matrix/${NAME}.txt ]]; then
                bsub -q research-rh6 -M 2000 -K -o $CD/annotationAllPc/matrix/logs/${NAME}.log -n 1 -R 'rusage[mem=2000]' computeMatrix reference-point --referencePoint TSS --regionsFileName $CD/annotationAllPc/posConAnnot.bed --sortRegions "no" --outFileNameMatrix $CD/annotationAllPc/matrix/${NAME}.txt --scoreFileName $i --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --binSize 10 --missingDataAsZero --outFileName /dev/null -p 1 &
        fi
done

wait

TIS="Gm12878 H1hesc K562 Hsmm"
for CELL in $TIS; do
        LABELS="H3k4me3"
        SECLABELS="$CELL"
        FILES="$CD/annotationAllPc/matrix/${CELL}H3k4me3.txt"
        for FILE in $(ls $CD/annotationAllPc/matrix/$CELL!(*H3k4me3).txt); do
                MOD=$(basename $FILE .txt | sed "s/$CELL//" | sed "s/StdSig.*//" | sed "s/UcdSig.*//");
                        LABELS="$LABELS $MOD"
                        SECLABELS="$SECLABELS $CELL"
                        FILES="$FILES $FILE"
        done
        FILEf=""
                Rscript $BIN/dtCompare.R --noHeat --files $FILES --labels $LABELS  --profileStdErr --heatH 14 --heatW 14 --profW 10 --profH 10 --sortByFirst --title $CELL --profileSplitLabels --outFile $CD/annotationAllPc/$CELL
done;



# Calculate H3K27me3 coverage
mkdir -p $CD/Coverage_H3K27me3
# Annotate pcRNA 1kb promoters
awk 'BEGIN{OFS=FS="\t"}{if($6=="+"){print $1,$2-1000,$2+1000,$4,$5,"."}if($6=="-"){print $1,$3-1000,$3+1000,$4,$5,"."}}' $BASEDIR/nc/posConNCwithCodingPartner.bed | sort -k1,1 -k2,2n -k3,3n > $CD/Coverage_H3K27me3/pcRNApromoters.bed


for i in $CD/data/*H3k27me3*; do
	BWNAME=$(basename $i .bw)
	mkdir -p $CD/data/bedgraphs/
	if [[ ! -f $CD/data/bedgraphs/$BWNAME.bedg ]]; then
		$BIGWIG_TO_BEDGRAPH $i $CD/data/bedgraphs/$BWNAME.bedg
	fi
	$MAP -a $CD/Coverage_H3K27me3/pcRNApromoters.bed -b $CD/data/bedgraphs/$BWNAME.bedg -c 4 -o mean -null 0 > $CD/Coverage_H3K27me3/${BWNAME}.map
done

for i in $CD/data/*H3k4me3*; do
        BWNAME=$(basename $i .bw)
        mkdir -p $CD/data/bedgraphs/
        if [[ ! -f $CD/data/bedgraphs/$BWNAME.bedg ]]; then
                $BIGWIG_TO_BEDGRAPH $i $CD/data/bedgraphs/$BWNAME.bedg
        fi
        $MAP -a $CD/Coverage_H3K27me3/pcRNApromoters.bed -b $CD/data/bedgraphs/$BWNAME.bedg -c 4 -o mean -null 0 > $CD/Coverage_H3K27me3/${BWNAME}.map
done

Rscript -e "source(\"$BASEDIR/../scripts/other_analysis/Coverage_H3K27me3.R\")"

touch $BASEDIR/.histones


