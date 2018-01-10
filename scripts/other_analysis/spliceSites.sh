#!/bin/bash

set -xe -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh

PC=$BASEDIR/../results/posConNCwithCodingPartner.bed
DIR=$BASEDIR/downstream/spliceSites
mkdir -p $DIR


# Print 5p exon
awk 'BEGIN{OFS=FS="\t"}{
	split($11,SIZES,",");
	if($6=="+"){
		print $1,$2,$2+SIZES[1],".",$5,$6
	}
	else if($6=="-"){
		print $1,$3-SIZES[$10],$3,".",$5,$6
	}
}' $PC | sort -k1,1 -k2,2n -k3,3n -k6,6 | uniq | awk '$3-$2>10' > $DIR/pcRNAs_5p_exon.bed

# Print 3p exon
awk 'BEGIN{OFS=FS="\t"}{
	split($11,SIZES,",");
	if($6=="+"){
		print $1,$3-SIZES[$10],$3,".",$5,$6
	}
	else if($6=="-"){
		print $1,$2,$2+SIZES[1],".",$5,$6
	}
}' $PC | sort -k1,1 -k2,2n -k3,3n -k6,6 | uniq | awk '$3-$2>10' > $DIR/pcRNAs_3p_exon.bed


awk 'BEGIN{OFS=FS="\t"}$10>2{
	split($11,SIZES,",");
	split($12,STARTS,",");
	START=$2+STARTS[2];
	ENDa=$2+STARTS[$10-1]+SIZES[$10-1];
	NEX=$10-2;
	NSIZES=""
	NSTARTS=""
	for(i=2;i<$10;i++){
		NSIZES=NSIZES""SIZES[i]","
		NSTARTS=NSTARTS""STARTS[i]-STARTS[2]","
	}
	print $1,START,ENDa,".",$5,$6,START,START,$9,NEX,NSIZES,NSTARTS;
}' $PC | $B12ToB6 | sort -k1,1 -k2,2n -k3,3n -k6,6 | uniq | awk '$3-$2>10' > $DIR/pcRNAs_internal_exons.bed




(cat $DIR/pcRNAs_5p_exon.bed
echo '#5p'
cat $DIR/pcRNAs_3p_exon.bed
echo '#3p'
cat $DIR/pcRNAs_internal_exons.bed
echo '#Internal')> $DIR/pcRNAs_all_exons_computeMatrixFormat.bed


mkdir -p $DIR/logs
bsub -q highpri -oo $DIR/logs/matrix.log -M 10000 -n 4 -R 'rusage[mem=10000]' computeMatrix scale-regions -b 1000 -a 1000 -m 10000 --regionsFileName $DIR/pcRNAs_all_exons_computeMatrixFormat.bed --sortRegions "no" --outFileNameMatrix $DIR/pcRNAs_spliceConsMatrix.txt --scoreFileName $DATA/hg38.phastCons100way.bw --missingDataAsZero --outFileName $DIR/matrix_out.gz -p 4
rm -f $DIR/matrix_out.gz

Rscript $BIN/dtCompare.R --files $DIR/pcRNAs_spliceConsMatrix.txt --labels Conservation --profileStdErr --profW 10 --profH 10 --outFile $DIR/spliceConservationPlot
convert -density 300 $DIR/spliceConservationPlot_profile.pdf -quality 100 $DIR/spliceConservationPlot_profile.png
convert -density 300 $DIR/spliceConservationPlot_heatmap.pdf -quality 100 $DIR/spliceConservationPlot_heatmap.png



